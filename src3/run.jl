# src3/run.jl — Threaded driver and result struct for the per-photon MC.

"""
    MCResult

Outcome of one MC run for one EffectiveSource.

Fields:
  * `name`           — source name (e.g. "CB_Bi214")
  * `isotope`        — `:Bi214`, `:Tl208`, `:Tl208c`
  * `counts`         — Dict{Symbol,Int}, one entry per outcome category
  * `n_total`        — total photons sampled
  * `runtime_s`      — wall-clock seconds
  * `γ_per_yr_total` — `eff.total_per_yr` (γ/yr at ICV inner exit surface)
  * `f_SS_in_ROI`    — `counts[:SS_in_ROI] / n_total`
  * `bg_per_yr`      — `γ_per_yr_total × f_SS_in_ROI` (events/yr in ROI)
  * `r_comp`         — companion-reach probability used (0 for non-Tl208)
"""
struct MCResult
    name::String
    isotope::Symbol
    counts::Dict{Symbol, Int}
    n_total::Int
    runtime_s::Float64
    γ_per_yr_total::Float64
    f_SS_in_ROI::Float64
    bg_per_yr::Float64
    r_comp::Float64
    histograms::Union{HistogramSet, Nothing}
    rej_hist::Union{RejectionHistograms, Nothing}
end

const _MC_OUTCOMES = (:escaped, :MS_rejected, :skin_vetoed,
                      :outside_FV, :SS_outside_ROI, :SS_in_ROI,
                      :companion_vetoed)

"""
    run_mc(det, eff, comp_eff, xcom, params, n_samples; mc_seed=1234,
            use_stack_tracker=false, ...) -> MCResult

Run the per-photon MC with `n_samples` photons sampled from `eff`,
parallelised across `Threads.nthreads()` threads. Each thread is seeded
with `mc_seed + thread_id` for reproducibility.

If `eff.isotope === :Tl208` and `comp_eff` is provided, the companion
veto is applied: each `:SS_in_ROI` candidate is converted to
`:companion_vetoed` with probability
`r_comp · P(companion produces visible deposit)`, where
`r_comp = companion_reach_prob(comp_eff)`. For Bi-214 sources, pass
`comp_eff = nothing`.

Tracker selection (kwarg `use_stack_tracker`):
  * `false` (default): legacy `track_one_photon!` path. `with_histograms`
    fills the control `HistogramSet`; `early_skin_reject` and
    `early_fv_reject` toggle the legacy in-tracking reject paths.
  * `true`: new pipeline. Per event:
      copy!(snapshot, rng);  fv = fast_veto(...; rej_hist)
      if fv === :pass
          copy!(rng, snapshot)
          status = track_photon_stack(rng, ..., stack)
          clusters = build_clusters(rng, stack, params)
          outcome  = classify_event(status, stack, clusters, params)
      else
          outcome = (fv === :vetoed_skin) ? :skin_vetoed : :outside_FV
      end
    `with_histograms` is silently ignored under this path —
    control-histogram re-sourcing from the stack is a later step.
    `early_skin_reject` / `early_fv_reject` are also ignored (the new
    pipeline's veto/reject is in `fast_veto` + `classify_event`).
"""
function run_mc(det::LXeDetector, eff::EffectiveSource,
                comp_eff::Union{EffectiveSource, Nothing},
                xcom::XCOMTable, params::MCParams,
                n_samples::Integer; mc_seed::Integer=1234,
                verbose::Bool=false,
                with_histograms::Bool=false,
                early_skin_reject::Bool=true,
                early_fv_reject::Bool=true,
                with_rejection_histograms::Bool=true,
                use_stack_tracker::Bool=false)::MCResult
    n_threads  = Threads.nthreads()
    base       = div(n_samples, n_threads)
    rem        = n_samples - base * n_threads

    apply_veto = comp_eff !== nothing && eff.isotope === :Tl208
    r_comp     = apply_veto ? companion_reach_prob(comp_eff) : 0.0

    # Histograms are not yet wired to the new tracker; force off if requested.
    fill_hists = with_histograms && !use_stack_tracker

    thread_counts = [Dict{Symbol,Int}(o => 0 for o in _MC_OUTCOMES)
                     for _ in 1:n_threads]
    thread_hists  = fill_hists ?
                    [HistogramSet() for _ in 1:n_threads] :
                    Vector{HistogramSet}()
    thread_scratch = [PhotonScratch() for _ in 1:n_threads]
    thread_rej_hist = with_rejection_histograms ?
                      [RejectionHistograms(r2_max_cm2 = det.R_ICV_inner^2,
                                            z_min_cm   = det.z_LXe_bottom,
                                            z_max_cm   = det.z_gate)
                       for _ in 1:n_threads] :
                      Vector{RejectionHistograms}()
    # New-tracker per-thread state (allocated only when needed).
    thread_stack     = use_stack_tracker ?
                       [PhotonStack() for _ in 1:n_threads] :
                       Vector{PhotonStack}()
    thread_snapshot  = use_stack_tracker ?
                       [MersenneTwister(0) for _ in 1:n_threads] :
                       Vector{MersenneTwister}()

    base_thread1 = base + (1 <= rem ? 1 : 0)
    report_every = max(1, base_thread1 ÷ 10)

    t0 = time()
    Threads.@threads for tid in 1:n_threads
        n_local      = base + (tid <= rem ? 1 : 0)
        rng          = MersenneTwister(mc_seed + tid)
        local_counts  = thread_counts[tid]
        local_hist    = fill_hists ? thread_hists[tid] : nothing
        local_scratch = thread_scratch[tid]
        local_rej     = with_rejection_histograms ? thread_rej_hist[tid] : nothing
        local_stack   = use_stack_tracker ? thread_stack[tid]    : PhotonStack()
        local_snap    = use_stack_tracker ? thread_snapshot[tid] : MersenneTwister(0)
        for i in 1:n_local
            if use_stack_tracker
                # New pipeline: fast_veto → (track + classify) on pass.
                copy!(local_snap, rng)
                fv = fast_veto(rng, det, eff, xcom, params; rej_hist=local_rej)
                if fv === :pass
                    copy!(rng, local_snap)
                    empty!(local_stack)
                    status   = track_photon_stack(rng, det, eff, xcom,
                                                  params, local_stack)
                    clusters = build_clusters(rng, local_stack, params)
                    outcome  = classify_event(status, local_stack,
                                              clusters, params)
                elseif fv === :vetoed_skin
                    outcome = :skin_vetoed
                else  # :rejected_fv
                    outcome = :outside_FV
                end
            else
                # Legacy pipeline.
                outcome = track_one_photon!(rng, det, eff, xcom, params;
                                             hist=local_hist,
                                             scratch=local_scratch,
                                             rej_hist=local_rej,
                                             early_skin_reject=early_skin_reject,
                                             early_fv_reject=early_fv_reject)
            end
            if apply_veto && outcome === :SS_in_ROI
                if rand(rng) < r_comp
                    if companion_visible!(rng, det, comp_eff, xcom, params)
                        outcome = :companion_vetoed
                    end
                end
            end
            local_counts[outcome] += 1
            if verbose && tid == 1 && (i % report_every == 0)
                est_done = i * n_threads
                elapsed  = time() - t0
                @printf("    progress %5.1f %%  (%d / %d)  elapsed %.1f s\r",
                        100.0 * est_done / n_samples,
                        est_done, n_samples, elapsed)
            end
        end
    end
    runtime = time() - t0
    if verbose
        @printf("    progress 100.0 %%  (%d / %d)  elapsed %.1f s\n",
                n_samples, n_samples, runtime)
    end

    counts = Dict{Symbol,Int}(
        o => sum(thread_counts[i][o] for i in 1:n_threads)
        for o in _MC_OUTCOMES
    )
    n_total  = sum(values(counts))
    f_ss_roi = counts[:SS_in_ROI] / n_total
    bg       = f_ss_roi * eff.total_per_yr

    merged_hist = if fill_hists
        h = HistogramSet()
        for th in thread_hists
            merge_histograms!(h, th)
        end
        h
    else
        nothing
    end

    merged_rej = if with_rejection_histograms
        rh = RejectionHistograms(r2_max_cm2 = det.R_ICV_inner^2,
                                  z_min_cm   = det.z_LXe_bottom,
                                  z_max_cm   = det.z_gate)
        for tr in thread_rej_hist
            merge_rejection_histograms!(rh, tr)
        end
        rh
    else
        nothing
    end

    MCResult(eff.name, eff.isotope, counts, n_total, runtime,
             eff.total_per_yr, f_ss_roi, bg, r_comp, merged_hist, merged_rej)
end

"""
    run_mc_all(det, effs, xcom, params, n_samples; mc_seed=1234) -> Vector{MCResult}

Run the MC for the 6 *main* effective sources (Bi-214 and Tl-208 main,
in regions CB / CTH / CBH) and skip the Tl-208 companion sources
(those are paired automatically with their main partner via the name
suffix). Returns one `MCResult` per main source.
"""
function run_mc_all(det::LXeDetector, effs::Vector{EffectiveSource},
                    xcom::XCOMTable, params::MCParams,
                    n_samples::Integer;
                    mc_seed::Integer=1234,
                    verbose::Bool=false,
                    with_histograms::Bool=false,
                    early_skin_reject::Bool=true,
                    early_fv_reject::Bool=true,
                    with_rejection_histograms::Bool=true,
                    use_stack_tracker::Bool=false)::Vector{MCResult}
    by_name = Dict(e.name => e for e in effs)
    results = MCResult[]
    main_names = ["CB_Bi214", "CTH_Bi214", "CBH_Bi214",
                  "CB_Tl208", "CTH_Tl208", "CBH_Tl208"]
    for (i, mname) in enumerate(main_names)
        haskey(by_name, mname) || continue
        eff = by_name[mname]
        comp_eff = if eff.isotope === :Tl208
            cname = replace(mname, "_Tl208" => "_Tl208c")
            get(by_name, cname, nothing)
        else
            nothing
        end
        seed = mc_seed + (i - 1) * 1_000_000
        verbose && @printf("\n── Running %s (%d / %d) ──\n",
                            mname, i, length(main_names))
        push!(results, run_mc(det, eff, comp_eff, xcom, params, n_samples;
                              mc_seed=seed, verbose=verbose,
                              with_histograms=with_histograms,
                              early_skin_reject=early_skin_reject,
                              early_fv_reject=early_fv_reject,
                              with_rejection_histograms=with_rejection_histograms,
                              use_stack_tracker=use_stack_tracker))
    end
    results
end
