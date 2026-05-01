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
  * `stack_hists`    — StackHistogramSet diagnostic accumulator
                       (chain depth, first-interaction type, per-region
                       deposit counts, inclusive edep, Δz, E_first).
                       `nothing` if `with_stack_histograms=false`.
  * `cluster_hists`  — ClusterHistogramSet diagnostic accumulator
                       (per-cluster energy, multiplicity, pair distances,
                       r²/D vs z heatmaps).
                       `nothing` if `with_cluster_histograms=false`.
  * `rej_hist`       — RejectionHistograms (skin/FV early-reject diagnostics).
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
    stack_hists::Union{StackHistogramSet, Nothing}
    cluster_hists::Union{ClusterHistogramSet, Nothing}
    rej_hist::Union{RejectionHistograms, Nothing}
end

const _MC_OUTCOMES = (:escaped, :MS_rejected, :skin_vetoed,
                      :outside_FV, :SS_outside_ROI, :SS_in_ROI,
                      :companion_vetoed)

"""
    run_mc(det, eff, comp_eff, xcom, params, n_samples; mc_seed=1234,
            with_rejection_histograms=true, verbose=false) -> MCResult

Run the per-photon MC with `n_samples` photons sampled from `eff`,
parallelised across `Threads.nthreads()` threads. Each thread is seeded
with `mc_seed + thread_id` for reproducibility.

Per-event pipeline (see `src3/tracker.jl` for the algorithm):

    copy!(snapshot, rng)
    fv = fast_veto(rng, det, eff, xcom, params; rej_hist)
    if fv === :pass
        copy!(rng, snapshot)                    # restore for the same event
        empty!(stack)
        status = track_photon_stack(rng, det, eff, xcom, params, stack)
        clusters = build_clusters(rng, stack, params)
        outcome = classify_event(status, stack, clusters, params)
    elseif fv === :vetoed_skin: outcome = :skin_vetoed
    else (:rejected_fv):       outcome = :outside_FV

If `eff.isotope === :Tl208` and `comp_eff` is provided, the companion
veto is applied: each `:SS_in_ROI` candidate is converted to
`:companion_vetoed` with probability
`r_comp · P(companion produces visible deposit)`, where
`r_comp = companion_reach_prob(comp_eff)`. For Bi-214 sources, pass
`comp_eff = nothing`.

The kwarg `with_histograms` is reserved for stack-based control
histograms (currently a no-op; `MCResult.histograms` is always `nothing`).
"""
function run_mc(det::LXeDetector, eff::EffectiveSource,
                comp_eff::Union{EffectiveSource, Nothing},
                xcom::XCOMTable, params::MCParams,
                n_samples::Integer; mc_seed::Integer=1234,
                verbose::Bool=false,
                with_stack_histograms::Bool=true,
                with_cluster_histograms::Bool=true,
                with_rejection_histograms::Bool=true,
                sample_stack::Int=0,
                sample_stack_path::Union{AbstractString, Nothing}=nothing
                )::MCResult
    n_threads  = Threads.nthreads()
    base       = div(n_samples, n_threads)
    rem        = n_samples - base * n_threads

    apply_veto = comp_eff !== nothing && eff.isotope === :Tl208
    r_comp     = apply_veto ? companion_reach_prob(comp_eff) : 0.0

    thread_counts = [Dict{Symbol,Int}(o => 0 for o in _MC_OUTCOMES)
                     for _ in 1:n_threads]
    thread_rej_hist = with_rejection_histograms ?
                      [RejectionHistograms(r2_max_cm2 = det.R_ICV_inner^2,
                                            z_min_cm   = det.z_LXe_bottom,
                                            z_max_cm   = det.z_gate)
                       for _ in 1:n_threads] :
                      Vector{RejectionHistograms}()
    thread_stack_hist = with_stack_histograms ?
                        [StackHistogramSet() for _ in 1:n_threads] :
                        Vector{StackHistogramSet}()
    thread_cluster_hist = with_cluster_histograms ?
                          [ClusterHistogramSet(
                              r2_max_cm2 = det.R_ICV_inner^2,
                              z_min_cm   = det.z_LXe_bottom,
                              z_max_cm   = det.z_gate)
                           for _ in 1:n_threads] :
                          Vector{ClusterHistogramSet}()
    thread_stack    = [PhotonStack()         for _ in 1:n_threads]
    thread_snapshot = [MersenneTwister(0)    for _ in 1:n_threads]

    # --sample-stack: per-thread temp file paths (concatenated post-loop).
    do_sample = sample_stack != 0 && sample_stack_path !== nothing
    sample_tmp_paths = if do_sample
        mkpath(dirname(sample_stack_path))
        ["$(sample_stack_path).tid$tid.tmp" for tid in 1:n_threads]
    else
        String[]
    end

    base_thread1 = base + (1 <= rem ? 1 : 0)
    report_every = max(1, base_thread1 ÷ 10)

    t0 = time()
    Threads.@threads for tid in 1:n_threads
        n_local      = base + (tid <= rem ? 1 : 0)
        rng          = MersenneTwister(mc_seed + tid)
        local_counts = thread_counts[tid]
        local_rej    = with_rejection_histograms ? thread_rej_hist[tid] : nothing
        local_stack_hist   = with_stack_histograms   ? thread_stack_hist[tid]   : nothing
        local_cluster_hist = with_cluster_histograms ? thread_cluster_hist[tid] : nothing
        local_stack  = thread_stack[tid]
        local_snap   = thread_snapshot[tid]

        # Per-thread stack-sample writer.
        # `sample_stack`:
        #   = -1   write every event with non-empty stack.
        #   > 0    write the FIRST max(1, cld(sample_stack, n_threads))
        #          events per thread that have a non-empty stack —
        #          gives ~sample_stack total events across all threads.
        local_sample_io = nothing
        local_sample_target = 0
        local_sample_written = 0
        if do_sample
            local_sample_io = open(sample_tmp_paths[tid], "w")
            local_sample_target = sample_stack == -1 ?
                                   typemax(Int) :
                                   max(1, cld(sample_stack, n_threads))
        end
        for i in 1:n_local
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
                local_stack_hist   !== nothing &&
                    update_stack_histograms!(local_stack_hist, local_stack, params)
                local_cluster_hist !== nothing &&
                    update_cluster_histograms!(local_cluster_hist, clusters, params)
            elseif fv === :vetoed_skin
                outcome = :skin_vetoed
            else  # :rejected_fv
                outcome = :outside_FV
            end
            if apply_veto && outcome === :SS_in_ROI
                if rand(rng) < r_comp
                    if companion_visible!(rng, det, comp_eff, xcom, params)
                        outcome = :companion_vetoed
                    end
                end
            end
            local_counts[outcome] += 1

            # Optional stack-sample dump: write the first `local_sample_target`
            # events per thread with non-empty stacks.
            if local_sample_io !== nothing && length(local_stack) > 0 &&
               local_sample_written < local_sample_target
                @inbounds for r in local_stack.rows
                    println(local_sample_io,
                            "$(tid),$(i),$(r.ng),$(r.nm),",
                            "$(r.parent_region),$(r.region),$(r.interaction),",
                            "$(r.x),$(r.y),$(r.z),$(r.epre),$(r.edep)")
                end
                local_sample_written += 1
            end

            if verbose && tid == 1 && (i % report_every == 0)
                est_done = i * n_threads
                elapsed  = time() - t0
                @printf("    progress %5.1f %%  (%d / %d)  elapsed %.1f s\r",
                        100.0 * est_done / n_samples,
                        est_done, n_samples, elapsed)
            end
        end
        local_sample_io !== nothing && close(local_sample_io)
    end
    runtime = time() - t0
    if verbose
        @printf("    progress 100.0 %%  (%d / %d)  elapsed %.1f s\n",
                n_samples, n_samples, runtime)
    end

    # Concatenate per-thread sample-stack temps into the final CSV.
    if do_sample
        open(sample_stack_path, "w") do f
            println(f, "tid,event_idx,ng,nm,parent_region,region,",
                       "interaction,x,y,z,epre,edep")
            for tmp in sample_tmp_paths
                isfile(tmp) || continue
                for line in eachline(tmp)
                    println(f, line)
                end
                rm(tmp)
            end
        end
    end

    counts = Dict{Symbol,Int}(
        o => sum(thread_counts[i][o] for i in 1:n_threads)
        for o in _MC_OUTCOMES
    )
    n_total  = sum(values(counts))
    f_ss_roi = counts[:SS_in_ROI] / n_total
    bg       = f_ss_roi * eff.total_per_yr

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

    merged_stack_hist = if with_stack_histograms
        sh = StackHistogramSet()
        for th in thread_stack_hist
            merge_stack_histograms!(sh, th)
        end
        sh
    else
        nothing
    end

    merged_cluster_hist = if with_cluster_histograms
        ch = ClusterHistogramSet(r2_max_cm2 = det.R_ICV_inner^2,
                                  z_min_cm   = det.z_LXe_bottom,
                                  z_max_cm   = det.z_gate)
        for tc in thread_cluster_hist
            merge_cluster_histograms!(ch, tc)
        end
        ch
    else
        nothing
    end

    MCResult(eff.name, eff.isotope, counts, n_total, runtime,
             eff.total_per_yr, f_ss_roi, bg, r_comp,
             merged_stack_hist, merged_cluster_hist, merged_rej)
end

"""
    run_mc_all(det, effs, xcom, params, n_samples; mc_seed=1234) -> Vector{MCResult}

Run the MC for the 6 *main* effective sources (Bi-214 and Tl-208 main,
in regions CB / CTH / CBH). Tl-208 sources are automatically paired with
their `*_Tl208c` companion source for the cascade-companion veto.
"""
function run_mc_all(det::LXeDetector, effs::Vector{EffectiveSource},
                    xcom::XCOMTable, params::MCParams,
                    n_samples::Integer;
                    mc_seed::Integer=1234,
                    verbose::Bool=false,
                    with_stack_histograms::Bool=true,
                    with_cluster_histograms::Bool=true,
                    with_rejection_histograms::Bool=true,
                    sample_stack::Int=0,
                    sample_stack_dir::Union{AbstractString, Nothing}=nothing
                    )::Vector{MCResult}
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
        sample_path = if sample_stack != 0 && sample_stack_dir !== nothing
            joinpath(sample_stack_dir, "hist_$(eff.name)", "stack_sample.csv")
        else
            nothing
        end
        push!(results, run_mc(det, eff, comp_eff, xcom, params, n_samples;
                              mc_seed=seed, verbose=verbose,
                              with_stack_histograms=with_stack_histograms,
                              with_cluster_histograms=with_cluster_histograms,
                              with_rejection_histograms=with_rejection_histograms,
                              sample_stack=sample_stack,
                              sample_stack_path=sample_path))
    end
    results
end
