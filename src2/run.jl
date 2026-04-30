# src2/run.jl — Threaded driver and result struct for the per-photon MC.

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
end

const _MC_OUTCOMES = (:escaped, :MS_rejected, :skin_vetoed,
                      :SS_outside_FV, :SS_outside_ROI, :SS_in_ROI,
                      :companion_vetoed)

"""
    run_mc(det, eff, comp_eff, xcom, params, n_samples; mc_seed=1234) -> MCResult

Run the per-photon MC with `n_samples` photons sampled from `eff`,
parallelised across `Threads.nthreads()` threads. Each thread is seeded
with `mc_seed + thread_id` for reproducibility.

If `eff.isotope === :Tl208` and `comp_eff` is provided, the companion
veto is applied: each `:SS_in_ROI` candidate is converted to
`:companion_vetoed` with probability
`r_comp · P(companion produces visible deposit)`, where
`r_comp = companion_reach_prob(comp_eff)`. For Bi-214 sources, pass
`comp_eff = nothing`.
"""
function run_mc(det::LXeDetector, eff::EffectiveSource,
                comp_eff::Union{EffectiveSource, Nothing},
                xcom::XCOMTable, params::MCParams,
                n_samples::Integer; mc_seed::Integer=1234)::MCResult
    n_threads  = Threads.nthreads()
    base       = div(n_samples, n_threads)
    rem        = n_samples - base * n_threads

    apply_veto = comp_eff !== nothing && eff.isotope === :Tl208
    r_comp     = apply_veto ? companion_reach_prob(comp_eff) : 0.0

    thread_counts = [Dict{Symbol,Int}(o => 0 for o in _MC_OUTCOMES)
                     for _ in 1:n_threads]

    t0 = time()
    Threads.@threads for tid in 1:n_threads
        n_local      = base + (tid <= rem ? 1 : 0)
        rng          = MersenneTwister(mc_seed + tid)
        local_counts = thread_counts[tid]
        for _ in 1:n_local
            outcome = track_one_photon!(rng, det, eff, xcom, params)
            if apply_veto && outcome === :SS_in_ROI
                if rand(rng) < r_comp
                    if companion_visible!(rng, det, comp_eff, xcom, params)
                        outcome = :companion_vetoed
                    end
                end
            end
            local_counts[outcome] += 1
        end
    end
    runtime = time() - t0

    counts = Dict{Symbol,Int}(
        o => sum(thread_counts[i][o] for i in 1:n_threads)
        for o in _MC_OUTCOMES
    )
    n_total  = sum(values(counts))
    f_ss_roi = counts[:SS_in_ROI] / n_total
    bg       = f_ss_roi * eff.total_per_yr

    MCResult(eff.name, eff.isotope, counts, n_total, runtime,
             eff.total_per_yr, f_ss_roi, bg, r_comp)
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
                    mc_seed::Integer=1234)::Vector{MCResult}
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
        push!(results, run_mc(det, eff, comp_eff, xcom, params, n_samples;
                              mc_seed=seed))
    end
    results
end
