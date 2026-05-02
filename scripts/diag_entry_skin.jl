# scripts/diag_entry_skin.jl
#
# Diagnostic for the per-source per-outcome breakdown — answers two questions:
#
#   1. Of escaped events, what's the entry-direction distribution?
#      Entry direction is parameterised by u = cos θ from the inward
#      surface normal at the entry point (u=1 head-on inward, u=0 tangential).
#      Hypothesis: escapes are concentrated at small u (grazing chord through LXe).
#
#   2. Of skin-vetoed events, what energy actually deposits in :Skin?
#      Skin veto fires when the MAIN γ deposits > E_skin_veto in :Skin.
#      Hypothesis (to verify): the depositing energy belongs to the main γ,
#      not the cascade companion (which is sampled separately and only
#      gates :SS_in_ROI events).
#
# Usage:
#   julia --project=. -t 1 scripts/diag_entry_skin.jl                # default CB_Tl208 N=50000
#   julia --project=. -t 1 scripts/diag_entry_skin.jl --source CB_Bi214 --n-samples 100000

using Random, Statistics, Printf
using ArgParse

include(joinpath(@__DIR__, "..", "src3", "XLZD3.jl"))
using .XLZD3

const PROJECT_ROOT = abspath(joinpath(@__DIR__, ".."))
const TI_PATH      = joinpath(PROJECT_ROOT, "data", "nist_ti.csv")
const LXE_NIST     = joinpath(PROJECT_ROOT, "data", "nist_lxe.csv")
const XCOM_PATH    = joinpath(PROJECT_ROOT, "data", "nist.csv")
const LXE_CSV      = joinpath(PROJECT_ROOT, "data", "lxe_detector.csv")
const GEOM_CSV     = joinpath(PROJECT_ROOT, "data", "lz_cryo_geometry.csv")
const EXTRAS_CSV   = joinpath(PROJECT_ROOT, "data", "lz_cryo_extras.csv")
const SURFACES_CSV = joinpath(PROJECT_ROOT, "data", "lz_cryo_surface_sources.csv")

function parse_cli()
    s = ArgParseSettings(description="diag: entry direction + skin breakdown")
    @add_arg_table! s begin
        "--source"
            default = "CB_Tl208"
            help    = "Source name to analyse"
        "--n-samples"
            arg_type = Int
            default  = 50_000
        "--seed"
            arg_type = Int
            default  = 1234
        "--n-print"
            arg_type = Int
            default  = 5
            help     = "How many full skin-vetoed event tracks to print"
    end
    parse_args(s)
end

# Inward normal at the source entry point (depends on source region).
function inward_normal(eff_region::Symbol, x, y, z, det)
    if eff_region === :barrel
        r = sqrt(x*x + y*y)
        return (-x/r, -y/r, 0.0)
    elseif eff_region === :endcap_top
        # Approximate: inward normal at top dome points downward.
        return (0.0, 0.0, -1.0)
    elseif eff_region === :endcap_bottom
        return (0.0, 0.0, +1.0)
    else
        return (0.0, 0.0, 0.0)
    end
end

function main()
    args = parse_cli()

    mat_LXe = load_material("LXe", 2.953, LXE_NIST)
    mat_Ti  = load_material("Ti",  4.510, TI_PATH)
    det     = build_lxe_detector(LXE_CSV, mat_LXe)
    cryo    = build_cryostat(GEOM_CSV, EXTRAS_CSV, SURFACES_CSV)
    indiv   = build_individual_sources(cryo, mat_Ti)
    effs    = build_effective_sources(indiv, cryo, mat_Ti)
    xcom    = load_xcom(XCOM_PATH)
    params  = MCParams()
    by_name = Dict(e.name => e for e in effs)

    haskey(by_name, args["source"]) || error("unknown source: $(args["source"])")
    eff = by_name[args["source"]]

    rng   = MersenneTwister(args["seed"])
    snap  = MersenneTwister(0)
    stack = PhotonStack()
    N     = args["n-samples"]

    u_escape = Float64[]; u_skin = Float64[]
    u_outFV  = Float64[]; u_other = Float64[]

    for _ in 1:N
        # Sample direction → compute u, then restore rng.
        copy!(snap, rng)
        x, y, z, dx, dy, dz = sample_entry(rng, det, eff)
        nx, ny, nz = inward_normal(eff.region, x, y, z, det)
        u = dx*nx + dy*ny + dz*nz

        copy!(rng, snap)
        fv = fast_veto(rng, det, eff, xcom, params)

        if fv === :vetoed_skin
            push!(u_skin, u)
        elseif fv === :rejected_fv
            push!(u_outFV, u)
        else  # :pass
            copy!(rng, snap)
            empty!(stack)
            status = track_photon_stack(rng, det, eff, xcom, params, stack)
            if status === :escaped
                push!(u_escape, u)
            else
                push!(u_other, u)
            end
        end
    end

    pct_lt(v, t) = isempty(v) ? 0.0 : 100.0 * sum(x -> x < t, v) / length(v)
    pct_gt(v, t) = isempty(v) ? 0.0 : 100.0 * sum(x -> x > t, v) / length(v)

    println("\n── Source: $(args["source"]), N = $N ──")
    @printf("  N escaped:      %d  (%.1f %%)\n", length(u_escape), 100*length(u_escape)/N)
    @printf("  N skin-vetoed:  %d  (%.1f %%)\n", length(u_skin),   100*length(u_skin)/N)
    @printf("  N outFV-vetoed: %d  (%.1f %%)\n", length(u_outFV),  100*length(u_outFV)/N)
    @printf("  N other:        %d  (%.1f %%)\n", length(u_other),  100*length(u_other)/N)

    println("\n  u = cos θ_inward (1 = head-on, 0 = tangential)")
    @printf("  %-12s  %6s  %6s  %14s  %14s\n",
            "category", "mean", "median", "%(u<0.3)tan", "%(u>0.7)head-on")
    for (label, v) in (("escaped",     u_escape),
                       ("skin-vetoed", u_skin),
                       ("outFV-vet.",  u_outFV),
                       ("other",       u_other))
        isempty(v) && continue
        @printf("  %-12s  %6.3f  %6.3f  %14.1f  %14.1f\n",
                label, mean(v), median(v), pct_lt(v, 0.3), pct_gt(v, 0.7))
    end

    # Track the first few skin-vetoed events with FULL track_photon_stack
    # (instead of bailing in fast_veto) so we can see where the energy
    # actually goes through the photon's cascade.
    println("\n── First $(args["n-print"]) skin-vetoed events: full cascade ──")
    println("  fast_veto fires on the FIRST :Skin deposit > E_skin_veto from")
    println("  the MAIN γ. The companion γ is a separate sampling and only")
    println("  veto's :SS_in_ROI events at the END — it is NOT what triggers")
    println("  these skin rejections.")

    rng2 = MersenneTwister(args["seed"])
    snap2 = MersenneTwister(0)
    n_done = 0
    i = 0
    while n_done < args["n-print"] && i < 2000
        i += 1
        copy!(snap2, rng2)
        fv = fast_veto(rng2, det, eff, xcom, params)
        if fv === :vetoed_skin
            copy!(rng2, snap2)
            empty!(stack)
            track_photon_stack(rng2, det, eff, xcom, params, stack)
            n_done += 1
            skin_E = sum((r.edep for r in stack.rows if r.region === :Skin); init=0.0)
            tpc_E  = sum((r.edep for r in stack.rows if r.region === :TPC);  init=0.0)
            @printf("\n  Event #%d (event idx %d in stream):\n", n_done, i)
            @printf("    stack rows  = %d\n", length(stack))
            @printf("    Σ edep Skin = %7.4f MeV (%.1f keV; threshold = %.1f keV)\n",
                    skin_E, skin_E*1000, params.E_skin_veto_keV)
            @printf("    Σ edep TPC  = %7.4f MeV\n", tpc_E)
            @printf("    %-3s %-3s %-8s %-12s %-8s   pos (cm)\n",
                    "ng", "nm", "region", "interaction", "edep")
            for r in stack.rows
                @printf("    %-3d %-3d %-8s %-12s %.4f   (% 6.1f, % 6.1f, % 6.1f)\n",
                        r.ng, r.nm, r.region, r.interaction, r.edep, r.x, r.y, r.z)
            end
        end
    end
end

main()
