# test/runtests3.jl — Aggregates src3/ tests, each in its own Julia process.
#
# Subprocess-per-file is required because:
#   1) `using XLZD` (registered package, src/) and `using .XLZD3` clash on
#      shared exports like `load_xcom`, so they cannot share a Main namespace.
#   2) Individual src3/ test files declare module-level constants
#      (e.g. `const ti_path = ...`) which collide if multiple test files
#      are `include`d into the same session.
#
# Run from project root:
#     julia --project=. -t 4 test/runtests3.jl

const TEST_FILES = [
    "test_geometry3.jl",
    "test_physics3.jl",
    "test_material3.jl",
    "test_pobjects3.jl",
    "test_pannulardisk3.jl",
    "test_field_cage3.jl",
    "test_sources3.jl",
    "test_transport3.jl",
    "test_effective_sources3.jl",
    "test_lxe_detector3.jl",
    "test_sampling3.jl",
    "test_stack3.jl",
    "test_build_clusters3.jl",
    "test_select3.jl",
    "test_classify_event3.jl",
    "test_tracker3.jl",
    "test_histograms3.jl",
    "test_stack_histograms3.jl",
    "test_cluster_histograms3.jl",
    "test_mc3.jl",
]

const TEST_DIR = @__DIR__
const PROJECT_ROOT = abspath(joinpath(TEST_DIR, ".."))
const NTHREADS = get(ENV, "JULIA_NUM_THREADS", "1")

failures = String[]
total_t  = 0.0

for f in TEST_FILES
    println("\n==================== $f ====================")
    t0 = time()
    cmd = `$(Base.julia_cmd()) --project=$PROJECT_ROOT -t $NTHREADS $(joinpath(TEST_DIR, f))`
    ok = success(run(ignorestatus(cmd)))
    dt = time() - t0
    global total_t += dt
    if !ok
        push!(failures, f)
        println("  ❌ $f FAILED  ($(round(dt, digits=1)) s)")
    else
        println("  ✓  $f passed   ($(round(dt, digits=1)) s)")
    end
end

println("\n──────── runtests3.jl summary ────────")
println("  files run     : $(length(TEST_FILES))")
println("  total time    : $(round(total_t, digits=1)) s")
if isempty(failures)
    println("  status        : ALL GREEN")
else
    println("  status        : $(length(failures)) FAILURES")
    foreach(f -> println("    - $f"), failures)
    exit(1)
end
