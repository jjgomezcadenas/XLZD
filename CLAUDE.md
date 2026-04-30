# CLAUDE.md — XLZD Field-Cage Ring Background MC

## What this project is

A Julia Monte Carlo that estimates the contribution of Bi-214 (²³⁸U-late chain, 2.448 MeV γ-line) decays in the copper field-cage rings of an XLZD-like detector to the 0νββ background budget at Q_ββ = 2458 keV in ¹³⁶Xe.

The MC implements full Compton tracking with Klein-Nishina kinematics, an XLZD-style single-site/multi-site classification using Δz clustering, and the experiment's ROI energy window. The output is the events-per-mBq calibration that feeds the comparison with plated ²²²Rn-daughter contributions.

Multithreaded with `Threads.@threads`. Runs ~10⁸ samples in under a minute on 8 cores.

## How to run

```bash
# default run (10⁸ samples, 2 μBq/kg copper, full XLZD geometry)
julia --project=. -t 8 scripts/run_mc.jl

# with overrides
julia --project=. -t 8 scripts/run_mc.jl --n-samples 1e7 --act-cu 1.0e-6 --output output/run1/
```

Tests:

```bash
julia --project=. -t 4 test/runtests.jl
```

Outputs (in `output/`): `results.h5`, `trajectories.json`, `viewer.html`, `summary.txt`, four PNG plots.

## File tree

```
.
├── CLAUDE.md                # this file
├── README.md                # human-facing
├── Project.toml             # Julia deps
├── src/
│   ├── XLZD.jl              # module root
│   ├── geometry.jl          # rings, masses, cylinder intersection
│   ├── physics.jl           # cross sections, Klein-Nishina, rotation
│   ├── trajectory.jl        # trajectory recording for viewer
│   ├── mc.jl                # MC loop, threading
│   └── output.jl            # HDF5, JSON, plots
├── scripts/run_mc.jl        # entry point
├── viewer/viewer.html       # Three.js viewer template
├── data/nist.csv            # NIST XCOM table (provided)
├── output/                  # gitignored, runtime artifacts
├── test/                    # @testset, geometry/physics/MC checks
└── planning/
    └── xlzd_code.md              # full implementation specification (read for non-trivial changes)
```

## Where the spec lives

**`planning/xlzd_code.md` is the authoritative implementation specification.** It contains:

- Physics derivations and parameter choices (NEXT→XLZD scaling, the 0.97/2.95 events ratio reasoning, the calibration to PTFE events/mBq, why we skip pair production, why z-only clustering)
- The complete `Params` struct definition with all defaults
- The MC algorithm pseudocode
- The HDF5 output schema
- The Three.js viewer scene + interaction spec
- The trajectory JSON schema
- The build order and tests

**Read `planning/code.md` before:**
- Adding or removing fields from `Params`
- Changing the MC algorithm (Klein-Nishina, SS/MS rules, ROI definition)
- Changing the HDF5 layout or trajectory JSON schema
- Adding new outcome categories
- Modifying the cross-section interpolation
- Changing the threading strategy

## Conventions

- **`Params` is immutable** (`@kwdef struct`). Add new fields with appropriate `geom_*`, `act_*`, `phys_*`, `cut_*`, `fv_*`, `mc_*`, `hist_*`, or `out_*` prefix. Document the field inline with its unit and meaning.
- **`Geometry` and `XCOMTable` are derived once at startup** from `Params` and the data file. They are read-only during the MC. Never store them inside `Params`.
- **`track_one_photon!` takes `Geometry`, `XCOMTable`, and `Params` as separate arguments**, not bundled. This makes the dependency explicit at the function signature.
- **All units in source: cm, g, MeV, Bq, seconds.** Conversions to keV happen at I/O boundaries (HDF5, plots, summary text).
- **Per-thread RNGs and per-thread histograms; merge after the parallel block.** Never share mutable state across threads.
- **Tests must pass before declaring a change done.** `julia --project=. -t 4 test/runtests.jl`.
- **Reproducibility:** the master seed is `params.mc_seed`; per-thread seeds are derived deterministically as `mc_seed + thread_id`.

## Quick numerical reference

For reasoning about results without re-reading the full spec:

- LXe cylinder: R = 149 cm, L = 297 cm → 61.6 t LXe
- 119 copper rings × ~10 kg = ~1184 kg Cu total
- At 2 μBq/kg: total Cu activity ≈ 2.4 mBq Bi-214
- 2.448 MeV γ in Xe: μ_lin ≈ 0.099 cm⁻¹ (Compton + photoelectric only), mfp ≈ 10 cm
- Photoelectric branching at 2.448 MeV: ~2.6%; Compton ~97.4%
- Bi-214 → 2.448 MeV branching: 1.55%
- ROI: Q_ββ ± 1 FWHM = 2458 ± 37.6 keV (XLZD σ_E/E = 0.65%)
- SS/MS threshold: |Δz| < 3 mm
- FV box (display-only, NOT a cut): z ∈ [50, 250] cm, r² ≤ 100² cm²
- Self-consistency target: ~0.4 events from FCR in the FV at default activity

The MC counts SS-in-ROI events across the **full LXe cylinder**; the FV-restricted count is a post-processing integral on the heatmap.