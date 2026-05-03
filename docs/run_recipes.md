# XLZD MC — sources, run commands, and outputs

Quick reference for `scripts/run_mc3.jl`, `py/show_summary.py` and
`py/plot_histograms.py`.

## Conventions

**Energy resolution** is expressed two ways across experiments. We
follow the LZ + nEXO convention everywhere:

| Convention                   | Symbol     | Quoted by                  |
|------------------------------|------------|----------------------------|
| Standard deviation (rms)     | `σ_rms`    | LZ, nEXO, this MC          |
| Full width at half maximum   | `FWHM`     | NEXT, KamLAND-Zen, …       |

For a Gaussian peak: `FWHM = 2√(2 ln 2) · σ_rms ≈ 2.3548 · σ_rms`.

The CLI flag `--sigma-e` is the resolution as `σ_rms / E` (LZ default
`0.007` = 0.7 %).

**ROI default** is `Q_ββ ± 1 FWHM` at the supplied `σ_rms`, computed
automatically. At default `σ_rms / E = 0.7 %`:

```
σ_rms at Q_ββ = 0.007 × 2458 = 17.21 keV
1 FWHM        = 2.355 × 17.21 = 40.52 keV
ROI window    = [Q_ββ − 40.52, Q_ββ + 40.52] keV = [2417.5, 2498.5] keV
```

To override the auto-default, pass `--roi-halfwidth N` (e.g. `17.2` for
±1σ, the previous default). The actual cut is always recorded in
`summary.csv` (`ROI_halfwidth_keV` column).

## Sources

The MC propagates γ from **18 main effective sources** (6 cryostat +
12 field cage), plus a Tl-208 cascade companion attached to each Tl-208
main source for the companion-veto step. Activities and masses are from
the LZ bb0nu paper Table I; geometry is summarized in
`design/lz_fieldcage_model.md` and `docs/LZ_detector_summary.md`.

### Cryostat (1–6)

| index | name        | isotope | E (MeV) | region        | role                    |
|-------|-------------|---------|---------|---------------|-------------------------|
| 1     | CB_Bi214    | Bi-214  | 2.448   | barrel        | main γ                  |
| 2     | CTH_Bi214   | Bi-214  | 2.448   | top endcap    | main γ                  |
| 3     | CBH_Bi214   | Bi-214  | 2.448   | bottom endcap | main γ                  |
| 4     | CB_Tl208    | Tl-208  | 2.615   | barrel        | main γ + companion veto |
| 5     | CTH_Tl208   | Tl-208  | 2.615   | top endcap    | main γ + companion veto |
| 6     | CBH_Tl208   | Tl-208  | 2.615   | bottom endcap | main γ + companion veto |

### Field cage (7–18)

Each field-cage component is one source per isotope; FC components sit
*inside* the LXe with no Ti to traverse. The "barrel" sources are
smeared cylindrical shells over the full FC axial extent
(z ∈ [−13.75, 145.6] cm); the two grid holders are annular slabs
(R_in = 72.8, R_out = 80.3 cm) with a single LXe-facing emission face.

| index | name        | component               | R (cm)       | z geometry                    | mass (kg) |
|-------|-------------|-------------------------|--------------|-------------------------------|-----------|
| 7     | FCRN_Bi214  | Field-cage rings (Ti)   | 74.3         | shell, full barrel            | 93.0      |
| 8     | FCRS_Bi214  | Field-cage resistors    | 74.3         | shell, full barrel            | 0.06      |
| 9     | FCSE_Bi214  | TPC sensors             | 74.3         | shell, full barrel            | 5.02      |
| 10    | FCPT_Bi214  | PTFE walls              | 72.8         | shell, full barrel            | 184.0     |
| 11    | FCTG_Bi214  | Top grids+holders (anode+gate) | 72.8–80.3 | annular slab at z=145.6, opens up   | 44.55 |
| 12    | FCBG_Bi214  | Cathode grid+holder     | 72.8–80.3    | annular slab at z=0, opens down     | 22.28 |
| 13–18 | FCRN_Tl208 … FCBG_Tl208 | (same components, Tl-208 isotope) | … | … | … |

The bottom-shield grid (deepest in LXe) is intentionally dropped — the
γ paths to the FV are blocked by the cathode + RFR, so its
contribution is negligible in this model. The cathode mass is taken as
¼ of the bb0nu paper's lumped 89.1 kg "field grids and holders" (4
grids of comparable mass).

Specific activities, masses, and the per-component μ·t at 2.448 MeV
are tabulated in the test output of `test/test_field_cage3.jl`.

### Cascade companions (auto-attached)

For every Tl-208 source the matching `*_Tl208c` companion source (lumped
at 583 keV, 99% of Tl-208 decays) is built but not user-selectable. The
driver attaches it automatically when a Tl-208 main source is run.

The driver runs **exactly one** main source per process. Production
jobs are per-source by design (e.g. 1e9 events × 18 sources is not a
single-process workload). Specify the source via `--source-name` or
`--source-index`; both can be passed if they agree.

## CLI flags (`scripts/run_mc3.jl`)

| flag                    | default          | meaning                                                         |
|-------------------------|------------------|-----------------------------------------------------------------|
| `--n-samples N`         | `1_000_000`      | photons sampled                                                 |
| `--seed N`              | `1234`           | master RNG seed (per-thread = seed + tid)                       |
| `--run-name STR`        | `last_run_v2`    | run dir under `output/<STR>/`                                   |
| `--source-name NAME`    | (empty)          | canonical source name (e.g. `CB_Tl208`)                         |
| `--source-index N`      | `0`              | 1..18 alternative (1–6 cryostat, 7–12 FC Bi-214, 13–18 FC Tl-208) |
| `--sigma-e σ`           | `0.007`          | `σ_rms / E` (LZ convention)                                     |
| `--roi-halfwidth keV`   | auto             | default = 1 FWHM at the supplied σ; pass to override            |
| `--fv-z-min cm`         | `26.0`           | FV box, lower z                                                 |
| `--fv-z-max cm`         | `96.0`           | FV box, upper z                                                 |
| `--fv-r-max cm`         | `39.0`           | FV box, max radial (`fv_r2_max_cm2 = r_max²`)                   |
| `--histos true|false`   | `true`           | write histogram CSVs (12 per source) and summary.csv            |
| `--sample-stack N`      | `0`              | stack-row CSV: 0=off, −1=all events, n>0 ≈ n events             |

Pass `-t N` to `julia` for `Threads.nthreads()`.

### Source-selection validation

| `--source-index` | `--source-name` | action                                                |
|------------------|-----------------|-------------------------------------------------------|
| given            | given           | accept iff `names[index] == name`; abort on mismatch  |
| given            | (empty)         | accept `names[index]`; abort if out of [1, 6]         |
| 0                | given           | accept `name`; abort if not in known list             |
| 0                | (empty)         | abort: must specify one                               |

## Common recipes

### Single-source production run (1e9 events, default cuts)
```bash
julia --project=. -t 8 scripts/run_mc3.jl \
    --n-samples 1000000000 --source-name CB_Tl208 \
    --run-name FV_default
```
Produces `output/FV_default/CB_Tl208/`. Repeat for the other 17 sources
(each in its own job).

### Field-cage source — annular grid holder (FCBG)
```bash
julia --project=. -t 8 scripts/run_mc3.jl \
    --n-samples 1000000000 --source-name FCBG_Bi214 \
    --run-name FV_default
```
FCBG sits in the LXe at z=0 (cathode plane) — emits over an annular
slab spanning R ∈ [72.8, 80.3]. Self-shielding is non-negligible
(μ·t ≈ 0.24 axial), already folded into the source's dN/du.

### Loose FV + ±1σ ROI (LZ-paper-like)
```bash
julia --project=. -t 8 scripts/run_mc3.jl \
    --n-samples 1000000000 --source-name CB_Tl208 \
    --run-name LZ_paper \
    --roi-halfwidth 17.2 \
    --fv-z-min 10 --fv-z-max 110 --fv-r-max 60
```

### Quick smoke (single source, no histos, very fast)
```bash
julia --project=. -t 4 scripts/run_mc3.jl \
    --n-samples 100000 --source-index 4 \
    --run-name smoke --histos false
```
Only writes `output/smoke/CB_Tl208/summary.csv`.

### Stack-content sampling (debug photon cascades)
```bash
julia --project=. -t 4 scripts/run_mc3.jl \
    --n-samples 100000 --source-name CB_Tl208 \
    --sample-stack 50 --run-name cascade_dump
```
Writes 50 sample stacks to
`output/cascade_dump/CB_Tl208/stack_sample.csv`. `-1` dumps every
event with a non-empty stack (large).

## Output layout

```
output/<run-name>/                          # set by --run-name
├── run_report.txt                          # written by py/show_summary.py
├── run_report.csv                          # written by py/show_summary.py
└── <source>/                               # one dir per source-job
    ├── summary.csv                         # single row, 24 columns
    ├── cut1_h_u_sampled.csv                ┐ Cut histos (8 files;
    ├── cut2_first_interaction_r_z.csv      │  drives cuts_acceptance.png
    ├── cut3_dz_inclusive.csv               │  + cuts_classification.png).
    ├── cut3_n_visible.csv                  │
    ├── cut3_E_total.csv                    │
    ├── cut4_ss_ec_pre_roi.csv              │
    ├── cut4_ss_es_pre_roi.csv              │
    ├── cut4_ss_r_z.csv                     ┘
    ├── diag_interaction_type_freq.csv      ┐ Diagnostic histos (4 files;
    ├── diag_path_length_LXe.csv            │  drives diagnostics.png).
    ├── diag_region_interaction.csv         │
    ├── diag_cluster_Ec.csv                 ┘
    ├── stack_sample.csv                    # only if --sample-stack ≠ 0
    ├── cuts_acceptance.png                 # written by py/plot_histograms.py
    ├── cuts_classification.png
    └── diagnostics.png
```

### `summary.csv` columns

`source, isotope, n_total, gamma_per_yr_total, f_SS_in_ROI, bg_per_yr,
r_comp, runtime_s, n_escaped, n_MS, n_skin_vetoed, n_outside_FV,
n_SS_outside_ROI, n_SS_in_ROI, n_companion_vetoed, Q_betabeta_keV,
sigma_E_over_E, ROI_halfwidth_keV, fv_z_min_cm, fv_z_max_cm,
fv_r_max_cm, E_visible_keV, E_skin_veto_keV, R_LXe_outer_cm`.

The 7 `n_*` outcome counts sum to `n_total`. The cut/geometry columns
record the *actual* settings used (so a later compare can tell which
σ_rms, ROI, FV produced any row).

### Cut-flow histogram CSVs (8)

| file                              | shape       | content                                                          |
|-----------------------------------|-------------|------------------------------------------------------------------|
| `cut1_h_u_sampled.csv`            | 100         | `u = cos θ_inward` of every sampled photon (pre-propagation).    |
| `cut2_first_interaction_r_z.csv`  | 100×100     | (r, z) of every event's first interaction in an LXe region.      |
| `cut3_dz_inclusive.csv`           | 200         | every consecutive-visible-deposit `|Δz|`, [0, 5] cm (zoomed).    |
| `cut3_n_visible.csv`              | 21 (int)    | n_visible per cut-2-passing event.                               |
| `cut3_E_total.csv`                | 270         | Σ visible-cluster energy per cut-2-passing event, [0, 2.7] MeV.  |
| `cut4_ss_ec_pre_roi.csv`          | 270         | true SS-cluster energy `ec` before ROI cut.                      |
| `cut4_ss_es_pre_roi.csv`          | 270         | smeared SS-cluster energy `es` before ROI cut (cut acts here).   |
| `cut4_ss_r_z.csv`                 | 100×100     | (r, z) of the SS cluster centroid before ROI cut.                |

### Diagnostic histogram CSVs (4)

| file                                | shape   | content                                                             |
|-------------------------------------|---------|---------------------------------------------------------------------|
| `diag_interaction_type_freq.csv`    | 4       | per-stack-row PHOTO/COMPTON/PAIR/BELOW_THRESH (∝ cross sections).   |
| `diag_path_length_LXe.csv`          | 100     | per-photon Σ distance through LXe (TPC+Skin+Inert), [0, 500] cm.    |
| `diag_region_interaction.csv`       | 3×4     | (TPC/Skin/Inert) × (PHOTO/COMPTON/PAIR/BELOW_THRESH) counts.        |
| `diag_cluster_Ec.csv`               | 270     | per-cluster energy `ec` across all clusters, all events.            |

## Aggregating a run — `py/show_summary.py`

After running one or more sources under the same `--run-name`,
aggregate them into a single text + CSV report:

```bash
conda activate itacatf
python py/show_summary.py output/<run-name>/
```

Output:
- stdout — three blocks (per-source results table, analysis cuts,
  per-source funnel).
- `output/<run-name>/run_report.txt` — same text, easy cut-and-paste.
- `output/<run-name>/run_report.csv` — aggregated multi-row CSV
  (one row per source, same columns as per-source `summary.csv`).

Optional flags:
- `--source NAME` — restrict to one source.
- `--no-files` — stdout only (skip the `.txt` / `.csv` writes).

## Plotting histograms — `py/plot_histograms.py`

Three combined PNG panels per source:

```bash
conda activate itacatf
python py/plot_histograms.py output/<run-name>/<source>/
```

Outputs:
- `cuts_acceptance.png` — Cut 1 (`h_u_sampled` with geometric-acceptance
  marker for barrel sources) + Cut 2 (`first_interaction_r_z` with FV
  box overlay).
- `cuts_classification.png` — 2×3 grid: Cut 3 row (`dz_inclusive` with
  3 mm marker, `n_visible`, `E_total`) and Cut 4 row (`ss_ec`, `ss_es`
  with ROI band, `ss_r_z` with FV box).
- `diagnostics.png` — 2×2: interaction-type frequency, path length,
  region × interaction matrix, per-cluster `Ec`.

Cut overlays (FV box, ROI window, 3 mm marker, `u_min` geometric
acceptance) are read from the per-source `summary.csv` in the same
directory.

Options:
- `--family acceptance|classification|diagnostics|all` — render a
  single panel.
- `--separate` — also save one PNG per CSV.
- `--show` — open in a window.
- `--output DIR` — write PNGs to DIR (default: input dir).

### Loop over all sources under a run

```bash
conda activate itacatf
for src in output/run_1e9/*/ ; do
    [ -f "$src/summary.csv" ] && python py/plot_histograms.py "$src"
done
```

## Field cage cross-check

Reference numbers from the LZ bb0nu paper (FV [26, 96] × r ≤ 39 cm),
reproduced from the spec in `design/lz_fieldcage_model.md` §3:

| component       | spec γ/yr | bb0nu γ/yr |
|-----------------|-----------|------------|
| Field-cage rings (FCRN)            | 0.50 | 0.30 |
| Field-cage resistors+sensors (FCRS+FCSE) | 1.70 | 1.39 |
| PTFE walls (FCPT)                  | 0.14 | 0.14 |
| Top grid holders (FCTG)            | 0.03 | —    |
| Bottom holder (full mass, both grids) | 0.38 | —    |
| **Holders total (paper)**          | **0.41** | **0.23** |
| **Field-cage total (paper)**       | **2.75** | **2.06** |

Our model (cathode-only FCBG = ¼ of the lumped 89.1 kg, vs the spec's
"full bottom holder" using ½) trims the FCBG mass from 44.55 to 22.28
kg. So our predicted FCBG should land near `0.38 / 2 ≈ 0.19` γ/yr.

A 1e6-sample smoke pass at default cuts (σ_rms/E = 0.7%, ROI = ±1 FWHM,
default FV) reproduced these orders of magnitude:

| component | smoke γ/yr | spec γ/yr | ratio |
|-----------|-----------|-----------|-------|
| FCRN_Bi214 | 0.90 | 0.50 | 1.8× |
| FCRS_Bi214 | 2.12 | (in 1.70 lump) | — |
| FCSE_Bi214 | 0.79 | (in 1.70 lump) | — |
| FCRS+FCSE | 2.91 | 1.70 | 1.7× |
| FCPT_Bi214 | 0.090 | 0.14 | 0.6× |
| FCTG_Bi214 | 0.032 | 0.03 | 1.1× |
| FCBG_Bi214 | 0.181 | 0.19 | 0.9× |
| **Bi-214 total** | **4.11** | **2.75 (full mass)** | **1.5×** |

A 1e9-sample production pass tightens the per-source uncertainty to ≪1%
and resolves the Tl-208 contributions (zero at 1e6 because the 2.615
MeV line sits above ROI; only Compton-edge fluctuations land in the
ROI window, requiring much higher statistics to count).
