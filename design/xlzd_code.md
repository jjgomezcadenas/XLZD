# code.md — XLZD Field-Cage Ring Background Monte Carlo (Build Specification)

> **Status**: original implementation specification.
> **Audience**: the agent implementing this project from scratch, and anyone making non-trivial changes later.
> **Read order**: this is the authoritative source for design decisions, parameter choices, and the algorithmic spec. Read it in full before writing or refactoring code.
> **Project root**: the surrounding repository (`Project.toml`, `src/`, `scripts/`, etc.).



## 1. Purpose

Estimate the contribution of Bi-214 (²³⁸U-late chain, 2.448 MeV γ-line) decays in the copper field-cage rings of an XLZD-like detector to the 0νββ background budget at Q_ββ = 2458 keV in ¹³⁶Xe. The output is the number of single-site (SS) events landing in a fiducial volume cut and a region-of-interest (ROI) energy window, plus 2D spatial heatmaps and a 3D interactive visualization.

The MC implements full Compton tracking with Klein-Nishina kinematics, an XLZD-style SS/MS classification using Δz clustering, and the experiment's ROI energy window. The single-number output (events per emitted photon, scaled by Bi-214 → 2.448 MeV branching ratio and total ring activity) is the calibration that feeds the comparison with plated ²²²Rn-daughter contributions to the FCR background.

The code must be **Julia, multithreaded with Threads.@threads**, and run ~10⁸ samples in under a minute on 8 cores.

## 2. Physical Setup

A cylindrical liquid-xenon volume of radius R_lxe and length L_lxe is surrounded by a stack of copper field-cage rings just outside the cylinder. The rings emit isotropic 2.448 MeV gammas distributed uniformly throughout their copper bulk by mass. Each photon is tracked through the LXe with full Compton kinematics until it either:

- escapes the cylinder, or
- is photoabsorbed (full remaining energy deposited locally), or
- is dropped early because its current energy fell below the tracking cutoff (deposit remaining energy at current position and stop).

At each interaction site, the deposited energy is checked against the SS clustering criterion (z-only) and either added to the cluster or used to reject the event as multi-site (MS).

Pair production is **included in the total cross section** for correct mean free path and interaction branching ratios. At 2.448 MeV, pair production is ~12% of the total cross section in Xe. When pair production is sampled as the interaction type, the event is **immediately classified as MS-rejected** without tracking the annihilation gammas: the two back-to-back 511 keV γ travel mfp ~6 cm in LXe and virtually always produce spatially separated deposits that fail the SS z-cut. This gives the correct branching (photo ~2.4%, Compton ~85%, pair ~12%) and correct μ_total (mfp ~8.7 cm rather than ~10.1 cm without pair), while avoiding the complexity of tracking secondary gammas.

## 3. Geometry — XLZD Defaults (NEXT-style ring scaling)

| Parameter | Symbol | Default | Notes |
|---|---|---:|---|
| LXe cylinder inner radius | `R_lxe` | 149 cm | Half of 298 cm diameter |
| LXe cylinder length | `L_lxe` | 297 cm | Field-cage length |
| Ring axial pitch | `pitch` | 2.5 cm | 25 mm |
| Ring cross-section area | `A_cs` | 1.2 cm² | NEXT geometry: rounded rectangle 10 mm × 8 mm with 1 mm corners |
| Ring radial thickness | `t_radial` | 1.0 cm | 10 mm |
| Ring axial width | `w_axial` | 1.2 cm² / t_radial = 1.2 cm | Set by A_cs / t_radial to preserve area |
| Ring mean radius | `R_ring_mean` | R_lxe + t_radial/2 = 149.5 cm | |
| Cu density | `ρ_Cu` | 8.96 g/cm³ | |
| LXe density | `ρ_LXe` | 2.953 g/cm³ | Saturated liquid ~165 K |
| Cu specific activity (²³⁸U-late) | `act_Cu` | 2 μBq/kg | XLZD-projected target; user parameter |

### Ring axial positions

```
N_rings = floor(Int, L_lxe / pitch)             # = 118 for defaults? No, use ceiling
N_rings = floor(Int, L_lxe / pitch) + 1          # default: 119
z_ring[i] = (i - 1) * pitch + pitch/2  for i = 1, …, N_rings
```

This places the first ring at z = pitch/2 = 1.25 cm and the last at z = (N_rings − 1) × pitch + pitch/2 ≈ L_lxe − pitch/2. Centered between ends with a half-pitch margin.

Confirm at runtime: `z_ring[end] ≤ L_lxe`. If not, decrement N_rings.

### Sanity checks (print at startup)

For defaults (R_lxe=149, L_lxe=297, pitch=2.5, A_cs=1.2):

- N_rings = 119
- Ring volume = 2π × 149.5 × 1.2 = 1127 cm³
- Mass per ring = 1127 × 8.96 / 1000 = 10.10 kg
- Total Cu mass = 119 × 10.10 = 1202 kg
- Total Cu activity at 2 μBq/kg = 2.40 mBq
- LXe mass = π × 149² × 297 × 2.953 / 1000 = 61.6 t

The exact numbers depend on the precise definition of mean radius. Allow ~1% tolerance vs the table above.

## 4. Physics

### 4.1 Cross sections

Pre-tabulate σ_photo(E) and σ_Compton(E) for xenon (Z=54) from NIST XCOM data, supplied by the user as `data/nist.csv`.

**XCOM column mapping** (the data file uses the standard XCOM output format, all units cm²/g):

| Column index | Column name | Used as |
|---:|---|---|
| 1 | Photon Energy (MeV) | energy grid |
| 2 | Coherent Scatter (Rayleigh) | **ignored** (elastic, deposits no energy) |
| 3 | **Incoherent Scatter** | **`σ_Compton(E)`** ← this is the Compton/Klein-Nishina cross section |
| 4 | **Photoelectric Absorption** | **`σ_photo(E)`** |
| 5 | Nuclear Pair Production | **ignored** (per user instruction) |
| 6 | Electron Pair Production | **ignored** (per user instruction) |
| 7 | Total w/ Coherent | informational |
| 8 | Total w/o Coherent | informational |

The three cross sections used in the MC are:
- `σ_photo(E)` = column 4 (photoelectric absorption)
- `σ_Compton(E)` = column 3 (incoherent scatter; this is the Compton cross section, NOT the coherent one)
- `σ_pair(E)` = column 5 + column 6 (nuclear + electronic pair production; zero below 1.022 MeV)
- `σ_total(E) = σ_photo(E) + σ_Compton(E) + σ_pair(E)` for sampling step length and interaction branching

Convert from cm²/g to linear attenuation cm⁻¹ by multiplying by ρ_LXe = 2.953 g/cm³.

**Energy range and the K-edge:**

The XCOM table spans 30 keV to 3.0 MeV. It includes a duplicate-energy row at 34.56 keV showing the discontinuity at the Xe K-edge:

```
3.456E-02 6.119E-01 1.008E-01 5.417E+00 ... ← just below K-edge
3.456E-02 6.119E-01 1.008E-01 3.244E+01 ... ← just above K-edge
```

Photoelectric cross section jumps by ~6× across the K-edge.

**The MC never queries cross sections below the photon tracking cutoff (default 40 keV)**, which is above the K-edge. The interpolator should:
- Use log-log linear interpolation between adjacent grid points.
- For the K-edge: use the **above-edge** value (the second 34.56 keV row) for E ≥ 34.56 keV. The below-edge row exists only for documentation; the MC should drop it during table loading.
- Error out if asked for E < 30 keV (should never happen given the tracking cutoff).

**Interpolation to 2.448 MeV** (the Bi-214 line, the photon's initial energy):

Bracketed by 2.044 and 3.000 MeV rows. Log-log linear interpolation gives:

- σ_photo(2.448) ≈ 8.75 × 10⁻⁴ cm²/g
- σ_Compton(2.448) ≈ 3.27 × 10⁻² cm²/g
- σ_pair(2.448) ≈ 4.5 × 10⁻³ cm²/g
- σ_total(2.448) ≈ 3.81 × 10⁻² cm²/g
- μ_lin(2.448) = σ_total × ρ_LXe ≈ 0.113 cm⁻¹ → mfp ≈ 8.9 cm

Interaction-type branching at 2.448 MeV:
- Photoelectric: ~2.3%
- Compton: ~85.8%
- Pair production: ~11.8%

When pair production is sampled, the event is immediately MS-rejected (see Section 2). This gives correct physics: shorter mean free path than without pair, and ~12% of interactions produce guaranteed MS events.

### 4.2 Klein-Nishina sampling

Sample the scattered photon energy E' from the Compton differential cross section. Use a standard rejection algorithm (e.g., Kahn's method for κ = E/m_e c² ≥ ~0.1).

Once E' is sampled, the polar scattering angle θ is determined by:

cos θ = 1 − m_e c² × (1/E' − 1/E)

Azimuthal angle φ is uniform in [0, 2π). Construct the new direction by rotating the previous direction by (θ, φ).

Use a numerically stable rotation (build orthogonal basis perpendicular to current direction, then rotate within that basis). Avoid Euler-angle representations.

### 4.3 Bi-214 → 2.448 MeV branching

Bi-214 has multiple γ-lines. The 2.448 MeV line (one of the highest-energy in the chain) has a branching ratio of **1.55%** per Bi-214 decay (ENSDF). The MC simulates **only this line**. The final per-Bq event rate is multiplied by 0.0155 to convert from "events per emitted γ" to "events per Bi-214 decay."

The 2.448 MeV line was chosen because it sits closest to Q_ββ = 2458 keV. Other Bi-214 γ at lower energy contribute to the ROI only via Compton up-scattering (negligible) or via random coincidence (not modeled here).

## 5. Event Selection (XLZD-style SS/MS + ROI; FV applied post-MC)

### 5.1 SS/MS classification (z-only clustering)

XLZD has good z-resolution (drift time) but coarse (x, y) resolution. SS/MS classification is therefore based on **|Δz|**, not 3D distance.

For each emitted photon, track the cluster z-position (`z_cluster`) and accumulated energy (`E_cluster`). The cluster is initialized on the **first visible deposit** (E_dep > E_visible_threshold). The cluster's (x, y) is also recorded at the first visible deposit, but **only for histogramming and post-processing**, not for SS/MS classification.

For each subsequent interaction at z_new with deposit E_dep:
- If E_dep < E_visible_threshold: deposit is invisible. Photon continues with reduced energy E'. **Do not** update cluster.
- If E_dep ≥ E_visible_threshold:
  - If |z_new − z_cluster| < Δz_threshold: add E_dep to E_cluster (same site).
  - Else: event is **MS-rejected**. Stop tracking the photon.

### 5.2 ROI energy window

ROI = [Q_ββ − ROI_halfwidth, Q_ββ + ROI_halfwidth]

Default ROI_halfwidth = 1 FWHM = 2.355 × σ_E = 2.355 × 0.0065 × Q_ββ ≈ 37.6 keV.

The window is symmetric in energy around Q_ββ = 2458 keV. The total `E_cluster` is checked against this window only **at the end** of photon tracking — i.e., we require the event to be SS first, then check ROI.

### 5.3 Fiducial volume — display-only and post-processed

**The FV is NOT applied as an event-selection cut during the MC.** The MC counts all SS-in-ROI events across the **full LXe cylinder**, and the (z, r²) heatmap covers the full cylinder. This preserves flexibility — the user can read off events for any chosen FV from the heatmap.

The fiducial volume parameters (`fv_z_min_cm`, `fv_z_max_cm`, `fv_r2_max_cm2`) serve two purposes:

1. **Heatmap overlay**: drawn as a rectangle on the (z, r²) heatmap so the reader sees XLZD's fiducial region.
2. **Post-processing integral**: the writer integrates the signal heatmap inside the FV box and reports `events_in_default_FV` in the HDF5 derived/ group and in `summary.txt`. This is the number directly comparable to XLZD's quoted ~0.4 events.

Cluster position used for histogramming is the position of the first visible deposit (where `z_cluster`, `x_cluster`, `y_cluster` were set).

## 6. Monte Carlo Algorithm

### 6.1 Per-photon procedure

```
For each emitted photon:
  
  # 1. Sample emission point
  Sample i ~ Uniform({1, …, N_rings}).  z_emit_center = z_ring[i].
  Sample φ ~ Uniform[0, 2π).
  Sample Δr ~ Uniform[-t_radial/2, +t_radial/2].
  Sample Δz_ring ~ Uniform[-w_axial/2, +w_axial/2].
  
  x_emit = (R_ring_mean + Δr) * cos(φ)
  y_emit = (R_ring_mean + Δr) * sin(φ)
  z_emit = z_emit_center + Δz_ring
  
  # 2. Sample isotropic direction
  u = 2 * rand() - 1            # cos θ_iso ∈ [-1, 1]
  v = 2π * rand()                # azimuthal
  s = sqrt(1 - u²)
  dx = s * cos(v); dy = s * sin(v); dz = u
  
  # 3. Backward-rejection (free 2× speedup)
  dot_radial = x_emit * dx + y_emit * dy
  if dot_radial > 0:
      → outcome = "outward". Continue to next photon.
  
  # 4. Initialize photon state
  x, y, z = x_emit, y_emit, z_emit
  E = 2.448  # MeV
  E_cluster = 0; z_cluster = NaN; x_cluster = NaN; y_cluster = NaN
  cluster_started = false
  outcome = "in_progress"
  n_interactions = 0
  trajectory = (record_this_photon ? [(x, y, z, E, "emit", 0.0)] : nothing)
  
  # 5. Photon-tracking loop
  while outcome == "in_progress":
      
      μ_total = σ_photo(E) + σ_Compton(E)  (linear, in cm⁻¹)
      
      # Sample step length
      d_int = -log(rand()) / μ_total
      
      # Compute distance to cylinder boundary along (dx, dy, dz)
      d_boundary = path_to_cylinder_exit(x, y, z, dx, dy, dz, R_lxe, L_lxe)
      
      if d_boundary < 0:    # already outside cylinder (shouldn't happen, sanity check)
          outcome = "error"
          break
      
      if d_int >= d_boundary:
          # Photon escapes before next interaction
          outcome = (cluster_started ? finalize_outcome(E_cluster) : "escaped")
          break
      
      # Move to interaction point
      x += d_int * dx; y += d_int * dy; z += d_int * dz
      n_interactions += 1
      
      # Sample interaction type
      if rand() < σ_photo(E) / μ_total:
          interaction_type = :photoelectric
      else:
          interaction_type = :compton
      
      if interaction_type == :photoelectric:
          E_dep = E
          if record_this_photon: push_trajectory((x,y,z,E,:photo,E_dep))
          
          if E_dep >= E_visible_threshold:
              if not cluster_started:
                  z_cluster = z; x_cluster = x; y_cluster = y
                  E_cluster = E_dep
                  cluster_started = true
              else:
                  if abs(z - z_cluster) < Δz_threshold:
                      E_cluster += E_dep
                  else:
                      outcome = "MS_rejected"
                      break
          # Photon fully absorbed
          if outcome == "in_progress":
              outcome = finalize_outcome(E_cluster)
          break
      
      else:  # Compton
          E_scatt, theta = sample_klein_nishina(E)   # E_scatt < E
          E_dep = E - E_scatt
          if record_this_photon: push_trajectory((x,y,z,E_scatt,:compton,E_dep))
          
          if E_dep >= E_visible_threshold:
              if not cluster_started:
                  z_cluster = z; x_cluster = x; y_cluster = y
                  E_cluster = E_dep
                  cluster_started = true
              else:
                  if abs(z - z_cluster) < Δz_threshold:
                      E_cluster += E_dep
                  else:
                      outcome = "MS_rejected"
                      break
          # else: invisible deposit, photon continues
          
          # Update direction (rotate by theta around current direction)
          φ_az = 2π * rand()
          (dx, dy, dz) = rotate_direction(dx, dy, dz, theta, φ_az)
          E = E_scatt
          
          # Energy cutoff
          if E < E_tracking_cutoff:
              # Treat remaining E as deposited at current location
              E_dep_final = E
              if E_dep_final >= E_visible_threshold:
                  if not cluster_started:
                      z_cluster = z; x_cluster = x; y_cluster = y
                      E_cluster = E_dep_final
                      cluster_started = true
                  else:
                      if abs(z - z_cluster) < Δz_threshold:
                          E_cluster += E_dep_final
                      else:
                          outcome = "MS_rejected"
                          break
              if record_this_photon: push_trajectory((x,y,z,0.0,:cutoff,E_dep_final))
              if outcome == "in_progress":
                  outcome = finalize_outcome(E_cluster)
              break
          
          # Otherwise loop continues with new E and direction
  
  # End of tracking. Record outcome.
  increment counter[outcome]
  if outcome == "SS_in_ROI": fill_2d_histogram(z_cluster, x_cluster² + y_cluster²)
  fill_diagnostic_histograms(...)
  if record_this_photon: save_trajectory_to_json(trajectory, outcome, ...)
```

### 6.2 `finalize_outcome` function

```
function finalize_outcome(E_cluster):
    if E_cluster == 0:
        return "escaped"           # nothing visible deposited
    
    # ROI check (NO FV cut here — FV is post-processed)
    in_ROI = abs(E_cluster - Q_ββ) ≤ ROI_halfwidth
    
    if in_ROI:
        return "SS_in_ROI"          # the headline signal — counted across full LXe volume
    else:
        return "SS_outside_ROI"
```

The cluster's (x, y, z) is recorded for histogramming regardless of outcome. The post-processing step in `write_outputs` integrates `H_signal` (which contains all `SS_in_ROI` events) inside the default FV box and reports the FV-restricted count separately.

### 6.3 Cylinder intersection

For a ray (x, y, z) + t × (dx, dy, dz), the cylinder R_lxe radius / [0, L_lxe] axial:

Solve quadratic (x + t·dx)² + (y + t·dy)² = R_lxe² for t. Get roots t₁ ≤ t₂.

If discriminant < 0: ray misses cylinder → return d_boundary = 0 (will be caught as "escape immediately").

Otherwise:
- t_radial_exit = max(t₁, t₂) clipped to ≥ 0
- t_axial_exit_zmin = (0 - z) / dz   if dz < 0; else +∞
- t_axial_exit_zmax = (L_lxe - z) / dz   if dz > 0; else +∞
- d_boundary = min(t_radial_exit, t_axial_exit_zmin, t_axial_exit_zmax)

This is the path length until the photon leaves the cylinder via either the side wall, top, or bottom.

If the photon starts on or just inside the boundary, d_boundary should be a small positive number — handle with care (use a tolerance ~10⁻⁶ cm).

### 6.4 Direction rotation

Given current direction d = (dx, dy, dz) and a scattering polar angle θ relative to d, plus azimuth φ, compute the new direction d' = R(d, θ, φ) · d.

Standard implementation: build an orthogonal basis (u, v, d) where u, v are perpendicular to d (any choice), then:

```
d'_x = sin(θ)cos(φ) · u_x + sin(θ)sin(φ) · v_x + cos(θ) · d_x
d'_y = sin(θ)cos(φ) · u_y + sin(θ)sin(φ) · v_y + cos(θ) · d_y
d'_z = sin(θ)cos(φ) · u_z + sin(θ)sin(φ) · v_z + cos(θ) · d_z
```

Building (u, v): if |d_z| < 0.9, set u = normalize(cross(ẑ, d)); else use x̂. Then v = cross(d, u).

Verify d'·d ≈ cos(θ) and |d'| ≈ 1 in unit tests.

## 7. Parameters — the `Params` struct

All user-tunable inputs go into a single flat struct, `Params`. Field names use prefixes (`geom_`, `act_`, `phys_`, `cut_`, `mc_`, `hist_`, `out_`) to indicate which concern each parameter belongs to without nested deref. Use `Base.@kwdef` so partial-override construction works (e.g., `Params(act_Cu=1e-6, mc_N_samples=10^7)`).

The struct is **immutable** (`struct`, not `mutable struct`). It is passed by reference into geometry construction, the MC, and the output writer. It is never modified after construction.

### 7.1 The struct (literal Julia code to use)

```julia
"""
    Params

All user-tunable inputs to the XLZD field-cage-ring background MC.
A single instance is constructed at startup (from CLI args or defaults)
and passed to `build_geometry`, `run_mc`, and `write_outputs`. It is
immutable.

Fields are prefix-grouped:
  geom_*  — detector geometry
  act_*   — source activity
  phys_*  — physics constants and external data paths
  cut_*   — event-selection cuts (SS/MS, ROI, FV)
  mc_*    — Monte Carlo control (N samples, seed, trajectory recording)
  hist_*  — output histogram binning
  out_*   — output paths
"""
Base.@kwdef struct Params
    # ─── Geometry (geom_*) ──────────────────────────────────────────────
    geom_R_lxe::Float64        = 149.0    # LXe cylinder inner radius (cm)
    geom_L_lxe::Float64        = 297.0    # LXe cylinder length (cm)
    geom_pitch::Float64        = 2.5      # ring axial pitch (cm)
    geom_A_cs::Float64         = 1.2      # ring cross-section area (cm²)
    geom_t_radial::Float64     = 1.0      # ring radial thickness (cm)
    geom_ρ_Cu::Float64         = 8.96     # copper density (g/cm³)
    geom_ρ_LXe::Float64        = 2.953    # LXe density at saturation, ~165 K (g/cm³)

    # ─── Source activity (act_*) ────────────────────────────────────────
    act_Cu::Float64            = 2.0e-6   # Cu specific activity, Bq/kg (= 2 μBq/kg)
    act_BR_2448::Float64       = 0.0155   # Bi-214 → 2.448 MeV branching (ENSDF)

    # ─── Physics constants and tables (phys_*) ──────────────────────────
    phys_E_initial_MeV::Float64    = 2.448   # Bi-214 γ line, photon initial energy
    phys_Q_betabeta_keV::Float64   = 2458.0  # ¹³⁶Xe Q value
    phys_σ_E_over_E::Float64       = 0.0065  # XLZD energy resolution (1σ relative)
    phys_xcom_data_path::String    = "data/nist.csv"

    # ─── Event-selection cuts (cut_*) ───────────────────────────────────
    # These genuinely gate the MC outcome (SS/MS, ROI, energy thresholds).
    cut_ROI_halfwidth_keV::Float64       = 37.6      # ±1 FWHM around Q_ββ by default
    cut_Δz_threshold_mm::Float64         = 3.0       # SS/MS threshold in z
    cut_E_visible_threshold_keV::Float64 = 5.0       # below this, deposit invisible
    cut_E_tracking_cutoff_keV::Float64   = 40.0      # below this, photon stops, deposits remaining E

    # ─── Fiducial volume (fv_*) — display-only, NOT applied as cut ──────
    # The MC counts SS-in-ROI events across the full LXe cylinder.
    # These parameters are used only for (a) drawing the FV box on heatmaps
    # and (b) the post-processing integral inside the FV (reported in summary.txt
    # and HDF5 derived/ group). Direct comparison with XLZD's quoted ~0.4 events.
    fv_z_min_cm::Float64             = 50.0
    fv_z_max_cm::Float64             = 250.0
    fv_r2_max_cm2::Float64           = 10000.0   # = 100² cm²

    # ─── Monte Carlo control (mc_*) ─────────────────────────────────────
    mc_N_samples::Int                = 100_000_000   # total emissions across all threads
    mc_seed::Int                     = 1234          # master RNG seed; per-thread seeds derived from this
    mc_n_traj_per_outcome::Int       = 100           # stratified trajectory recording per outcome category

    # ─── Histogram binning (hist_*) ─────────────────────────────────────
    # Heatmap covers the FULL LXe cylinder so the FV box can be drawn as overlay.
    hist_z_min_cm::Float64         = 0.0
    hist_z_max_cm::Float64         = 297.0           # = geom_L_lxe by default
    hist_r2_min_cm2::Float64       = 0.0
    hist_r2_max_cm2::Float64       = 22201.0         # = geom_R_lxe² by default
    hist_n_z_bins::Int             = 100
    hist_n_r2_bins::Int            = 100
    hist_n_E_bins::Int             = 250
    hist_E_max_keV::Float64        = 2500.0

    # ─── Output (out_*) ─────────────────────────────────────────────────
    out_dir::String                = "output/"
    out_viewer_template::String    = "viewer/viewer.html"
end
```

### 7.2 Default values summary (for reference and CLI mapping)

| Field | Default | Description |
|---|---:|---|
| `geom_R_lxe` | 149.0 | LXe cylinder radius (cm) |
| `geom_L_lxe` | 297.0 | LXe cylinder length (cm) |
| `geom_pitch` | 2.5 | Ring axial pitch (cm) |
| `geom_A_cs` | 1.2 | Ring cross-section area (cm²) |
| `geom_t_radial` | 1.0 | Ring radial thickness (cm) |
| `geom_ρ_Cu` | 8.96 | Copper density (g/cm³) |
| `geom_ρ_LXe` | 2.953 | LXe density (g/cm³) |
| `act_Cu` | 2.0e-6 | Cu specific activity, Bq/kg (= 2 μBq/kg) |
| `act_BR_2448` | 0.0155 | Bi-214 → 2.448 MeV branching |
| `phys_E_initial_MeV` | 2.448 | Photon initial energy (MeV) |
| `phys_Q_betabeta_keV` | 2458.0 | ¹³⁶Xe Q value (keV) |
| `phys_σ_E_over_E` | 0.0065 | Energy resolution (1σ) |
| `phys_xcom_data_path` | "data/nist.csv" | Path to NIST XCOM table |
| `cut_ROI_halfwidth_keV` | 37.6 | ±1 FWHM around Q_ββ |
| `cut_Δz_threshold_mm` | 3.0 | SS/MS threshold |
| `cut_E_visible_threshold_keV` | 5.0 | Visible-deposit threshold |
| `cut_E_tracking_cutoff_keV` | 40.0 | Photon-tracking energy cutoff |
| `fv_z_min_cm` | 50.0 | FV box, z lower (display + post-processing only) |
| `fv_z_max_cm` | 250.0 | FV box, z upper |
| `fv_r2_max_cm2` | 10000.0 | FV box, r² upper (= 100² cm²) |
| `mc_N_samples` | 100_000_000 | Total emissions |
| `mc_seed` | 1234 | Master RNG seed |
| `mc_n_traj_per_outcome` | 100 | Trajectories per outcome category |
| `hist_z_min_cm` | 0.0 | Heatmap z range, lower |
| `hist_z_max_cm` | 297.0 | Heatmap z range, upper (= L_lxe) |
| `hist_r2_min_cm2` | 0.0 | Heatmap r² range, lower |
| `hist_r2_max_cm2` | 22201.0 | Heatmap r² range, upper (= R_lxe²) |
| `hist_n_z_bins` | 100 | Heatmap z bins |
| `hist_n_r2_bins` | 100 | Heatmap r² bins |
| `hist_n_E_bins` | 250 | 1D energy spectrum bins |
| `hist_E_max_keV` | 2500.0 | 1D spectrum upper limit |
| `out_dir` | "output/" | Output directory |
| `out_viewer_template` | "viewer/viewer.html" | Viewer template path |

### 7.3 CLI argument mapping

`scripts/run_mc.jl` parses common overrides from the command line. Every CLI flag maps to one `Params` field:

| CLI flag | Maps to field |
|---|---|
| `--n-samples` | `mc_N_samples` |
| `--act-cu` | `act_Cu` |
| `--seed` | `mc_seed` |
| `--output` | `out_dir` |
| `--xcom` | `phys_xcom_data_path` |
| `--n-traj` | `mc_n_traj_per_outcome` |
| `--roi-halfwidth` | `cut_ROI_halfwidth_keV` |
| `--dz-threshold` | `cut_Δz_threshold_mm` |
| `--e-visible` | `cut_E_visible_threshold_keV` |
| `--e-cutoff` | `cut_E_tracking_cutoff_keV` |

Other parameters can be set by editing defaults at the top of `scripts/run_mc.jl` or programmatically in a Julia REPL/notebook.


## 8. Threading and Performance

### 8.1 Threading

Use `Threads.@threads` static partition. Each thread:
- Has its own `Random.MersenneTwister(seed + thread_id)` (or `TaskLocalRNG()` seeded per-thread)
- Has its own per-thread histograms and counters
- After all threads finish: master sums histograms and merges trajectory lists

```julia
using Random
using Base.Threads

n_t = Threads.nthreads()
samples_per_thread = div(N_samples, n_t)
remainder = N_samples - samples_per_thread * n_t

# Per-thread storage
thread_results = [ThreadResult() for _ in 1:n_t]

@threads for tid in 1:n_t
    n_local = samples_per_thread + (tid <= remainder ? 1 : 0)
    rng = MersenneTwister(seed + tid)
    for _ in 1:n_local
        run_one_photon!(thread_results[tid], rng, params, xcom_tables)
    end
end

# Merge
merged = merge_results(thread_results)
```

### 8.2 Performance targets

- 10⁸ samples in ≤ 60 seconds on 8 threads
- Memory ≤ 2 GB peak
- Per-photon work: ~5–10 interactions on average × ~10 floating-point ops per interaction = ~100 ns per photon. 10⁸ × 100 ns / 8 threads ≈ 1.25 s. Reality will be slower due to RNG calls, branch mispredictions, etc. — target is 1 minute.

### 8.3 Trajectory recording strategy

Stratified sampling for the visualization JSON. For each outcome category, save the first `mc_n_traj_per_outcome` photons that produce that outcome. Allocate a per-thread, per-outcome buffer; merge after all threads finish; cap down to global count.

Outcomes to record:
- `SS_in_ROI` (the headline signal)
- `SS_outside_ROI`
- `MS_rejected`
- `escaped`
- `outward` (sample only ~5, not 100, since they're trivial)

Default total: ~400 trajectories saved → ~1 MB JSON.

Each trajectory is a record with:
- `emission`: x, y, z, ring_idx
- `direction`: dx, dy, dz (initial)
- `outcome`: string
- `interactions`: list of (x, y, z, E_after_interaction, type, E_deposited)
- `cluster`: {z, x, y, E_total} (or null if no cluster started)
- `n_interactions`: count

## 9. Output

### 9.1 HDF5 output file (`output/results.h5`)

```
/geometry/
    R_lxe, L_lxe, pitch, A_cs, t_radial
    n_rings, ring_z_positions[]
    ρ_Cu, ρ_LXe
    total_Cu_mass_kg, total_Cu_activity_Bq
    LXe_mass_kg
/mc_params/
    n_emitted, n_threads, seed
    runtime_seconds
/physics/
    E_gamma_initial_MeV
    Q_ββ_keV
    σ_E_over_E, ROI_halfwidth_keV
    Δz_threshold_mm, E_visible_threshold_keV, E_tracking_cutoff_keV
    BR_2448
    xcom_table (energy_MeV, σ_photo_cm2g, σ_Compton_cm2g, σ_total_cm2g)
/fv/
    z_min_cm, z_max_cm, r2_max_cm2     # display-only FV box, NOT a cut
/counts/
    n_outward                  # rejected at emission (radial dot > 0)
    n_escaped                  # entered LXe but escaped without visible deposit
    n_MS_rejected              # multi-site (Δz > threshold for some deposit)
    n_SS_outside_ROI           # single-site, total cluster E outside ROI window
    n_SS_in_ROI                # single-site, in ROI — the headline signal (full LXe)
    (sum should equal n_emitted)
/histograms_2D/
    bin_edges_z[]              # cm — covers full LXe cylinder (0 to L_lxe)
    bin_edges_r2[]             # cm² — covers full LXe cylinder (0 to R_lxe²)
    H_first_interaction[]      # all events with at least one interaction
    H_cluster_position[]       # all events with cluster started (any outcome)
    H_signal[]                 # only SS_in_ROI events (the result heatmap)
    H_total_energy[]           # weighted by deposited E (sum, not count)
/histograms_1D/
    bin_edges_E_keV[]          # 0 to hist_E_max_keV
    E_total_cluster_all_SS[]   # spectrum of all SS events (any cluster position)
    n_interactions_per_photon[]
    path_length_in_LXe_cm[]
/derived/
    P_signal_per_emission_full_LXe       # n_SS_in_ROI / n_emitted (whole volume)
    events_per_decay_full_LXe            # P_signal_full × BR_2448
    events_per_Bq_per_year_full_LXe      # events_per_decay × 3.156e7 s/yr
    events_at_default_activity_full_LXe  # × total_Cu_activity → calibration in full LXe

    # Post-processing FV integral (from H_signal summed inside the FV box)
    n_signal_in_FV_box                   # integer count, sum of H_signal inside FV
    P_signal_per_emission_in_FV          # n_signal_in_FV_box / n_emitted
    events_per_decay_in_FV
    events_at_default_activity_in_FV     # the number directly comparable to XLZD ~0.4
```

### 9.2 Trajectory JSON (`output/trajectories.json`)

For 3D visualization (see Section 10).

### 9.3 Summary text (`output/summary.txt`)

Plain-text printout of:
- All input parameters
- Geometry checks (masses, activities)
- MC counts and percentages
- Calibration result (events at default activity)
- Comparison to XLZD's quoted ~0.4 events from FCR (sanity check)

### 9.4 Default plots

- `output/heatmap_signal.png`: H_signal as a 2D heatmap (z vs r²) with FV box overlaid
- `output/heatmap_first_interaction.png`: H_first_interaction (shows where photons first hit)
- `output/spectrum_SS.png`: 1D energy spectrum of all SS events with ROI window highlighted
- `output/diagnostic.png`: 4-panel grid showing counts breakdown, n_interactions distribution, path-length distribution, etc.

Use Plots.jl with the `gr()` backend for static PNG output.

## 10. 3D Interactive Visualization (Three.js)

A self-contained `viewer.html` file that loads `trajectories.json` and renders:

### 10.1 Scene elements

- **LXe cylinder**: semi-transparent blue (opacity ~0.15), wireframe edges visible
- **Field-cage rings**: thin gold/copper-colored tori at each z_ring position, radial extent [R_lxe, R_lxe + t_radial]
- **FV box** (display-only reference, not a cut): green wireframe cylinder of radius √(fv_r2_max_cm2) ≈ 100 cm, z ∈ [fv_z_min_cm, fv_z_max_cm]
- **Coordinate axes** with labels (z = drift direction)

### 10.2 Trajectory rendering

For each trajectory in the JSON:
- **Polyline** from emission through each interaction to final point
- **Spheres** at each interaction point, color-coded by interaction type and sized by log(E_deposited)
- **Trajectory color** by outcome:
  - SS_in_ROI: bright green (the headline signal)
  - SS_outside_ROI: yellow
  - MS_rejected: red
  - escaped: gray
  - outward: very faint gray (hidden by default)

### 10.3 UI controls

- **Outcome filter checkboxes**: toggle each outcome category on/off
- **Opacity sliders**: cylinder, rings, FV box (independent)
- **Trajectory count slider**: show 1 to all loaded trajectories
- **Camera presets**: top view (down z), side view (across z), perspective, FV closeup
- **Click-to-select trajectory**: highlights selected trajectory in white, shows panel with: outcome, n_interactions, total deposited energy, cluster position, energy at each interaction
- **Reset view** button

### 10.4 Implementation

Single `viewer.html` file, ~300 lines including embedded CSS and JS. Uses Three.js loaded from CDN (`https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js`) and OrbitControls.

Loads `trajectories.json` from the same directory via `fetch()`.

## 11. Code Organization

This is a standard Julia Pkg project (`Project.toml` at the root, `src/` for the module, `scripts/` for entry points, `test/` for tests, `data/` for static inputs, `viewer/` for the JS viewer template, `output/` for runtime artifacts).

### 11.1 Directory layout

```
XLZD/
├── Project.toml                # Pkg manifest — declares deps, makes XLZD a proper module
├── Manifest.toml               # auto-generated by Pkg
├── CLAUDE.md                   # this file
├── README.md                   # how to run, brief
├── .gitignore                  # output/, Manifest.toml, *.h5, *.json
├── src/
│   ├── XLZD.jl                 # module root: declares `module XLZD`, includes others
│   ├── geometry.jl             # detector geometry, ring sampling, cylinder intersection
│   ├── physics.jl              # cross sections, Klein-Nishina, direction rotation
│   ├── trajectory.jl           # per-photon trajectory recording for the 3D viewer
│   ├── mc.jl                   # main MC loop, threading, single-photon tracker
│   └── output.jl               # HDF5, JSON, plots, summary.txt
├── scripts/
│   └── run_mc.jl               # entry point: parses CLI args, builds Params, runs MC
├── viewer/
│   └── viewer.html             # self-contained Three.js viewer template
├── data/
│   └── nist.csv                # NIST XCOM table for xenon (provided)
├── output/                     # gitignored, created at runtime
│   ├── results.h5              # main analysis data
│   ├── trajectories.json       # for viewer
│   ├── viewer.html             # copied from viewer/ at runtime
│   ├── summary.txt
│   ├── heatmap_signal.png
│   ├── heatmap_first_interaction.png
│   ├── spectrum_SS.png
│   └── diagnostic.png
└── test/
    ├── runtests.jl             # @testset entry point, includes the rest
    ├── test_geometry.jl        # mass, ring count, cylinder intersection
    ├── test_physics.jl         # Klein-Nishina, rotation, cross-section interp
    └── test_mc.jl              # short MC run, sanity checks on counts
```

Run with:

```bash
julia --project=. -t 8 scripts/run_mc.jl --n-samples 1e8 --act-cu 2.0e-6 --output output/
```

Test with:

```bash
julia --project=. -t 4 test/runtests.jl
```

### 11.2 File responsibilities

#### `src/XLZD.jl` — Module root

- `module XLZD ... end`
- `using` external packages: `Random, Base.Threads, HDF5, JSON3, Plots, StatsBase, DelimitedFiles, Printf`
- `include`s the other source files in dependency order: `geometry.jl`, `physics.jl`, `trajectory.jl`, `mc.jl`, `output.jl`
- `export`s the public API: `Params`, `Geometry`, `XCOMTable`, `Result`, `build_geometry`, `load_xcom`, `run_mc`, `write_outputs`
- Contains no logic — just assembly.

#### `src/geometry.jl` — Detector geometry

**Structs:**
- `Geometry` — derived from `Params` at startup. Holds: `R_lxe, L_lxe, pitch, A_cs, t_radial, w_axial, R_ring_mean, ρ_Cu, ρ_LXe, n_rings, ring_z_positions::Vector{Float64}, total_Cu_mass_kg, total_Cu_activity_Bq_at_default, LXe_mass_kg`. Immutable.

**Functions:**
- `build_geometry(params::Params) -> Geometry` — derives `n_rings`, `ring_z_positions`, ring volume/mass, total Cu mass, LXe mass. Validates that all rings fit inside `L_lxe`. Prints sanity-check summary to stdout.
- `sample_emission_point(rng, geom::Geometry) -> (x, y, z, ring_idx)` — samples one emission point uniformly in the toroidal ring volumes (uniform pick of ring, uniform azimuth, uniform offset within the cross-section rectangle).
- `path_to_cylinder_exit(x, y, z, dx, dy, dz, R_lxe, L_lxe) -> Float64` — ray–cylinder intersection. Returns the path length until the photon exits via side wall, top, or bottom (whichever comes first). Returns 0.0 if the ray misses or starts outside.

#### `src/physics.jl` — Cross sections, Klein-Nishina, rotation

**Structs:**
- `XCOMTable` — derived from the NIST CSV at startup. Holds: `energy_MeV::Vector{Float64}`, `σ_photo::Vector{Float64}`, `σ_Compton::Vector{Float64}`, `log_E::Vector{Float64}` (precomputed for log-log interp), `log_σ_photo::Vector{Float64}`, `log_σ_Compton::Vector{Float64}`. Immutable.

**Functions:**
- `load_xcom(path::String) -> XCOMTable` — reads `data/nist.csv`. Drops the below-K-edge duplicate row at 34.56 keV (uses only the above-edge value). Sets up log-log interpolation grids. Errors during load if the table doesn't bracket `phys_E_initial_MeV`.
- `σ_photo(table::XCOMTable, E_MeV::Float64) -> Float64` — log-log interpolated photoelectric cross section in cm²/g. Errors if `E_MeV < 30 keV`.
- `σ_Compton(table::XCOMTable, E_MeV::Float64) -> Float64` — log-log interpolated Compton (incoherent) cross section in cm²/g.
- `μ_total_lin(table, E_MeV, ρ_LXe) -> Float64` — total linear attenuation in cm⁻¹.
- `sample_klein_nishina(rng, E_MeV) -> (E_scattered_MeV, cos_theta)` — samples scattered photon energy from Klein-Nishina via Kahn's rejection method. Returns `cos(θ)` derived from the Compton kinematic formula `cos θ = 1 − m_e c² × (1/E' − 1/E)`.
- `rotate_direction(dx, dy, dz, cos_theta, phi) -> (dx', dy', dz')` — numerically stable rotation of a unit vector by polar angle `θ` (specified by `cos θ`) around its current direction with azimuth `phi`. Builds an orthogonal basis perpendicular to the input and rotates within that basis.

#### `src/trajectory.jl` — Trajectory recording for the 3D viewer

**Structs:**
- `InteractionPoint` — `(x, y, z, E_after, type::Symbol, E_dep)` where `type ∈ (:emit, :compton, :photo, :cutoff)`.
- `Trajectory` — holds emission point, initial direction, `Vector{InteractionPoint}`, outcome `Symbol`, cluster info `(z, x, y, E_total)` (or `nothing` if no cluster started), `n_interactions::Int`.
- `TrajectoryBuffer` — per-thread storage. Fixed-capacity buffer per outcome category: `Dict{Symbol, Vector{Trajectory}}`.

**Functions:**
- `should_record(buffer::TrajectoryBuffer, outcome::Symbol, params::Params) -> Bool` — returns true if this thread still has space in its buffer for `outcome`. Called early in `track_one_photon!` so the MC knows whether to allocate the (otherwise expensive) trajectory record. Note: outcome is only known at the end, so in practice the MC always allocates a small trajectory record while tracking; on completion, it's discarded if the buffer is full for its outcome.
- `commit!(buffer::TrajectoryBuffer, traj::Trajectory)` — stores the completed trajectory in the slot for its outcome, if there's room.
- `merge_buffers(buffers::Vector{TrajectoryBuffer}, params::Params) -> Vector{Trajectory}` — merges per-thread buffers into the final list, capped at `mc_n_traj_per_outcome` per outcome globally.
- `to_json(trajectories, geom::Geometry, params::Params) -> String` — serializes scene metadata + trajectories into the JSON format consumed by `viewer.html`.

#### `src/mc.jl` — Main MC loop and threading

**Structs:**
- `ThreadResult` — per-thread accumulator. Holds: counts `Dict{Symbol, Int}`, 2D histograms (`H_first_interaction`, `H_cluster_position`, `H_signal`, `H_total_energy`), 1D histograms (`E_total_cluster_all_SS`, `n_interactions_per_photon`, `path_length_in_LXe`), `TrajectoryBuffer`, runtime accumulator.
- `Result` — merged across all threads. Holds `params::Params`, `geometry::Geometry`, summed counts, summed histograms, merged trajectories, total runtime.

**Functions:**
- `track_one_photon!(state::ThreadResult, rng, geom::Geometry, xcom::XCOMTable, params::Params) -> Symbol` — the core single-photon tracker. Implements emission, isotropic direction sampling, backward rejection (`dot_radial > 0`), step-length sampling, cylinder intersection, photoelectric vs Compton branching, Klein-Nishina + direction rotation, z-only SS/MS clustering, tracking cutoff, FV/ROI checks. Updates `state` in place. Returns the outcome symbol.
- `run_mc(params::Params) -> Result` — top-level driver. Calls `build_geometry(params)` and `load_xcom(params.phys_xcom_data_path)` once. Allocates per-thread `ThreadResult`s with seeded `MersenneTwister(params.mc_seed + tid)`. Dispatches `Threads.@threads` over the partitioned sample count. After threads finish, calls `merge_results`. Returns the merged `Result`.
- `merge_results(thread_results::Vector{ThreadResult}, geom::Geometry, params::Params, runtime_seconds::Float64) -> Result` — sums histograms element-wise, sums counts, merges trajectory buffers via `merge_buffers`.

**Note on signature design:** `track_one_photon!` takes `Geometry`, `XCOMTable`, and `Params` as **separate arguments** rather than bundling derived data into `Params`. This makes the dependency explicit at the function signature: `Params` is purely user inputs, `Geometry` and `XCOMTable` are derived once, and `ThreadResult` is per-thread mutable state. Tests for individual physics functions can construct minimal versions of any of these without needing the others.

#### `src/output.jl` — Persistence and plots

**Functions:**
- `integrate_in_FV_box(H::Matrix, bin_edges_z, bin_edges_r2, params::Params) -> Int` — sums `H[i, j]` over all bins whose centers fall inside the FV box (`fv_z_min_cm ≤ z_center ≤ fv_z_max_cm` and `r²_center ≤ fv_r2_max_cm2`). Used to compute `n_signal_in_FV_box` from `H_signal`.
- `write_h5(result::Result, path::String)` — writes the HDF5 file with the structure specified in Section 9.1. Computes the post-processing FV integrals and writes them under `/derived/`.
- `write_summary(result::Result, path::String)` — human-readable `summary.txt`: parameters, geometry checks, outcome counts and percentages, calibration in full LXe AND in FV box, comparison vs XLZD.
- `write_trajectories_json(trajectories, geom, params, path::String)` — writes `trajectories.json`.
- `copy_viewer(template_path::String, output_path::String)` — copies the static `viewer/viewer.html` to `output/viewer.html`.
- `plot_heatmap_signal(result, path)` — `heatmap_signal.png`. **Heatmap covers full LXe cylinder** (z ∈ [0, L_lxe], r² ∈ [0, R_lxe²]). Overlay: FV box drawn as a green rectangle. Annotation: total events full-LXe and total inside-FV.
- `plot_heatmap_first_interaction(result, path)` — `heatmap_first_interaction.png`, also full cylinder, also with FV overlay.
- `plot_spectrum_SS(result, path)` — 1D energy spectrum of all SS events with ROI window highlighted.
- `plot_diagnostic(result, path)` — 4-panel diagnostic.
- `write_outputs(result::Result, output_dir::String)` — top-level: `mkpath(output_dir)`, then calls all of the above in sequence.

#### `scripts/run_mc.jl` — Entry point

Tiny — argument parsing and the `run_mc` → `write_outputs` call. Roughly:

```julia
using XLZD, ArgParse

function parse_cli()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--n-samples"; arg_type = Int;     default = 100_000_000
        "--act-cu";    arg_type = Float64; default = 2.0e-6
        "--seed";      arg_type = Int;     default = 1234
        "--output";    default = "output/"
        "--xcom";      default = "data/nist.csv"
        "--n-traj";    arg_type = Int;     default = 100
        # ... rest of mappings from Section 7.3 ...
    end
    return parse_args(s)
end

function main()
    args = parse_cli()
    params = Params(
        mc_N_samples           = args["n-samples"],
        act_Cu                 = args["act-cu"],
        mc_seed                = args["seed"],
        out_dir                = args["output"],
        phys_xcom_data_path    = args["xcom"],
        mc_n_traj_per_outcome  = args["n-traj"],
        # ... etc
    )
    result = XLZD.run_mc(params)
    XLZD.write_outputs(result, params.out_dir)
end

main()
```

#### `viewer/viewer.html` — Static template

Single self-contained HTML file with embedded CSS and JavaScript. Loads Three.js + OrbitControls from CDN. Fetches `./trajectories.json` from the same directory at startup. Renders the scene per Section 10. No build step. The script `copy_viewer` copies this file to `output/viewer.html` so the user opens it next to the JSON.

#### `data/nist.csv` — NIST XCOM table

Static input. Already provided by user. Format described in Section 4.1.

#### `test/runtests.jl` and `test/test_*.jl`

Tests as described in Section 13. `runtests.jl` does:

```julia
using Test, XLZD

@testset "XLZD" begin
    include("test_geometry.jl")
    include("test_physics.jl")
    include("test_mc.jl")
end
```

### 11.3 Internal call graph

```
scripts/run_mc.jl
  └─ XLZD.run_mc(params)                                     [src/mc.jl]
       ├─ build_geometry(params)                             [src/geometry.jl]
       ├─ load_xcom(params.phys_xcom_data_path)              [src/physics.jl]
       └─ @threads { track_one_photon!(...) }                [src/mc.jl]
            ├─ sample_emission_point(rng, geom)              [src/geometry.jl]
            ├─ path_to_cylinder_exit(...)                    [src/geometry.jl]
            ├─ σ_photo(xcom, E), σ_Compton(xcom, E)          [src/physics.jl]
            ├─ sample_klein_nishina(rng, E)                  [src/physics.jl]
            ├─ rotate_direction(...)                         [src/physics.jl]
            └─ commit!(buffer, traj)                         [src/trajectory.jl]
       merge_results(...) → Result                           [src/mc.jl]
  └─ XLZD.write_outputs(result, output_dir)                  [src/output.jl]
       ├─ write_h5
       ├─ write_summary
       ├─ write_trajectories_json
       ├─ copy_viewer
       └─ plot_*
```

## 12. Dependencies (Project.toml)

```toml
[deps]
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
ArgParse = "c7e460c6-2fb9-53a9-8c5b-16f535851c63"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
```

## 13. Tests (must pass before declaring done)

### 13.1 Geometry

- `test_geometry.jl`: with defaults, total Cu mass within 1% of 1202 kg, N_rings = 119, LXe mass within 1% of 61.6 t
- Ring positions monotonic, all within [0, L_lxe]

### 13.2 Physics

- `test_physics.jl`:
  - Klein-Nishina: scattered E' < E always; cos(θ) consistent with E' relation; for many samples, E_dep distribution has Compton edge at the right place
  - Rotation: |new_direction| ≈ 1; new_direction · old_direction ≈ cos(θ)
  - Cross-section interpolation: σ_total(2.448 MeV) within 2% of NIST value (~0.0388 cm²/g)
  - Photoelectric branching at 2.448 MeV: ~2–3% (within 0.5% absolute)

### 13.3 MC sanity checks

- `test_mc.jl`: with N_samples = 10⁵ (short run):
  - Sum of all outcome counts = N_emitted (n_outward + n_escaped + n_MS_rejected + n_SS_outside_ROI + n_SS_in_ROI)
  - n_outward ≈ 0.5 × N_emitted (within 1%, since emission is isotropic and rings are just outside cylinder)
  - Mean n_interactions for non-rejected, non-outward, non-escaped photons: 3–8
  - Energy spectrum of all SS events shows a clear peak at 2.448 MeV (full-energy peak from photoelectric or full Compton chains terminating in photo)
  - n_SS_in_ROI > 0 and is a small fraction (~10⁻³ to 10⁻⁴ of N_emitted)
  - Heatmaps not all-zero, not all-NaN
  - Post-processing FV integral: `n_signal_in_FV_box ≤ n_SS_in_ROI` (FV is a subset of full LXe)

### 13.4 Self-consistency check vs XLZD

After running with defaults (1e8 samples, act_Cu = 2 μBq/kg, full XLZD geometry):

- The post-processing FV integral `events_at_default_activity_in_FV` should be **on the order of 0.4 events** for the full XLZD exposure (assumes the FV box defaults match XLZD's reported FV).
- The full-LXe number `events_at_default_activity_full_LXe` will be larger than the FV-restricted one by ~3–10× (most events are near the ring wall).
- If the FV-restricted result is wildly different from 0.4 (factor 5+), something is wrong. If it's within a factor 2, that's already plausibly fine — the prescription (assume same self-shielding as PMTs) was a rough approximation by the user.

This is a calibration check, not a hard test. Print the comparison clearly in summary.txt.

## 14. Implementation Notes for Claude Code

Build the project in the order below. At each step, the corresponding tests in `test/` should pass before moving on.

1. **Project skeleton.** Create `Project.toml` (`Pkg.activate(".")`, `Pkg.add` each dep listed in Section 12). Create empty source files matching Section 11.1. Create `.gitignore` (output/, Manifest.toml, *.h5, *.json under output/). Create `README.md` with the run command from Section 11.1.

2. **`src/geometry.jl`.** Implement `Params` (Section 7.1) in `src/XLZD.jl` first (or in a small `params.jl` if you prefer; CLAUDE.md doesn't mandate a separate file for it). Then implement `Geometry`, `build_geometry`, `sample_emission_point`, `path_to_cylinder_exit`. Write `test/test_geometry.jl` and verify: total Cu mass ~1202 kg with defaults, n_rings = 119, ring positions monotonic and within [0, L_lxe], cylinder intersection returns sensible values for a few hand-checked rays.

3. **`src/physics.jl`.** Implement `XCOMTable`, `load_xcom`, `σ_photo`, `σ_Compton`, `μ_total_lin`. Test: `σ_total(2.448 MeV) ≈ 0.0336 cm²/g` (Compton + photoelectric only, since pair is excluded), photoelectric branching ~2.6%. Then implement `sample_klein_nishina` and `rotate_direction`. Test: KN scattered E always less than incident E; angular distribution matches differential cross section by χ² test on a histogram; rotation preserves unit norm and `cos(θ)` between input and output.

4. **`src/trajectory.jl`.** Implement the structs and per-thread buffer. Trivially testable on its own.

5. **Single-threaded MC.** Implement `track_one_photon!` in `src/mc.jl`. Run by hand on a few photons (10–20) with verbose tracing — print every interaction. Manually verify: energies are monotone non-increasing, interaction types follow branching ratios at each E, SS/MS classification is correct, FV/ROI cuts fire correctly. Write `test/test_mc.jl` for a short run (10⁵ samples) and check outcome counts add up to N_emitted, n_outward ≈ 0.5 × N, energy spectrum has a peak at 2.448 MeV.

6. **Threading.** Add `Threads.@threads` partitioning, per-thread storage, `merge_results`. Verify that single-threaded and 8-threaded runs give statistically identical histograms.

7. **`src/output.jl`.** HDF5 writer first (this is the deliverable); then the plots; then the JSON; then the summary. Each has straightforward implementation.

8. **`viewer/viewer.html`.** Write the static template last. Test by manually constructing a small `trajectories.json` (5–10 trajectories) and opening the HTML file in a browser. Verify all UI controls work.

9. **Sanity checks (Section 13.4).** Run with full defaults, compare events output to XLZD's quoted ~0.4 events. If far off, investigate.

10. **Documentation.** Inline comments for non-obvious physics: "Kahn's algorithm for Klein-Nishina, valid for κ = E/m_e c² ≥ 0.1"; reference NIST XCOM for cross-section data source; reference ENSDF for branching ratio. Brief module-level docstring on each `src/*.jl` file describing its responsibility per Section 11.2.

## 15. Context (for Claude Code: where this fits)

This MC is part of a paper-in-preparation about the ITACA detector, but this specific calculation is auxiliary — it's a cross-check on the XLZD background budget. The user has reasoned (in conversation) that:

1. XLZD's quoted ~0.4 events from "field cage rings" implicitly assumes 75% radiopurity reduction from LZ-as-built.
2. After this reduction, the irreducible cathode-region background is dominated by PTFE walls + Cu rings (resistors and other reducible components are essentially zeroed in the projection).
3. Plated ²²²Rn daughters at the cathode/walls add to this same surface budget linearly. Even at XLZD's design 0.1 μBq/kg target, plated Rn contributes ~0.79 events using a "PTFE-anchored" calibration (0.22 events/mBq).
4. This MC will refine that number by computing the actual γ-attenuation-corrected events per mBq of FCR Bi-214 in the FV cut, including realistic SS/MS classification and ROI energy cuts.

The user already walked through:
- NEXT field-cage geometry (Cu rings, 177.9 kg total)
- Scaled to XLZD (1184 kg Cu, ~10 kg per ring × 119 rings)
- 2 μBq/kg copper specific activity → 2.4 mBq total → 0.4 events under user's "events scale linearly with activity" assumption

This MC will compute the actual events per mBq directly from physics, and the resulting calibration will be applied to plated-Rn-derived activities at various ²²²Rn concentrations (0.1, 0.5, 0.9, 5 μBq/kg) to estimate plated-Rn contribution to the 0νββ background.

## 16. Open Questions / Optional Extensions (NOT required for v1)

- Add 511 keV annihilation γ tracking from pair production (currently skipped).
- Add other Bi-214 γ lines (1.764 MeV, etc.) to study contributions outside the ROI.
- Add PTFE walls as a second source population (similar geometry, distributed mass at R = R_lxe).
- Implement a more realistic detector response (energy resolution Gaussian smearing per deposit, position smearing).
- Compare against a Geant4 reference simulation if available.

These are noted for future work; not in scope for v1.

---

**End of CLAUDE.md.** The script + viewer + tests should be deliverable as a complete project that runs `julia --project=. -t 8 scripts/run_mc.jl` and produces all outputs in `output/`.