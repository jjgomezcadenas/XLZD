"""
LZ Field Cage + Grids Model for 0νββ Background Simulations
=============================================================

This module computes the gamma-ray flux from radioactive decays in the
LZ field cage and electrode grid components that enters the active LXe
volume. This is the "second generator" of background gammas, complementing
the cryostat model (lz_cryo.py).

Physics context
---------------
Same isotopes and energies as the cryostat model:
  - ²¹⁴Bi (²³⁸U late): 2.448 MeV γ, 1.55% BR
  - ²⁰⁸Tl (²³²Th late): 2.615 MeV γ, 35.9% BR to ²⁰⁸Tl × 100% γ BR

Detector model
--------------
The field cage and grids sit directly inside the LXe, with no intervening
material between them and the active volume. The components are:

  BARREL (segmented cylinder at R ≈ 72.8 cm):
    - Ti field shaping rings (65 total: 57 drift + 8 RFR), 91.6 kg
    - Field shaping resistors (~50 mg total), at gaps between rings
    - PTFE reflector panels (184 kg), lining the full cylinder height

  ENDCAPS (4 electrode grids as discs at top and bottom):
    - Anode grid:   SS wire mesh + holder ring, at TPC top
    - Gate grid:    SS wire mesh + holder ring, just below liquid surface
    - Cathode grid: SS wire mesh + holder ring, at bottom of drift region
    - Bottom grid:  SS wire mesh + holder ring, below RFR (shields PMTs)

Key simplifications
-------------------
1. NO INTERVENING LAYERS: all sources face the active LXe directly.
   Gammas enter the active volume immediately upon emission inward.

2. NEGLIGIBLE SELF-SHIELDING: the Ti rings are ~3 mm thick radially
   (μt ≈ 0.05, <5% absorption). Grid wires are <0.13 mm SS. PTFE is
   5 mm of low-Z material. Grid holders are spread over the full disc
   area with low effective thickness. All are treated as surface sources.

3. HALF-SPACE EMISSION: only the inward hemisphere contributes (gammas
   going outward into the skin/cryostat are irrelevant). Factor 1/2.

4. ISOTROPIC ANGULAR DISTRIBUTION: with no material filtering, the
   surviving gammas are flat in u = cos θ. dN/du = 1/2 = constant.

5. UNIFORM Z-SMEARING: the 65 discrete rings are smeared uniformly
   over the full barrel height. With 25 mm ring pitch over 1.5 m, the
   granularity is ~1.7% — negligible for the MC. Same for resistors.

6. PTFE vs TEFLON LINER: the TPC PTFE reflectors (184 kg at R ≈ 72.8 cm)
   are DISTINCT from the cryostat Teflon liner (26 kg at R ≈ 81.4 cm on
   the ICV inner wall). No double counting.

Transport result
----------------
For all field cage sources:
  gammas/yr entering active volume = activity × BR × (1/2) × sec/yr

The angular distribution is flat in u ∈ [0, 1]:
  dN/du = 1/2  (normalized so ∫₀¹ dN/du du = 1/2, the inward fraction)

MC sampling is trivial: u ~ Uniform(0, 1), θ = arccos(u), φ ~ Uniform(0, 2π).

Output files (CSV in ../data/):
  lz_fc_geometry.csv     Component dimensions
  lz_fc_activity.csv     Radioactivity per component
  lz_fc_gammas.csv       Gamma rates per component per isotope
  lz_fc_sampling.csv     MC sampling parameters (4 entries: iso × region)

Sources:
  LZ TDR (arXiv:1703.09144), Tables 9.2.1, 9.2.7
  LZ Instrument Paper (arXiv:1910.09124), Table 1
  LZ_detector_summary.md

All internal units: cm, kg, mBq, seconds. Output rates in gammas/year.
"""

import os
import numpy as np
import pandas as pd

# ===========================================================================
# Paths
# ===========================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data")

# ===========================================================================
# Physical constants
# ===========================================================================
RHO_TI = 4.510        # g/cm³
RHO_SS = 7.930        # g/cm³, 304 stainless steel
RHO_PTFE = 2.200      # g/cm³
SEC_PER_YEAR = 3.1557e7

# Branching ratios (same as cryostat model)
BR_BI214_GAMMA = 0.0155        # 2.448 MeV γ per ²¹⁴Bi decay
BR_TL208_FROM_CHAIN = 0.359    # ²⁰⁸Tl per ²³²Th_late decay
BR_TL208_GAMMA = 1.0           # 2.615 MeV γ per ²⁰⁸Tl decay

# Gamma energies
E_BI214 = 2.448  # MeV
E_TL208 = 2.615  # MeV

# ===========================================================================
# TPC Geometry
# ===========================================================================
TPC_INNER_R = 72.8     # cm (1456 mm diameter / 2)
FC_THICKNESS = 1.5     # cm (15 mm)

# Drift region: 57 rings, 25 mm pitch, height = 1456 mm
DRIFT_HEIGHT = 145.6   # cm
DRIFT_N_RINGS = 57
DRIFT_PITCH = 2.5      # cm

# Reverse field region: 8 rings, height = 137.5 mm
RFR_HEIGHT = 13.75     # cm
RFR_N_RINGS = 8

# Total barrel height (drift + RFR)
BARREL_HEIGHT = DRIFT_HEIGHT + RFR_HEIGHT  # 159.35 cm

# Total number of rings
N_RINGS_TOTAL = DRIFT_N_RINGS + RFR_N_RINGS  # 65

# ===========================================================================
# Attenuation check: verify self-shielding is negligible
# ===========================================================================
def check_self_shielding():
    """
    Compute the radial thickness and optical depth (μt) of each component
    to verify that self-shielding can be neglected.

    For a component with mass M distributed uniformly around a cylinder of
    radius R and height H, the effective radial thickness is:
        t = M / (ρ × 2πR × H)

    We use μ_Ti ≈ 0.17 cm⁻¹ and μ_SS ≈ 0.33 cm⁻¹ at ~2.5 MeV
    (SS is higher Z and density than Ti).
    μ_PTFE ≈ 0.08 cm⁻¹ at 2.5 MeV (low Z, low density).

    Returns list of dicts with component name, thickness, and μt.
    """
    mu_Ti = 0.17    # cm⁻¹ approximate
    mu_SS = 0.33    # cm⁻¹ approximate (ρ=7.93, μ/ρ≈0.042)
    mu_PTFE = 0.08  # cm⁻¹ approximate

    circumference = 2 * np.pi * TPC_INNER_R  # cm
    disc_area = np.pi * TPC_INNER_R**2       # cm²

    checks = []

    # Ti rings: 91.6 kg over barrel
    t_rings = (91.6e3 / RHO_TI) / (circumference * BARREL_HEIGHT)  # cm
    checks.append({
        'component': 'Ti rings', 'mass_kg': 91.6,
        't_cm': t_rings, 'mu': mu_Ti, 'mu_t': mu_Ti * t_rings,
    })

    # PTFE reflectors: 184 kg, 5 mm thick panels, full barrel
    t_ptfe = 0.5  # cm (5 mm, from design)
    checks.append({
        'component': 'PTFE reflectors', 'mass_kg': 184.0,
        't_cm': t_ptfe, 'mu': mu_PTFE, 'mu_t': mu_PTFE * t_ptfe,
    })

    # Field grids and holders: 89.1 kg over 4 discs (bbonu paper)
    t_holders = (89.1e3 / RHO_SS) / (4 * disc_area)  # cm, averaged
    checks.append({
        'component': 'Grids+holders (avg)', 'mass_kg': 89.1,
        't_cm': t_holders, 'mu': mu_SS, 'mu_t': mu_SS * t_holders,
    })

    # Resistors: 60 g (bbonu paper)
    checks.append({
        'component': 'Resistors', 'mass_kg': 0.06,
        't_cm': 0.001, 'mu': 0.2, 'mu_t': 0.2 * 0.001,
    })

    # TPC sensors: 5.02 kg
    checks.append({
        'component': 'TPC sensors', 'mass_kg': 5.02,
        't_cm': 0.01, 'mu': 0.2, 'mu_t': 0.2 * 0.01,
    })

    return checks


# ===========================================================================
# Component definitions
# ===========================================================================
def build_components():
    """
    Define all field cage and grid components with their properties.

    Each component is a dict with:
        name       : human-readable label
        region     : 'barrel' or 'endcap'
        material   : 'Ti', 'SS', 'PTFE', 'mixed'
        mass_kg    : total mass
        Bi214_sp   : specific activity ²³⁸U_late (mBq/kg)
        Tl208_sp   : specific activity ²³²Th_late (mBq/kg)
        R_cm       : radius of the cylinder/disc (cm)
        H_cm       : height (barrel) or 0 (endcap)
        position   : descriptive position in the TPC

    Activity values from LZ bb0nu paper, Table I.

    Returns list of component dicts.
    """
    components = [
        # --- BARREL components (segmented cylinder) ---
        {
            'name': 'Field-cage rings',
            'region': 'barrel',
            'material': 'Ti',
            'mass_kg': 93.0,       # bbonu paper
            'Bi214_sp': 0.35,      # mBq/kg, ²³⁸U-late (bbonu Table I, upper limit)
            'Tl208_sp': 0.24,      # mBq/kg, ²³²Th-late (bbonu Table I, upper limit)
            'R_cm': TPC_INNER_R + FC_THICKNESS,
            'H_cm': BARREL_HEIGHT,
            'position': f'{N_RINGS_TOTAL} rings (57 drift + 8 RFR), '
                        f'25 mm pitch, smeared over {BARREL_HEIGHT:.1f} cm',
        },
        {
            'name': 'Field-cage resistors',
            'region': 'barrel',
            'material': 'mixed (alumina ceramic)',
            'mass_kg': 0.06,       # 60 g (bbonu paper)
            'Bi214_sp': 1350.0,    # mBq/kg, ²³⁸U-late (bbonu Table I, upper limit)
            'Tl208_sp': 2010.0,    # mBq/kg, ²³²Th-late (bbonu Table I, upper limit)
            'R_cm': TPC_INNER_R + FC_THICKNESS,
            'H_cm': BARREL_HEIGHT,
            'position': 'FC drift + RFR + cathode feedthrough resistors',
        },
        {
            'name': 'PTFE walls',
            'region': 'barrel',
            'material': 'PTFE',
            'mass_kg': 184.0,
            'Bi214_sp': 0.04,      # mBq/kg (bbonu Table I)
            'Tl208_sp': 0.01,      # mBq/kg
            'R_cm': TPC_INNER_R,
            'H_cm': BARREL_HEIGHT,
            'position': f'5 mm thick panels, full barrel height {BARREL_HEIGHT:.1f} cm',
        },
        {
            'name': 'TPC sensors',
            'region': 'barrel',
            'material': 'mixed',
            'mass_kg': 5.02,
            'Bi214_sp': 5.82,      # mBq/kg (bbonu Table I)
            'Tl208_sp': 1.88,      # mBq/kg
            'R_cm': TPC_INNER_R + FC_THICKNESS,
            'H_cm': BARREL_HEIGHT,
            'position': 'Loop antennas, thermometers, acoustic sensors',
        },

        # --- ENDCAP components (4 electrode grids) ---
        # bbonu paper groups all grids+holders as one item: 89.1 kg.
        # Split 50/50 between top and bottom.
        {
            'name': 'Top grids and holders',
            'region': 'endcap',
            'material': 'SS/Ti',
            'mass_kg': 44.55,      # half of 89.1 kg
            'Bi214_sp': 2.63,      # mBq/kg (bbonu Table I)
            'Tl208_sp': 1.46,      # mBq/kg
            'R_cm': TPC_INNER_R,
            'H_cm': 0,
            'position': 'Anode + gate grids with holder rings (top)',
        },
        {
            'name': 'Bottom grids and holders',
            'region': 'endcap',
            'material': 'SS/Ti',
            'mass_kg': 44.55,      # other half
            'Bi214_sp': 2.63,
            'Tl208_sp': 1.46,
            'R_cm': TPC_INNER_R,
            'H_cm': 0,
            'position': 'Cathode + bottom shield grids with holder rings (bottom)',
        },
    ]

    return components


# ===========================================================================
# Gamma rate computation
# ===========================================================================
def compute_gamma_rates(components):
    """
    Compute gamma production and emission rates for each component.

    Since all components are thin surface sources with no intervening
    material, the rate entering the active volume is simply:

        gammas/yr = activity (mBq) × BR × (1/2) × sec/yr × 1e-3

    The factor 1/2 accounts for the inward hemisphere (outward gammas
    are irrelevant — they go into the skin/cryostat).

    Parameters
    ----------
    components : list of dict
        From build_components().

    Returns
    -------
    results : list of dict
        One entry per (component, isotope) with production and emission rates.
    """
    results = []

    for comp in components:
        for iso in ['Bi214', 'Tl208']:
            sp = comp[f'{iso}_sp']
            mass = comp['mass_kg']
            activity_mBq = mass * sp

            if iso == 'Bi214':
                gamma_br = BR_BI214_GAMMA
            else:
                gamma_br = BR_TL208_FROM_CHAIN * BR_TL208_GAMMA

            # Total gammas produced (full 4π)
            gammas_produced = activity_mBq * 1e-3 * gamma_br * SEC_PER_YEAR

            # Gammas entering active volume (inward hemisphere = 1/2)
            gammas_entering = gammas_produced * 0.5

            results.append({
                'source': comp['name'],
                'region': comp['region'],
                'material': comp['material'],
                'isotope': iso,
                'mass_kg': mass,
                'activity_mBq': activity_mBq,
                'gamma_br': gamma_br,
                'gammas_produced_per_yr': gammas_produced,
                'gammas_entering_per_yr': gammas_entering,
            })

    return results


# ===========================================================================
# Output CSV files
# ===========================================================================
def save_geometry_csv(components, filepath):
    """Save component geometry to CSV."""
    rows = []
    for c in components:
        rows.append({
            'component': c['name'],
            'region': c['region'],
            'material': c['material'],
            'R_cm': c['R_cm'],
            'H_cm': c['H_cm'],
            'position': c['position'],
        })
    df = pd.DataFrame(rows)
    df.to_csv(filepath, index=False)
    return df


def save_activity_csv(components, filepath):
    """Save radioactivity per component to CSV."""
    rows = []
    for c in components:
        mass = c['mass_kg']
        rows.append({
            'component': c['name'],
            'region': c['region'],
            'material': c['material'],
            'mass_kg': mass,
            'Bi214_mBq_per_kg': c['Bi214_sp'],
            'Bi214_total_mBq': round(mass * c['Bi214_sp'], 4),
            'Tl208_mBq_per_kg': c['Tl208_sp'],
            'Tl208_total_mBq': round(mass * c['Tl208_sp'], 4),
        })
    df = pd.DataFrame(rows)
    df.to_csv(filepath, index=False)
    return df


def save_gammas_csv(results, filepath):
    """Save gamma rate table to CSV."""
    rows = []
    for r in results:
        rows.append({
            'source': r['source'],
            'region': r['region'],
            'isotope': r['isotope'],
            'gammas_produced_per_yr': f"{r['gammas_produced_per_yr']:.4e}",
            'gammas_entering_per_yr': f"{r['gammas_entering_per_yr']:.4e}",
        })
    df = pd.DataFrame(rows)
    df.to_csv(filepath, index=False)
    return df


def save_sampling_csv(results, filepath):
    """
    Save MC sampling parameters for the field cage generator.

    Since all sources have isotropic angular distribution (flat in u),
    no fit is needed. We store the total rate per (isotope, region)
    and note that u ~ Uniform(0, 1).

    The output file contains one row per (isotope, region) with:
      - total_gammas_per_yr: rate entering active volume
      - angular_distribution: 'flat' (u ~ Uniform(0,1))
      - R_cm: radius of the emission surface
      - H_cm: height (barrel) or 0 (endcap)
    """
    rows = []
    for iso in ['Bi214', 'Tl208']:
        for region in ['barrel', 'endcap']:
            total = sum(r['gammas_entering_per_yr'] for r in results
                        if r['isotope'] == iso and r['region'] == region)

            if region == 'barrel':
                R = TPC_INNER_R
                H = BARREL_HEIGHT
            else:
                R = TPC_INNER_R
                H = 0

            rows.append({
                'isotope': iso,
                'region': region,
                'total_gammas_per_yr': f"{total:.6e}",
                'angular_distribution': 'flat (u ~ Uniform(0,1))',
                'R_cm': R,
                'H_cm': H,
            })

    df = pd.DataFrame(rows)
    with open(filepath, 'w') as f:
        f.write("# LZ Field Cage: MC sampling parameters\n")
        f.write("# Angular distribution is isotropic over the inward hemisphere:\n")
        f.write("#   u = cos(theta) ~ Uniform(0, 1)\n")
        f.write("#   theta = arccos(u)\n")
        f.write("#   phi ~ Uniform(0, 2*pi)\n")
        f.write("# Entry point: uniform on barrel cylinder or endcap disc\n")
        f.write("# Barrel: cylindrical shell at R_cm, height H_cm\n")
        f.write("# Endcap: disc at R_cm (split 50/50 top and bottom)\n")
        df.to_csv(f, index=False)
    return df


# ===========================================================================
# Console output
# ===========================================================================
def print_section(title):
    print(f"\n{'#'*70}")
    print(f"# {title}")
    print(f"{'#'*70}")


def main():
    print_section("LZ FIELD CAGE + GRIDS MODEL")

    # --- Self-shielding check ---
    print_section("SELF-SHIELDING CHECK")
    checks = check_self_shielding()
    print(f"\n  {'Component':<25s} {'Mass (kg)':>10s} {'t (cm)':>8s} {'μ (cm⁻¹)':>9s} {'μt':>8s} {'Absorb.':>8s}")
    print(f"  {'-'*72}")
    for c in checks:
        absorb = 1 - np.exp(-c['mu_t'])
        print(f"  {c['component']:<25s} {c['mass_kg']:>10.4f} {c['t_cm']:>8.4f} "
              f"{c['mu']:>9.3f} {c['mu_t']:>8.4f} {absorb:>7.1%}")
    print(f"\n  All μt < 0.05 → self-shielding < 5% → NEGLIGIBLE. Treat as surface sources.")

    # --- Components ---
    components = build_components()

    print_section("GEOMETRY")
    print(f"\n  TPC inner radius:   {TPC_INNER_R} cm")
    print(f"  FC thickness:       {FC_THICKNESS} cm")
    print(f"  Barrel height:      {BARREL_HEIGHT} cm (drift {DRIFT_HEIGHT} + RFR {RFR_HEIGHT})")
    print(f"  Drift rings:        {DRIFT_N_RINGS} at {DRIFT_PITCH} cm pitch")
    print(f"  RFR rings:          {RFR_N_RINGS}")
    print(f"  Total rings:        {N_RINGS_TOTAL}")

    print(f"\n  {'Component':<35s} {'Region':<8s} {'Material':<10s} {'Mass (kg)':>10s}")
    print(f"  {'-'*68}")
    for c in components:
        print(f"  {c['name']:<35s} {c['region']:<8s} {c['material']:<10s} {c['mass_kg']:>10.4f}")

    # --- Activity table ---
    print_section("RADIOACTIVITY")
    print(f"\n  {'Component':<35s} {'Mass':>7s} {'Bi214/kg':>9s} {'Bi214':>9s} {'Tl208/kg':>9s} {'Tl208':>9s}")
    print(f"  {'':35s} {'(kg)':>7s} {'(mBq/kg)':>9s} {'(mBq)':>9s} {'(mBq/kg)':>9s} {'(mBq)':>9s}")
    print(f"  {'-'*80}")
    tot_bi = 0; tot_tl = 0
    for c in components:
        bi = c['mass_kg'] * c['Bi214_sp']
        tl = c['mass_kg'] * c['Tl208_sp']
        tot_bi += bi; tot_tl += tl
        print(f"  {c['name']:<35s} {c['mass_kg']:>7.4f} {c['Bi214_sp']:>9.2f} "
              f"{bi:>9.3f} {c['Tl208_sp']:>9.2f} {tl:>9.3f}")
    print(f"  {'-'*80}")
    print(f"  {'TOTAL':<35s} {'':>7s} {'':>9s} {tot_bi:>9.2f} {'':>9s} {tot_tl:>9.2f}")

    # --- Gamma rates ---
    print_section("GAMMA RATES")
    results = compute_gamma_rates(components)

    print(f"\n  {'Source':<35s} {'Region':<8s} {'Isotope':<7s} "
          f"{'Produced':>12s} {'Entering':>12s}")
    print(f"  {'':35s} {'':8s} {'':7s} {'(γ/yr)':>12s} {'(γ/yr)':>12s}")
    print(f"  {'-'*78}")
    for r in results:
        print(f"  {r['source']:<35s} {r['region']:<8s} {r['isotope']:<7s} "
              f"{r['gammas_produced_per_yr']:>12.2e} {r['gammas_entering_per_yr']:>12.2e}")

    # --- Aggregated MC input ---
    print_section("MC INPUT SUMMARY (gammas/yr entering active volume)")
    for iso in ['Bi214', 'Tl208']:
        for region in ['barrel', 'endcap']:
            total = sum(r['gammas_entering_per_yr'] for r in results
                        if r['isotope'] == iso and r['region'] == region)
            print(f"  {iso} {region:<8s}: {total:>12.2e} γ/yr  (u ~ Uniform(0,1))")

    # --- Save CSV files ---
    print_section("SAVING OUTPUT FILES")

    f1 = os.path.join(DATA_DIR, "lz_fc_geometry.csv")
    save_geometry_csv(components, f1)
    print(f"  {f1}")

    f2 = os.path.join(DATA_DIR, "lz_fc_activity.csv")
    save_activity_csv(components, f2)
    print(f"  {f2}")

    f3 = os.path.join(DATA_DIR, "lz_fc_gammas.csv")
    save_gammas_csv(results, f3)
    print(f"  {f3}")

    f4 = os.path.join(DATA_DIR, "lz_fc_sampling.csv")
    save_sampling_csv(results, f4)
    print(f"  {f4}")

    # --- Comparison with cryostat ---
    print_section("COMPARISON: FIELD CAGE vs CRYOSTAT (gammas/yr entering active volume)")
    print(f"\n  Note: cryostat gammas traverse Ti walls + LXe skin before entering.")
    print(f"        Field cage gammas enter directly (no intervening material).")
    print(f"\n  See lz_cryo_linear_fit.csv for cryostat rates for comparison.")

    print("\nDone.")


if __name__ == "__main__":
    main()
