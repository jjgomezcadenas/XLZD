"""
LZ PMT Model for 0νββ Background Simulations
===============================================

This module computes the gamma-ray flux from radioactive decays in the
LZ PMT systems that enters the active LXe volume. This is the "third
generator" of background gammas, complementing the cryostat (lz_cryo.py)
and field cage (lz_fieldcage.py) models.

PMT groups
----------
Three distinct groups at different positions:

  1. TOP PMT ARRAY (endcap, top):
     253 × R11410-22 (3-inch), in gas phase, downward-looking.
     Plus 253 bases and the upper Ti support structure.
     No attenuation — gammas enter the LXe volume directly.

  2. BOTTOM PMT ARRAY (endcap, bottom):
     241 × R11410-22 (3-inch), in liquid, upward-looking.
     Plus 241 bases and the lower Ti support structure.
     No attenuation — gammas face the RFR LXe directly (the RFR is
     part of the detector volume and handled by the MC propagation).

  3. SKIN PMTs (mixed positions):
     a) 93 × R8520-406 (1-inch), side skin, on FC outer wall (barrel).
        Plus bases. No attenuation (same position as FC components).
     b) 20 × R8778 (2-inch), lower ring at ICV wall (barrel).
        Must traverse ~6 cm LXe skin → shaped angular distribution.
     c) 18 × R8778 (2-inch), dome below bottom array (endcap, bottom).
        Must traverse ~6 cm LXe dome → shaped angular distribution.

Activity units
--------------
CAUTION: R11410 PMT activities are quoted PER PMT (mBq/PMT), not per kg.
All other components (bases, structures, cabling, skin PMTs) are per kg.
The code handles both cases explicitly.

Transport
---------
  - 9 of 11 sources: no attenuation. Factor 1/2 (inward hemisphere).
    Angular distribution flat in u = cos θ ∈ [0, 1].
  - 2 sources (R8778 lower ring + dome): traverse ~6 cm LXe.
    Uses angular_spectrum_surface() from lz_cryo.py for the shaped
    distribution through one LXe slab.

Output files (CSV in ../data/):
  lz_pmt_activity.csv    Radioactivity per component
  lz_pmt_gammas.csv      Gamma rates per component per isotope
  lz_pmt_sampling.csv    MC sampling parameters

Sources:
  LZ Instrument Paper (arXiv:1910.09124), Sections 2.1-2.2
  LZ Radioactivity Paper (arXiv:2006.02506), Table 7
  LZ TDR (arXiv:1703.09144), Table 9.2.7
  LZ_detector_summary.md
"""

import os
import sys
import numpy as np
import pandas as pd

# Import transport functions from the cryostat model
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from lz_cryo import (
    angular_spectrum_surface,
    integrate_spectrum,
    load_attenuation_coefficients,
    U_BINS,
)

DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data")

# ===========================================================================
# Physical constants
# ===========================================================================
SEC_PER_YEAR = 3.1557e7

# Branching ratios (same as cryostat and field cage models)
BR_BI214_GAMMA = 0.0155
BR_TL208_FROM_CHAIN = 0.359
BR_TL208_GAMMA = 1.0

# Gamma energies
E_BI214 = 2.448  # MeV
E_TL208 = 2.615  # MeV

# LXe skin thickness for R8778 attenuation
SKIN_THICKNESS_CM = 6.0

# ===========================================================================
# Geometry reference positions
# ===========================================================================
TPC_INNER_R = 72.8     # cm
FC_OUTER_R = 74.3      # cm (TPC_INNER_R + 1.5 cm FC thickness)
ICV_INNER_R = 81.4     # cm (approximate, from cryostat model)

# ===========================================================================
# Component definitions
# ===========================================================================
def build_components():
    """
    Define all PMT-related components with their properties.

    Each component is a dict with:
        name         : human-readable label
        region       : 'barrel', 'endcap_top', or 'endcap_bottom'
        activity_type: 'per_pmt' or 'per_kg'
        count        : number of units (for per_pmt sources)
        mass_kg      : total mass (for per_kg sources, or total PMT mass)
        Bi214_val    : specific activity — mBq/PMT or mBq/kg
        Tl208_val    : specific activity — mBq/PMT or mBq/kg
        R_cm         : radial position
        attenuation  : 'none' or 'lxe_skin'
        note         : description

    Activity values from LZ bb0nu paper, Table I.
    All activities are now per kg (the paper reports per-kg values).
    The bbonu paper groups PMTs and bases as aggregate items.

    Component mapping from bbonu paper:
      TPC PMTs:           91.9 kg, 3.22/1.61 mBq/kg — split 253/241 top/bottom
      TPC PMT bases:      2.80 kg, 75.9/33.1 mBq/kg — split 253/241 top/bottom
      TPC PMT structures: 166 kg,  1.60/1.06 mBq/kg — split 50/50 top/bottom
      TPC PMT cables:     88.7 kg, 4.31/0.82 mBq/kg — barrel
      Skin PMTs+bases:    8.59 kg, 46.0/14.9 mBq/kg — split side/lower/dome
    """
    # TPC PMT mass split by count (253 top, 241 bottom out of 494)
    pmt_total_kg = 91.9
    pmt_top_kg = pmt_total_kg * 253 / 494
    pmt_bot_kg = pmt_total_kg * 241 / 494

    # Base mass split
    base_total_kg = 2.80
    base_top_kg = base_total_kg * 253 / 494
    base_bot_kg = base_total_kg * 241 / 494

    # Structure mass split 50/50
    struct_total_kg = 166.0
    struct_top_kg = struct_total_kg / 2.0
    struct_bot_kg = struct_total_kg / 2.0

    # Skin PMTs+bases: 8.59 kg total
    # Split: R8520 side (93 PMTs) ~60%, R8778 lower (20) ~21%, R8778 dome (18) ~19%
    skin_total_kg = 8.59
    skin_side_kg = skin_total_kg * 93 / 131
    skin_lower_kg = skin_total_kg * 20 / 131
    skin_dome_kg = skin_total_kg * 18 / 131

    components = [
        # --- TOP PMT ARRAY (endcap, top) ---
        {
            'name': 'Top TPC PMTs',
            'region': 'endcap_top',
            'activity_type': 'per_kg',
            'count': 253,
            'mass_kg': pmt_top_kg,
            'Bi214_val': 3.22,       # mBq/kg (bbonu Table I)
            'Tl208_val': 1.61,
            'R_cm': TPC_INNER_R,
            'attenuation': 'none',
            'note': 'R11410-22 in gas, downward-looking.',
        },
        {
            'name': 'Top PMT bases',
            'region': 'endcap_top',
            'activity_type': 'per_kg',
            'count': 253,
            'mass_kg': base_top_kg,
            'Bi214_val': 75.9,
            'Tl208_val': 33.1,
            'R_cm': TPC_INNER_R,
            'attenuation': 'none',
            'note': 'Kapton PCB + components.',
        },
        {
            'name': 'Top PMT structure',
            'region': 'endcap_top',
            'activity_type': 'per_kg',
            'count': 1,
            'mass_kg': struct_top_kg,
            'Bi214_val': 1.60,
            'Tl208_val': 1.06,
            'R_cm': TPC_INNER_R,
            'attenuation': 'none',
            'note': 'Ti support structure for top array.',
        },

        # --- BOTTOM PMT ARRAY (endcap, bottom) ---
        {
            'name': 'Bottom TPC PMTs',
            'region': 'endcap_bottom',
            'activity_type': 'per_kg',
            'count': 241,
            'mass_kg': pmt_bot_kg,
            'Bi214_val': 3.22,
            'Tl208_val': 1.61,
            'R_cm': TPC_INNER_R,
            'attenuation': 'none',
            'note': 'R11410-22 in liquid, upward-looking.',
        },
        {
            'name': 'Bottom PMT bases',
            'region': 'endcap_bottom',
            'activity_type': 'per_kg',
            'count': 241,
            'mass_kg': base_bot_kg,
            'Bi214_val': 75.9,
            'Tl208_val': 33.1,
            'R_cm': TPC_INNER_R,
            'attenuation': 'none',
            'note': 'Kapton PCB + components.',
        },
        {
            'name': 'Bottom PMT structure',
            'region': 'endcap_bottom',
            'activity_type': 'per_kg',
            'count': 1,
            'mass_kg': struct_bot_kg,
            'Bi214_val': 1.60,
            'Tl208_val': 1.06,
            'R_cm': TPC_INNER_R,
            'attenuation': 'none',
            'note': 'Ti support structure for bottom array.',
        },

        # --- PMT CABLING (barrel) ---
        {
            'name': 'TPC PMT cables',
            'region': 'barrel',
            'activity_type': 'per_kg',
            'count': 1,
            'mass_kg': 88.7,
            'Bi214_val': 4.31,
            'Tl208_val': 0.82,
            'R_cm': FC_OUTER_R,
            'attenuation': 'none',
            'note': 'Cables, 88.7 kg (bbonu paper).',
        },

        # --- SKIN PMTs+bases: side (barrel, no attenuation) ---
        {
            'name': 'Skin PMTs side',
            'region': 'barrel',
            'activity_type': 'per_kg',
            'count': 93,
            'mass_kg': skin_side_kg,
            'Bi214_val': 46.0,
            'Tl208_val': 14.9,
            'R_cm': FC_OUTER_R,
            'attenuation': 'none',
            'note': 'R8520 side skin PMTs+bases (bbonu grouped).',
        },

        # --- SKIN PMTs: R8778 lower ring (barrel, LXe skin attenuation) ---
        {
            'name': 'Skin PMTs lower ring',
            'region': 'barrel',
            'activity_type': 'per_kg',
            'count': 20,
            'mass_kg': skin_lower_kg,
            'Bi214_val': 46.0,
            'Tl208_val': 14.9,
            'R_cm': ICV_INNER_R,
            'attenuation': 'lxe_skin',
            'note': 'R8778 lower ring, traverse LXe skin.',
        },

        # --- SKIN PMTs: dome (endcap bottom, LXe attenuation) ---
        {
            'name': 'Skin PMTs dome',
            'region': 'endcap_bottom',
            'activity_type': 'per_kg',
            'count': 18,
            'mass_kg': skin_dome_kg,
            'Bi214_val': 46.0,
            'Tl208_val': 14.9,
            'R_cm': ICV_INNER_R,
            'attenuation': 'lxe_skin',
            'note': 'R8778 dome PMTs, traverse LXe dome.',
        },
    ]

    return components


# ===========================================================================
# Activity computation
# ===========================================================================
def compute_total_activity(comp, isotope):
    """
    Compute total activity in mBq for a component and isotope.

    Handles both per-PMT and per-kg activity units:
      - per_pmt: total = value × count
      - per_kg:  total = value × mass_kg

    Parameters
    ----------
    comp : dict
        Component definition.
    isotope : str
        'Bi214' or 'Tl208'.

    Returns
    -------
    activity_mBq : float
    """
    val = comp[f'{isotope}_val']
    if comp['activity_type'] == 'per_pmt':
        return val * comp['count']
    else:
        return val * comp['mass_kg']


# ===========================================================================
# Gamma rate computation
# ===========================================================================
def compute_gamma_rates(components, mu):
    """
    Compute gamma production and emission rates for each component.

    For sources with no attenuation:
        gammas/yr entering = activity × BR × (1/2) × sec/yr

    For sources with LXe skin attenuation:
        Uses angular_spectrum_surface() with one LXe slab to compute
        the shaped angular distribution and escape fraction.

    Parameters
    ----------
    components : list of dict
    mu : dict
        Attenuation coefficients from load_attenuation_coefficients().

    Returns
    -------
    results : list of dict
        One entry per (component, isotope).
    """
    results = []

    for comp in components:
        for iso in ['Bi214', 'Tl208']:
            activity_mBq = compute_total_activity(comp, iso)

            if iso == 'Bi214':
                gamma_br = BR_BI214_GAMMA
            else:
                gamma_br = BR_TL208_FROM_CHAIN * BR_TL208_GAMMA

            gammas_produced = activity_mBq * 1e-3 * gamma_br * SEC_PER_YEAR

            if comp['attenuation'] == 'none':
                # No material to traverse. Factor 1/2 for inward hemisphere.
                gammas_entering = gammas_produced * 0.5
                frac_entering = 0.5
                mean_u = 0.5  # isotropic: <u> = <cos θ> = 1/2
                angular_type = 'flat'
                dNdu_exiting = None
            else:
                # LXe skin attenuation: surface source behind one LXe slab.
                mu_key = f'LXe_{iso}'
                slab = [{"label": "LXe_skin", "mu": mu[mu_key],
                         "t_cm": SKIN_THICKNESS_CM}]
                dNdu = angular_spectrum_surface(slab, U_BINS)
                frac_entering, mean_u = integrate_spectrum(dNdu, U_BINS)
                gammas_entering = gammas_produced * frac_entering
                angular_type = 'shaped'
                dNdu_exiting = dNdu

            results.append({
                'source': comp['name'],
                'region': comp['region'],
                'activity_type': comp['activity_type'],
                'isotope': iso,
                'activity_mBq': activity_mBq,
                'gamma_br': gamma_br,
                'gammas_produced_per_yr': gammas_produced,
                'frac_entering': frac_entering,
                'gammas_entering_per_yr': gammas_entering,
                'mean_u': mean_u,
                'angular_type': angular_type,
                'dNdu_exiting': dNdu_exiting,
            })

    return results


# ===========================================================================
# Output CSV files
# ===========================================================================
def save_activity_csv(components, filepath):
    """Save radioactivity per component to CSV."""
    rows = []
    for c in components:
        bi_total = compute_total_activity(c, 'Bi214')
        tl_total = compute_total_activity(c, 'Tl208')
        rows.append({
            'component': c['name'],
            'region': c['region'],
            'activity_type': c['activity_type'],
            'count': c['count'],
            'mass_kg': round(c['mass_kg'], 3),
            'Bi214_val': c['Bi214_val'],
            'Bi214_unit': 'mBq/PMT' if c['activity_type'] == 'per_pmt' else 'mBq/kg',
            'Bi214_total_mBq': round(bi_total, 2),
            'Tl208_val': c['Tl208_val'],
            'Tl208_unit': 'mBq/PMT' if c['activity_type'] == 'per_pmt' else 'mBq/kg',
            'Tl208_total_mBq': round(tl_total, 2),
            'attenuation': c['attenuation'],
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
            'angular_type': r['angular_type'],
            'gammas_produced_per_yr': f"{r['gammas_produced_per_yr']:.4e}",
            'gammas_entering_per_yr': f"{r['gammas_entering_per_yr']:.4e}",
            'mean_cos_theta': f"{r['mean_u']:.4f}",
        })
    df = pd.DataFrame(rows)
    df.to_csv(filepath, index=False)
    return df


def save_sampling_csv(results, filepath):
    """
    Save MC sampling parameters.

    Aggregates rates by (isotope, region). For each combination, notes
    whether the angular distribution is flat or has a shaped contribution.
    Since the shaped sources (R8778) are minor, we note their fraction
    but provide a single effective angular type for the aggregate.
    """
    rows = []
    for iso in ['Bi214', 'Tl208']:
        for region in ['barrel', 'endcap_top', 'endcap_bottom']:
            matching = [r for r in results
                        if r['isotope'] == iso and r['region'] == region]
            if not matching:
                continue

            total_entering = sum(r['gammas_entering_per_yr'] for r in matching)
            flat_entering = sum(r['gammas_entering_per_yr'] for r in matching
                               if r['angular_type'] == 'flat')
            shaped_entering = sum(r['gammas_entering_per_yr'] for r in matching
                                 if r['angular_type'] == 'shaped')
            flat_frac = flat_entering / total_entering if total_entering > 0 else 1.0

            if flat_frac > 0.99:
                ang_type = 'flat (u ~ Uniform(0,1))'
            elif flat_frac < 0.01:
                ang_type = 'shaped (LXe skin attenuation)'
            else:
                ang_type = f'mixed (flat {flat_frac:.0%}, shaped {1-flat_frac:.0%})'

            rows.append({
                'isotope': iso,
                'region': region,
                'total_gammas_per_yr': f"{total_entering:.6e}",
                'flat_gammas_per_yr': f"{flat_entering:.6e}",
                'shaped_gammas_per_yr': f"{shaped_entering:.6e}",
                'angular_distribution': ang_type,
            })

    df = pd.DataFrame(rows)
    with open(filepath, 'w') as f:
        f.write("# LZ PMTs: MC sampling parameters\n")
        f.write("# Flat sources: u = cos(theta) ~ Uniform(0, 1), factor 1/2 already applied\n")
        f.write("# Shaped sources: R8778 PMTs behind ~6 cm LXe skin\n")
        f.write("# Entry point: uniform on barrel cylinder or endcap disc\n")
        f.write("#   endcap_top: disc at top of TPC\n")
        f.write("#   endcap_bottom: disc at bottom of TPC\n")
        f.write("#   barrel: cylindrical shell\n")
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
    print_section("LZ PMT MODEL — Third Gamma Generator")

    # --- Load attenuation coefficients (for R8778 transport) ---
    mu = load_attenuation_coefficients()
    print(f"\n  μ_LXe at 2.448 MeV: {mu['LXe_Bi214']:.4f} cm⁻¹")
    print(f"  μ_LXe at 2.615 MeV: {mu['LXe_Tl208']:.4f} cm⁻¹")
    print(f"  LXe skin thickness:  {SKIN_THICKNESS_CM} cm")
    print(f"  μt (Bi214): {mu['LXe_Bi214']*SKIN_THICKNESS_CM:.3f}")
    print(f"  μt (Tl208): {mu['LXe_Tl208']*SKIN_THICKNESS_CM:.3f}")

    # --- Build components ---
    components = build_components()

    # --- Activity table ---
    print_section("RADIOACTIVITY")
    print(f"\n  {'Component':<30s} {'Region':<14s} {'Type':<8s} {'N':>4s} "
          f"{'Mass':>7s} {'Bi214':>10s} {'Tl208':>10s}")
    print(f"  {'':30s} {'':14s} {'':8s} {'':>4s} "
          f"{'(kg)':>7s} {'(mBq)':>10s} {'(mBq)':>10s}")
    print(f"  {'-'*90}")

    tot_bi = 0; tot_tl = 0
    for c in components:
        bi = compute_total_activity(c, 'Bi214')
        tl = compute_total_activity(c, 'Tl208')
        tot_bi += bi; tot_tl += tl
        print(f"  {c['name']:<30s} {c['region']:<14s} {c['activity_type']:<8s} "
              f"{c['count']:>4d} {c['mass_kg']:>7.2f} {bi:>10.1f} {tl:>10.1f}")
    print(f"  {'-'*90}")
    print(f"  {'TOTAL':<30s} {'':14s} {'':8s} {'':>4s} {'':>7s} "
          f"{tot_bi:>10.1f} {tot_tl:>10.1f}")

    # --- Gamma rates ---
    print_section("GAMMA RATES")
    results = compute_gamma_rates(components, mu)

    print(f"\n  {'Source':<30s} {'Region':<14s} {'Iso':<7s} {'Type':<7s} "
          f"{'Produced':>12s} {'Entering':>12s} {'<cosθ>':>8s}")
    print(f"  {'-'*95}")
    for r in results:
        print(f"  {r['source']:<30s} {r['region']:<14s} {r['isotope']:<7s} "
              f"{r['angular_type']:<7s} {r['gammas_produced_per_yr']:>12.2e} "
              f"{r['gammas_entering_per_yr']:>12.2e} {r['mean_u']:>8.3f}")

    # --- Aggregated MC input ---
    print_section("MC INPUT SUMMARY (gammas/yr entering active volume)")
    for iso in ['Bi214', 'Tl208']:
        for region in ['endcap_top', 'endcap_bottom', 'barrel']:
            matching = [r for r in results
                        if r['isotope'] == iso and r['region'] == region]
            if not matching:
                continue
            total = sum(r['gammas_entering_per_yr'] for r in matching)
            flat = sum(r['gammas_entering_per_yr'] for r in matching
                       if r['angular_type'] == 'flat')
            shaped = total - flat
            note = ""
            if shaped > 0:
                note = f"  (flat: {flat:.2e}, shaped: {shaped:.2e})"
            print(f"  {iso} {region:<14s}: {total:>12.2e} γ/yr{note}")

    # --- Comparison with other generators ---
    print_section("COMPARISON: ALL THREE GENERATORS (gammas/yr entering active volume)")
    pmt_bi = sum(r['gammas_entering_per_yr'] for r in results if r['isotope'] == 'Bi214')
    pmt_tl = sum(r['gammas_entering_per_yr'] for r in results if r['isotope'] == 'Tl208')
    print(f"\n  {'Generator':<20s} {'²¹⁴Bi (γ/yr)':>14s} {'²⁰⁸Tl (γ/yr)':>14s}")
    print(f"  {'-'*50}")
    print(f"  {'Cryostat':20s} {'5.9e+04':>14s} {'7.1e+05':>14s}  (from lz_cryo)")
    print(f"  {'Field cage':20s} {'7.1e+03':>14s} {'3.3e+05':>14s}  (from lz_fieldcage)")
    print(f"  {'PMTs':20s} {pmt_bi:>14.2e} {pmt_tl:>14.2e}  (this model)")

    # --- Save CSV files ---
    print_section("SAVING OUTPUT FILES")

    f1 = os.path.join(DATA_DIR, "lz_pmt_activity.csv")
    save_activity_csv(components, f1)
    print(f"  {f1}")

    f2 = os.path.join(DATA_DIR, "lz_pmt_gammas.csv")
    save_gammas_csv(results, f2)
    print(f"  {f2}")

    f3 = os.path.join(DATA_DIR, "lz_pmt_sampling.csv")
    save_sampling_csv(results, f3)
    print(f"  {f3}")

    print("\nDone.")


if __name__ == "__main__":
    main()
