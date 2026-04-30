"""
LZ Parameter Consistency Check
===============================

Compares our model parameters against the LZ TDR Table 3.1.1 reference
values and the bb0nu paper Table I. Reads geometry from Julia source and
CSV output files, computes derived quantities, and flags discrepancies.

Usage:
  python py/lz_check_params.py

Reads from:
  src/geometry.jl       — Params defaults (parsed)
  data/lz_cryo_*.csv    — cryostat model outputs
  data/lz_fc_*.csv      — field cage model outputs
  data/lz_pmt_*.csv     — PMT model outputs

Outputs:
  Console comparison table
  data/lz_param_check.csv
"""

import os
import re
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")
DATA_DIR = os.path.join(PROJECT_DIR, "data")
SRC_DIR = os.path.join(PROJECT_DIR, "src")

RHO_LXE = 2.953  # g/cm³


# ===========================================================================
# Parse Julia Params defaults from geometry.jl
# ===========================================================================
def parse_params_defaults():
    """
    Extract default parameter values from the Params struct in geometry.jl.
    Returns dict of {field_name: value}.
    """
    path = os.path.join(SRC_DIR, "geometry.jl")
    params = {}
    with open(path) as f:
        in_struct = False
        for line in f:
            stripped = line.strip()
            if 'struct Params' in stripped:
                in_struct = True
                continue
            if in_struct and stripped == 'end':
                break
            if not in_struct:
                continue
            # Match lines like: field_name::Type = value  # comment
            m = re.match(r'(\w+)::\w+\s*=\s*([^#]+)', stripped)
            if m:
                name = m.group(1)
                val_str = m.group(2).strip().rstrip(',')
                # Skip string values
                if val_str.startswith('"'):
                    params[name] = val_str.strip('"')
                else:
                    try:
                        # Handle underscores in Julia numbers
                        val_str = val_str.replace('_', '')
                        params[name] = float(val_str)
                    except ValueError:
                        params[name] = val_str
    return params


# ===========================================================================
# Read CSV helpers
# ===========================================================================
def read_csv_column_sum(filepath, column):
    """Sum a numeric column from a CSV file (skip comments)."""
    total = 0.0
    header = None
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            if header is None:
                header = line.strip().split(',')
                col_idx = header.index(column)
                continue
            parts = line.strip().split(',')
            try:
                total += float(parts[col_idx])
            except (ValueError, IndexError):
                pass
    return total


def read_csv_rows(filepath):
    """Read CSV into list of dicts (skip comments)."""
    rows = []
    header = None
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            if header is None:
                header = line.strip().split(',')
                continue
            parts = line.strip().split(',')
            if len(parts) == len(header):
                rows.append(dict(zip(header, parts)))
    return rows


# ===========================================================================
# Reference values from TDR Table 3.1.1 and bb0nu paper
# ===========================================================================
TDR_REFERENCE = {
    'TPC active mass (kg)': 7000.0,
    'Skin mass side+dome (kg)': 2000.0,
    'Total LXe in cryostat (kg)': 9600.0,
    'TPC PMTs (top+bottom)': 494,
    'Top Skin PMTs (R8520)': 93,
    'Bottom Skin PMTs (R8778)': 38,
    'EL region gate-anode (mm)': 13.0,
    'Drift region cathode-gate (mm)': 1456.0,
    'RFR sub-cathode (mm)': 137.5,
    'TPC inner diameter (mm)': 1456.0,
    'Field cage thickness (mm)': 15.0,
    'Skin thickness at surface (mm)': 40.0,
    'Skin thickness at cathode (mm)': 80.0,
    'Drift field stages': 57,
    'RFR stages': 7,
    'Operating pressure (bar)': 1.8,
    'Temperature (K)': 175.8,
}

BBONU_REFERENCE = {
    'Ti cryostat vessel mass (kg)': 2590.0,
    'Cryostat insulation mass (kg)': 13.8,
    'Field-cage rings mass (kg)': 93.0,
    'Field-cage resistors mass (kg)': 0.06,
    'PTFE walls mass (kg)': 184.0,
    'TPC sensors mass (kg)': 5.02,
    'Field grids+holders mass (kg)': 89.1,
    'TPC PMTs mass (kg)': 91.9,
    'TPC PMT bases mass (kg)': 2.80,
    'TPC PMT structures mass (kg)': 166.0,
    'TPC PMT cables mass (kg)': 88.7,
    'Skin PMTs+bases mass (kg)': 8.59,
    'FV z_min (cm)': 26.0,
    'FV z_max (cm)': 96.0,
    'FV r_max (cm)': 39.0,
    'FV mass (kg)': 967.0,
    'Energy resolution sigma/E': 0.01,
    'SS/MS z threshold (mm)': 3.0,
}


# ===========================================================================
# Compute model values
# ===========================================================================
def compute_model_values(params):
    """Compute all model quantities from parsed Params."""
    model = {}

    R_tpc = params['geom_R_lxe']
    L_tpc = params['geom_L_lxe']
    R_skin_in = params['geom_R_skin_inner']
    R_skin_out = params['geom_R_skin_outer']
    rho = params['geom_ρ_LXe'] if 'geom_ρ_LXe' in params else RHO_LXE

    # TPC
    tpc_vol = np.pi * R_tpc**2 * L_tpc
    tpc_mass = tpc_vol * rho / 1e6  # tonnes → kg
    model['TPC active mass (kg)'] = tpc_vol * rho / 1e3

    # Skin (cylindrical shell)
    skin_vol = np.pi * (R_skin_out**2 - R_skin_in**2) * L_tpc
    model['Skin mass side+dome (kg)'] = skin_vol * rho / 1e3

    # Total (TPC + skin only, no dome/RFR)
    full_vol = np.pi * R_skin_out**2 * L_tpc
    model['Total LXe in cryostat (kg)'] = full_vol * rho / 1e3

    # Dimensions
    model['TPC inner diameter (mm)'] = R_tpc * 20.0
    model['Drift region cathode-gate (mm)'] = L_tpc * 10.0
    model['Field cage thickness (mm)'] = (R_skin_in - R_tpc) * 10.0
    model['Skin thickness avg (mm)'] = (R_skin_out - R_skin_in) * 10.0

    # FV
    fv_r = np.sqrt(params['fv_r2_max_cm2'])
    fv_z_min = params['fv_z_min_cm']
    fv_z_max = params['fv_z_max_cm']
    fv_vol = np.pi * fv_r**2 * (fv_z_max - fv_z_min)
    model['FV z_min (cm)'] = fv_z_min
    model['FV z_max (cm)'] = fv_z_max
    model['FV r_max (cm)'] = fv_r
    model['FV mass (kg)'] = fv_vol * rho / 1e3

    # Cuts
    model['SS/MS z threshold (mm)'] = params['cut_Δz_threshold_mm']
    model['Energy resolution sigma/E'] = params['phys_σ_E_over_E']
    model['ROI halfwidth (keV)'] = params['cut_ROI_halfwidth_keV']
    model['E visible threshold (keV)'] = params['cut_E_visible_threshold_keV']
    model['Companion veto (keV)'] = params['cut_companion_veto_keV']
    model['Skin veto (keV)'] = params['cut_skin_veto_keV']

    return model


def get_component_masses():
    """Read component masses from CSV files."""
    masses = {}

    # Cryostat
    try:
        rows = read_csv_rows(os.path.join(DATA_DIR, "lz_cryo_activity.csv"))
        ti_mass = sum(float(r['mass_kg']) for r in rows
                      if 'shell' in r['source'].lower() or 'heads' in r['source'].lower())
        mli_mass = sum(float(r['mass_kg']) for r in rows
                       if 'insulation' in r['source'].lower())
        masses['Ti cryostat vessel mass (kg)'] = ti_mass
        masses['Cryostat insulation mass (kg)'] = mli_mass
    except Exception:
        pass

    # Field cage
    try:
        rows = read_csv_rows(os.path.join(DATA_DIR, "lz_fc_activity.csv"))
        for r in rows:
            name = r['component'].lower()
            m = float(r['mass_kg'])
            if 'ring' in name and 'grid' not in name:
                masses['Field-cage rings mass (kg)'] = m
            elif 'resistor' in name:
                masses['Field-cage resistors mass (kg)'] = m
            elif 'ptfe' in name:
                masses['PTFE walls mass (kg)'] = m
            elif 'sensor' in name:
                masses['TPC sensors mass (kg)'] = m
        # Grids total
        grid_mass = sum(float(r['mass_kg']) for r in rows if 'grid' in r['component'].lower())
        masses['Field grids+holders mass (kg)'] = grid_mass
    except Exception:
        pass

    # PMTs
    try:
        rows = read_csv_rows(os.path.join(DATA_DIR, "lz_pmt_activity.csv"))
        pmt_mass = sum(float(r['mass_kg']) for r in rows if 'TPC PMTs' in r['component'])
        base_mass = sum(float(r['mass_kg']) for r in rows if 'bases' in r['component'].lower())
        struct_mass = sum(float(r['mass_kg']) for r in rows if 'structure' in r['component'].lower())
        cable_mass = sum(float(r['mass_kg']) for r in rows if 'cable' in r['component'].lower())
        skin_mass = sum(float(r['mass_kg']) for r in rows if 'Skin' in r['component'])

        n_tpc_pmts = sum(int(r['count']) for r in rows if 'TPC PMTs' in r['component'])
        n_skin_side = sum(int(r['count']) for r in rows if 'side' in r['component'].lower())
        n_skin_lower = sum(int(r['count']) for r in rows if 'lower' in r['component'].lower())
        n_skin_dome = sum(int(r['count']) for r in rows if 'dome' in r['component'].lower())

        masses['TPC PMTs mass (kg)'] = pmt_mass
        masses['TPC PMT bases mass (kg)'] = base_mass
        masses['TPC PMT structures mass (kg)'] = struct_mass
        masses['TPC PMT cables mass (kg)'] = cable_mass
        masses['Skin PMTs+bases mass (kg)'] = skin_mass
        masses['TPC PMTs (top+bottom)'] = n_tpc_pmts
        masses['Top Skin PMTs (R8520)'] = n_skin_side
        masses['Bottom Skin PMTs (R8778)'] = n_skin_lower + n_skin_dome
    except Exception:
        pass

    return masses


# ===========================================================================
# Comparison
# ===========================================================================
def compare(label, reference, model, tolerance=0.05):
    """
    Compare reference values against model values.
    Returns list of (parameter, reference, model, status) tuples.
    """
    results = []
    for param, ref_val in reference.items():
        mod_val = model.get(param, None)
        if mod_val is None:
            status = "NOT IN MODEL"
        elif isinstance(ref_val, (int, float)) and isinstance(mod_val, (int, float)):
            if ref_val == 0:
                status = "OK" if mod_val == 0 else "MISMATCH"
            else:
                diff = abs(mod_val - ref_val) / abs(ref_val)
                if diff < 0.001:
                    status = "OK"
                elif diff < tolerance:
                    status = f"~OK ({diff:.1%})"
                else:
                    status = f"** DIFF {diff:.1%} **"
        else:
            status = "OK" if mod_val == ref_val else "MISMATCH"
        results.append((param, ref_val, mod_val, status))
    return results


# ===========================================================================
# Main
# ===========================================================================
def main():
    print("=" * 78)
    print("LZ PARAMETER CONSISTENCY CHECK")
    print("=" * 78)

    # Parse Julia defaults
    params = parse_params_defaults()
    model = compute_model_values(params)

    # Add component masses from CSVs
    comp_masses = get_component_masses()
    model.update(comp_masses)

    # --- TDR Table 3.1.1 comparison ---
    print("\n--- vs TDR Table 3.1.1 ---")
    print(f"{'Parameter':<40s} {'TDR':>12s} {'Model':>12s} {'Status':>15s}")
    print("-" * 82)

    tdr_results = compare("TDR", TDR_REFERENCE, model)
    for param, ref, mod, status in tdr_results:
        ref_str = f"{ref}" if isinstance(ref, int) else f"{ref:.1f}"
        mod_str = "—" if mod is None else (f"{mod}" if isinstance(mod, int) else f"{mod:.1f}")
        print(f"{param:<40s} {ref_str:>12s} {mod_str:>12s} {status:>15s}")

    # --- bb0nu paper comparison ---
    print(f"\n--- vs bb0nu Paper Table I ---")
    print(f"{'Parameter':<40s} {'Paper':>12s} {'Model':>12s} {'Status':>15s}")
    print("-" * 82)

    bbonu_results = compare("bbonu", BBONU_REFERENCE, model)
    for param, ref, mod, status in bbonu_results:
        ref_str = f"{ref}" if isinstance(ref, int) else f"{ref:.2f}"
        mod_str = "—" if mod is None else (f"{mod}" if isinstance(mod, int) else f"{mod:.2f}")
        print(f"{param:<40s} {ref_str:>12s} {mod_str:>12s} {status:>15s}")

    # --- Summary ---
    all_results = tdr_results + bbonu_results
    n_ok = sum(1 for _, _, _, s in all_results if s.startswith("OK") or s.startswith("~OK"))
    n_diff = sum(1 for _, _, _, s in all_results if s.startswith("**"))
    n_missing = sum(1 for _, _, _, s in all_results if s == "NOT IN MODEL")

    print(f"\n--- Summary ---")
    print(f"  OK / ~OK:      {n_ok}")
    print(f"  Discrepancy:   {n_diff}")
    print(f"  Not in model:  {n_missing}")
    print(f"  Total checked: {len(all_results)}")

    # --- Save CSV ---
    csv_path = os.path.join(DATA_DIR, "lz_param_check.csv")
    with open(csv_path, 'w') as f:
        f.write("source,parameter,reference,model,status\n")
        for param, ref, mod, status in tdr_results:
            mod_str = "" if mod is None else f"{mod}"
            f.write(f"TDR,{param},{ref},{mod_str},{status}\n")
        for param, ref, mod, status in bbonu_results:
            mod_str = "" if mod is None else f"{mod}"
            f.write(f"bbonu,{param},{ref},{mod_str},{status}\n")
    print(f"\n  Saved: {csv_path}")


if __name__ == "__main__":
    main()
