"""
Tests for the Python gamma generator scripts.

Verifies that the three generators (cryostat, field cage, PMTs) produce
consistent output files with correct activity budgets and valid PDFs.

Run with: python -m pytest py/test_generators.py -v
Or simply: python py/test_generators.py
"""

import os
import sys
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data")
sys.path.insert(0, SCRIPT_DIR)


def test_cryo_total_ti_activity():
    """Ti activity = 2590 kg × specific activity, bbonu paper."""
    import lz_cryo as lc
    total_mass = lc.CRYO_TOTAL_TI_MASS
    assert total_mass == 2590.0

    expected_bi = total_mass * lc.TI_BI214_SP
    expected_tl = total_mass * lc.TI_TL208_SP

    assert abs(expected_bi - 2590 * 0.08) < 0.01
    assert abs(expected_tl - 2590 * 0.22) < 0.01


def test_cryo_mass_scale():
    """Scaled masses sum to 2590 kg for Ti sectors."""
    import lz_cryo as lc
    mu = lc.load_attenuation_coefficients()
    ocv = lc.compute_vessel(**lc.OCV_CFG)
    icv = lc.compute_vessel(**lc.ICV_CFG)
    sectors = lc.build_source_sectors(ocv, icv, mu)

    # Sum Ti sector masses (exclude MLI)
    ti_mass = sum(s['mass_kg'] for s in sectors if 'shell' in s['name'] or 'heads' in s['name'])
    assert abs(ti_mass - 2590.0) < 1.0, f"Ti mass = {ti_mass}, expected 2590"


def test_cryo_geometry_unchanged():
    """ICV and OCV wall thicknesses solved from TDR masses, not bbonu."""
    import lz_cryo as lc
    ocv = lc.compute_vessel(**lc.OCV_CFG)
    icv = lc.compute_vessel(**lc.ICV_CFG)

    assert abs(ocv['m_total'] - 1115.0) < 0.1
    assert abs(icv['m_total'] - 950.0) < 0.1


def test_cryo_no_seals_no_liner():
    """Seals and Teflon liner are no longer modeled."""
    import lz_cryo as lc
    assert 'Seals' not in lc.ANCILLARY
    assert 'Teflon liner' not in lc.ANCILLARY


def test_cryo_mli_bbonu_values():
    """MLI uses bbonu paper values."""
    import lz_cryo as lc
    mli = lc.ANCILLARY['MLI']
    assert mli['mass'] == 13.8
    assert mli['Bi214_sp'] == 11.1
    assert mli['Tl208_sp'] == 7.79


def test_cryo_linear_fit_csv():
    """Linear fit CSV exists and has valid content."""
    path = os.path.join(DATA_DIR, "lz_cryo_linear_fit.csv")
    assert os.path.isfile(path), f"Missing: {path}"

    with open(path) as f:
        lines = [l for l in f if not l.startswith('#')]

    header = lines[0].strip().split(',')
    assert 'isotope' in header
    assert 'a' in header
    assert 'b' in header

    # Should have 4 data rows
    data_lines = [l for l in lines[1:] if l.strip()]
    assert len(data_lines) == 4


def test_cryo_angular_pdfs():
    """Angular PDF files exist and integrate to ~1."""
    for iso in ['Bi214', 'Tl208']:
        for region in ['barrel', 'endcap']:
            path = os.path.join(DATA_DIR, f"lz_cryo_angular_{iso}_{region}.csv")
            assert os.path.isfile(path), f"Missing: {path}"

            u, pdf = [], []
            with open(path) as f:
                for line in f:
                    if line.startswith('#') or line.startswith('u'):
                        continue
                    parts = line.strip().split(',')
                    u.append(float(parts[0]))
                    pdf.append(float(parts[1]))

            u = np.array(u)
            pdf = np.array(pdf)
            integral = np.trapz(pdf, u)
            assert abs(integral - 1.0) < 0.05, \
                f"{iso}_{region}: PDF integral = {integral}, expected ~1"


def test_fc_bbonu_values():
    """Field cage uses bbonu paper values."""
    import lz_fieldcage as lf
    comps = lf.build_components()

    rings = [c for c in comps if 'rings' in c['name'].lower()]
    assert len(rings) == 1
    assert rings[0]['mass_kg'] == 93.0
    assert rings[0]['Bi214_sp'] == 0.35

    resistors = [c for c in comps if 'resistor' in c['name'].lower()]
    assert len(resistors) == 1
    assert resistors[0]['mass_kg'] == 0.06
    assert resistors[0]['Bi214_sp'] == 1350.0


def test_fc_grids_split():
    """Grids+holders split 50/50 top/bottom, total 89.1 kg."""
    import lz_fieldcage as lf
    comps = lf.build_components()

    grids = [c for c in comps if 'grids' in c['name'].lower()]
    assert len(grids) == 2
    total_grid_mass = sum(c['mass_kg'] for c in grids)
    assert abs(total_grid_mass - 89.1) < 0.1


def test_fc_sampling_csv():
    """FC sampling CSV exists with 4 rows."""
    path = os.path.join(DATA_DIR, "lz_fc_sampling.csv")
    assert os.path.isfile(path)

    with open(path) as f:
        lines = [l for l in f if not l.startswith('#')]
    data = [l for l in lines[1:] if l.strip()]
    assert len(data) == 4


def test_pmt_bbonu_values():
    """PMTs use bbonu paper per-kg activities."""
    import lz_pmts as lp
    comps = lp.build_components()

    tpc_pmts = [c for c in comps if 'TPC PMTs' in c['name']]
    assert len(tpc_pmts) == 2  # top + bottom
    for c in tpc_pmts:
        assert c['activity_type'] == 'per_kg'
        assert c['Bi214_val'] == 3.22
        assert c['Tl208_val'] == 1.61


def test_pmt_total_masses():
    """PMT component masses match bbonu paper."""
    import lz_pmts as lp
    comps = lp.build_components()

    # Total TPC PMT mass should be 91.9 kg
    pmt_mass = sum(c['mass_kg'] for c in comps if 'TPC PMTs' in c['name'])
    assert abs(pmt_mass - 91.9) < 0.1

    # Total structure mass should be 166 kg
    struct_mass = sum(c['mass_kg'] for c in comps if 'structure' in c['name'])
    assert abs(struct_mass - 166.0) < 0.1

    # Cable mass
    cable = [c for c in comps if 'cable' in c['name'].lower()]
    assert len(cable) == 1
    assert cable[0]['mass_kg'] == 88.7

    # Skin total
    skin = [c for c in comps if 'Skin' in c['name']]
    skin_mass = sum(c['mass_kg'] for c in skin)
    assert abs(skin_mass - 8.59) < 0.1


def test_pmt_sampling_csv():
    """PMT sampling CSV exists."""
    path = os.path.join(DATA_DIR, "lz_pmt_sampling.csv")
    assert os.path.isfile(path)


def test_all_csvs_parseable():
    """All output CSVs are readable with consistent format."""
    csv_files = [
        "lz_cryo_geometry.csv", "lz_cryo_activity.csv", "lz_cryo_gammas.csv",
        "lz_cryo_linear_fit.csv",
        "lz_cryo_angular_Bi214_barrel.csv", "lz_cryo_angular_Bi214_endcap.csv",
        "lz_cryo_angular_Tl208_barrel.csv", "lz_cryo_angular_Tl208_endcap.csv",
        "lz_fc_geometry.csv", "lz_fc_activity.csv", "lz_fc_gammas.csv",
        "lz_fc_sampling.csv",
        "lz_pmt_activity.csv", "lz_pmt_gammas.csv", "lz_pmt_sampling.csv",
    ]

    for fname in csv_files:
        path = os.path.join(DATA_DIR, fname)
        assert os.path.isfile(path), f"Missing: {fname}"
        # Try reading: skip comments, check header exists
        with open(path) as f:
            lines = [l for l in f if not l.startswith('#')]
        assert len(lines) >= 2, f"{fname}: too few lines (need header + data)"
        header = lines[0].strip()
        assert ',' in header, f"{fname}: no commas in header"


def test_branching_ratios():
    """Verify BR constants match bbonu paper."""
    import lz_cryo as lc
    assert lc.BR_BI214_GAMMA == 0.0155
    assert lc.BR_TL208_FROM_CHAIN == 0.359
    assert lc.BR_TL208_GAMMA == 1.0


if __name__ == "__main__":
    import traceback
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    passed = 0
    failed = 0
    for test in tests:
        try:
            test()
            print(f"  PASS  {test.__name__}")
            passed += 1
        except Exception as e:
            print(f"  FAIL  {test.__name__}: {e}")
            traceback.print_exc()
            failed += 1

    print(f"\n{passed} passed, {failed} failed out of {passed + failed} tests")
    if failed > 0:
        sys.exit(1)
