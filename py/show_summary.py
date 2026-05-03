"""Render a run-level report from the per-source `summary.csv` files
written by `scripts/run_mc3.jl`.

Reads `output/<run-name>/<source>/summary.csv` for every source dir
under the run dir and emits three blocks (mirroring what the Julia
driver prints at the end of a single-source run, but aggregated):

    1. Per-source results table  (source / iso / γ/yr / f_SS_ROI /
                                   bg γ/yr / runtime), with TOTAL.
    2. Analysis cuts             (read once from any row; identical
                                   across rows of the same run).
    3. Per-source funnel          (N0..N6 with cum= and acc= columns).

Output: stdout (default), `<run-dir>/run_report.txt`, and an aggregated
`<run-dir>/run_report.csv` (one row per source, same columns as the
per-source summary.csv).

Usage
-----
    python py/show_summary.py output/<run-name>/
    python py/show_summary.py output/<run-name>/ --source CB_Tl208
    python py/show_summary.py output/<run-name>/ --no-files   # stdout only
"""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

import pandas as pd


# Source ordering for the report. Sources not in this list are appended
# alphabetically at the end so future additions still print.
_PREFERRED_ORDER = [
    # Cryostat (6)
    "CB_Bi214", "CTH_Bi214", "CBH_Bi214",
    "CB_Tl208", "CTH_Tl208", "CBH_Tl208",
    # Field cage Bi-214 (6): rings, resistors, sensors, PTFE, top grids,
    # bottom grid (cathode-only)
    "FCRN_Bi214", "FCRS_Bi214", "FCSE_Bi214",
    "FCPT_Bi214", "FCTG_Bi214", "FCBG_Bi214",
    # Field cage Tl-208 (6)
    "FCRN_Tl208", "FCRS_Tl208", "FCSE_Tl208",
    "FCPT_Tl208", "FCTG_Tl208", "FCBG_Tl208",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "run_dir",
        type=Path,
        help="Run directory: output/<run-name>/. The script globs "
              "**/summary.csv under it.",
    )
    p.add_argument(
        "--source",
        type=str,
        default=None,
        help="Optional source filter: only show rows whose `source` "
             "column matches this value.",
    )
    p.add_argument(
        "--no-files",
        action="store_true",
        help="Skip writing run_report.txt and run_report.csv; print "
             "to stdout only.",
    )
    return p.parse_args()


def _load_summaries(run_dir: Path, source_filter: str | None) -> pd.DataFrame:
    csvs = sorted(run_dir.glob("*/summary.csv"))
    if not csvs:
        sys.exit(f"no */summary.csv found under {run_dir}")
    frames = []
    for csv in csvs:
        try:
            frames.append(pd.read_csv(csv))
        except Exception as e:
            print(f"  ! skipping {csv} ({e})", file=sys.stderr)
    if not frames:
        sys.exit(f"no readable summary.csv under {run_dir}")
    df = pd.concat(frames, ignore_index=True)
    if source_filter is not None:
        df = df[df["source"] == source_filter]
        if df.empty:
            sys.exit(f"no source matched --source {source_filter!r}")
    # Stable ordering: preferred list first, then anything new alphabetically.
    df["__rank"] = df["source"].apply(
        lambda s: _PREFERRED_ORDER.index(s) if s in _PREFERRED_ORDER
                   else 1000 + ord(s[0])
    )
    df = df.sort_values("__rank").drop(columns="__rank").reset_index(drop=True)
    return df


def _render_results_table(df: pd.DataFrame, out: io.StringIO) -> None:
    out.write("── Per-source results ──\n")
    out.write(f"  {'source':<12} {'iso':<7} "
              f"{'γ/yr_in':>12} {'f_SS_ROI':>10} "
              f"{'bg γ/yr':>12} {'runtime (s)':>12}\n")
    out.write("  " + "─" * 77 + "\n")
    total_bg = 0.0
    for _, r in df.iterrows():
        out.write(f"  {r['source']:<12} {r['isotope']:<7} "
                  f"{r['gamma_per_yr_total']:>12.3e} "
                  f"{r['f_SS_in_ROI']:>10.3e} "
                  f"{r['bg_per_yr']:>12.3e} "
                  f"{r['runtime_s']:>12.2f}\n")
        total_bg += float(r["bg_per_yr"])
    out.write("  " + "─" * 77 + "\n")
    out.write(f"  TOTAL background : {total_bg:.4e} events/yr\n\n")


def _render_cuts(df: pd.DataFrame, out: io.StringIO) -> None:
    """Print analysis-cuts block. Reads once from the first row; assumes
    cuts are identical across rows of the same run (they are, by
    construction in run_mc3.jl)."""
    r = df.iloc[0]
    Q   = float(r["Q_betabeta_keV"])
    HW  = float(r["ROI_halfwidth_keV"])
    σ   = float(r["sigma_E_over_E"])
    z_lo, z_hi = float(r["fv_z_min_cm"]), float(r["fv_z_max_cm"])
    r_max = float(r["fv_r_max_cm"])
    Evis  = float(r["E_visible_keV"])
    Eskv  = float(r["E_skin_veto_keV"])
    Rlxe  = float(r["R_LXe_outer_cm"])
    out.write("── Analysis cuts ──\n")
    out.write(f"  Q_ββ                : {Q:.1f} keV\n")
    out.write(f"  σ_E / E             : {σ:.4f}\n")
    out.write(f"  ROI window          : Q_ββ ± {HW:.2f} keV  →  "
              f"[{Q-HW:.2f}, {Q+HW:.2f}] keV\n")
    out.write(f"  FV box              : z ∈ [{z_lo:.1f}, {z_hi:.1f}] cm,  "
              f"r ≤ {r_max:.1f} cm\n")
    out.write(f"  Visible threshold   : {Evis:.1f} keV  (per-cluster)\n")
    out.write(f"  Skin-veto threshold : {Eskv:.1f} keV  (cumulative on :Skin)\n")
    out.write(f"  R_LXe outer         : {Rlxe:.1f} cm\n\n")


def _render_funnel(df: pd.DataFrame, out: io.StringIO) -> None:
    out.write("── Per-source funnel ──\n")
    rows = (
        ("Into detector",          "n_escaped"),
        ("Pass skin veto",         "n_skin_vetoed"),
        ("Inside FV",              "n_outside_FV"),
        ("Accepted as SS",         "n_MS"),
        ("In ROI (pre-companion)", "n_SS_outside_ROI"),
        ("After companion veto",   "n_companion_vetoed"),
    )
    for _, r in df.iterrows():
        n0 = int(r["n_total"])
        out.write(f"\n  ── {r['source']} ──\n")
        out.write(f"    {'Total events run':<24} N  = {n0:>12d}  "
                  f"cum=1.000e+00\n")
        prev = n0
        cum = n0
        for i, (label, drop_col) in enumerate(rows, start=1):
            cum -= int(r[drop_col])
            cum_frac = cum / n0 if n0 > 0 else 0.0
            acc = (cum / prev) if prev > 0 else None
            acc_str = f"{acc:.3e}" if acc is not None else "  ---  "
            out.write(f"    {label:<24} N{i:<1d} = {cum:>12d}  "
                      f"cum={cum_frac:.3e}  acc={acc_str}\n")
            prev = cum
    out.write("\n")


def _render_full_report(df: pd.DataFrame) -> str:
    buf = io.StringIO()
    _render_results_table(df, buf)
    _render_cuts(df, buf)
    _render_funnel(df, buf)
    return buf.getvalue()


def main() -> None:
    args = parse_args()
    if not args.run_dir.is_dir():
        sys.exit(f"not a directory: {args.run_dir}")
    df = _load_summaries(args.run_dir, args.source)

    text = _render_full_report(df)
    print(text)

    if not args.no_files:
        txt_path = args.run_dir / "run_report.txt"
        csv_path = args.run_dir / "run_report.csv"
        txt_path.write_text(text)
        df.to_csv(csv_path, index=False)
        print(f"  → wrote {txt_path}")
        print(f"  → wrote {csv_path}")


if __name__ == "__main__":
    main()
