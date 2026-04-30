"""
Run all Bi-214 sources and collect results.

Runs each of the 11 Bi-214 source indices sequentially via Julia,
parses the summary.txt output, and produces a comparison table
against the LZ bb0nu paper.

Usage:
  python py/run_bi214.py
  python py/run_bi214.py --n-samples 10000000   # quick test with 10M
"""

import os
import sys
import subprocess
import re

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")
DATA_DIR = os.path.join(PROJECT_DIR, "data")
OUTPUT_BASE = os.path.join(PROJECT_DIR, "output", "bi214_run")  # Julia writes to output/<--output>/LABEL/

# Bi-214 sources: (index, label, paper_counts_per_1000d)
# Paper counts from bb0nu Table I, ²³⁸U column
BI214_SOURCES = [
    (1,  "cryo_barrel_Bi214",     None),    # Ti vessel barrel (part of 1.30)
    (2,  "cryo_head_top_Bi214",   None),    # Ti vessel top head (part of 1.30)
    (3,  "cryo_head_bot_Bi214",   None),    # Ti vessel bot head + extra (part of 1.30)
    (7,  "fc_rings_Bi214",        0.82),    # Field-cage rings
    (9,  "fc_ptfe_Bi214",         0.39),    # PTFE walls
    (11, "fc_ressens_Bi214",      3.82),    # Resistors + sensors (2.63 + 1.19)
    (13, "fc_holder_top_Bi214",   None),    # Grid holders top (part of 0.62)
    (15, "fc_holder_bot_Bi214",   None),    # Grid holders bot (part of 0.62)
    (17, "pmt_top_Bi214",         None),    # PMT top (part of 9.31)
    (19, "pmt_bottom_Bi214",      None),    # PMT bottom (part of 9.31)
    # (21, "pmt_barrel_Bi214")    # EXCLUDED — needs proper cable/skin modeling
]

# Paper grouped totals for comparison
PAPER_GROUPS = {
    "Cryostat (Ti+MLI)": 2.20,      # counts/1000d
    "FC rings":          0.82,
    "PTFE walls":        0.39,
    "FC res+sens":       3.82,
    "FC holders":        0.62,
    "PMTs (top+bot)":    7.12,       # PMTs+bases+structures only (excl cables+skin)
}
# Not modeled: PMT barrel (cables+skin) = 2.19 cts/1000d = 0.80 evt/yr

# Default FV and cuts
FV_ARGS = ["--fv-z-min", "26.0", "--fv-z-max", "96.0", "--fv-r-max", "39.0"]
CUT_ARGS = ["--e-visible", "10.0"]


def run_source(index, label, n_samples=100_000_000):
    """Run one source and return the parsed results."""
    outdir = os.path.join(OUTPUT_BASE, label)
    cmd = [
        "julia", "--project=.", "-t", "8",
        "scripts/run_mc.jl",
        "--source-index", str(index),
        "--n-samples", str(n_samples),
        "--output", os.path.join("bi214_run", label),
    ] + FV_ARGS + CUT_ARGS

    print(f"  Running source {index}: {label} ...", end=" ", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=PROJECT_DIR)

    if result.returncode != 0:
        print(f"FAILED (exit {result.returncode})")
        print("STDERR:", result.stderr[-500:] if result.stderr else "no stderr")
        print("STDOUT:", result.stdout[-500:] if result.stdout else "no stdout")
        return None

    # Find summary.txt (Julia creates OUTDIR/LABEL/summary.txt)
    summary_path = None
    for root, dirs, files in os.walk(OUTPUT_BASE):
        if "summary.txt" in files and label in root:
            summary_path = os.path.join(root, "summary.txt")
            break

    if summary_path is None:
        print("FAILED (no summary.txt)")
        return None

    # Parse summary
    with open(summary_path) as f:
        txt = f.read()

    f_roi = None
    bg_rate = None
    counts = {}

    m = re.search(r"f\(SS-in-ROI\)\s*=\s*([0-9.e+-]+)", txt)
    if m:
        f_roi = float(m.group(1))

    m = re.search(r"Background rate\s*=\s*([0-9.e+-]+)", txt)
    if m:
        bg_rate = float(m.group(1))

    for m in re.finditer(r"(\w+)\s+(\d+)\s+\(\s*[\d.]+%\)", txt):
        counts[m.group(1)] = int(m.group(2))

    print(f"done  bg={bg_rate:.4e} events/yr" if bg_rate else "done (no result)")

    return {
        "index": index,
        "label": label,
        "f_roi": f_roi,
        "bg_rate": bg_rate,
        "counts": counts,
    }


def main():
    n_samples = 100_000_000
    if len(sys.argv) > 1 and sys.argv[1] == "--n-samples":
        n_samples = int(float(sys.argv[2]))

    os.makedirs(OUTPUT_BASE, exist_ok=True)

    print(f"=" * 70)
    print(f"Running all 11 Bi-214 sources ({n_samples:.0e} events each)")
    print(f"FV: z=[26, 96] cm, r<=39 cm")
    print(f"=" * 70)

    results = []
    for index, label, paper_cts in BI214_SOURCES:
        r = run_source(index, label, n_samples)
        if r:
            r["paper_cts"] = paper_cts
            results.append(r)

    # Print results table
    print(f"\n{'=' * 80}")
    print(f"RESULTS: Bi-214 (²³⁸U-late), {n_samples:.0e} events per source")
    print(f"{'=' * 80}")
    print(f"{'Index':>5s}  {'Label':<25s}  {'f_SS_ROI':>12s}  {'events/yr':>12s}  {'cts/1000d':>10s}")
    print(f"{'-' * 70}")

    total_bg = 0.0
    group_totals = {
        "Cryostat (Ti+MLI)": 0.0,
        "FC rings": 0.0,
        "PTFE walls": 0.0,
        "FC res+sens": 0.0,
        "FC holders": 0.0,
        "PMTs (top+bot)": 0.0,
    }

    for r in results:
        bg = r["bg_rate"] if r["bg_rate"] else 0.0
        cts_1000d = bg * 1000 / 365.25
        total_bg += bg
        print(f"{r['index']:>5d}  {r['label']:<25s}  {r['f_roi']:>12.4e}  {bg:>12.4e}  {cts_1000d:>10.2f}")

        # Assign to groups
        lbl = r["label"]
        if "cryo" in lbl:
            group_totals["Cryostat (Ti+MLI)"] += bg
        elif "rings" in lbl:
            group_totals["FC rings"] += bg
        elif "ptfe" in lbl:
            group_totals["PTFE walls"] += bg
        elif "ressens" in lbl:
            group_totals["FC res+sens"] += bg
        elif "holder" in lbl:
            group_totals["FC holders"] += bg
        elif "pmt" in lbl:
            group_totals["PMTs (top+bot)"] += bg

    print(f"{'-' * 70}")
    print(f"{'':>5s}  {'TOTAL':<25s}  {'':>12s}  {total_bg:>12.4e}  {total_bg*1000/365.25:>10.2f}")

    # Grouped comparison with paper
    print(f"\n{'=' * 70}")
    print(f"COMPARISON WITH LZ bb0nu PAPER")
    print(f"{'=' * 70}")
    print(f"{'Group':<25s}  {'Our (evt/yr)':>12s}  {'Paper (cts/1000d)':>18s}  {'Paper (evt/yr)':>14s}  {'Ratio':>7s}")
    print(f"{'-' * 80}")

    total_paper = 0.0
    for group, paper_cts in PAPER_GROUPS.items():
        our = group_totals[group]
        paper_evts = paper_cts * 365.25 / 1000
        total_paper += paper_evts
        ratio = our / paper_evts if paper_evts > 0 else float('inf')
        print(f"{group:<25s}  {our:>12.4f}  {paper_cts:>18.2f}  {paper_evts:>14.4f}  {ratio:>7.2f}×")

    print(f"{'-' * 80}")
    print(f"{'TOTAL':<25s}  {total_bg:>12.4f}  {sum(PAPER_GROUPS.values()):>18.2f}  {total_paper:>14.4f}  {total_bg/total_paper:>7.2f}×")

    # Save CSV
    csv_path = os.path.join(DATA_DIR, "bi214_results.csv")
    with open(csv_path, "w") as f:
        f.write("index,label,f_SS_ROI,events_per_yr,counts_per_1000d\n")
        for r in results:
            bg = r["bg_rate"] if r["bg_rate"] else 0.0
            f.write(f"{r['index']},{r['label']},{r['f_roi']:.6e},{bg:.6e},{bg*1000/365.25:.4f}\n")
        f.write(f"0,TOTAL,0,{total_bg:.6e},{total_bg*1000/365.25:.4f}\n")
    print(f"\nSaved: {csv_path}")


if __name__ == "__main__":
    main()
