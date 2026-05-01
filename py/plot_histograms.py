"""Plot the six MC control histograms produced by `scripts/run_mc2.jl`.

Reads six CSVs from a `hist_<source>` directory:
    ssms_counts.csv      delta_z.csv         e_first.csv
    e_cluster.csv        n_clusters.csv      n_extra_clusters.csv

Renders one combined 2x3 PNG (default) and optionally six individual PNGs.

Usage
-----
    python py/plot_histograms.py output/<run>/hist_<source>/
    python py/plot_histograms.py output/<run>/hist_CBH_Tl208/ --separate
    python py/plot_histograms.py output/<run>/hist_CBH_Tl208/ --log-y --show
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


CSV_NAMES: tuple[str, ...] = (
    "ssms_counts.csv",
    "delta_z.csv",
    "e_first.csv",
    "e_cluster.csv",
    "n_clusters.csv",
    "n_extra_clusters.csv",
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing the six histogram CSVs (hist_<source>/).",
    )
    p.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help="Where to write PNGs (default: same as input_dir).",
    )
    p.add_argument(
        "--separate",
        action="store_true",
        help="Also save one PNG per histogram, in addition to the combined panel.",
    )
    p.add_argument(
        "--show",
        action="store_true",
        help="Open the combined plot in a window.",
    )
    p.add_argument(
        "--log-y",
        action="store_true",
        help="Log scale on the y axis for histogram panels (Δz / E_first / E_cluster).",
    )
    p.add_argument(
        "--title",
        type=str,
        default=None,
        help="Title for the combined plot. Default: directory name.",
    )
    return p.parse_args()


def _read(input_dir: Path, name: str) -> pd.DataFrame:
    path = input_dir / name
    if not path.exists():
        sys.exit(f"missing CSV: {path}")
    return pd.read_csv(path)


def plot_ssms(ax: plt.Axes, df: pd.DataFrame) -> None:
    counts = df.set_index("bucket")["count"]
    ax.bar(counts.index, counts.values, color=["#3b6", "#c54", "#888"])
    ax.set_ylabel("events")
    ax.set_title("SS / MS / no_cluster")
    for x, y in zip(counts.index, counts.values):
        ax.text(x, y, f"{y:,}", ha="center", va="bottom", fontsize=8)


def plot_delta_z(ax: plt.Axes, df: pd.DataFrame, log_y: bool = False) -> None:
    centres = 0.5 * (df["bin_left_cm"] + df["bin_right_cm"])
    width = df["bin_right_cm"] - df["bin_left_cm"]
    ax.bar(centres, df["count"], width=width, color="#46a", edgecolor="none")
    ax.set_xlabel("Δz from first interaction (cm)")
    ax.set_ylabel("entries")
    ax.set_title("Δz to subsequent interactions")
    if log_y:
        ax.set_yscale("log")


def plot_e_hist(
    ax: plt.Axes,
    df: pd.DataFrame,
    title: str,
    log_y: bool = False,
) -> None:
    centres = 0.5 * (df["bin_left_MeV"] + df["bin_right_MeV"])
    width = df["bin_right_MeV"] - df["bin_left_MeV"]
    ax.bar(centres, df["count"], width=width, color="#a64", edgecolor="none")
    ax.set_xlabel("E (MeV)")
    ax.set_ylabel("entries")
    ax.set_title(title)
    if log_y:
        ax.set_yscale("log")


def plot_n(ax: plt.Axes, df: pd.DataFrame, title: str) -> None:
    nz = df[df["count"] > 0]
    ax.bar(df["n"], df["count"], color="#684", edgecolor="none")
    ax.set_xlabel("n")
    ax.set_ylabel("events")
    ax.set_title(title)
    ax.set_xticks(df["n"][:: max(1, len(df) // 10)])
    if not nz.empty:
        for x, y in zip(nz["n"], nz["count"]):
            ax.text(x, y, f"{y}", ha="center", va="bottom", fontsize=7)


def render_combined(
    dfs: dict[str, pd.DataFrame],
    title: str,
    log_y: bool,
) -> plt.Figure:
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle(title)
    plot_ssms(axes[0, 0], dfs["ssms_counts.csv"])
    plot_delta_z(axes[0, 1], dfs["delta_z.csv"], log_y=log_y)
    plot_e_hist(axes[0, 2], dfs["e_first.csv"], "E of first interaction", log_y=log_y)
    plot_e_hist(axes[1, 0], dfs["e_cluster.csv"], "E of clusters", log_y=log_y)
    plot_n(axes[1, 1], dfs["n_clusters.csv"], "N clusters per photon")
    plot_n(axes[1, 2], dfs["n_extra_clusters.csv"], "N extra clusters (≥2nd)")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def render_individual(
    dfs: dict[str, pd.DataFrame],
    log_y: bool,
) -> dict[str, plt.Figure]:
    figs: dict[str, plt.Figure] = {}

    fig, ax = plt.subplots(figsize=(6, 4))
    plot_ssms(ax, dfs["ssms_counts.csv"])
    fig.tight_layout()
    figs["ssms_counts"] = fig

    fig, ax = plt.subplots(figsize=(7, 4))
    plot_delta_z(ax, dfs["delta_z.csv"], log_y=log_y)
    fig.tight_layout()
    figs["delta_z"] = fig

    fig, ax = plt.subplots(figsize=(7, 4))
    plot_e_hist(ax, dfs["e_first.csv"], "E of first interaction", log_y=log_y)
    fig.tight_layout()
    figs["e_first"] = fig

    fig, ax = plt.subplots(figsize=(7, 4))
    plot_e_hist(ax, dfs["e_cluster.csv"], "E of clusters", log_y=log_y)
    fig.tight_layout()
    figs["e_cluster"] = fig

    fig, ax = plt.subplots(figsize=(7, 4))
    plot_n(ax, dfs["n_clusters.csv"], "N clusters per photon")
    fig.tight_layout()
    figs["n_clusters"] = fig

    fig, ax = plt.subplots(figsize=(7, 4))
    plot_n(ax, dfs["n_extra_clusters.csv"], "N extra clusters (≥2nd)")
    fig.tight_layout()
    figs["n_extra_clusters"] = fig

    return figs


def main() -> None:
    args = parse_args()
    if not args.input_dir.is_dir():
        sys.exit(f"not a directory: {args.input_dir}")
    out = args.output if args.output else args.input_dir
    out.mkdir(parents=True, exist_ok=True)
    title = args.title if args.title else args.input_dir.name

    dfs = {name: _read(args.input_dir, name) for name in CSV_NAMES}

    fig = render_combined(dfs, title=title, log_y=args.log_y)
    combined_path = out / "histograms.png"
    fig.savefig(combined_path, dpi=150)
    print(f"  wrote {combined_path}")

    if args.separate:
        for name, f in render_individual(dfs, log_y=args.log_y).items():
            p = out / f"{name}.png"
            f.savefig(p, dpi=150)
            print(f"  wrote {p}")

    if args.show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()
