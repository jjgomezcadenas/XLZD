"""Plot the diagnostic histograms produced by `scripts/run_mc3.jl`.

Reads the per-source `hist_<source>/` directory written by run_mc3 and
renders three combined PNG panels:

    stack.png      — StackHistogramSet (10 plots: ng_max, first_interaction,
                     n_photo / n_compton / n_pair / n_below_thresh,
                     inclusive_edep, E_first, delta_z, region_interaction)
    cluster.png    — ClusterHistogramSet (11 plots: Ec / Emax / Emin / Einc,
                     closest/furthest D3, closest/furthest dz,
                     N_clusters, r²-vs-z heatmap, D-vs-z heatmap)
    rejection.png  — RejectionHistograms (4 plots: skin/fv × E and r²-vs-z)

Each family is rendered only when its CSVs are present.

Usage
-----
    python py/plot_histograms.py output/<run>/hist_<source>/
    python py/plot_histograms.py output/<run>/hist_CB_Bi214/ --separate
    python py/plot_histograms.py output/<run>/hist_CB_Bi214/ --log-y --show
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# CSV file inventory by family
# ---------------------------------------------------------------------------
STACK_CSVS = (
    "stack_ng_max.csv",
    "stack_first_interaction.csv",
    "stack_n_photo.csv",
    "stack_n_compton.csv",
    "stack_n_pair.csv",
    "stack_n_below_thresh.csv",
    "stack_inclusive_edep.csv",
    "stack_E_first.csv",
    "stack_delta_z.csv",
    "stack_region_interaction.csv",
)

CLUSTER_CSVS = (
    "cluster_Ec.csv",
    "cluster_Emax.csv",
    "cluster_Emin.csv",
    "cluster_Einc.csv",
    "cluster_closest_D3.csv",
    "cluster_furthest_D3.csv",
    "cluster_closest_dz.csv",
    "cluster_furthest_dz.csv",
    "cluster_N_clusters.csv",
    "cluster_r2_vs_z.csv",
    "cluster_D_vs_z.csv",
)

REJECTION_CSVS = (
    "rejected_skin_E.csv",
    "rejected_skin_r2z.csv",
    "rejected_fv_E.csv",
    "rejected_fv_r2z.csv",
)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing histogram CSVs (hist_<source>/).",
    )
    p.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help="Where to write PNGs (default: same as input_dir).",
    )
    p.add_argument(
        "--family",
        choices=("stack", "cluster", "rejection", "all"),
        default="all",
        help="Which histogram family to plot (default: all detected).",
    )
    p.add_argument(
        "--separate",
        action="store_true",
        help="Also save one PNG per histogram, in addition to the combined panels.",
    )
    p.add_argument(
        "--show",
        action="store_true",
        help="Open the combined plots in a window.",
    )
    p.add_argument(
        "--log-y",
        action="store_true",
        help="Log scale on the y axis for 1D histograms.",
    )
    p.add_argument(
        "--title",
        type=str,
        default=None,
        help="Title prefix for the combined plots. Default: directory name.",
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# IO helpers
# ---------------------------------------------------------------------------
def _read(input_dir: Path, name: str) -> pd.DataFrame | None:
    path = input_dir / name
    if not path.exists():
        return None
    return pd.read_csv(path)


def _present(input_dir: Path, names: tuple[str, ...]) -> bool:
    return all((input_dir / n).exists() for n in names)


# ---------------------------------------------------------------------------
# Generic plot primitives
# ---------------------------------------------------------------------------
def _plot_int_bar(ax: plt.Axes, df: pd.DataFrame, *, title: str, xlabel: str) -> None:
    ax.bar(df.iloc[:, 0], df["count"], color="#684", edgecolor="none")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("events")
    ax.set_title(title)
    nonzero = df[df["count"] > 0]
    if not nonzero.empty and len(nonzero) <= 12:
        for x, y in zip(nonzero.iloc[:, 0], nonzero["count"]):
            ax.text(x, y, f"{y}", ha="center", va="bottom", fontsize=7)


def _plot_categorical_bar(
    ax: plt.Axes, df: pd.DataFrame, *, title: str, xlabel: str
) -> None:
    cats = df.iloc[:, 0].astype(str).tolist()
    ax.bar(cats, df["count"], color=["#3b6", "#46a", "#a64", "#888"])
    ax.set_xlabel(xlabel)
    ax.set_ylabel("events")
    ax.set_title(title)
    for i, y in enumerate(df["count"]):
        ax.text(i, y, f"{y:,}", ha="center", va="bottom", fontsize=7)
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")


def _plot_1d(
    ax: plt.Axes,
    df: pd.DataFrame,
    *,
    title: str,
    xlabel: str,
    log_y: bool = False,
    color: str = "#46a",
) -> None:
    left = df.iloc[:, 0]
    right = df.iloc[:, 1]
    centres = 0.5 * (left + right)
    width = right - left
    ax.bar(centres, df["count"], width=width, color=color, edgecolor="none")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("entries")
    ax.set_title(title)
    if log_y:
        # Avoid log(0): set y_min to 0.5 so empty bars stay visible.
        ax.set_yscale("log")
        if df["count"].max() > 0:
            ax.set_ylim(bottom=0.5)


def _plot_2d_heatmap(
    ax: plt.Axes,
    df: pd.DataFrame,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    """Render a 2D histogram from a CSV with columns:
        x_left, x_right, y_left, y_right, count.
    """
    x_lefts = df.iloc[:, 0]
    y_lefts = df.iloc[:, 2]
    x_edges = np.unique(
        np.concatenate([df.iloc[:, 0].values, df.iloc[:, 1].values])
    )
    y_edges = np.unique(
        np.concatenate([df.iloc[:, 2].values, df.iloc[:, 3].values])
    )
    nx = len(x_edges) - 1
    ny = len(y_edges) - 1
    M = np.zeros((nx, ny), dtype=float)
    # Map (x_left, y_left) → bin indices via searchsorted.
    xi = np.searchsorted(x_edges[:-1], x_lefts.values, side="right") - 1
    yi = np.searchsorted(y_edges[:-1], y_lefts.values, side="right") - 1
    for i, j, c in zip(xi, yi, df["count"].values):
        if 0 <= i < nx and 0 <= j < ny:
            M[i, j] = c
    im = ax.pcolormesh(x_edges, y_edges, M.T, cmap="viridis", shading="flat")
    plt.colorbar(im, ax=ax, label="count")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


def _plot_region_interaction_table(ax: plt.Axes, df: pd.DataFrame) -> None:
    """3×4 region × interaction-type counts — render as a heatmap with text."""
    pivot = df.pivot(index="region", columns="interaction", values="count")
    # Force ordering.
    region_order = [r for r in ("TPC", "Skin", "Inert") if r in pivot.index]
    inter_order = [
        c for c in ("PHOTO", "COMPTON", "PAIR", "BELOW_THRESH") if c in pivot.columns
    ]
    pivot = pivot.loc[region_order, inter_order]
    M = pivot.values.astype(float)
    im = ax.imshow(M, cmap="Blues", aspect="auto")
    ax.set_xticks(range(len(inter_order)), inter_order, rotation=20, ha="right")
    ax.set_yticks(range(len(region_order)), region_order)
    ax.set_title("region × interaction (deposit counts)")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{int(M[i, j]):,}", ha="center", va="center",
                    color="white" if M[i, j] > 0.5 * M.max() else "black",
                    fontsize=8)
    plt.colorbar(im, ax=ax, label="count")


# ---------------------------------------------------------------------------
# Family renderers
# ---------------------------------------------------------------------------
def render_stack_panel(input_dir: Path, title: str, log_y: bool) -> plt.Figure:
    fig, axes = plt.subplots(3, 4, figsize=(20, 11))
    fig.suptitle(f"{title} — stack histograms")
    flat = axes.flatten()

    df = _read(input_dir, "stack_ng_max.csv")
    _plot_int_bar(flat[0], df, title="ng_max (chain depth)", xlabel="ng_max")

    df = _read(input_dir, "stack_first_interaction.csv")
    _plot_categorical_bar(flat[1], df,
                          title="first interaction type",
                          xlabel="interaction")

    df = _read(input_dir, "stack_n_photo.csv")
    _plot_int_bar(flat[2], df, title="n_photo per event", xlabel="n")

    df = _read(input_dir, "stack_n_compton.csv")
    _plot_int_bar(flat[3], df, title="n_compton per event", xlabel="n")

    df = _read(input_dir, "stack_n_pair.csv")
    _plot_int_bar(flat[4], df, title="n_pair per event", xlabel="n")

    df = _read(input_dir, "stack_n_below_thresh.csv")
    _plot_int_bar(flat[5], df, title="n_below_thresh per event", xlabel="n")

    df = _read(input_dir, "stack_inclusive_edep.csv")
    _plot_1d(flat[6], df, title="inclusive Σ edep",
             xlabel="E (MeV)", log_y=log_y, color="#a64")

    df = _read(input_dir, "stack_E_first.csv")
    _plot_1d(flat[7], df, title="E of first :TPC deposit",
             xlabel="E (MeV)", log_y=log_y, color="#a64")

    df = _read(input_dir, "stack_delta_z.csv")
    _plot_1d(flat[8], df, title="Δz from first :TPC deposit",
             xlabel="Δz (cm)", log_y=log_y, color="#46a")

    df = _read(input_dir, "stack_region_interaction.csv")
    _plot_region_interaction_table(flat[9], df)

    # Hide unused panels.
    for k in (10, 11):
        flat[k].axis("off")

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def render_cluster_panel(input_dir: Path, title: str, log_y: bool) -> plt.Figure:
    fig, axes = plt.subplots(3, 4, figsize=(20, 11))
    fig.suptitle(f"{title} — cluster histograms")
    flat = axes.flatten()

    df = _read(input_dir, "cluster_Ec.csv")
    _plot_1d(flat[0], df, title="cluster energy Ec",
             xlabel="E (MeV)", log_y=log_y, color="#a64")

    df = _read(input_dir, "cluster_Emax.csv")
    _plot_1d(flat[1], df, title="per-event Emax",
             xlabel="E (MeV)", log_y=log_y, color="#a64")

    df = _read(input_dir, "cluster_Emin.csv")
    _plot_1d(flat[2], df, title="per-event Emin",
             xlabel="E (MeV)", log_y=log_y, color="#a64")

    df = _read(input_dir, "cluster_Einc.csv")
    _plot_1d(flat[3], df, title="per-event Σ Ec",
             xlabel="E (MeV)", log_y=log_y, color="#a64")

    df = _read(input_dir, "cluster_closest_D3.csv")
    _plot_1d(flat[4], df, title="closest pair distance (3D)",
             xlabel="D (cm)", log_y=log_y, color="#46a")

    df = _read(input_dir, "cluster_furthest_D3.csv")
    _plot_1d(flat[5], df, title="furthest pair distance (3D)",
             xlabel="D (cm)", log_y=log_y, color="#46a")

    df = _read(input_dir, "cluster_closest_dz.csv")
    _plot_1d(flat[6], df, title="closest pair |Δz|",
             xlabel="|Δz| (cm)", log_y=log_y, color="#46a")

    df = _read(input_dir, "cluster_furthest_dz.csv")
    _plot_1d(flat[7], df, title="furthest pair |Δz|",
             xlabel="|Δz| (cm)", log_y=log_y, color="#46a")

    df = _read(input_dir, "cluster_N_clusters.csv")
    _plot_int_bar(flat[8], df, title="N clusters per photon", xlabel="n")

    df = _read(input_dir, "cluster_r2_vs_z.csv")
    _plot_2d_heatmap(flat[9], df,
                     title="r² vs z (cluster centroid)",
                     xlabel="r² (cm²)", ylabel="z (cm)")

    df = _read(input_dir, "cluster_D_vs_z.csv")
    _plot_2d_heatmap(flat[10], df,
                     title="D vs z (cluster centroid)",
                     xlabel="D (cm)", ylabel="z (cm)")

    flat[11].axis("off")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def render_rejection_panel(input_dir: Path, title: str, log_y: bool) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"{title} — fast-veto rejection histograms")

    df = _read(input_dir, "rejected_skin_E.csv")
    _plot_1d(axes[0, 0], df, title="rejected_skin: triggering E",
             xlabel="E (MeV)", log_y=log_y, color="#c54")

    df = _read(input_dir, "rejected_skin_r2z.csv")
    _plot_2d_heatmap(axes[0, 1], df,
                     title="rejected_skin: r² vs z",
                     xlabel="r² (cm²)", ylabel="z (cm)")

    df = _read(input_dir, "rejected_fv_E.csv")
    _plot_1d(axes[1, 0], df, title="rejected_fv: triggering E",
             xlabel="E (MeV)", log_y=log_y, color="#c54")

    df = _read(input_dir, "rejected_fv_r2z.csv")
    _plot_2d_heatmap(axes[1, 1], df,
                     title="rejected_fv: r² vs z",
                     xlabel="r² (cm²)", ylabel="z (cm)")

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


# ---------------------------------------------------------------------------
# Per-histogram (--separate) renderer
# ---------------------------------------------------------------------------
def _save_one(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=150)
    print(f"  wrote {path}")


def render_individual(input_dir: Path, out: Path, log_y: bool) -> None:
    """One PNG per CSV; auto-dispatch by filename suffix / column shape."""
    for csv in (*STACK_CSVS, *CLUSTER_CSVS, *REJECTION_CSVS):
        df = _read(input_dir, csv)
        if df is None:
            continue
        stem = csv.removesuffix(".csv")
        cols = list(df.columns)

        fig, ax = plt.subplots(figsize=(7, 4.5))
        if cols == ["interaction", "count"]:
            _plot_categorical_bar(ax, df, title=stem, xlabel="interaction")
        elif cols == ["region", "interaction", "count"]:
            _plot_region_interaction_table(ax, df)
        elif cols[-1] == "count" and len(cols) == 2:
            _plot_int_bar(ax, df, title=stem, xlabel=cols[0])
        elif cols[-1] == "count" and len(cols) == 3:
            _plot_1d(ax, df, title=stem, xlabel=cols[0], log_y=log_y)
        elif cols[-1] == "count" and len(cols) == 5:
            fig.set_size_inches(8, 6)
            _plot_2d_heatmap(ax, df, title=stem, xlabel=cols[0], ylabel=cols[2])
        else:
            ax.text(0.5, 0.5, f"unknown shape:\n{cols}",
                    ha="center", va="center")
        fig.tight_layout()
        _save_one(fig, out / f"{stem}.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()
    if not args.input_dir.is_dir():
        sys.exit(f"not a directory: {args.input_dir}")
    out = args.output if args.output else args.input_dir
    out.mkdir(parents=True, exist_ok=True)
    title = args.title if args.title else args.input_dir.name

    families = {
        "stack":     (STACK_CSVS,     "stack.png",     render_stack_panel),
        "cluster":   (CLUSTER_CSVS,   "cluster.png",   render_cluster_panel),
        "rejection": (REJECTION_CSVS, "rejection.png", render_rejection_panel),
    }
    todo = (
        ["stack", "cluster", "rejection"] if args.family == "all" else [args.family]
    )

    for fam in todo:
        csvs, png_name, renderer = families[fam]
        if not _present(args.input_dir, csvs):
            print(f"  skip {fam} (CSVs missing in {args.input_dir})")
            continue
        fig = renderer(args.input_dir, title=title, log_y=args.log_y)
        _save_one(fig, out / png_name)

    if args.separate:
        render_individual(args.input_dir, out, log_y=args.log_y)

    if args.show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()
