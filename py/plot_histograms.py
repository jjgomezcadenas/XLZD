"""Plot the diagnostic histograms produced by `scripts/run_mc3.jl`.

Reads the per-source `hist_<source>/` directory written by run_mc3 and
renders three combined PNG panels:

    stack.png      — StackHistogramSet (9 plots: first_interaction,
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

# Detection list for the cluster.png panel. Some CSVs (Emax, Emin, Einc,
# D_vs_z) are still produced by the MC but intentionally NOT plotted in
# the combined panel — they are kept on disk for ad-hoc analysis.
CLUSTER_CSVS = (
    "cluster_Ec.csv",
    "cluster_closest_D3.csv",
    "cluster_furthest_D3.csv",
    "cluster_closest_dz.csv",
    "cluster_furthest_dz.csv",
    "cluster_N_clusters.csv",
    "cluster_r2_vs_z.csv",
    "cluster_Ec_vs_dz.csv",
)

REJECTION_CSVS = (
    "rejected_skin_E.csv",
    "rejected_skin_r2z.csv",
    "rejected_fv_E.csv",
    "rejected_fv_r2z.csv",
)

# SS pre-ROI spectrum: the cluster energy the ROI cut would act on, for
# events that survived skin / FV / SC. Two CSVs (true ec, smeared es) →
# two panels in SS_energy_precut.png with the ROI window overlaid.
SS_PRECUT_CSVS = (
    "cluster_ss_ec_pre_roi.csv",
    "cluster_ss_es_pre_roi.csv",
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
        choices=("stack", "cluster", "rejection", "ss_precut", "all"),
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
def _plot_int_bar(
    ax: plt.Axes, df: pd.DataFrame, *,
    title: str, xlabel: str,
    log_y: bool = False,
    normalize: bool = False,
    show_values: bool = True,
    xlim_max: int | None = None,
) -> None:
    counts = df["count"].astype(float)
    total = counts.sum()
    values = counts / total if (normalize and total > 0) else counts
    ylabel = "fraction" if normalize else "events"
    ax.bar(df.iloc[:, 0], values, color="#684", edgecolor="none")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if show_values:
        nonzero_mask = counts > 0
        nz_x = df.iloc[:, 0][nonzero_mask]
        nz_y = values[nonzero_mask]
        if len(nz_x) <= 12:
            for x, y in zip(nz_x, nz_y):
                ax.text(x, y, f"{y:.2e}" if normalize else f"{int(y)}",
                        ha="center", va="bottom", fontsize=7)
    if log_y and (values.max() > 0):
        ax.set_yscale("log")
        ax.set_ylim(bottom=values[values > 0].min() / 2 if normalize else 0.5)
    if xlim_max is not None:
        ax.set_xlim(-0.5, xlim_max + 0.5)


def _plot_categorical_bar(
    ax: plt.Axes, df: pd.DataFrame, *,
    title: str, xlabel: str,
    log_y: bool = False,
    normalize: bool = False,
) -> None:
    counts = df["count"].astype(float)
    total = counts.sum()
    values = counts / total if (normalize and total > 0) else counts
    ylabel = "fraction" if normalize else "events"
    cats = df.iloc[:, 0].astype(str).tolist()
    ax.bar(cats, values, color=["#3b6", "#46a", "#a64", "#888"])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    # Place text JUST INSIDE the top of each bar so the value sits
    # at the same visual height regardless of bar magnitude (important
    # on log scale where centre-of-bar lands far below short bars).
    for i, y in enumerate(values):
        if y > 0:
            ax.text(i, y,
                    f"{y:.2e}" if normalize else f"{int(y):,}",
                    ha="center", va="top", fontsize=11, color="white",
                    weight="bold")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", fontsize=8)
    if log_y and values.max() > 0:
        ax.set_yscale("log")
        ax.set_ylim(bottom=values[values > 0].min() / 2 if normalize else 0.5)


def _plot_1d(
    ax: plt.Axes,
    df: pd.DataFrame,
    *,
    title: str,
    xlabel: str,
    log_y: bool = False,
    normalize: bool = False,
    color: str = "#46a",
) -> None:
    left = df.iloc[:, 0]
    right = df.iloc[:, 1]
    counts = df["count"].astype(float)
    total = counts.sum()
    values = counts / total if (normalize and total > 0) else counts
    ylabel = "fraction" if normalize else "entries"
    centres = 0.5 * (left + right)
    width = right - left
    ax.bar(centres, values, width=width, color=color, edgecolor="none")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if log_y and values.max() > 0:
        ax.set_yscale("log")
        ax.set_ylim(bottom=values[values > 0].min() / 2 if normalize else 0.5)


def _r2_to_r(df: pd.DataFrame) -> pd.DataFrame:
    """Convert a 2D-heatmap DataFrame whose first two columns are r²
    bin edges (cm²) into one with r bin edges (cm). Bin widths in r
    become non-uniform (smaller at large r), but pcolormesh handles
    non-uniform edges fine."""
    out = df.copy()
    out.iloc[:, 0] = np.sqrt(out.iloc[:, 0].values)
    out.iloc[:, 1] = np.sqrt(out.iloc[:, 1].values)
    return out


def _plot_2d_heatmap(
    ax: plt.Axes,
    df: pd.DataFrame,
    *,
    title: str,
    rlabel: str,   # axis label for the r²/D dimension (CSV columns 0,1)
    zlabel: str,   # axis label for the z dimension     (CSV columns 2,3)
) -> None:
    """Render a 2D histogram from a CSV with columns:
        rL, rR, zL, zR, count.

    Rendered with z on X axis and r² (or D) on Y axis (physics convention).
    """
    r_edges = np.unique(
        np.concatenate([df.iloc[:, 0].values, df.iloc[:, 1].values])
    )
    z_edges = np.unique(
        np.concatenate([df.iloc[:, 2].values, df.iloc[:, 3].values])
    )
    nr = len(r_edges) - 1
    nz = len(z_edges) - 1
    M = np.zeros((nr, nz), dtype=float)
    ri = np.searchsorted(r_edges[:-1], df.iloc[:, 0].values, side="right") - 1
    zi = np.searchsorted(z_edges[:-1], df.iloc[:, 2].values, side="right") - 1
    for i, j, c in zip(ri, zi, df["count"].values):
        if 0 <= i < nr and 0 <= j < nz:
            M[i, j] = c
    # Place z on X axis, r²/D on Y axis: pcolormesh(x_edges, y_edges, C[ny,nx]).
    im = ax.pcolormesh(z_edges, r_edges, M, cmap="viridis", shading="flat")
    plt.colorbar(im, ax=ax, label="count")
    ax.set_xlabel(zlabel)
    ax.set_ylabel(rlabel)
    ax.set_title(title)
    # Span the full bin range (full LXe region), not just where data falls.
    ax.set_xlim(z_edges[0], z_edges[-1])
    ax.set_ylim(r_edges[0], r_edges[-1])


def _plot_region_interaction_table(ax: plt.Axes, df: pd.DataFrame) -> None:
    """3×4 region × interaction-type counts — render as a heatmap of
    fractions of total deposits with the value as text in each cell."""
    pivot = df.pivot(index="region", columns="interaction", values="count")
    # Force ordering.
    region_order = [r for r in ("TPC", "Skin", "Inert") if r in pivot.index]
    inter_order = [
        c for c in ("PHOTO", "COMPTON", "PAIR", "BELOW_THRESH") if c in pivot.columns
    ]
    pivot = pivot.loc[region_order, inter_order]
    M = pivot.values.astype(float)
    total = M.sum()
    F = M / total if total > 0 else M
    im = ax.imshow(F, cmap="Blues", aspect="auto")
    ax.set_xticks(range(len(inter_order)), inter_order, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(region_order)), region_order)
    ax.set_title("region × interaction (fraction of all deposits)")
    fmax = F.max() if F.size else 0.0
    for i in range(F.shape[0]):
        for j in range(F.shape[1]):
            ax.text(j, i, f"{F[i, j]:.2e}", ha="center", va="center",
                    color="white" if F[i, j] > 0.5 * fmax else "black",
                    fontsize=12, weight="bold")
    plt.colorbar(im, ax=ax, label="fraction")


# ---------------------------------------------------------------------------
# Family renderers
# ---------------------------------------------------------------------------
def render_stack_panel(input_dir: Path, title: str, log_y: bool) -> plt.Figure:
    # 9 plots in a 3×3 grid (ng_max dropped per user preference).
    fig, axes = plt.subplots(3, 3, figsize=(20, 14), constrained_layout=True)
    fig.suptitle(f"{title} — stack histograms", fontsize=14)
    flat = axes.flatten()

    df = _read(input_dir, "stack_first_interaction.csv")
    _plot_categorical_bar(flat[0], df,
                          title="first interaction type",
                          xlabel="interaction", log_y=True,
                          normalize=True)

    df = _read(input_dir, "stack_n_photo.csv")
    _plot_int_bar(flat[1], df, title="n_photo per event", xlabel="n",
                  log_y=True, normalize=True,
                  show_values=False, xlim_max=10)

    df = _read(input_dir, "stack_n_compton.csv")
    _plot_int_bar(flat[2], df, title="n_compton per event", xlabel="n",
                  log_y=True, normalize=True, show_values=False)

    df = _read(input_dir, "stack_n_pair.csv")
    _plot_int_bar(flat[3], df, title="n_pair per event", xlabel="n",
                  log_y=True, normalize=True,
                  show_values=False, xlim_max=5)

    df = _read(input_dir, "stack_n_below_thresh.csv")
    _plot_int_bar(flat[4], df, title="n_below_thresh per event", xlabel="n",
                  log_y=True, normalize=True, show_values=False)

    df = _read(input_dir, "stack_inclusive_edep.csv")
    _plot_1d(flat[5], df, title="inclusive Σ edep",
             xlabel="E (MeV)", log_y=True, normalize=True, color="#a64")

    df = _read(input_dir, "stack_E_first.csv")
    _plot_1d(flat[6], df, title="E of first :TPC deposit",
             xlabel="E (MeV)", log_y=True, normalize=True, color="#a64")

    df = _read(input_dir, "stack_delta_z.csv")
    _plot_1d(flat[7], df, title="Δz from first :TPC deposit",
             xlabel="Δz (cm)", log_y=True, normalize=True, color="#46a")

    df = _read(input_dir, "stack_region_interaction.csv")
    _plot_region_interaction_table(flat[8], df)

    return fig


def render_cluster_panel(input_dir: Path, title: str, log_y: bool) -> plt.Figure:
    # 11 plots in a 3×4 grid (1 cell unused).
    # 8 plots in a 2×4 grid (Emax / Emin / Einc / D_vs_z dropped per
    # user preference; CSVs are still on disk for ad-hoc analysis).
    fig, axes = plt.subplots(2, 4, figsize=(22, 11), constrained_layout=True)
    fig.suptitle(f"{title} — cluster histograms", fontsize=14)
    flat = axes.flatten()

    df = _read(input_dir, "cluster_Ec.csv")
    _plot_1d(flat[0], df, title="cluster energy Ec",
             xlabel="E (MeV)", log_y=True, color="#a64")

    df = _read(input_dir, "cluster_N_clusters.csv")
    _plot_int_bar(flat[1], df, title="N clusters per photon", xlabel="n",
                  log_y=True)

    df = _read(input_dir, "cluster_closest_D3.csv")
    _plot_1d(flat[2], df, title="closest pair distance (3D)",
             xlabel="D (cm)", log_y=True, color="#46a")

    df = _read(input_dir, "cluster_furthest_D3.csv")
    _plot_1d(flat[3], df, title="furthest pair distance (3D)",
             xlabel="D (cm)", log_y=True, color="#46a")

    df = _read(input_dir, "cluster_closest_dz.csv")
    _plot_1d(flat[4], df, title="closest pair |Δz|",
             xlabel="|Δz| (cm)", log_y=True, color="#46a")

    df = _read(input_dir, "cluster_furthest_dz.csv")
    _plot_1d(flat[5], df, title="furthest pair |Δz|",
             xlabel="|Δz| (cm)", log_y=True, color="#46a")

    df = _read(input_dir, "cluster_r2_vs_z.csv")
    _plot_2d_heatmap(flat[6], _r2_to_r(df),
                     title="r vs z (cluster centroid)",
                     rlabel="r (cm)", zlabel="z (cm)")

    df = _read(input_dir, "cluster_Ec_vs_dz.csv")
    _plot_2d_heatmap(flat[7], df,
                     title="Ec vs |Δz to nearest cluster| (per cluster)",
                     rlabel="Ec (MeV)", zlabel="|Δz| (cm)")
    return fig


def render_rejection_panel(input_dir: Path, title: str, log_y: bool) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    fig.suptitle(f"{title} — fast-veto rejection histograms", fontsize=14)

    df = _read(input_dir, "rejected_skin_E.csv")
    _plot_1d(axes[0, 0], df, title="rejected_skin: triggering E",
             xlabel="E (MeV)", log_y=True, color="#c54")

    df = _read(input_dir, "rejected_skin_r2z.csv")
    _plot_2d_heatmap(axes[0, 1], _r2_to_r(df),
                     title="rejected_skin: r vs z",
                     rlabel="r (cm)", zlabel="z (cm)")

    df = _read(input_dir, "rejected_fv_E.csv")
    _plot_1d(axes[1, 0], df, title="rejected_fv: triggering E",
             xlabel="E (MeV)", log_y=True, color="#c54")

    df = _read(input_dir, "rejected_fv_r2z.csv")
    _plot_2d_heatmap(axes[1, 1], _r2_to_r(df),
                     title="rejected_fv: r vs z",
                     rlabel="r (cm)", zlabel="z (cm)")
    return fig


# ---------------------------------------------------------------------------
# SS pre-ROI energy spectrum (true ec + smeared es) with ROI window overlay
# ---------------------------------------------------------------------------
def _read_cuts_from_summary(input_dir: Path) -> tuple[float, float]:
    """Return (Q_ββ in MeV, ROI half-width in MeV) read from
    `<parent>/summary.csv`. Falls back to MCParams defaults
    (Q_ββ = 2.458 MeV, halfwidth = 0.0172 MeV) if the file is missing or
    the columns aren't present (older runs)."""
    Q_DEFAULT_MeV  = 2.458
    HW_DEFAULT_MeV = 0.0172
    summary_path = input_dir.parent / "summary.csv"
    if not summary_path.is_file():
        return Q_DEFAULT_MeV, HW_DEFAULT_MeV
    try:
        s = pd.read_csv(summary_path)
        if "Q_betabeta_keV" in s.columns and "ROI_halfwidth_keV" in s.columns:
            return float(s["Q_betabeta_keV"].iloc[0]) / 1000.0, \
                   float(s["ROI_halfwidth_keV"].iloc[0]) / 1000.0
    except Exception:
        pass
    return Q_DEFAULT_MeV, HW_DEFAULT_MeV


def _plot_1d_with_roi(
    ax: plt.Axes,
    df: pd.DataFrame,
    *,
    title: str,
    Q_MeV: float,
    HW_MeV: float,
    color: str,
) -> None:
    _plot_1d(ax, df, title=title, xlabel="E (MeV)", log_y=True, color=color)
    # ROI window: shaded band + dashed line at Q_ββ.
    ax.axvspan(Q_MeV - HW_MeV, Q_MeV + HW_MeV, color="#fc6", alpha=0.35,
               label=f"ROI: Q_ββ ± {HW_MeV*1000:.1f} keV")
    ax.axvline(Q_MeV, color="#a40", linestyle="--", linewidth=1.2,
               label=f"Q_ββ = {Q_MeV*1000:.0f} keV")
    ax.legend(loc="upper left", fontsize=8)


def render_ss_precut_panel(input_dir: Path, title: str, log_y: bool) -> plt.Figure:
    Q_MeV, HW_MeV = _read_cuts_from_summary(input_dir)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
    fig.suptitle(f"{title} — SS energy pre-ROI cut", fontsize=14)

    df = _read(input_dir, "cluster_ss_ec_pre_roi.csv")
    _plot_1d_with_roi(axes[0], df,
                       title="SS true energy ec (no smearing)",
                       Q_MeV=Q_MeV, HW_MeV=HW_MeV, color="#46a")

    df = _read(input_dir, "cluster_ss_es_pre_roi.csv")
    _plot_1d_with_roi(axes[1], df,
                       title="SS smeared energy es (ROI cut acts on this)",
                       Q_MeV=Q_MeV, HW_MeV=HW_MeV, color="#a64")
    return fig


# ---------------------------------------------------------------------------
# Per-histogram (--separate) renderer
# ---------------------------------------------------------------------------
def _save_one(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=150)
    print(f"  wrote {path}")


def render_individual(input_dir: Path, out: Path, log_y: bool) -> None:
    """One PNG per CSV; auto-dispatch by filename suffix / column shape."""
    for csv in (*STACK_CSVS, *CLUSTER_CSVS, *REJECTION_CSVS, *SS_PRECUT_CSVS):
        df = _read(input_dir, csv)
        if df is None:
            continue
        stem = csv.removesuffix(".csv")
        cols = list(df.columns)

        fig, ax = plt.subplots(figsize=(7, 4.5))
        if cols == ["interaction", "count"]:
            _plot_categorical_bar(ax, df, title=stem, xlabel="interaction", log_y=True)
        elif cols == ["region", "interaction", "count"]:
            _plot_region_interaction_table(ax, df)
        elif cols[-1] == "count" and len(cols) == 2:
            _plot_int_bar(ax, df, title=stem, xlabel=cols[0])
        elif cols[-1] == "count" and len(cols) == 3:
            _plot_1d(ax, df, title=stem, xlabel=cols[0], log_y=True)
        elif cols[-1] == "count" and len(cols) == 5:
            fig.set_size_inches(8, 6)
            # Convert r²-vs-z CSVs to r-vs-z on the fly.
            is_r2 = "r2" in cols[0] or "r²" in cols[0]
            if is_r2:
                _plot_2d_heatmap(ax, _r2_to_r(df),
                                 title=stem.replace("r2", "r"),
                                 rlabel="r (cm)", zlabel=cols[2])
            else:
                _plot_2d_heatmap(ax, df, title=stem,
                                 rlabel=cols[0], zlabel=cols[2])
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
        "stack":     (STACK_CSVS,     "stack.png",             render_stack_panel),
        "cluster":   (CLUSTER_CSVS,   "cluster.png",           render_cluster_panel),
        "rejection": (REJECTION_CSVS, "rejection.png",         render_rejection_panel),
        "ss_precut": (SS_PRECUT_CSVS, "SS_energy_precut.png",  render_ss_precut_panel),
    }
    todo = (
        ["stack", "cluster", "rejection", "ss_precut"]
        if args.family == "all" else [args.family]
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
