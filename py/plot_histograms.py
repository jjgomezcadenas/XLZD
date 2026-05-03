"""Plot the cut-flow + diagnostic histograms produced by `scripts/run_mc3.jl`.

Reads a per-source `output/<run-name>/<source>/` directory and renders
three combined PNG panels organised by analysis stage:

    cuts_acceptance.png       — cut 1 (h_u_sampled, geom-acceptance line)
                                + cut 2 (first-interaction r vs z, FV box).
    cuts_classification.png   — cut 3 (Δz inclusive with 3 mm marker,
                                n_visible, E_total) on top row;
                                cut 4 (ss_ec / ss_es with ROI band, ss_r_z
                                with FV box) on bottom row.
    diagnostics.png           — interaction_type_freq, path_length_LXe,
                                region_interaction matrix, per-cluster Ec.

Each family is rendered only when its CSVs are present. Cut overlays
(FV box, ROI band, 3 mm Δz line, geometric-acceptance u_min) are read
from the per-source `summary.csv` in the same directory.

Usage
-----
    python py/plot_histograms.py output/<run>/<source>/
    python py/plot_histograms.py output/<run>/CB_Bi214/ --separate
    python py/plot_histograms.py output/<run>/CB_Bi214/ --family classification
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
ACCEPTANCE_CSVS = (
    "cut1_h_u_sampled.csv",
    "cut2_first_interaction_r_z.csv",
)

CLASSIFICATION_CSVS = (
    "cut3_dz_inclusive.csv",
    "cut3_n_visible.csv",
    "cut3_E_total.csv",
    "cut4_ss_ec_pre_roi.csv",
    "cut4_ss_es_pre_roi.csv",
    "cut4_ss_r_z.csv",
)

DIAGNOSTIC_CSVS = (
    "diag_interaction_type_freq.csv",
    "diag_path_length_LXe.csv",
    "diag_region_interaction.csv",
    "diag_cluster_Ec.csv",
)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "input_dir",
        type=Path,
        help="Per-source directory: output/<run-name>/<source>/.",
    )
    p.add_argument(
        "--output", "-o",
        type=Path,
        default=None,
        help="Where to write PNGs (default: same as input_dir).",
    )
    p.add_argument(
        "--family",
        choices=("acceptance", "classification", "diagnostics", "all"),
        default="all",
        help="Which panel to plot (default: all detected).",
    )
    p.add_argument(
        "--separate",
        action="store_true",
        help="Also save one PNG per histogram, in addition to the panels.",
    )
    p.add_argument(
        "--show",
        action="store_true",
        help="Open the plots in a window after saving.",
    )
    p.add_argument(
        "--title",
        type=str,
        default=None,
        help="Title prefix for the panels. Default: input directory name.",
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
# Cuts inventory: read run-wide cuts from <run>/summary.csv (one level up
# from the per-source histogram folder). Returns a dict with the values
# the plotter needs for overlays. Falls back to MCParams defaults if the
# file is missing or columns aren't present (older runs).
# ---------------------------------------------------------------------------
_DEFAULT_CUTS: dict[str, float] = {
    "Q_betabeta_keV":     2458.0,
    "ROI_halfwidth_keV":  17.2,
    "fv_z_min_cm":        26.0,
    "fv_z_max_cm":        96.0,
    "fv_r_max_cm":        39.0,
    "Δz_threshold_mm":    3.0,
    "R_LXe_outer_cm":     80.0,
}


def _read_cuts(input_dir: Path) -> dict[str, float]:
    """Read run-wide cut config from `input_dir/summary.csv`. Each
    per-source dir holds its own single-row summary.csv; the cut columns
    are identical across sources, so any one row is sufficient. Falls
    back to MCParams defaults if the file is missing or columns aren't
    present (older runs)."""
    cuts = dict(_DEFAULT_CUTS)
    summary_path = input_dir / "summary.csv"
    if not summary_path.is_file():
        return cuts
    try:
        s = pd.read_csv(summary_path)
        for k in cuts:
            if k in s.columns:
                cuts[k] = float(s[k].iloc[0])
    except Exception:
        pass
    return cuts


def _u_min_geometric(cuts: dict[str, float], source_name: str) -> float | None:
    """Geometric-acceptance threshold for cut-1 markers.
    Definition (barrel only): u_min = sqrt(1 − (r_FV/R_LXe)²) — below this
    u, the photon's straight-line chord cannot reach the FV cylinder.
    Returns None for endcap sources (no single-line marker is meaningful)."""
    if not source_name.startswith("CB_"):
        return None
    R = cuts["R_LXe_outer_cm"]
    r_fv = cuts["fv_r_max_cm"]
    if R <= 0 or r_fv <= 0 or r_fv >= R:
        return None
    return float(np.sqrt(1.0 - (r_fv / R) ** 2))


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


def _plot_2d_heatmap(
    ax: plt.Axes,
    df: pd.DataFrame,
    *,
    title: str,
    rlabel: str,
    zlabel: str,
) -> None:
    """Render a 2D histogram from a CSV with columns: rL, rR, zL, zR, count.
    z on X axis, r (or D) on Y axis (physics convention)."""
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
    im = ax.pcolormesh(z_edges, r_edges, M, cmap="viridis", shading="flat")
    plt.colorbar(im, ax=ax, label="count")
    ax.set_xlabel(zlabel)
    ax.set_ylabel(rlabel)
    ax.set_title(title)
    ax.set_xlim(z_edges[0], z_edges[-1])
    ax.set_ylim(r_edges[0], r_edges[-1])


def _draw_fv_box(ax: plt.Axes, cuts: dict[str, float], color: str = "#fff") -> None:
    """Draw the FV box on a (z, r) heatmap. z on X, r on Y."""
    z_lo, z_hi = cuts["fv_z_min_cm"], cuts["fv_z_max_cm"]
    r_max = cuts["fv_r_max_cm"]
    ax.plot([z_lo, z_hi, z_hi, z_lo, z_lo],
            [0,    0,    r_max, r_max, 0],
            color=color, linewidth=1.5, linestyle="--",
            label=f"FV: z∈[{z_lo:.0f},{z_hi:.0f}] r≤{r_max:.0f} cm")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.7)


def _plot_region_interaction_table(ax: plt.Axes, df: pd.DataFrame) -> None:
    """3×4 region × interaction-type counts as a heatmap of fractions."""
    pivot = df.pivot(index="region", columns="interaction", values="count")
    region_order = [r for r in ("TPC", "Skin", "Inert") if r in pivot.index]
    inter_order = [
        c for c in ("PHOTO", "COMPTON", "PAIR", "BELOW_THRESH") if c in pivot.columns
    ]
    pivot = pivot.loc[region_order, inter_order]
    M = pivot.values.astype(float)
    total = M.sum()
    F = M / total if total > 0 else M
    im = ax.imshow(F, cmap="Blues", aspect="auto")
    ax.set_xticks(range(len(inter_order)), inter_order, rotation=30, ha="right",
                   fontsize=8)
    ax.set_yticks(range(len(region_order)), region_order)
    ax.set_title("region × interaction (fraction of all deposits)")
    fmax = F.max() if F.size else 0.0
    for i in range(F.shape[0]):
        for j in range(F.shape[1]):
            ax.text(j, i, f"{F[i, j]:.2e}", ha="center", va="center",
                    color="white" if F[i, j] > 0.5 * fmax else "black",
                    fontsize=11, weight="bold")
    plt.colorbar(im, ax=ax, label="fraction")


# ---------------------------------------------------------------------------
# Family renderers
# ---------------------------------------------------------------------------
def render_acceptance_panel(
    input_dir: Path, title: str, source_name: str
) -> plt.Figure:
    """Cut 1 (h_u_sampled with geom-acceptance line) + Cut 2
    (first-interaction r vs z with FV box)."""
    cuts = _read_cuts(input_dir)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
    fig.suptitle(f"{title} — cuts 1 + 2 (acceptance)", fontsize=14)

    # Cut 1
    df = _read(input_dir, "cut1_h_u_sampled.csv")
    _plot_1d(axes[0], df,
              title="Cut 1: dN/du sampled (cos θ_inward)",
              xlabel="u",
              log_y=False, color="#46a")
    u_min = _u_min_geometric(cuts, source_name)
    if u_min is not None:
        axes[0].axvline(u_min, color="#a40", linestyle="--", linewidth=1.2,
                        label=f"u_min geom (FV cyl) = {u_min:.3f}")
        axes[0].legend(loc="upper left", fontsize=8)
    elif source_name.startswith(("CTH_", "CBH_")):
        axes[0].text(0.02, 0.98,
                     "endcap source: no single-line geom-acceptance",
                     transform=axes[0].transAxes,
                     ha="left", va="top", fontsize=8, color="#a40")

    # Cut 2
    df = _read(input_dir, "cut2_first_interaction_r_z.csv")
    _plot_2d_heatmap(axes[1], df,
                      title="Cut 2: first interaction (r, z), FV box overlaid",
                      rlabel="r (cm)", zlabel="z (cm)")
    _draw_fv_box(axes[1], cuts)
    return fig


def render_classification_panel(input_dir: Path, title: str) -> plt.Figure:
    """Cut 3 (Δz inclusive, n_visible, E_total) on top row;
    Cut 4 (ss_ec, ss_es, ss_r_z) on bottom row."""
    cuts = _read_cuts(input_dir)
    fig, axes = plt.subplots(2, 3, figsize=(20, 11), constrained_layout=True)
    fig.suptitle(f"{title} — cuts 3 + 4 (classification)", fontsize=14)

    # Cut 3 — top row
    df = _read(input_dir, "cut3_dz_inclusive.csv")
    _plot_1d(axes[0, 0], df,
              title="Cut 3: Δz inclusive (consecutive visible deposits)",
              xlabel="|Δz| (cm)", log_y=True, color="#46a")
    dz_thr_cm = cuts["Δz_threshold_mm"] / 10.0
    axes[0, 0].axvline(dz_thr_cm, color="#a40", linestyle="--", linewidth=1.2,
                        label=f"SS/MS = {cuts['Δz_threshold_mm']:.1f} mm")
    axes[0, 0].legend(loc="upper right", fontsize=8)

    df = _read(input_dir, "cut3_n_visible.csv")
    _plot_int_bar(axes[0, 1], df,
                   title="Cut 3: n_visible per event",
                   xlabel="n",
                   log_y=True, normalize=False, show_values=False,
                   xlim_max=15)

    df = _read(input_dir, "cut3_E_total.csv")
    _plot_1d(axes[0, 2], df,
              title="Cut 3: E_total (Σ visible cluster energies)",
              xlabel="E (MeV)", log_y=True, color="#a64")

    # Cut 4 — bottom row
    Q_MeV = cuts["Q_betabeta_keV"] / 1000.0
    HW_MeV = cuts["ROI_halfwidth_keV"] / 1000.0

    df = _read(input_dir, "cut4_ss_ec_pre_roi.csv")
    _plot_1d(axes[1, 0], df,
              title="Cut 4: SS true energy ec (pre-ROI)",
              xlabel="E (MeV)", log_y=True, color="#46a")
    _draw_roi_band(axes[1, 0], Q_MeV, HW_MeV)

    df = _read(input_dir, "cut4_ss_es_pre_roi.csv")
    _plot_1d(axes[1, 1], df,
              title="Cut 4: SS smeared energy es (pre-ROI; cut acts here)",
              xlabel="E (MeV)", log_y=True, color="#a64")
    _draw_roi_band(axes[1, 1], Q_MeV, HW_MeV)

    df = _read(input_dir, "cut4_ss_r_z.csv")
    _plot_2d_heatmap(axes[1, 2], df,
                      title="Cut 4: SS cluster (r, z), FV box overlaid",
                      rlabel="r (cm)", zlabel="z (cm)")
    _draw_fv_box(axes[1, 2], cuts)
    return fig


def _draw_roi_band(ax: plt.Axes, Q_MeV: float, HW_MeV: float) -> None:
    ax.axvspan(Q_MeV - HW_MeV, Q_MeV + HW_MeV, color="#fc6", alpha=0.35,
               label=f"ROI ± {HW_MeV*1000:.1f} keV")
    ax.axvline(Q_MeV, color="#a40", linestyle="--", linewidth=1.2,
               label=f"Q_ββ = {Q_MeV*1000:.0f} keV")
    ax.legend(loc="upper left", fontsize=8)


def render_diagnostic_panel(input_dir: Path, title: str) -> plt.Figure:
    """4 diagnostic plots in a 2×2 grid."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), constrained_layout=True)
    fig.suptitle(f"{title} — diagnostics", fontsize=14)

    df = _read(input_dir, "diag_interaction_type_freq.csv")
    _plot_categorical_bar(axes[0, 0], df,
                           title="interaction type frequency (per stack row)",
                           xlabel="interaction", log_y=True, normalize=True)

    df = _read(input_dir, "diag_path_length_LXe.csv")
    _plot_1d(axes[0, 1], df,
              title="per-photon path length in LXe",
              xlabel="path length (cm)", log_y=True, color="#684")

    df = _read(input_dir, "diag_region_interaction.csv")
    _plot_region_interaction_table(axes[1, 0], df)

    df = _read(input_dir, "diag_cluster_Ec.csv")
    _plot_1d(axes[1, 1], df,
              title="per-cluster energy Ec",
              xlabel="E (MeV)", log_y=True, color="#a64")

    return fig


# ---------------------------------------------------------------------------
# Per-histogram (--separate) renderer — auto-dispatch by CSV column shape
# ---------------------------------------------------------------------------
def _save_one(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=150)
    print(f"  wrote {path}")


def render_individual(input_dir: Path, out: Path) -> None:
    for csv in (*ACCEPTANCE_CSVS, *CLASSIFICATION_CSVS, *DIAGNOSTIC_CSVS):
        df = _read(input_dir, csv)
        if df is None:
            continue
        stem = csv.removesuffix(".csv")
        cols = list(df.columns)

        fig, ax = plt.subplots(figsize=(7, 4.5))
        if cols == ["interaction", "count"]:
            _plot_categorical_bar(ax, df, title=stem, xlabel="interaction",
                                   log_y=True, normalize=True)
        elif cols == ["region", "interaction", "count"]:
            _plot_region_interaction_table(ax, df)
        elif cols[-1] == "count" and len(cols) == 2:
            _plot_int_bar(ax, df, title=stem, xlabel=cols[0], log_y=True)
        elif cols[-1] == "count" and len(cols) == 3:
            _plot_1d(ax, df, title=stem, xlabel=cols[0], log_y=True)
        elif cols[-1] == "count" and len(cols) == 5:
            fig.set_size_inches(8, 6)
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
    # Source name = the directory name with the leading "hist_" stripped.
    # Per-source dir is named <source> (e.g. CB_Bi214) directly under
    # the run dir. No prefix to strip.
    source_name = args.input_dir.name

    families = {
        "acceptance":     (ACCEPTANCE_CSVS,     "cuts_acceptance.png",
                            lambda d, t: render_acceptance_panel(d, t, source_name)),
        "classification": (CLASSIFICATION_CSVS, "cuts_classification.png",
                            render_classification_panel),
        "diagnostics":    (DIAGNOSTIC_CSVS,     "diagnostics.png",
                            render_diagnostic_panel),
    }
    todo = (
        ["acceptance", "classification", "diagnostics"]
        if args.family == "all" else [args.family]
    )

    for fam in todo:
        csvs, png_name, renderer = families[fam]
        if not _present(args.input_dir, csvs):
            print(f"  skip {fam} (CSVs missing in {args.input_dir})")
            continue
        fig = renderer(args.input_dir, title)
        _save_one(fig, out / png_name)

    if args.separate:
        render_individual(args.input_dir, out)

    if args.show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()
