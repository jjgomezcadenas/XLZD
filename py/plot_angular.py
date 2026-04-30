"""
Plot angular distributions of gammas exiting the LXe skin.

Reads the CSV output files from lz_cryo.py and produces plots of
dN/du × total_rate (gammas/yr per unit u) for each isotope, broken
down by source sector and by region (barrel/endcap).

Produces 4 plots:
  1. Bi214 barrel  — all barrel source sectors
  2. Bi214 endcap  — all endcap source sectors
  3. Tl208 barrel  — all barrel source sectors
  4. Tl208 endcap  — all endcap source sectors

Each plot also shows the total (sum over sectors) as a thick black line.

Usage:
  python plot_angular.py

Reads from ../data/lz_cryo_gammas.csv and ../data/lz_cryo_angular_*.csv
Saves plots to ../data/lz_cryo_angular_*.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data")

# Energy labels for plot titles
ENERGY = {'Bi214': '2.448 MeV', 'Tl208': '2.615 MeV'}
ISOTOPE_LABEL = {'Bi214': '²¹⁴Bi', 'Tl208': '²⁰⁸Tl'}


def load_gammas():
    """Load the gamma rates table."""
    path = os.path.join(DATA_DIR, "lz_cryo_gammas.csv")
    return pd.read_csv(path)


def load_angular_pdf(iso, region):
    """
    Load an aggregated angular PDF and its total rate.

    Returns (u_bins, pdf, total_rate).
    """
    path = os.path.join(DATA_DIR, f"lz_cryo_angular_{iso}_{region}.csv")
    with open(path, 'r') as f:
        header_lines = []
        for line in f:
            if line.startswith('#'):
                header_lines.append(line.strip())
            else:
                break

    total_rate = None
    for line in header_lines:
        if 'total_gammas_per_yr' in line:
            total_rate = float(line.split('=')[1].strip())

    df = pd.read_csv(path, comment='#')
    return df['u'].values, df['dNdu_pdf'].values, total_rate


def reconstruct_sector_spectra(gammas_df, iso, region):
    """
    Reconstruct per-sector absolute angular spectra by running the
    transport calculation from lz_cryo. We import the module to reuse
    its functions rather than duplicating code.

    Returns list of (source_name, u_bins, dNdu_absolute) tuples.
    """
    import sys
    sys.path.insert(0, SCRIPT_DIR)
    import lz_cryo as lc

    mu = lc.load_attenuation_coefficients()
    ocv = lc.compute_vessel(**lc.OCV_CFG)
    icv = lc.compute_vessel(**lc.ICV_CFG)
    sectors = lc.build_source_sectors(ocv, icv, mu)
    u_bins = lc.U_BINS

    spectra = []
    for sec in sectors:
        if sec['region'] != region:
            continue

        sp = sec[f'{iso}_sp']
        mass = sec['mass_kg']
        activity_mBq = mass * sp

        if iso == 'Bi214':
            gamma_br = lc.BR_BI214_GAMMA
        else:
            gamma_br = lc.BR_TL208_FROM_CHAIN * lc.BR_TL208_GAMMA

        gammas_per_yr = activity_mBq * 1e-3 * gamma_br * lc.SEC_PER_YEAR

        layers_info = sec['layers_fn'](iso)
        slabs_thru_skin = layers_info['to_skin'] + [layers_info['skin']]

        if sec['src_type'] == 'volume':
            dNdu = lc.angular_spectrum_volume(
                layers_info['src_mu'], layers_info['src_t'],
                slabs_thru_skin, u_bins)
        else:
            dNdu = lc.angular_spectrum_surface(slabs_thru_skin, u_bins)

        # Absolute spectrum: gammas/yr per unit u
        dNdu_abs = dNdu * gammas_per_yr
        spectra.append((sec['name'], u_bins, dNdu_abs))

    return spectra


def plot_one(iso, region, ax):
    """
    Plot all sector angular spectra for one (isotope, region) combination.
    """
    spectra = reconstruct_sector_spectra(None, iso, region)

    total = np.zeros_like(spectra[0][1])
    for name, u, dNdu in spectra:
        ax.plot(u, dNdu, label=name, linewidth=1.2)
        total += dNdu

    ax.plot(spectra[0][1], total, 'k-', linewidth=2.0, label='Total')

    ax.set_xlabel(r'$u = \cos\theta$', fontsize=12)
    ax.set_ylabel(r'$dN/du$ (gammas/yr per unit $u$)', fontsize=11)
    ax.set_title(f'{ISOTOPE_LABEL[iso]} ({ENERGY[iso]}) — {region}', fontsize=13)
    ax.legend(fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    # Set y-axis lower limit to avoid empty log plots
    yvals = total[total > 0]
    if len(yvals) > 0:
        ax.set_ylim(bottom=yvals.min() * 0.1)


def main():
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    configs = [
        ('Bi214', 'barrel',  axes[0, 0]),
        ('Bi214', 'endcap',  axes[0, 1]),
        ('Tl208', 'barrel',  axes[1, 0]),
        ('Tl208', 'endcap',  axes[1, 1]),
    ]

    for iso, region, ax in configs:
        plot_one(iso, region, ax)

    fig.suptitle('LZ Cryostat: Angular Distribution of Gammas Exiting LXe Skin',
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    outpath = os.path.join(DATA_DIR, "lz_cryo_angular_spectra.png")
    fig.savefig(outpath, dpi=150)
    print(f"Saved: {outpath}")
    plt.show()


if __name__ == "__main__":
    main()
