"""
LZ Cryostat Model for 0νββ Background Simulations
===================================================

This module computes the gamma-ray flux from radioactive decays in the
LZ cryostat that reaches the active LXe volume, for use in 0νββ background
studies.

Physics context
---------------
The dominant radiogenic backgrounds at Q_ββ = 2458 keV are:
  - ²¹⁴Bi (from ²³⁸U late chain): emits a 2.448 MeV γ with 1.55% BR
  - ²⁰⁸Tl (from ²³²Th late chain): emits a 2.615 MeV γ with ~100% BR,
    but ²⁰⁸Tl itself is produced in only 35.9% of ²¹²Bi decays

The measured specific activities are for the full late chains (²³⁸U_late
and ²³²Th_late). To get the relevant gamma production rates:
  - ²¹⁴Bi γ rate = ²³⁸U_late activity × 0.0155 (BR of 2.448 MeV line)
  - ²⁰⁸Tl γ rate = ²³²Th_late activity × 0.359 (BR to ²⁰⁸Tl) × 1.0 (γ BR)

Detector model
--------------
Two concentric cylinders (OCV and ICV) with ancillary layers:

  From outside in (at mid-height):
    OCV wall        Ti cylinder, R=91.5 cm, H=304 cm, t solved for mass
    vacuum gap      ~8 cm, transparent to gammas
    MLI             thin surface source on ICV outer wall
    ICV wall        Ti cylinder, R=83.0 cm, H=259 cm, t solved for mass
    Teflon liner    thin surface source on ICV inner wall
    LXe skin        ~6 cm liquid xenon, attenuator
    Field cage      outer wall of active volume

  Wall thicknesses are solved so that the geometric mass matches the
  TDR values (950 kg ICV, 1115 kg OCV). Head thicknesses are fixed
  to nominal values slightly above ASME code minimums.

  An additional 345 kg of Ti (support structure, legs, bolts) exists
  in the TDR mass budget but is far from the fiducial volume and is
  NOT included in this model. This is a conservative omission.

Gamma transport
---------------
Each source sector (e.g., "ICV shell", "OCV heads", "MLI") has an
associated layer stack: the ordered list of material layers between
the source and the exit of the LXe skin. Gammas are emitted
isotropically; only the inward hemisphere (factor 1/2) contributes.

For a gamma traveling at polar angle θ from the inward normal
(u = cos θ), the path through layer i of thickness tᵢ is tᵢ/u.
The survival probability for a single ray direction is:

    P(u) = exp( -Σᵢ μᵢ tᵢ / u )

For volume sources (Ti walls), gammas are born uniformly at depth
x ∈ [0, t_src] from the inner face. The depth integral is analytic:

    ∫₀ᵗ exp(-μ x/u) dx = (u/μ)(1 - exp(-μt/u))

The angular spectrum of gammas exiting the last layer is:

    dN/du(u) = (rate/2t) × (u/μ_src)(1 - exp(-μ_src t_src/u))
                × exp(-Σ_slabs μᵢ tᵢ / u)

For surface sources, the depth integral is absent:

    dN/du(u) = (rate/2) × exp(-Σ_slabs μᵢ tᵢ / u)

The integral over u ∈ [0,1] gives the total gamma rate at each stage.
The full dN/du distribution is stored for MC sampling.

We compute three stages:
  1. Gammas PRODUCED (activity × BR × time)
  2. Gammas ENTERING the LXe skin (after traversing all layers except skin)
  3. Gammas EXITING the skin (after traversing all layers including skin)

The dN/du at stage 3 is the MC input: sample u from this PDF, then
  θ = arccos(u), φ = uniform[0, 2π], entry point uniform on skin surface.

Output files (all CSV in ../data/):
  lz_cryo_geometry.csv       Model element dimensions and masses
  lz_cryo_activity.csv       Radioactivity per source sector
  lz_cryo_gammas.csv         Gamma rates at each stage per sector
  lz_cryo_angular_*.csv      dN/du PDFs for MC sampling (4 files)

Sources:
  LZ TDR (arXiv:1703.09144), LZ Titanium Paper (arXiv:1702.02646),
  LZ Instrument Paper (arXiv:1910.09124), LZ Radioactivity Paper
  (arXiv:2006.02506), LZ_detector_summary.md.

All internal units: cm, kg, mBq, seconds. Output rates in gammas/year.
"""

import os
import numpy as np
import pandas as pd
# No scipy needed — we use np.trapz on uniform bins

# ===========================================================================
# Paths
# ===========================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data")

# ===========================================================================
# Physical constants
# ===========================================================================
RHO_TI = 4.510       # g/cm³, Grade 1 titanium
RHO_LXE = 2.953      # g/cm³, liquid xenon at ~170 K
SEC_PER_YEAR = 3.1557e7  # seconds in a Julian year

# Branching ratios for gamma production from chain activities:
#   ²¹⁴Bi 2.448 MeV γ: 1.55% of ²¹⁴Bi decays
#   ²⁰⁸Tl 2.615 MeV γ: 35.9% of ²¹²Bi decays go to ²⁰⁸Tl,
#     and ~100% of ²⁰⁸Tl decays emit the 2.615 MeV γ.
#   So from ²³²Th_late chain activity: γ rate = activity × 0.359
BR_BI214_GAMMA = 0.0155    # 2.448 MeV γ per ²¹⁴Bi decay
BR_TL208_FROM_CHAIN = 0.359  # ²⁰⁸Tl produced per ²³²Th_late decay
BR_TL208_GAMMA = 1.0      # 2.615 MeV γ per ²⁰⁸Tl decay

# Gamma energies (MeV)
E_BI214 = 2.448
E_TL208 = 2.615


# ===========================================================================
# NIST attenuation data
# ===========================================================================
def load_nist_mu_rho(filepath):
    """
    Load a NIST XCOM attenuation table and return (energies_MeV, mu_over_rho).

    Reads the 'Tot. w/ Coherent' column, which is the total mass attenuation
    coefficient including coherent scattering, in cm²/g.

    The NIST files have a 2-line header, then whitespace-separated columns.
    For Ti: columns are (Energy, Tot_w_Coherent).
    For LXe: columns are (Energy, Coherent, Incoherent, Photoelectric,
             Nuclear_PP, Electron_PP, Tot_w_Coherent, Tot_wo_Coherent).

    Parameters
    ----------
    filepath : str
        Path to the NIST CSV file.

    Returns
    -------
    energies : ndarray
        Photon energies in MeV.
    mu_over_rho : ndarray
        Total mass attenuation coefficient (with coherent) in cm²/g.
    """
    energies = []
    mu_rho = []
    with open(filepath, 'r') as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        try:
            e = float(parts[0])
            # Tot. w/ Coherent is the last column for Ti (2-col file)
            # or the second-to-last column for LXe (8-col file)
            if len(parts) == 2:
                val = float(parts[1])
            elif len(parts) >= 7:
                val = float(parts[6])  # "Tot. w/ Coherent" is column index 6
            else:
                continue
            energies.append(e)
            mu_rho.append(val)
        except (ValueError, IndexError):
            continue

    return np.array(energies), np.array(mu_rho)


def interpolate_mu(energies, mu_rho, E_target, rho):
    """
    Interpolate mass attenuation coefficient and convert to linear.

    Uses log-log linear interpolation, standard for attenuation data.

    Parameters
    ----------
    energies : ndarray
        NIST table energies (MeV).
    mu_rho : ndarray
        Mass attenuation coefficients (cm²/g).
    E_target : float
        Target energy (MeV).
    rho : float
        Material density (g/cm³).

    Returns
    -------
    mu_linear : float
        Linear attenuation coefficient (cm⁻¹) at E_target.
    """
    log_e = np.log(energies)
    log_mu = np.log(mu_rho)
    log_target = np.log(E_target)
    mu_rho_interp = np.exp(np.interp(log_target, log_e, log_mu))
    return mu_rho_interp * rho


def load_attenuation_coefficients():
    """
    Load NIST tables for Ti and LXe, compute μ_linear at 2.448 and 2.615 MeV.

    Returns
    -------
    mu : dict
        Keys like 'Ti_Bi214', 'Ti_Tl208', 'LXe_Bi214', 'LXe_Tl208'.
        Values are linear attenuation coefficients in cm⁻¹.
    """
    ti_e, ti_mu = load_nist_mu_rho(os.path.join(DATA_DIR, "nist_ti.csv"))
    lxe_e, lxe_mu = load_nist_mu_rho(os.path.join(DATA_DIR, "nist_lxe.csv"))

    mu = {
        'Ti_Bi214':  interpolate_mu(ti_e, ti_mu, E_BI214, RHO_TI),
        'Ti_Tl208':  interpolate_mu(ti_e, ti_mu, E_TL208, RHO_TI),
        'LXe_Bi214': interpolate_mu(lxe_e, lxe_mu, E_BI214, RHO_LXE),
        'LXe_Tl208': interpolate_mu(lxe_e, lxe_mu, E_TL208, RHO_LXE),
    }
    return mu


# ===========================================================================
# Geometry helpers
# ===========================================================================
def ellipsoidal_head_area(R, aspect_ratio):
    """
    Inside surface area of one ellipsoidal head (half oblate spheroid).

    For an n:1 ellipsoidal head on a cylinder of inside radius R,
    the head depth is R/n. The surface area is half of the full
    oblate spheroid with semi-axes a=R, b=R/n.

    Parameters
    ----------
    R : float
        Inside radius of the cylinder (cm).
    aspect_ratio : float
        n in n:1 (e.g., 2 for 2:1 head where depth = R/2).

    Returns
    -------
    area : float
        Surface area in cm².
    """
    a = R
    b = R / aspect_ratio
    e = np.sqrt(1.0 - (b / a) ** 2)
    S_full = 2 * np.pi * a**2 + np.pi * (b**2 / e) * np.log((1 + e) / (1 - e))
    return S_full / 2.0


def ellipsoidal_head_depth(R, aspect_ratio):
    """Depth (height) of an n:1 ellipsoidal head: R/n, in cm."""
    return R / aspect_ratio


def cylinder_lateral_area(R, H):
    """Lateral surface area of a cylinder: 2πRH, in cm²."""
    return 2 * np.pi * R * H


def mass_from_area(area_cm2, thickness_mm, rho=RHO_TI):
    """
    Mass of a shell element.

    Parameters
    ----------
    area_cm2 : float
        Surface area (cm²).
    thickness_mm : float
        Wall thickness (mm).
    rho : float
        Material density (g/cm³). Default: titanium.

    Returns
    -------
    mass_kg : float
    """
    return area_cm2 * (thickness_mm / 10.0) * rho / 1000.0


# ===========================================================================
# Vessel computation
# ===========================================================================
def compute_vessel(label, R, H_total, top_ar, bot_ar,
                   t_top_mm, t_bot_mm, M_target):
    """
    Compute vessel mass breakdown with fixed geometry and head thicknesses,
    solving for the cylindrical wall thickness to match M_target.

    The vessel is modeled as a cylinder of inside radius R and total height
    H_total, with ellipsoidal heads at top and bottom. The head depths are
    subtracted from H_total to get the cylindrical section height. Head
    masses are computed from their area and fixed thickness. The wall
    thickness is then solved analytically:

        t_wall = (M_target - M_heads) × 1000 / (A_cyl × ρ_Ti)

    This absorbs the real-world mass from flanges, ports, and stiffening
    rings into a slightly thicker effective wall — the correct approach for
    a model that needs to scale to larger detectors.

    Parameters
    ----------
    label : str
        Vessel name (e.g., "ICV", "OCV").
    R : float
        Inside radius (cm).
    H_total : float
        Total vessel height including heads (cm).
    top_ar, bot_ar : float
        Aspect ratios of top and bottom ellipsoidal heads (e.g., 2 for 2:1).
    t_top_mm, t_bot_mm : float
        Nominal thicknesses of top and bottom heads (mm).
    M_target : float
        Target total vessel mass (kg) from TDR.

    Returns
    -------
    v : dict
        Complete vessel description including solved wall thickness,
        all areas, masses, and dimensions.
    """
    h_top = ellipsoidal_head_depth(R, top_ar)
    h_bot = ellipsoidal_head_depth(R, bot_ar)
    H_cyl = H_total - h_top - h_bot
    assert H_cyl > 0, f"{label}: H_cyl={H_cyl:.1f} cm is negative"

    A_top = ellipsoidal_head_area(R, top_ar)
    A_bot = ellipsoidal_head_area(R, bot_ar)
    A_cyl = cylinder_lateral_area(R, H_cyl)

    m_top = mass_from_area(A_top, t_top_mm)
    m_bot = mass_from_area(A_bot, t_bot_mm)
    m_heads = m_top + m_bot

    m_shell_needed = M_target - m_heads
    assert m_shell_needed > 0, f"{label}: heads ({m_heads:.1f} kg) exceed target ({M_target:.1f} kg)"

    # Solve: m_shell = A_cyl × (t_wall/10) × ρ / 1000
    t_wall_mm = m_shell_needed * 1000.0 / (A_cyl * RHO_TI) * 10.0

    return {
        'label': label, 'R': R, 'H_total': H_total, 'H_cyl': H_cyl,
        'h_top': h_top, 'h_bot': h_bot, 'top_ar': top_ar, 'bot_ar': bot_ar,
        'A_top': A_top, 'A_bot': A_bot, 'A_cyl': A_cyl,
        't_top_mm': t_top_mm, 't_bot_mm': t_bot_mm, 't_wall_mm': t_wall_mm,
        'm_top': m_top, 'm_bot': m_bot, 'm_heads': m_heads,
        'm_shell': m_shell_needed, 'm_total': M_target,
    }


# ===========================================================================
# Gamma transport through layer stacks
# ===========================================================================
def angular_spectrum_volume(mu_src, t_src_cm, slab_layers, u_bins):
    """
    Angular spectrum of gammas escaping inward from a uniform volume source
    through a sequence of slab layers.

    The source is distributed uniformly in a slab of thickness t_src with
    linear attenuation coefficient mu_src. Gammas are emitted isotropically;
    we consider only the inward hemisphere (factor 1/2 implicit in the
    normalization — the returned spectrum integrates to the escape fraction,
    i.e., the fraction of all emitted gammas that exit the last layer heading
    inward).

    For a gamma at angle θ from the inward normal (u = cos θ), born at
    depth x from the inner face of the source slab:

        P(x, u) = exp(-(μ_src × x + Σ μᵢ tᵢ) / u)

    Integrating over birth depth x analytically:

        ∫₀ᵗ exp(-μ_src x/u) dx = (u/μ_src)(1 - exp(-μ_src t/u))

    So the angular spectrum per unit u is:

        dN/du = (1 / 2t_src) × (u / μ_src) × (1 - exp(-μ_src t_src / u))
                × exp(-τ_slabs / u)

    where τ_slabs = Σ μᵢ tᵢ is the total optical depth of downstream slabs.

    Parameters
    ----------
    mu_src : float
        Linear attenuation coefficient of the source material (cm⁻¹).
    t_src_cm : float
        Thickness of the source slab (cm).
    slab_layers : list of dict
        Downstream layers, each {"mu": float, "t_cm": float, "label": str}.
        These are traversed after the gamma escapes the source slab.
    u_bins : ndarray
        Bin centers in u = cos θ ∈ (0, 1].

    Returns
    -------
    dNdu : ndarray
        Angular spectrum (same length as u_bins). Units: fraction per unit u.
        Multiply by the gamma production rate to get gammas/yr per unit u.
    """
    tau_slabs = sum(layer["mu"] * layer["t_cm"] for layer in slab_layers)
    dNdu = np.zeros_like(u_bins)

    for i, u in enumerate(u_bins):
        # Analytic integral over source depth x
        depth_integral = (u / mu_src) * (1.0 - np.exp(-mu_src * t_src_cm / u))
        # Attenuation through downstream slabs
        slab_atten = np.exp(-tau_slabs / u)
        # Factor 1/(2 t_src): 1/2 for inward hemisphere, 1/t_src for averaging
        # over birth depth
        dNdu[i] = (1.0 / (2.0 * t_src_cm)) * depth_integral * slab_atten

    return dNdu


def angular_spectrum_surface(slab_layers, u_bins):
    """
    Angular spectrum of gammas from a thin surface source through slab layers.

    The source is infinitely thin (no self-shielding). Gammas are emitted
    isotropically; we take the inward hemisphere. For direction u = cos θ:

        dN/du = (1/2) × exp(-τ_slabs / u)

    Parameters
    ----------
    slab_layers : list of dict
        Layers to traverse, each {"mu": float, "t_cm": float, "label": str}.
    u_bins : ndarray
        Bin centers in u = cos θ ∈ (0, 1].

    Returns
    -------
    dNdu : ndarray
        Angular spectrum (fraction per unit u).
    """
    tau_slabs = sum(layer["mu"] * layer["t_cm"] for layer in slab_layers)
    dNdu = np.zeros_like(u_bins)
    for i, u in enumerate(u_bins):
        dNdu[i] = 0.5 * np.exp(-tau_slabs / u)
    return dNdu


def integrate_spectrum(dNdu, u_bins):
    """
    Integrate dN/du over u to get total fraction, and compute ⟨cos θ⟩.

    Uses numpy trapezoid rule directly on the binned data. This is both
    fast and well-behaved, unlike scipy.quad on interpolated binned data
    which struggles near u → 0 where exp(-τ/u) drops steeply.

    With 100 uniform bins in u ∈ [0.005, 0.995], the trapezoid rule is
    accurate to better than 0.1% for the smooth spectra encountered here.

    Parameters
    ----------
    dNdu : ndarray
        Angular spectrum values at u_bins.
    u_bins : ndarray
        Bin centers.

    Returns
    -------
    total : float
        ∫ dN/du du — the total escape/transmission fraction.
    mean_u : float
        ⟨u⟩ = ⟨cos θ⟩ = ∫ u dN/du du / ∫ dN/du du.
    """
    total = np.trapz(dNdu, u_bins)
    moment = np.trapz(u_bins * dNdu, u_bins)
    mean_u = moment / total if total > 0 else 0.0
    return total, mean_u


# ===========================================================================
# Configuration: geometry, materials, radioactivity
# ===========================================================================

# --- OCV: outer cryostat vessel (straight cylinder) ---
# Wall thicknesses at TDR code minimum values (real plate thicknesses).
# This gives physical self-shielding consistent with the actual vessel.
# The extra mass (2590 - 1430 = 1160 kg from flanges, ports, support)
# is assigned to the bottom endcap source (heavily shielded by dome LXe).
OCV_CFG = {
    'label': 'OCV', 'R': 91.5, 'H_total': 304.0,
    'top_ar': 2, 'bot_ar': 2,
    't_top_mm': 9.0, 't_bot_mm': 15.0,
    'M_target': 778.6,   # geometric mass at real wall thickness (7 mm)
}

# --- ICV: inner cryostat vessel (modeled as cylinder at top diameter) ---
ICV_CFG = {
    'label': 'ICV', 'R': 83.0, 'H_total': 259.0,
    'top_ar': 2, 'bot_ar': 3,
    't_top_mm': 8.0, 't_bot_mm': 12.0,
    'M_target': 651.1,   # geometric mass at real wall thickness (9 mm)
}

# Extra Ti mass (flanges, ports, stiffeners, support) → bottom endcap
EXTRA_TI_MASS = 2590.0 - 778.6 - 651.1  # ≈ 1160 kg

# --- LXe skin ---
SKIN_THICKNESS_CM = 6.0    # average of 4 cm (top) to 8 cm (cathode)

# --- TPC geometry (for reference) ---
TPC_INNER_R = 72.8   # cm, 1456 mm diameter / 2
FC_THICKNESS = 1.5    # cm, 15 mm

# --- Specific activities (mBq/kg) ---
# From LZ bb0nu paper (arXiv:2004.XXXXX), Table I
# ²³⁸U-late = ²²⁶Ra and after;  ²³²Th-late = ²²⁴Ra and after
TI_BI214_SP = 0.08     # mBq/kg, Ti cryostat vessel ²³⁸U-late
TI_TL208_SP = 0.22     # mBq/kg, Ti cryostat vessel ²³²Th-late

# Cryostat insulation (MLI) — from bbonu paper Table I
ANCILLARY = {
    'MLI':  {'mass': 13.8, 'Bi214_sp': 11.1,  'Tl208_sp': 7.79},
}

# u bins for angular spectra
N_UBINS = 100
U_BINS = np.linspace(0.005, 0.995, N_UBINS)  # avoid u=0 (divergent path)


# ===========================================================================
# Source sector definitions
# ===========================================================================
def build_source_sectors(ocv, icv, mu):
    """
    Define all source sectors with their layer stacks for both isotopes.

    Each source sector is a dict with:
        name     : human-readable label
        region   : 'barrel' or 'endcap'
        src_type : 'volume' or 'surface' or 'circumference'
        mass_kg  : total mass of this source
        Bi214_sp : specific activity for ²¹⁴Bi chain (mBq/kg)
        Tl208_sp : specific activity for ²⁰⁸Tl chain (mBq/kg)
        layers_Bi214 : layer stack for 2.448 MeV gammas (list of dicts)
        layers_Tl208 : layer stack for 2.615 MeV gammas (list of dicts)

    Layer stacks are ordered from source outward (toward LXe). The first
    entry is the source itself (type "volume" or "surface"), followed by
    any intervening material slabs, ending with the LXe skin.

    For volume sources, the first layer's mu and t_cm define the self-
    shielding integral. For surface sources, all layers are slabs.

    The stacks are built for two stages:
        layers_to_skin  : layers up to (not including) the skin
        layers_thru_skin: layers up to and including the skin

    Parameters
    ----------
    ocv, icv : dict
        Vessel descriptions from compute_vessel().
    mu : dict
        Attenuation coefficients from load_attenuation_coefficients().

    Returns
    -------
    sectors : list of dict
        All source sectors with complete layer stacks.
    """
    t_ocv = ocv['t_wall_mm'] / 10.0  # cm
    t_icv = icv['t_wall_mm'] / 10.0
    t_ocv_top = ocv['t_top_mm'] / 10.0
    t_ocv_bot = ocv['t_bot_mm'] / 10.0
    t_icv_top = icv['t_top_mm'] / 10.0
    t_icv_bot = icv['t_bot_mm'] / 10.0
    t_skin = SKIN_THICKNESS_CM

    # Average head thickness for "endcap" sources (top+bottom combined)
    t_ocv_head_avg = (ocv['m_top'] * t_ocv_top + ocv['m_bot'] * t_ocv_bot) / ocv['m_heads']
    t_icv_head_avg = (icv['m_top'] * t_icv_top + icv['m_bot'] * t_icv_bot) / icv['m_heads']

    def make_layers(material, t_cm, isotope):
        """Helper to create a slab layer dict for a given isotope."""
        mu_key = f"{material}_{isotope}"
        return {"label": material, "mu": mu[mu_key], "t_cm": t_cm}

    sectors = []

    # --- OCV shell (barrel, volume source) ---
    # Path: OCV wall (self-shielding) → vacuum gap (transparent) → ICV wall
    # LXe skin is NOT included — it will be tracked in the Julia MC.
    for iso in ['Bi214', 'Tl208']:
        pass  # We build both isotope layer stacks in the sector dict

    # "skin" layer is now a zero-thickness dummy (no LXe attenuation)
    def _no_skin(iso):
        return {"label": "none", "mu": 0.0, "t_cm": 0.0}

    def ocv_shell_layers(iso):
        return {
            'src_mu': mu[f'Ti_{iso}'], 'src_t': t_ocv,
            'to_skin': [make_layers('Ti', t_icv, iso)],
            'skin':     _no_skin(iso),
        }

    def ocv_heads_layers(iso):
        return {
            'src_mu': mu[f'Ti_{iso}'], 'src_t': t_ocv_head_avg,
            'to_skin': [make_layers('Ti', t_icv_head_avg, iso)],
            'skin':     _no_skin(iso),
        }

    def icv_shell_layers(iso):
        return {
            'src_mu': mu[f'Ti_{iso}'], 'src_t': t_icv,
            'to_skin': [],
            'skin':     _no_skin(iso),
        }

    def icv_heads_layers(iso):
        return {
            'src_mu': mu[f'Ti_{iso}'], 'src_t': t_icv_head_avg,
            'to_skin': [],
            'skin':     _no_skin(iso),
        }

    def mli_layers(iso):
        """MLI sits on ICV outer surface → must traverse ICV wall only.
        LXe skin will be tracked in the Julia MC."""
        return {
            'to_skin': [make_layers('Ti', t_icv, iso)],
            'skin':     _no_skin(iso),
        }

    # Ti vessel sectors at real wall thicknesses.
    # Extra mass (flanges, support, ports) assigned to bottom endcap.
    sectors = [
        {
            'name': 'OCV shell', 'region': 'barrel', 'src_type': 'volume',
            'mass_kg': ocv['m_shell'],
            'Bi214_sp': TI_BI214_SP, 'Tl208_sp': TI_TL208_SP,
            'layers_fn': ocv_shell_layers,
        },
        {
            'name': 'OCV heads', 'region': 'endcap', 'src_type': 'volume',
            'mass_kg': ocv['m_heads'],
            'Bi214_sp': TI_BI214_SP, 'Tl208_sp': TI_TL208_SP,
            'layers_fn': ocv_heads_layers,
        },
        {
            'name': 'ICV shell', 'region': 'barrel', 'src_type': 'volume',
            'mass_kg': icv['m_shell'],
            'Bi214_sp': TI_BI214_SP, 'Tl208_sp': TI_TL208_SP,
            'layers_fn': icv_shell_layers,
        },
        {
            'name': 'ICV heads', 'region': 'endcap', 'src_type': 'volume',
            'mass_kg': icv['m_heads'],
            'Bi214_sp': TI_BI214_SP, 'Tl208_sp': TI_TL208_SP,
            'layers_fn': icv_heads_layers,
        },
        {
            # Extra Ti mass (flanges, ports, support structure)
            # Placed at bottom endcap — heavily shielded by dome LXe
            'name': 'Extra Ti (flanges etc)', 'region': 'endcap', 'src_type': 'surface',
            'mass_kg': EXTRA_TI_MASS,
            'Bi214_sp': TI_BI214_SP, 'Tl208_sp': TI_TL208_SP,
            'layers_fn': icv_heads_layers,  # same attenuation as ICV heads
        },
        {
            'name': 'Cryostat insulation', 'region': 'barrel', 'src_type': 'surface',
            'mass_kg': ANCILLARY['MLI']['mass'],
            'Bi214_sp': ANCILLARY['MLI']['Bi214_sp'],
            'Tl208_sp': ANCILLARY['MLI']['Tl208_sp'],
            'layers_fn': mli_layers,
        },
    ]

    return sectors


# ===========================================================================
# Gamma rate and transport computation
# ===========================================================================
def compute_gamma_tables(sectors, u_bins):
    """
    For each source sector and each isotope, compute:
      1. Gammas produced per year
      2. Angular spectrum and total rate entering the skin
      3. Angular spectrum and total rate exiting the skin

    The angular spectra at the skin exit are the MC inputs.

    Parameters
    ----------
    sectors : list of dict
        Source sectors from build_source_sectors().
    u_bins : ndarray
        Bin centers for angular spectra.

    Returns
    -------
    results : list of dict
        One entry per (sector, isotope) with all rates and spectra.
    """
    results = []

    for sec in sectors:
        for iso, iso_label in [('Bi214', 'Bi214'), ('Tl208', 'Tl208')]:
            # --- Activity and gamma production ---
            sp = sec[f'{iso}_sp']
            mass = sec['mass_kg']
            activity_mBq = mass * sp  # mBq

            if iso == 'Bi214':
                gamma_br = BR_BI214_GAMMA
            else:
                gamma_br = BR_TL208_FROM_CHAIN * BR_TL208_GAMMA

            # Gammas per second = activity (mBq) × 1e-3 (Bq/mBq) × BR
            gammas_per_sec = activity_mBq * 1e-3 * gamma_br
            gammas_per_yr = gammas_per_sec * SEC_PER_YEAR

            # --- Layer stacks ---
            layers_info = sec['layers_fn'](iso)

            # Layers up to skin entrance (not including skin)
            slabs_to_skin = layers_info['to_skin']
            # Layers up to skin exit (including skin)
            slabs_thru_skin = slabs_to_skin + [layers_info['skin']]

            # --- Angular spectra ---
            if sec['src_type'] == 'volume':
                mu_src = layers_info['src_mu']
                t_src = layers_info['src_t']
                dNdu_entering = angular_spectrum_volume(mu_src, t_src,
                                                       slabs_to_skin, u_bins)
                dNdu_exiting = angular_spectrum_volume(mu_src, t_src,
                                                      slabs_thru_skin, u_bins)
            else:  # surface or circumference
                dNdu_entering = angular_spectrum_surface(slabs_to_skin, u_bins)
                dNdu_exiting = angular_spectrum_surface(slabs_thru_skin, u_bins)

            # --- Integrate spectra ---
            frac_entering, mean_u_entering = integrate_spectrum(dNdu_entering, u_bins)
            frac_exiting, mean_u_exiting = integrate_spectrum(dNdu_exiting, u_bins)

            gammas_entering = gammas_per_yr * frac_entering
            gammas_exiting = gammas_per_yr * frac_exiting

            results.append({
                'source': sec['name'],
                'region': sec['region'],
                'src_type': sec['src_type'],
                'isotope': iso_label,
                'mass_kg': mass,
                'activity_mBq': activity_mBq,
                'gamma_br': gamma_br,
                'gammas_produced_per_yr': gammas_per_yr,
                'frac_entering_skin': frac_entering,
                'gammas_entering_skin_per_yr': gammas_entering,
                'mean_u_entering': mean_u_entering,
                'frac_exiting_skin': frac_exiting,
                'gammas_exiting_skin_per_yr': gammas_exiting,
                'mean_u_exiting': mean_u_exiting,
                'dNdu_exiting': dNdu_exiting,  # for MC sampling
            })

    return results


# ===========================================================================
# Output: tables and CSV files
# ===========================================================================
def save_geometry_csv(ocv, icv, filepath):
    """Save model geometry to CSV."""
    rows = []
    for v in [ocv, icv]:
        rows.append({
            'element': v['label'],
            'R_cm': v['R'],
            'H_total_cm': v['H_total'],
            'H_cyl_cm': v['H_cyl'],
            't_wall_mm': round(v['t_wall_mm'], 2),
            't_top_mm': v['t_top_mm'],
            't_bot_mm': v['t_bot_mm'],
            'top_ar': v['top_ar'],
            'bot_ar': v['bot_ar'],
            'mass_shell_kg': round(v['m_shell'], 1),
            'mass_heads_kg': round(v['m_heads'], 1),
            'mass_total_kg': round(v['m_total'], 1),
        })
    rows.append({
        'element': 'LXe_skin',
        'R_cm': round(icv['R'] - icv['t_wall_mm'] / 10.0, 1),
        'H_total_cm': icv['H_cyl'],
        'H_cyl_cm': icv['H_cyl'],
        't_wall_mm': SKIN_THICKNESS_CM * 10,
        't_top_mm': 0, 't_bot_mm': 0,
        'top_ar': 0, 'bot_ar': 0,
        'mass_shell_kg': 0, 'mass_heads_kg': 0, 'mass_total_kg': 0,
    })
    df = pd.DataFrame(rows)
    df.to_csv(filepath, index=False)
    return df


def save_activity_csv(sectors, filepath):
    """Save radioactivity assignment to CSV."""
    rows = []
    for sec in sectors:
        mass = sec['mass_kg']
        rows.append({
            'source': sec['name'],
            'region': sec['region'],
            'source_type': sec['src_type'],
            'mass_kg': round(mass, 1),
            'Bi214_mBq_per_kg': sec['Bi214_sp'],
            'Bi214_total_mBq': round(mass * sec['Bi214_sp'], 2),
            'Tl208_mBq_per_kg': sec['Tl208_sp'],
            'Tl208_total_mBq': round(mass * sec['Tl208_sp'], 2),
        })
    df = pd.DataFrame(rows)
    df.to_csv(filepath, index=False)
    return df


def save_gammas_csv(results, filepath):
    """Save gamma rate tables to CSV."""
    rows = []
    for r in results:
        rows.append({
            'source': r['source'],
            'region': r['region'],
            'isotope': r['isotope'],
            'gammas_produced_per_yr': f"{r['gammas_produced_per_yr']:.2e}",
            'gammas_entering_skin_per_yr': f"{r['gammas_entering_skin_per_yr']:.2e}",
            'gammas_exiting_skin_per_yr': f"{r['gammas_exiting_skin_per_yr']:.2e}",
            'mean_cos_theta_exiting': f"{r['mean_u_exiting']:.4f}",
        })
    df = pd.DataFrame(rows)
    df.to_csv(filepath, index=False)
    return df


def save_angular_csv(results, u_bins, iso, region, filepath):
    """
    Save the aggregated dN/du PDF for one (isotope, region) combination.

    Sums the dN/du × gammas_produced_per_yr over all source sectors matching
    the given isotope and region, then normalizes to a PDF (integral = 1).
    The total rate is stored in the CSV header comment.

    Parameters
    ----------
    results : list of dict
        From compute_gamma_tables().
    u_bins : ndarray
    iso : str
        'Bi214' or 'Tl208'.
    region : str
        'barrel' or 'endcap'.
    filepath : str
    """
    # Accumulate absolute dN/du (gammas/yr per unit u)
    dNdu_total = np.zeros_like(u_bins)
    for r in results:
        if r['isotope'] == iso and r['region'] == region:
            dNdu_total += r['dNdu_exiting'] * r['gammas_produced_per_yr']

    # Total rate
    total_rate = np.trapz(dNdu_total, u_bins)

    # Normalize to PDF
    if total_rate > 0:
        pdf = dNdu_total / total_rate  # pdf integrates to 1
    else:
        pdf = np.zeros_like(u_bins)

    # Save with total rate in header
    df = pd.DataFrame({'u': u_bins, 'dNdu_pdf': pdf})
    with open(filepath, 'w') as f:
        f.write(f"# total_gammas_per_yr = {total_rate:.6e}\n")
        f.write(f"# isotope = {iso}\n")
        f.write(f"# region = {region}\n")
        df.to_csv(f, index=False)

    return total_rate, pdf


# ===========================================================================
# Linear fit for MC sampling
# ===========================================================================
def fit_linear_pdf(u_bins, pdf, u_min=0.3, u_max=1.0):
    """
    Fit a linear model dN/du = a + b*u to the PDF in the range [u_min, u_max].

    Below u_min, the PDF is negligible (<2% of the integral) and can be
    safely ignored in MC sampling. Above u_min, the PDF is nearly linear,
    so this fit captures the shape accurately.

    The fit uses numpy's polyfit (degree 1) weighted by the PDF values
    to emphasize the region where most gammas are.

    Parameters
    ----------
    u_bins : ndarray
        Bin centers of the angular distribution.
    pdf : ndarray
        Normalized PDF values (integral = 1 over full range).
    u_min : float
        Lower cutoff for the fit (default 0.3).
    u_max : float
        Upper cutoff for the fit (default 1.0).

    Returns
    -------
    a : float
        Intercept of the linear fit (value at u=0).
    b : float
        Slope of the linear fit.
    norm : float
        Integral of (a + b*u) over [u_min, u_max], for normalization.
        To make a proper PDF: pdf(u) = (a + b*u) / norm.
    """
    mask = (u_bins >= u_min) & (u_bins <= u_max)
    u_fit = u_bins[mask]
    pdf_fit = pdf[mask]

    # polyfit returns [b, a] for degree 1 (highest power first)
    b, a = np.polyfit(u_fit, pdf_fit, 1)

    # Normalization: ∫_{u_min}^{u_max} (a + b*u) du
    norm = a * (u_max - u_min) + b / 2.0 * (u_max**2 - u_min**2)

    return a, b, norm


def sample_linear_pdf(a, b, u_min, u_max, n=1, rng=None):
    """
    Sample from the linear PDF p(u) = (a + b*u) / norm via inverse CDF.

    The CDF is: F(u) = [a*(u - u_min) + b/2*(u² - u_min²)] / norm

    Setting F(u) = r and solving the quadratic a*u + b/2*u² = C:
        u = (-a + sqrt(a² + 2*b*C)) / b    (for b ≠ 0)
    where C = r * norm + a*u_min + b/2*u_min²

    Parameters
    ----------
    a, b : float
        Linear PDF coefficients: p(u) ∝ (a + b*u).
    u_min, u_max : float
        Range of the distribution.
    n : int
        Number of samples.
    rng : numpy.random.Generator or None
        Random number generator. If None, uses default.

    Returns
    -------
    u_samples : ndarray
        Array of n sampled u values in [u_min, u_max].
    """
    if rng is None:
        rng = np.random.default_rng()

    norm = a * (u_max - u_min) + b / 2.0 * (u_max**2 - u_min**2)
    r = rng.uniform(0.0, 1.0, size=n)

    # Quadratic: b/2 * u² + a * u - (r*norm + a*u_min + b/2*u_min²) = 0
    # Using the standard form: A*u² + B*u + C = 0 with A=b/2, B=a
    offset = a * u_min + b / 2.0 * u_min**2
    C_vals = r * norm + offset  # = a*u + b/2*u² evaluated at the target

    if abs(b) < 1e-15:
        # Degenerate case: uniform PDF (b≈0), CDF is linear
        u_samples = u_min + r * (u_max - u_min)
    else:
        # Solve b/2 * u² + a * u = C  →  u = (-a + sqrt(a² + 2*b*C)) / b
        discriminant = a**2 + 2.0 * b * C_vals
        u_samples = (-a + np.sqrt(discriminant)) / b

    return u_samples


def save_linear_fit_csv(results, u_bins, filepath, u_min=0.3, u_max=1.0):
    """
    For each (isotope, region), fit a linear PDF to the angular distribution
    in [u_min, u_max] and save the coefficients to a CSV file.

    The output file contains one row per (isotope, region) with:
      - isotope, region: identifiers
      - a, b: linear PDF coefficients, p(u) = (a + b*u)/norm
      - u_min, u_max: fit range
      - norm: integral of (a + b*u) over [u_min, u_max]
      - total_gammas_per_yr: total gamma rate exiting skin for this combination
      - frac_in_fit_range: fraction of the full PDF integral within [u_min, u_max]

    To sample in MC:
      1. Draw r ~ Uniform(0,1)
      2. u = (-a + sqrt(a² + 2*b*(r*norm + a*u_min + b/2*u_min²))) / b
      3. θ = arccos(u), φ ~ Uniform(0, 2π)

    Parameters
    ----------
    results : list of dict
        From compute_gamma_tables().
    u_bins : ndarray
    filepath : str
    u_min, u_max : float
        Fit range.
    """
    rows = []

    for iso in ['Bi214', 'Tl208']:
        for region in ['barrel', 'endcap']:
            # Accumulate absolute dN/du
            dNdu_abs = np.zeros_like(u_bins)
            for r in results:
                if r['isotope'] == iso and r['region'] == region:
                    dNdu_abs += r['dNdu_exiting'] * r['gammas_produced_per_yr']

            # Total rate (all angles)
            total_rate = np.trapz(dNdu_abs, u_bins)

            if total_rate <= 0:
                rows.append({
                    'isotope': iso, 'region': region,
                    'a': 0, 'b': 0, 'u_min': u_min, 'u_max': u_max,
                    'norm': 0, 'total_gammas_per_yr': 0,
                    'frac_in_fit_range': 0,
                })
                continue

            # Normalize to PDF
            pdf = dNdu_abs / total_rate

            # Fraction in fit range
            mask = (u_bins >= u_min) & (u_bins <= u_max)
            frac = np.trapz(pdf[mask], u_bins[mask])

            # Linear fit
            a, b, norm = fit_linear_pdf(u_bins, pdf, u_min, u_max)

            rows.append({
                'isotope': iso,
                'region': region,
                'a': a,
                'b': b,
                'u_min': u_min,
                'u_max': u_max,
                'norm': norm,
                'total_gammas_per_yr': f"{total_rate:.6e}",
                'frac_in_fit_range': f"{frac:.4f}",
            })

    df = pd.DataFrame(rows)
    with open(filepath, 'w') as f:
        f.write("# LZ Cryostat: linear fit to angular PDF for MC sampling\n")
        f.write("# p(u) = (a + b*u) / norm  for u in [u_min, u_max]\n")
        f.write("#\n")
        f.write("# To sample u in MC:\n")
        f.write("#   1. Draw r ~ Uniform(0, 1)\n")
        f.write("#   2. C = r * norm + a * u_min + b/2 * u_min^2\n")
        f.write("#   3. u = (-a + sqrt(a^2 + 2*b*C)) / b\n")
        f.write("#   4. theta = arccos(u), phi ~ Uniform(0, 2*pi)\n")
        f.write("#\n")
        f.write("# total_gammas_per_yr: rate of gammas exiting skin (all angles)\n")
        f.write("# frac_in_fit_range: fraction of PDF integral within [u_min, u_max]\n")
        f.write("#   (gammas below u_min are negligible and safely ignored)\n")
        df.to_csv(f, index=False)
    return df


# ===========================================================================
# Console output
# ===========================================================================
def print_section(title):
    print(f"\n{'#'*70}")
    print(f"# {title}")
    print(f"{'#'*70}")


def print_vessel(v):
    """Print vessel geometry and mass breakdown."""
    label = v['label']
    code_min = {'ICV': 9.0, 'OCV': 7.0}
    print(f"\n  {label}: R = {v['R']:.1f} cm, H = {v['H_total']:.1f} cm")
    print(f"    Top head: {v['top_ar']}:1, depth {v['h_top']:.1f} cm, t = {v['t_top_mm']:.1f} mm")
    print(f"    Wall:     H_cyl = {v['H_cyl']:.1f} cm, t = {v['t_wall_mm']:.2f} mm [SOLVED]", end="")
    if label in code_min:
        print(f"  (code min {code_min[label]:.0f} mm, +{v['t_wall_mm']-code_min[label]:.2f} mm)")
    else:
        print()
    print(f"    Bot head: {v['bot_ar']}:1, depth {v['h_bot']:.1f} cm, t = {v['t_bot_mm']:.1f} mm")
    print(f"    Mass: shell {v['m_shell']:.1f} + heads {v['m_heads']:.1f} = {v['m_total']:.1f} kg")


def print_gamma_tables(results):
    """Print the three gamma rate tables to console."""
    print(f"\n  {'Source':<16s} {'Region':<8s} {'Isotope':<7s} "
          f"{'Produced':>12s} {'Enter skin':>12s} {'Exit skin':>12s} {'<cosθ>':>8s}")
    print(f"  {'':16s} {'':8s} {'':7s} "
          f"{'(γ/yr)':>12s} {'(γ/yr)':>12s} {'(γ/yr)':>12s} {'':>8s}")
    print(f"  {'-'*75}")
    for r in results:
        print(f"  {r['source']:<16s} {r['region']:<8s} {r['isotope']:<7s} "
              f"{r['gammas_produced_per_yr']:>12.2e} "
              f"{r['gammas_entering_skin_per_yr']:>12.2e} "
              f"{r['gammas_exiting_skin_per_yr']:>12.2e} "
              f"{r['mean_u_exiting']:>8.3f}")


def print_mc_summary(results, u_bins):
    """Print the four MC input summaries."""
    for iso in ['Bi214', 'Tl208']:
        for region in ['barrel', 'endcap']:
            dNdu_total = np.zeros_like(u_bins)
            for r in results:
                if r['isotope'] == iso and r['region'] == region:
                    dNdu_total += r['dNdu_exiting'] * r['gammas_produced_per_yr']
            total = np.trapz(dNdu_total, u_bins)
            if total > 0:
                mean_u = np.trapz(u_bins * dNdu_total, u_bins) / total
            else:
                mean_u = 0
            print(f"  {iso} {region:<8s}: {total:>12.2e} γ/yr, <cosθ> = {mean_u:.3f}")


# ===========================================================================
# Main
# ===========================================================================
def main():
    print_section("LZ CRYOSTAT MODEL — Gamma Transport to Active Volume")

    # --- Load attenuation coefficients ---
    mu = load_attenuation_coefficients()
    print("\n  Attenuation coefficients (cm⁻¹):")
    for key, val in mu.items():
        print(f"    μ_{key} = {val:.4f}")

    # --- Compute vessels ---
    print_section("VESSEL GEOMETRY AND MASSES")
    ocv = compute_vessel(**OCV_CFG)
    icv = compute_vessel(**ICV_CFG)
    print_vessel(ocv)
    print_vessel(icv)

    # Radial layout
    icv_inner_R = icv['R'] - icv['t_wall_mm'] / 10.0
    skin_outer_R = icv_inner_R
    fc_outer_R = TPC_INNER_R + FC_THICKNESS
    gap = ocv['R'] - icv['R']

    print(f"\n  Radial stack (mid-height, from axis):")
    print(f"    0 — {TPC_INNER_R} cm       active LXe")
    print(f"    {TPC_INNER_R} — {fc_outer_R} cm   field cage")
    print(f"    {fc_outer_R} — {skin_outer_R:.1f} cm  LXe skin ({skin_outer_R-fc_outer_R:.1f} cm)")
    print(f"    {skin_outer_R:.1f} — {icv['R']} cm    ICV wall ({icv['t_wall_mm']:.1f} mm)")
    print(f"    {icv['R']} — {ocv['R']} cm      vacuum gap ({gap:.1f} cm)")
    print(f"    {ocv['R']} — {ocv['R']+ocv['t_wall_mm']/10:.1f} cm   OCV wall ({ocv['t_wall_mm']:.1f} mm)")

    # --- Build source sectors ---
    sectors = build_source_sectors(ocv, icv, mu)

    # --- Radioactivity table ---
    print_section("RADIOACTIVITY")
    print(f"\n  {'Source':<16s} {'Region':<8s} {'Type':<10s} {'Mass':>7s} "
          f"{'Bi214/kg':>9s} {'Bi214':>9s} {'Tl208/kg':>9s} {'Tl208':>9s}")
    print(f"  {'-'*75}")
    tot_bi = 0; tot_tl = 0
    for s in sectors:
        bi = s['mass_kg'] * s['Bi214_sp']
        tl = s['mass_kg'] * s['Tl208_sp']
        tot_bi += bi; tot_tl += tl
        print(f"  {s['name']:<16s} {s['region']:<8s} {s['src_type']:<10s} {s['mass_kg']:>7.1f} "
              f"{s['Bi214_sp']:>9.2f} {bi:>9.1f} {s['Tl208_sp']:>9.2f} {tl:>9.1f}")
    print(f"  {'-'*75}")
    print(f"  {'TOTAL':<16s} {'':8s} {'':10s} {'':>7s} {'':>9s} {tot_bi:>9.1f} {'':>9s} {tot_tl:>9.1f}")

    # --- Gamma transport ---
    print_section("GAMMA TRANSPORT")
    results = compute_gamma_tables(sectors, U_BINS)
    print_gamma_tables(results)

    # --- MC input summary ---
    print_section("MC INPUT SUMMARY (gammas/yr exiting skin)")
    print_mc_summary(results, U_BINS)

    # --- Save CSV files ---
    print_section("SAVING OUTPUT FILES")

    f1 = os.path.join(DATA_DIR, "lz_cryo_geometry.csv")
    save_geometry_csv(ocv, icv, f1)
    print(f"  {f1}")

    f2 = os.path.join(DATA_DIR, "lz_cryo_activity.csv")
    save_activity_csv(sectors, f2)
    print(f"  {f2}")

    f3 = os.path.join(DATA_DIR, "lz_cryo_gammas.csv")
    save_gammas_csv(results, f3)
    print(f"  {f3}")

    for iso in ['Bi214', 'Tl208']:
        for region in ['barrel', 'endcap']:
            fname = f"lz_cryo_angular_{iso}_{region}.csv"
            fpath = os.path.join(DATA_DIR, fname)
            total, _ = save_angular_csv(results, U_BINS, iso, region, fpath)
            print(f"  {fpath}  ({total:.2e} γ/yr)")

    # --- Linear fit for MC sampling ---
    print_section("LINEAR FIT FOR MC SAMPLING (u > 0.3)")

    f_lin = os.path.join(DATA_DIR, "lz_cryo_linear_fit.csv")
    df_lin = save_linear_fit_csv(results, U_BINS, f_lin)
    print(f"  {f_lin}\n")
    print(f"  {'Isotope':<8s} {'Region':<8s} {'a':>10s} {'b':>10s} "
          f"{'norm':>10s} {'γ/yr':>12s} {'frac>0.3':>9s}")
    print(f"  {'-'*70}")
    for _, row in df_lin.iterrows():
        print(f"  {row['isotope']:<8s} {row['region']:<8s} "
              f"{float(row['a']):>10.6f} {float(row['b']):>10.6f} "
              f"{float(row['norm']):>10.6f} {row['total_gammas_per_yr']:>12s} "
              f"{row['frac_in_fit_range']:>9s}")

    print(f"\n  To sample u in MC:")
    print(f"    1. Draw r ~ Uniform(0, 1)")
    print(f"    2. C = r * norm + a * u_min + b/2 * u_min²")
    print(f"    3. u = (-a + sqrt(a² + 2*b*C)) / b")
    print(f"    4. θ = arccos(u), φ ~ Uniform(0, 2π)")

    print("\nDone.")


if __name__ == "__main__":
    main()
