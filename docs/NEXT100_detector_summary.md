# NEXT-100 Detector Summary for 0νββ Background Studies

Sources: NEXT-100 Detector Paper (arXiv:2505.17848),
NEXT-100 Cathode/EL Paper (arXiv:2311.03528, Mistry et al.),
NEXT-100 HV Paper (arXiv:2505.01002, Rogers et al.).

All radioactivities focus on **²¹⁴Bi** (from ²³⁸U chain, γ at 2.448 MeV)
and **²⁰⁸Tl** (from ²³²Th chain, γ at 2.615 MeV), the dominant radiogenic
backgrounds at Q_ββ = 2,458 keV.

---

## 1. Pressure Vessel and Shielding

### 1.1 Pressure Vessel

- **Material:** Stainless steel 316Ti (radiopure titanium-stabilized alloy)
- **Volume:** 1.7 m³
- **Maximum certified pressure:** 15 bar
- **Operating pressure:** 13.5 bar (design), currently 4.2 bar during commissioning

### 1.2 Inner Copper Shielding (ICS)

- **Material:** Ultra-pure electrolytic copper
- **Barrel:** 40 individual copper blocks, **120 mm thick**
- **Energy Plane (EP) copper plate:** 140 mm thick (single piece, with PMT coupling holes)
- **Tracking Plane (TP) copper plate:** 120 mm thick (single piece, with SiPM cable slots)
- **Total copper mass:** ~6,100 kg (barrel) + ~1,400 kg (EP plate) + ~1,500 kg (TP plate)
  (approximate, from Table 2 quantities)
- Purpose: blocks X-rays and secondary radiation from vessel walls

### 1.3 Outer Lead Shielding

- Lead bricks surrounding the pressure vessel
- Reduces cosmic-ray and natural radioactivity gammas from LSC cavern

---

## 2. TPC Dimensions

| Parameter                                | Value              |
|------------------------------------------|--------------------|
| Drift region length (cathode to gate)    | **1,187 mm**       |
| Buffer region length (cathode to EP)     | **241 mm**         |
| Total field cage length                  | **1,427 mm**       |
| Light tube (active volume) diameter      | **983 mm**         |
| Field cage outer diameter                | **1,105 mm**       |
| Copper ring internal diameter            | 1,014 mm           |
| Copper ring external diameter            | 1,038 mm           |
| EL gap                                   | 9.70 ± 0.15 mm     |
| Active Xe mass (at 13.5 bar)             | ~70.5 kg ¹³⁶Xe     |
| Drift field (design, 13.5 bar)           | 400 V/cm           |
| Cathode voltage (design)                 | −65 kV to −70 kV   |
| Operating pressure (design)              | 13.5 bar           |
| Energy resolution                        | ≤ 1% FWHM at Q_ββ  |
| Background target                        | 4×10⁻⁴ counts/(keV·kg·yr) after cuts |

### Current commissioning parameters (4.2 bar)

| Gas   | V_cathode (V) | V_gate (V) | Drift field (V/cm) |
|-------|--------------|------------|---------------------|
| Argon | 14,700       | 6,700      | 74                  |
| Xenon | 23,000       | 9,000      | 118                 |

---

## 3. Field Cage

### 3.1 Structure

- **52 copper field-shaping rings** total: 48 (drift) + 4 (buffer)
- Supported on **18 HDPE struts** running full detector length
- **PTFE reflector panels** (5 mm thick, TPB-coated) attached via dovetail joints
- **HDPE "Poly Wrap"** insulator: 2 layers of 6.35 mm (¼") HDPE around
  the field cage, separating it from the ICS

### 3.2 Copper Ring Specifications

- Internal diameter: **1,014 mm**, external diameter: **1,038 mm**
- Cross section: **10 mm × 8 mm**, rounded rectangle (1 mm corner radius)
- Cross-sectional area: 1.2 cm²
- Each ring made in **3 segments**, joined by 12 brass countersunk screws
- Ring pitch: **24 mm** (drift region), **38 mm** (buffer region)
- Rings provide the main structural support for the entire field cage

### 3.3 Resistor Chain

- **3 parallel arrays** of 100 MΩ film-chip resistors (Vishay Techno CHRV100MEDKR-ND)
  between each ring pair (drift region only)
- Mounted on **Cuflon boards** (PTFE + Cu, 28 × 6 × 2.7 mm)
- Boards on inner surface of rings, rotated 120° azimuthally
- Final ring to cathode: 150 MΩ resistor board + silicon bronze button
- Buffer region: **no resistors** (relies on HDPE resistivity)
- Total field cage resistance: **1.65 GΩ** (measured 1.68 GΩ)

### 3.4 Field Cage Masses

| Component            | Material | Mass (kg) |
|----------------------|----------|-----------|
| Field shaping rings  | Copper   | **177.9** |
| HDPE struts          | HDPE     | **54.0**  |
| PTFE reflectors      | PTFE     | **38.8**  |
| Resistors            | mixed    | **0.050** (50 mg) |
| **Total (excl. Poly Wrap)** |  | **270.7** |
| HDPE wrap            | HDPE     | ~50 kg (est.) |

### 3.5 Field Cage Radioactivity (from NEXT-HV paper, Table 1)

| Component        | Material | ²³²Th (ppt) | ²³⁸U (ppt) | Mass (kg) |
|------------------|----------|-------------|-------------|-----------|
| FC rings         | Copper   | < 0.8       | < 0.6       | 177.9     |
| PTFE reflectors  | PTFE     | 2.6 ± 0.1   | 4.9 ± 1.0   | 55.2*     |
| Poly Wrap+struts | HDPE     | 15 ± 5      | 3 ± 1       | 88.3      |
| FC resistors     | mixed    | (0.22±0.02)×10⁶ | (0.28±0.01)×10⁶ | 0.050 |

*Mass differs slightly between NEXT-HV (55.2 kg) and NEXT-100 papers (38.8 kg for reflectors alone).

Converting to mBq/kg (1 ppt ²³⁸U ≈ 0.0124 mBq/kg; 1 ppt ²³²Th ≈ 0.00406 mBq/kg):

| Component        | ²³⁸U (mBq/kg)  | ²³²Th (mBq/kg) | Total ²¹⁴Bi (mBq) | Total ²⁰⁸Tl (mBq)** |
|------------------|-----------------|-----------------|--------------------|--------------------|
| FC rings (Cu)    | < 0.007         | < 0.003         | < 1.3              | < 0.2              |
| PTFE reflectors  | 0.061 ± 0.012   | 0.011 ± 0.0004  | ~3.4               | ~0.6               |
| HDPE wrap+struts | 0.037 ± 0.012   | 0.061 ± 0.020   | ~3.3               | ~5.4               |
| FC resistors     | ~3,472          | ~893            | ~0.17              | ~0.04              |

**²⁰⁸Tl activity = ²³²Th_late × 0.359 (BR of ²⁰⁸Tl from ²¹²Bi). Values here are
for ²³²Th_late assuming secular equilibrium; effective ²⁰⁸Tl is 35.9% of these.

---

## 4. Electrode Grids (Cathode + EL Region)

### 4.1 Mesh Specifications

**All meshes:** photoetched from a **127 µm thick SS316Ti stainless steel sheet**,
hexagonal pattern. Manufactured by PCM Products, Inc. Diameter: ~1 m.

| Mesh     | Hex inner Ø | Optical transparency | Tension   | Position                 |
|----------|-------------|---------------------|-----------|--------------------------|
| Cathode  | 5 mm        | 95%                 | ~0.6 kN   | Between drift and buffer |
| Gate     | 2.5 mm      | 90%                 | 835–990 N | Entry to EL region       |
| Anode    | 2.5 mm      | 90%                 | ~0.9 kN   | Exit of EL region (ground) |

- Wire width: **127 µm** (same as sheet thickness — etched, not woven)
- EL gap: **9.70 ± 0.15 mm** (gate to anode)
- EL field: ~20 kV/cm (at operating pressure)
- Cathode local field: ~2 kV/cm
- Breakdown tested to 28 kV without structural damage

### 4.2 Frame Specifications

- **Material:** Silicon bronze (CuSn8P alloy)
  - Chosen over SS because no radiopure SS source at required scale could be found
  - Young's modulus: 115 GPa; yield strength: 205 MPa
- Each mesh has a **base ring** + **tensioning ring** (two-part frame)
- Cathode frame thickness: **13.5 mm**
- EL frames connected by **8 HDPE U-shaped brackets** (ridged inner surface
  to minimize leakage current)

### 4.3 Grid Masses (from Mistry et al., Tables 1–2)

| Component                       | Material  | Mass        |
|---------------------------------|-----------|-------------|
| **EL + cathode frames**         | Si-bronze | **26.1 kg** |
| EL/cathode meshes (3 total)     | SS316Ti   | **509 g**   |
| Dowel pins and screws (540+216) | CuSn8P    | 329 g       |
| EL brackets (8)                 | HDPE      | 324 g       |
| Gate screws (16)                | PEEK      | 2.7 g       |
| Anode screws (16)               | 18-8 SS   | 22.4 g      |

### 4.4 Grid Radioactivity

**Silicon bronze frames (selected batches):**

| Batch   | ²³²Th (ppt)   | ²³⁸U (ppt)    |
|---------|---------------|----------------|
| 2023-03 | 4.9 ± 1.1    | 14.6 ± 1.3    |
| 2023-08 | 4.9 ± 0.6    | 2.3 ± 0.7     |

Converting: ²³⁸U ~0.03–0.18 mBq/kg, ²³²Th ~0.02 mBq/kg.

**Other grid components:**

| Component       | ²³²Th (ppt)    | ²³⁸U (ppt)    | Mass   |
|-----------------|---------------|----------------|--------|
| SS316Ti mesh    | 202 ± 7       | 604 ± 38       | 509 g  |
| CuSn8P screws   | < 100         | < 100          | 329 g  |
| HDPE brackets   | 15 ± 5        | 3 ± 1          | 324 g  |
| PEEK screws     | < 4.1 (²⁰⁸Tl)| < 0.6 (²¹⁴Bi) | 2.7 g  |
| SS anode screws | 2630 ± 1490   | 480 ± 230      | 22.4 g |

**Total grid contribution to background:** 7.4×10⁻³ counts/yr (²⁰⁸Tl) and
8.8×10⁻² counts/yr (²¹⁴Bi) — well within the NEXT-100 specification.

---

## 5. Sensors

### 5.1 Energy Plane (EP) — PMTs

- **53 PMTs** installed (60 positions available; 7 sapphire windows failed QC)
- **Model:** Hamamatsu R11410-10, 3-inch (64 mm surface diameter)
- Bialkali photocathode, operated in **vacuum** behind the EP copper plate
- Gain: 5 × 10⁶
- Coupled to gas volume through **sapphire windows** with PEDOT coating
  (resistive, transparent, defines ground at EP, prevents charge buildup)
- Sapphire windows also coated with **TPB** wavelength shifter
- PMT-to-window optical coupling: NOA76 gel (Edmund Optics)
- PMT bases: Kapton PCB, covered with copper cap filled with radiopure epoxy
- Signal readout: pseudo-differential (last and third-to-last dynodes)
- 12 LEDs distributed on EP for calibration

### 5.2 Tracking Plane (TP) — SiPMs

- **3,584 SiPMs** (Hamamatsu S13372-1350TE)
- Effective photosensitive area: 1.3 × 1.3 mm²
- Mounted on **Kapton DICE-Boards** (flexible PCBs)
- In direct contact with xenon gas
- Provides x,y position via S2 light pattern

### 5.3 Sensor Radioactivity (from Table 2 of NEXT-100 paper)

| Component                           | Quantity     | ²¹⁴Bi (mBq) | ²⁰⁸Tl (mBq) |
|-------------------------------------|-------------|-------------|-------------|
| **PMTs (R11410-10)**                | 60 units    | 21          | 11          |
| **PMT bases**                       | 60 units    | 41          | 15          |
| Sapphire windows + PEDOT + Cu cans  | 60 units    | < 66        | < 22        |
| SS316Ti M18 screws (EP Cu plate)    | 36 kg       | 80          | 47          |
| Optical gel (NuSil)                 | 0.24 kg     | < 7.9       | < 3.2       |
| EP copper plate                     | ~1,400 kg   | 1.6         | 0.14        |
| Cathode mesh (SS316Ti)              | 0.17 kg     | 1.3         | 0.05        |
| **Kapton connectors FX11LA**        | 5.6 units   | **160**     | **120**     |
| PTFE masks (TP)                     | 5.6 units   | < 14        | < 4.2       |
| SiPMs (S13372-1350TE)              | 3,600 units | < 12        | < 2.8       |
| Kapton DICE-Boards                  | 5.6 units   | 3.9         | < 0.58      |
| Gate + anode meshes (SS316Ti)       | 0.34 kg     | 2.5         | 0.10        |
| TP copper plate                     | ~1,500 kg   | 1.8         | 0.15        |
| Si-bronze gate + anode rings        | 17 kg       | 1.1         | 0.28        |

**Note:** The Kapton cable connectors (Hirose FX11LA) are flagged as the only
component exceeding requirements by ~10×. They will be replaced in a future
detector intervention.

---

## 6. Barrel Region Components (from Table 2)

| Component                         | Mass (kg) | ²¹⁴Bi (mBq) | ²⁰⁸Tl (mBq) |
|-----------------------------------|-----------|-------------|-------------|
| Electrolytic Cu (ICS + FC rings)  | 6,100     | < 45        | < 7.1       |
| FC resistors (Vishay Techno)      | 160 units | 27          | 2.5         |
| SS316Ti M6×20 screws (ICS)       | 0.72 kg   | 8.9         | 0.74        |
| PTFE reflector panels             | 55 kg     | 2.8         | 0.29        |
| HDPE wrap                         | 50 kg     | 1.8         | 1.1         |
| HDPE struts                       | 39 kg     | 1.4         | 0.85        |
| Other                             | —         | < 1.1       | < 0.28      |
| **Barrel sub-total**              |           | **< 88**    | **< 13**    |

---

## 7. Full Radioactivity Budget Summary

### By detector region (from Table 2 of NEXT-100 paper)

| Region              | ²¹⁴Bi total (mBq) | ²⁰⁸Tl total (mBq) |
|---------------------|--------------------|---------------------|
| Energy end-cap      | < 222              | < 99                |
| Tracking end-cap    | < 192              | < 130               |
| Barrel region       | < 88               | < 13                |
| **Total (in-vessel)** | **< 500**        | **< 240**           |

### Dominant ²¹⁴Bi sources (ranked by total activity)

| Component                        | ²¹⁴Bi (mBq) |
|----------------------------------|-------------|
| Kapton connectors FX11LA (TP)    | 160         |
| SS screws for EP Cu plate        | 80          |
| Sapphire windows + coatings      | < 66        |
| Electrolytic copper (ICS+rings)  | < 45        |
| PMT bases                        | 41          |
| FC resistors                     | 27          |
| PMTs (R11410-10)                 | 21          |

### Dominant ²⁰⁸Tl sources (ranked by total activity)

| Component                        | ²⁰⁸Tl (mBq) |
|----------------------------------|-------------|
| Kapton connectors FX11LA (TP)    | 120         |
| SS screws for EP Cu plate        | 47          |
| Sapphire windows + coatings      | < 22        |
| PMT bases                        | 15          |
| PMTs (R11410-10)                 | 11          |
| Electrolytic copper (ICS+rings)  | < 7.1       |

---

## 8. Key Design Differences: NEXT-100 vs LZ

| Feature             | NEXT-100                        | LZ                              |
|---------------------|---------------------------------|----------------------------------|
| Medium              | HPGXe (13.5 bar gas)            | LXe (liquid)                     |
| Active mass         | ~70.5 kg ¹³⁶Xe                 | 7 tonnes nat. Xe                 |
| TPC diameter        | 983 mm (light tube)             | 1,456 mm                         |
| Drift length        | 1,187 mm                        | 1,456 mm                         |
| Drift field         | 400 V/cm                        | 300 V/cm                         |
| FC rings            | 52 × copper                     | 65 × titanium                    |
| FC ring mass        | 177.9 kg (Cu)                   | 91.6 kg (Ti)                     |
| FC ring radiopurity | < 0.6 ppt U, < 0.8 ppt Th      | ~7 ppt U, ~56 ppt Th             |
| Grids               | 3 photoetched SS hexagonal mesh | 4 woven SS square mesh            |
| Grid frames         | 26.1 kg silicon bronze           | 62.2 kg SS/Ti (holders)          |
| Grid mesh mass      | 509 g SS316Ti                   | 800 g SS304                      |
| Amplification       | Electroluminescence (gas)       | Electroluminescence (gas)        |
| Energy readout      | 53 PMTs (R11410-10)             | 494 PMTs (R11410-22) top+bottom  |
| Tracking readout    | 3,584 SiPMs                     | (from S2 pattern in top PMTs)    |
| Vessel              | SS 316Ti, 15 bar rated          | Ti Grade 1 (ICV+OCV), 4 bar     |
| Inner shielding     | 120 mm Cu                       | LXe skin (40–80 mm)             |
| Energy resolution   | ≤ 1% FWHM at Q_ββ              | ~3% FWHM at Q_ββ (estimated)    |
| Background target   | 4×10⁻⁴ cts/(keV·kg·yr)         | (DM-optimized, not 0νββ)        |
