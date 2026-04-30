# LZ Detector Summary for 0νββ Background Studies

Sources: LZ TDR (arXiv:1703.09144), LZ Instrument Paper (arXiv:1910.09124),
LZ Radioactivity Paper (arXiv:2006.02506), LZ Titanium Paper (arXiv:1702.02646).

**NOTE:** For the 0νββ background simulation, the masses and specific
activities from the **LZ bb0nu paper** (Table I) supersede the values
compiled here from the TDR and radioactivity paper. The bb0nu paper values
are used in the simulation code (`py/lz_cryo.py`, `py/lz_fieldcage.py`,
`py/lz_pmts.py`) and documented in the design files
(`design/lz_cryostat_model.md`, `design/lz_fieldcage_model.md`,
`design/lz_pmt_model.md`). This file retains the original TDR/radioactivity
paper compilation for reference.

All radioactivities are **²³⁸U_late** (²¹⁴Bi chain, relevant γ: 2.448 MeV at 1.55% BR)
and **²³²Th_late** (²⁰⁸Tl chain, relevant γ: 2.615 MeV at 99.8% BR from ²⁰⁸Tl;
²⁰⁸Tl itself is 35.9% BR from ²¹²Bi).

---

## 1. Cryostat

**Material:** ASTM Grade 1 titanium (TIMET HN3469), EBCH melt, 0% scrap.
A 5-tonne Ti Gr-1 slab from TIMET supplied all stock for vessels, flanges,
ports and welding wire. Fabricated by Loterios, Milan. ASME BPVC VIII Div. 1.

Engineering drawing: `docs/cryostat.png` (authoritative for geometry).

### 1.1 Inner Cryostat Vessel (ICV)

The ICV has a 2:1 ellipsoidal top head, a cylindrical body that tapers near
half-height to reduce passive LXe near the cathode, and a 3:1 ellipsoidal
dished bottom (chosen specifically to minimize LXe volume below the cathode).
Split near the top by a flange pair. A stiffening ring sits at the top of the
conical section to allow thinner walls.

| Quantity                                 | Value         | Source         |
|------------------------------------------|---------------|----------------|
| Top neck flange ID                       | **1,143 mm**  | drawing        |
| Upper region (top head) ID               | **1,700 mm**  | drawing        |
| Upper cylindrical body ID                | **1,662 mm**  | drawing        |
| Lower cylindrical body ID (after taper)  | **1,578 mm**  | drawing        |
| Dished-end ID near cathode               | **1,458 mm**  | drawing        |
| Body ID (TDR cylindrical-section range)  | 1.58–1.66 m   | TDR Table 1.2.1|
| Height                                   | **2.59 m**    | TDR Table 1.2.1|
| Bare vessel mass                         | **950 kg**    | TDR Table 5.4.1|

Wall thicknesses (TDR Table 5.4.1, ASME-code minimum as-built):

| Segment        | Wall thickness (mm) |
|----------------|---------------------|
| Top head       | 7                   |
| Upper wall     | 9                   |
| Conical section| 9                   |
| Lower wall     | 9                   |
| Dished end     | 11                  |

### 1.2 Outer Cryostat Vessel (OCV)

The OCV is split into three flanged segments (top head, middle cylinder,
bottom head) so each piece fits in the Yates shaft. Top head is 2:1
ellipsoidal; bottom head is 2:1 ellipsoidal with three welded support feet.
The ICV is suspended from the OCV top head by three tie-bar assemblies with
external leveling adjustment.

| Quantity                          | Value                                     | Source           |
|-----------------------------------|-------------------------------------------|------------------|
| Body OD                           | **1,842–1,844 mm**                        | drawing          |
| Top flange OD                     | **1,924 mm**                              | drawing          |
| Bottom OD (at dished end ring)    | **1,856 mm**                              | drawing          |
| Inside diameter                   | **1.83 m**                                | TDR Table 1.2.1  |
| Total height                      | **3.04 m** (= 599 + 1,142 + 1,298 mm)     | TDR + drawing    |
| Bare vessel mass                  | **1,115 kg**                              | TDR Table 5.4.1  |

| Segment    | Wall thickness (mm) |
|------------|---------------------|
| Top head   | 8                   |
| Side wall  | 7                   |
| Dished end | 14                  |

### 1.3 Full Cryostat Assembly

- Three Cryostat Support (CS) legs at the OCV base; CS shelves carry the lower
  acrylic vessels of the Outer Detector.
- Bare vessel walls (TDR Table 5.4.1): ICV 950 kg + OCV 1,115 kg = **2,065 kg**.
- **Full Ti cryostat assembly mass: 2,590 kg** (bb0nu Table I — used for activity
  normalization). The ~525 kg above the bare-wall total accounts for flanges,
  ports, tie-bar hardware, coldhead fins and other welded Ti elements.
- Cryostat Seals: 33.7 kg (TDR Table 9.2.7) †
- Cryostat Insulation (MLI): 13.8 kg (bb0nu Table I)
- Cryostat Teflon Liner: 26 kg (TDR Table 9.2.7) †

† Not listed as a separate row in bb0nu Table I; values retained from the
LZ TDR Materials Table (Table 9.2.7) for completeness. In bb0nu these
contributions are folded into the "Ti cryostat vessel", "Cryostat insulation"
or "Other components" rows.

### 1.4 Cryostat Radioactivity

For the 0νββ background simulation, the **bb0nu paper Table I** values are
authoritative. The seals and Teflon liner rows (italicized) are kept from
TDR Table 9.2.7 since bb0nu does not break them out separately.

| Component              | Mass (kg) | ²³⁸U_late (mBq/kg) | ²³²Th_late (mBq/kg) | Total ²³⁸U_l (mBq) | Total ²³²Th_l (mBq) | Source            |
|------------------------|-----------|---------------------|----------------------|---------------------|----------------------|-------------------|
| Ti Cryostat Vessel     | **2,590** | **0.08** (UL)       | **0.22** (UL)        | ~207                | ~570                 | bb0nu Table I     |
| Cryostat Insulation    | **13.8**  | **11.1**            | **7.79**             | ~153                | ~107                 | bb0nu Table I     |
| *Cryostat Seals*       | *33.7*    | *26.2*              | *4.24*               | *~883*              | *~143*               | *TDR Table 9.2.7* † |
| *Cryostat Teflon Liner*| *26*      | *0.02*              | *0.03*               | *~0.5*              | *~0.8*               | *TDR Table 9.2.7* † |

UL = upper limit (denoted "†" in bb0nu Table I).

---

## 2. TPC Dimensions

| Parameter                               | Value           |
|-----------------------------------------|-----------------|
| Drift region height (cathode to gate)   | 1,456 mm        |
| TPC inner diameter                      | 1,456 mm        |
| Field cage thickness                    | 15 mm           |
| Reverse field region height             | 137.5 mm        |
| Skin thickness at surface (at cathode)  | 40 mm (80 mm)   |
| Active LXe mass                         | 7 tonnes        |
| Fiducial LXe mass                       | 5.6 tonnes      |
| Total Xe inventory                      | 10 tonnes       |
| Drift field (nominal)                   | 300 V/cm        |
| Cathode voltage (nominal / design)      | −50 kV / −100 kV|
| Operating pressure                      | 1.8 bara        |
| Gas gap (liquid surface to anode)       | 8 mm            |
| Gate-anode separation                   | 13 mm           |
| Maximum drift time                      | ~800 µs         |

---

## 3. Electrode Grids

**Material:** 304 stainless steel, ultra-finish wire, woven into 90° meshes.
Glued onto holder rings with MasterBond EP29LPSP cryogenic epoxy.
Final cathode, gate, and anode grids passivated in citric acid.

### 3.1 Grid Parameters (from LZ Instrument Paper, Table 1)

| Grid     | Voltage (kV) | Wire Ø (µm) | Pitch (mm) | Num. wires | Position               |
|----------|-------------|-------------|------------|------------|------------------------|
| Anode    | +5.75       | 100         | 2.5        | 1,169      | In gas, 8 mm above surface |
| Gate     | −5.75       | 75          | 5.0        | 583        | In liquid, 5 mm below surface |
| Cathode  | −50.0       | 100         | 5.0        | 579        | Bottom of drift region |
| Bottom   | −1.5        | 75          | 5.0        | 565        | Below RFR, shields PMTs |

### 3.2 Grid Masses and Radioactivity (from Table 9.2.7)

The TDR groups grids into two categories:

| Component    | Mass (kg) | ²³⁸U_late (mBq/kg) | ²³²Th_late (mBq/kg) | Total ²³⁸U_l (mBq) | Total ²³²Th_l (mBq) |
|--------------|-----------|---------------------|----------------------|---------------------|----------------------|
| Grid Wires (SS)   | 0.8  | 0.27                | 0.49                 | ~0.2                | ~0.4                 |
| Grid Holders (SS/Ti) | 62.2 | 0.27             | 0.49                 | ~16.8               | ~30.5                |

The grid holders (the rings onto which the wire meshes are glued) are the dominant
mass. The wire mass is small. The radioactivity values (0.27 / 0.49 mBq/kg) correspond
to stainless steel (NIRONIT) used for fasteners and structural elements. The holder
rings themselves may include both Ti and SS components.

### 3.3 Notes on Grids

- All grids are woven meshes (not stretched single-direction wires). This gives
  azimuthally uniform load on the holder ring and better field uniformity.
- Each wire tensioned with 250 g weights, mesh glued to holder ring with
  cryogenic epoxy containing acrylic beads.
- A second metal ring captures the glued region.
- The gate grid was passivated in citric acid (2 hours at ~125°F) to suppress
  spurious electron emission from cathodic wire surfaces.
- Grid deflection: gate-anode gap decreases ~1.6 mm at center at operating voltage.

---

## 4. Field Cage and Field Shaping Rings

### 4.1 Structure

- **57 Ti field-shaping rings** in drift region, equally spaced
- **8 Ti rings** in Reverse Field Region (RFR)
- **Total: 65 rings**
- PTFE layers: 58 (drift, 25 mm tall each) + 8 (RFR, ~15 mm tall each)
- PTFE azimuthally segmented 24 times per layer
- Ring shape: **"I"-shaped** (drift region); **oval** (RFR)
- Resistor ladder (drift): pairs of **2 GΩ** per section
  (first step: 1 GΩ ∥ 2 GΩ for field tuning near cathode)
- Resistor ladder (RFR): pairs of **5 GΩ** per section

### 4.2 Mass and Radioactivity

| Component              | Mass (kg) | ²³⁸U_late (mBq/kg) | ²³²Th_late (mBq/kg) | Total ²³⁸U_l (mBq) | Total ²³²Th_l (mBq) |
|------------------------|-----------|---------------------|----------------------|---------------------|----------------------|
| Field Shaping Rings (Ti) | **91.6** | **0.09**            | **0.23**             | ~8.2                | ~21.1                |
| TPC PTFE (reflectors)  | 184       | 0.02                | 0.03                 | ~3.7                | ~5.5                 |

Average ring mass: ~1.4 kg each. Same TIMET HN3469 Ti as cryostat.

---

## 5. Field Shaping Resistors

**From Table 9.2.1 (TDR Materials Table), under "TPC Items":**

| Chain          | Activity (mBq/kg) |
|----------------|-------------------|
| ²³⁸U early     | 5,679             |
| ²³⁸U late      | **1.73**          |
| ²³²Th early    | 0.57              |
| ²³²Th late     | **0.57**          |

- Total mass: very small (surface-mount resistors, order of grams)
- Material: alumina ceramic substrate — extremely high ²³⁸U_early from the ceramic,
  but the late-chain activities (relevant for 0νββ) are modest.
- Glass-coated (no epoxy potting) to minimize mass.
- Smallest available ceramic mass for the voltage rating.
- The TDR explicitly identifies resistors as a **major background contributor for 0νββ**
  (Section 2.4, alongside PMTs and the xenon vessel), primarily due to the ²⁰⁸Tl
  2,615 keV line from ²³²Th and the ²¹⁴Bi 2,448 MeV line from ²³⁸U.
- RFR resistors are physically larger (higher voltage rating) than drift resistors
  and thus more radioactively challenging.

---

## 6. PMTs

### 6.1 PMT Arrays and Positioning

| Array                  | PMT Model        | Diameter | Count | Location                           |
|------------------------|------------------|----------|-------|------------------------------------|
| TPC Top                | R11410-22        | 3"       | 253   | In gas phase, downward-looking     |
| TPC Bottom             | R11410-22        | 3"       | 241   | In liquid, upward-looking          |
| Side Skin (upper)      | R8520-406        | 1"       | 93    | Outside field cage, below surface  |
| Side Skin (lower ring) | R8778            | 2"       | 20    | Bottom of vessel, upward-looking   |
| Dome Skin              | R8778            | 2"       | 18    | Below bottom array, horizontal     |
| Outer Detector         | R5912            | 8"       | 120   | In water, behind Tyvek curtain     |
| **Total**              |                  |          | **745** |                                  |

**TPC PMT layout:**
- Top array (253): hybrid pattern — hexagonal at center, circular at periphery.
  Final PMT row overhangs the field cage inner walls for optimal wall-event
  reconstruction. Located in gas phase.
- Bottom array (241): close-packed hexagonal. In liquid phase.
- Both arrays mounted on Ti support structures (same Ti batch as cryostat).

**PMT performance (R11410-22):**
- Average cold QE: 30.9%
- Nominal gain: 3.5 × 10⁶
- PDE for S1: ~12%

### 6.2 PMT Radioactivity (²¹⁴Bi and ²⁰⁸Tl chains only)

| Component              | Mass (kg) | ²³⁸U_late             | ²³²Th_late            | Total ²³⁸U_l (mBq) | Total ²³²Th_l (mBq) |
|------------------------|-----------|----------------------|----------------------|---------------------|----------------------|
| R11410 TPC PMTs (×494) | 91.9      | 0.9 ± 0.2 mBq/PMT   | 0.8 ± 0.2 mBq/PMT   | ~445                | ~395                 |
| R11410 PMT Bases       | 2.8       | 75.8 mBq/kg          | 27.9 mBq/kg          | ~212                | ~78                  |
| R8778 2" Skin PMTs     | 6.1       | 59.4 mBq/kg          | 16.9 mBq/kg          | ~362                | ~103                 |
| R8520 1" Skin PMTs     | 2.2       | 5.19 mBq/kg          | 4.75 mBq/kg          | ~11                 | ~10                  |
| R8520 Skin PMT Bases   | 0.2       | 108 mBq/kg           | 37.6 mBq/kg          | ~22                 | ~8                   |
| R5912 OD PMTs          | 205       | 470 mBq/kg           | 388 mBq/kg           | ~96,350             | ~79,540              |

The OD PMTs have very high intrinsic activity but are far from the fiducial volume
(shielded by cryostat + LXe + water gap). For 0νββ, the TPC R11410 PMTs and their
bases are the dominant PMT background source.

233 of the 494 R11410 tubes were individually screened by HPGe. Unscreened tubes
assigned the average value of 12.2 mBq/PMT for ⁴⁰K (not relevant for 0νββ).

---

## 7. Other Components

| Component                | Mass (kg) | ²³⁸U_late (mBq/kg) | ²³²Th_late (mBq/kg) | Total ²³⁸U_l (mBq) | Total ²³²Th_l (mBq) |
|--------------------------|-----------|---------------------|----------------------|---------------------|----------------------|
| Upper PMT Structure (Ti) | 40.5      | 0.23                | 0.38                 | ~9.3                | ~15.4                |
| Lower PMT Structure (Ti) | 69.9      | 0.13                | 0.24                 | ~9.1                | ~16.8                |
| PMT Cabling              | 104       | 1.47                | 3.15                 | ~153                | ~328                 |
| HV Components            | 138       | 2.00                | 0.60                 | ~276                | ~83                  |
| Conduits                 | 200       | 0.40                | 0.66                 | ~80                 | ~132                 |
| TPC Sensors              | 0.90      | 13.5                | 14.2                 | ~12                 | ~13                  |
| Xe Tubing                | 15.1      | 0.18                | 0.33                 | ~2.7                | ~5.0                 |

---

## 8. Total Activity Inventory (²¹⁴Bi and ²⁰⁸Tl chains, near TPC)

Summary of the dominant sources by total activity, excluding the distant OD:

| Component              | Total ²³⁸U_late (mBq) | Total ²³²Th_late (mBq) |
|------------------------|----------------------|----------------------|
| Cryostat Seals         | ~883                 | ~143                 |
| R11410 TPC PMTs        | ~445                 | ~395                 |
| Cryostat Insulation    | ~450                 | ~82                  |
| R8778 Skin PMTs        | ~362                 | ~103                 |
| HV Components          | ~276                 | ~83                  |
| Cryostat Vessel (Ti)   | < 217                | ~603                 |
| R11410 PMT Bases       | ~212                 | ~78                  |
| PMT Cabling            | ~153                 | ~328                 |
| Conduits               | ~80                  | ~132                 |
| Grid Holders           | ~17                  | ~31                  |
| Upper+Lower PMT Struct.| ~18                  | ~32                  |
| Field Shaping Rings    | ~8                   | ~21                  |
| TPC PTFE               | ~4                   | ~6                   |
| Grid Wires             | ~0.2                 | ~0.4                 |

**For 0νββ at Q_ββ = 2,458 keV:**
- The ²¹⁴Bi 2,448 MeV γ (1.55% BR) sits essentially at Q_ββ.
  Dominant sources: cryostat seals, cryostat insulation, PMTs, HV components.
- The ²⁰⁸Tl 2,615 MeV γ (effective BR = 35.9% from ²¹²Bi) is above Q_ββ
  but its Compton edge and tail reach into the ROI.
  Dominant sources: cryostat vessel (Ti), PMTs, PMT cabling.
- Self-shielding of the LXe (~10 cm attenuation length at 2.5 MeV) means that
  only the innermost components (PMTs, field cage, grids, PTFE) contribute
  significantly to the fiducial background. The cryostat, despite high total
  activity, is partially shielded by the skin LXe layer (4–8 cm).
