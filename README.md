# Soil Health Impacts on Water Storage and Hydraulic Conductivity

Calculations and visualizations of how bulk density and organic matter content impact soil water-holding properties and hydraulic conductivity for every USDA texture class. These relationships help quantify how soil health practices (i.e. decompaction; compost amendments; diverse, dense, deep-rooted vegatation) can improve resiliency to drought and storms by increasing water storage and infiltration.

Calculations are based on USDA [Rosetta v3](https://github.com/usda-ars-ussl/rosetta-soil) pedotransfer functions recalibrated in 2017 and on the meta-analysis of organic matter effects on water storage by [Minasny & McBratney (2018)](https://bsssjournals.onlinelibrary.wiley.com/doi/abs/10.1111/ejss.12475). Dynamic visualizations are built with the interactive HoloViews/Bokeh libraries for exploration in Jupyter Notebooks or exportable as HTML.

## Notebooks

| Notebook | Purpose |
| --- | --- |
| [1_rosetta_porosity_by_texture.ipynb](notebooks/1_rosetta_porosity_by_texture.ipynb) | **ROSETTA texture × bulk-density baseline** — total porosity, field capacity, wilting point, available & drainable water, plus Ksat and the Mualem–van Genuchten parameters for each texture class across bulk densities 0.8–1.9 g/cm³. Writes `rosetta_porosity_by_texture.csv`. |
| [2_organic_matter_water_holding.ipynb](notebooks/2_organic_matter_water_holding.ipynb) | **Organic-matter effects** — layers SOM/SOC sensitivity on top, two ways: a ROSETTA + Minasny–McBratney (2018) blend (anchored to a UNSODA-derived mineral-baseline OC), and Saxton–Rawls (2006). **Reads the CSV from Notebook 1.** |
| [3_rosetta_hydraulic_conductivity.ipynb](notebooks/3_rosetta_hydraulic_conductivity.ipynb) | **Hydraulic conductivity & infiltration** — saturated Ksat, unsaturated K(h) from the Mualem–van Genuchten parameters, and Green–Ampt infiltration f(t)/F(t), in stormwater units. **Reads the CSV from Notebook 1.** |

**Run order:** execute `1_rosetta_porosity_by_texture.ipynb` first (it writes
`rosetta_porosity_by_texture.csv`); then `2_organic_matter_water_holding.ipynb` and
`3_rosetta_hydraulic_conductivity.ipynb`, which both read that CSV (independently of each other).

## Notebook 1 — ROSETTA porosity by texture and bulk density

Uses the [Rosetta](https://github.com/usda-ars-ussl/rosetta-soil) pedotransfer functions to
estimate, for every USDA soil texture class (using representative *median* sand/silt/clay values)
and for bulk densities from 0.8 to 1.9 g/cm³ in 0.1 steps:

- **Total porosity** — saturated volumetric water content θₛ
- **Field-capacity porosity** — volumetric water content at 33 kPa (330 cm suction)
- **Permanent-wilting-point porosity** — volumetric water content at 1500 kPa (15000 cm suction)
- **Available water capacity** — field capacity minus wilting point
- **Drainable water** — saturation minus field capacity (shown in the band plots)
- **Saturated hydraulic conductivity** — Rosetta Ksat, as `ksat_cm_day` (cm/day) and `ksat_in_hr` (in/hr, stormwater units)
- **Mualem–van Genuchten conductivity parameters** — `theta_r, vg_alpha_1cm, vg_n, k0_cm_day, mualem_L`, plus K at field capacity and wilting point (`k_fc_cm_day`, `k_wp_cm_day`) — carried in the CSV and used by Notebook 3

Each output table also carries an **NRCS hydrologic soil group** column, inferred from
texture class (see caveat below). Output is a 144-row DataFrame (12 texture classes × 12 bulk
densities, 0.8–1.9 g/cm³), written to `notebooks/rosetta_porosity_by_texture.csv`. Interactive plots:
porosity/FC/WP vs. BD (§6), a stacked-bar partition (§7), and a transposed FAO-style band diagram
(§8). Hydraulic conductivity and infiltration are visualized in **Notebook 3**.

Notebook 1 also exports **per-bulk-density snapshot tables** to [`data/ROSETTA/`](data/ROSETTA)
(`rosetta_bd_0.8.csv` … `rosetta_bd_1.9.csv`) — one CSV per bulk-density step, rows by
**texture (HSG)** with columns *wilting point, field capacity, saturation, available water,
drainable water* (cm³/cm³, the same six fields as the home-page slider table). The writer
(`export_water_storage_tables` in `_helpers.py`) is shared so Notebook 2 can export its blended
outputs the same way.

### Key methodologies (baseline)

- **Rosetta version 3** (2017 recalibration), called with input columns
  `[sand %, silt %, clay %, bulk_density]` → model code 3. The package returns **linear**
  van Genuchten parameters `[θr, θs, α, n, Ksat, K0, L]` (α in 1/cm, n dimensionless —
  not log₁₀, despite what some documentation implies; confirmed empirically).
- **Total porosity** = θₛ (saturated water content).
- **Field capacity** = θ at 33 kPa (330 cm) and **permanent wilting point** = θ at 1500 kPa
  (15000 cm), computed from the van Genuchten (1980) retention curve:

  $$\theta(h) = \theta_r + \frac{\theta_s - \theta_r}{\left[1 + (\alpha\,h)^{n}\right]^{m}}, \qquad m = 1 - \tfrac{1}{n}$$

  The two suction set-points are constants at the top of the relevant cell — change them
  (e.g. to 10 kPa = 100 cm for field capacity in coarse soils) as needed.
- **Median texture values** — the standard central sand/silt/clay triplet for each of the
  12 USDA classes (e.g. sand 92/5/3, loam 40/40/20, clay 20/20/60). Each triplet sums to
  100 % and plots inside the correct region of the USDA texture triangle. Defined in the
  `TEXTURE_CLASSES` dict; swap for the Levi (2017) modified-centroid values if you prefer
  that convention.
- **Bulk-density sweep** — Rosetta was trained on realistic texture + density combinations,
  so the extremes (e.g. clay at 0.8 or sand at 1.9 g/cm³) are extrapolations: plausible-looking
  but outside the well-constrained range. (The sweep is capped at **1.9 g/cm³** — BD 2.0 is
  implausible for most textures.) The **`implausible_bd`** column flags cells where
  Rosetta's θₛ exceeds the BD-implied pore space (1 − BD/2.65) — i.e. physically impossible, hence
  extrapolation. The flag propagates through the visuals: **grey dashed tails** in the §6 line plots
  (clearest symptom: the spurious **high-BD Ksat upturn** for fine textures like silt — a
  neural-network artifact, not a real rise in conductivity), a red **⚠** on flagged bars in §7, and
  **greyed texture columns** in the §8 band diagram and in Notebook 2's BD-sensitive blend diagram.
- **Hydrologic soil group (HSG)** — each row is tagged with an NRCS HSG immediately to
  the right of `texture_class`, via the `HYDROLOGIC_SOIL_GROUP` mapping:

  | HSG | Texture classes | Runoff / infiltration |
  | --- | --- | --- |
  | A | sand, loamy sand, sandy loam | low runoff / high infiltration |
  | B | loam, silt loam | moderate |
  | B/D | silt | dual group — see note |
  | C | sandy clay loam | slow infiltration |
  | D | clay loam, silty clay loam, sandy clay, silty clay, clay | high runoff / very slow infiltration |

  A **dual group** such as **B/D** is the NRCS convention for soils that are group **D** in their
  natural (undrained, seasonally-high-water-table) state but the lettered group (**B**) once
  drained; silt is classified B/D here.

  **Caveat:** HSG is *not* strictly defined by texture. The official NRCS definition
  (National Engineering Handbook, Part 630, Ch. 7) is based on saturated hydraulic
  conductivity, depth to a water-impermeable layer, and depth to the seasonal high water
  table — so the same texture can fall in different groups depending on soil structure,
  rock fragments, and profile conditions. The texture-only mapping used here is the common
  approximation applied when Ksat is unknown (e.g. HYSOGs250m, SWAT); treat it as a
  first-order estimate, not an authoritative classification.

### Sanity checks

θₛ ≥ field capacity ≥ wilting point holds for every row, and per-class values track the
published [ROSETTA centroid reference](https://ncss-tech.github.io/AQP/aqp/water-retention-curves.html).

---

## Notebook 2 — Organic-matter effects on water holding

ROSETTA has no organic-matter input — it only "sees" SOM through bulk density. This notebook adds
the organic-matter dimension on top of the ROSETTA baseline, with **available water** (FC − WP,
plant-available) and **drainable water** (SAT − FC, the fast-draining pore space relevant to
**stormwater** storage/infiltration) given equal billing.

- **Section 1 — Mineral-baseline organic carbon (UNSODA 2.0).** Estimates the organic carbon
  implicit in ROSETTA's mineral baseline from the 367 UNSODA 2.0 samples reporting both BD and OM.
  The OC-vs-BD scatter overlays **OLS regressions fit separately for topsoil, subsoil, and all
  mineral data**, with the slopes / intercepts / R² / p reported in a parameter table. The
  **all-mineral fit** (OC% = −3.16·BD + 5.43, r ≈ −0.63) becomes the blend's **BD-dependent
  baseline** `OC_base(BD)` (≈ 2.9 % at BD 0.8, ≈ 1.0 % at the mean BD ≈ 1.4, floored to 0 by
  BD ≈ 1.72); `OC_BASELINE_PCT` ≈ 1 % is kept only as the reference value at the mean BD.
- **Section 2 — ROSETTA + Minasny & McBratney (2018) blend (recommended).** Keeps ROSETTA's
  texture + bulk-density skill for the mineral baseline, then adds M&M's empirical organic-carbon
  increments (Table 2 ΔWP / ΔAWC / ΔSAT slopes, by coarse/medium/fine texture group). We apply the
  WP, AWC and SAT slopes and derive FC = WP + AWC, so the blend reproduces M&M's headline AWC
  sensitivity exactly, applied *relative to* the Section 1 BD-dependent baseline OC rather than from
  0 % OC. The FAO-style diagram has **two sliders — mineral bulk density and organic matter
  (0–8 %)** — plus dedicated AWC- and drainable-vs-OC line plots.
- **Section 3 — Saxton & Rawls (2006).** An independent, self-contained PTF taking sand/clay/OM
  directly (calibrated for OM ≤ 8 %; higher OM flagged as extrapolation), with AWC and drainable
  line plots and an OM-slider diagram. **§3.1 validates** its OC sensitivity against M&M and shows
  it runs systematically lower — and negative for clays (a known Rawls-lineage behaviour) — which
  is why Section 2 builds the blend on the M&M increments rather than Saxton–Rawls.

### Key methodologies (organic matter)

- **M&M additive increments** (per +1 % organic carbon = +10 g C kg⁻¹, mm 100 mm⁻¹ ≡ vol %):
  coarse ΔWP/ΔAWC/ΔSAT = 0.86/1.94/4.59; medium 0.68/1.79/3.59; fine 0.54/1.41/3.23. Note
  1.16 mm 100 mm⁻¹ = 0.0116 cm³/cm³, so the per-1 %-OC AWC effect is ~0.01–0.02 cm³/cm³ — modest,
  largest in sands, smallest/negative in clays, while saturation (and drainable water) rises more.
- **Mineral-baseline OC (UNSODA)** — §1 reads `data_temp/unsoda_bd_om.csv` (367 UNSODA 2.0
  samples reporting both BD and OM): mineral all-horizon mean ≈ 0.9 % OC, topsoil median ≈ 1.1 %,
  with OC declining as BD rises (r ≈ −0.6). The **all-mineral OC~BD regression** sets the
  BD-dependent baseline `oc_baseline_for_bd(BD)` (constants `OC_BD_SLOPE`/`OC_BD_INTERCEPT` in
  `_helpers.py`, which §1 re-derives from UNSODA and asserts still match).
- **Baseline & sliders** — ROSETTA's prediction is a *nominal* baseline at the **BD-dependent**
  `OC_base(BD)` (not OC = 0); M&M increments are applied relative to it (values floored at ≥ 0), so
  the slider's 0 % end is a truly organic-free mineral soil (drier than ROSETTA), and the firebrick
  reference line marking `OC_base(BD)` moves with the BD slider. The
  BD and OM sliders are *independent what-if axes*; since organic matter physically lowers bulk
  density, the realistic region runs low-BD ↔ high-OM, and the extreme corner double-counts
  porosity. The modifier is linear (M&M found diminishing returns), so it may overstate gains at
  high OM; OM ≈ OC / 0.58 (van Bemmelen).

### Data & reproducibility (UNSODA)

The §1 mineral-baseline estimate uses **UNSODA 2.0** (Nemes et al., 2001), distributed as a
Microsoft Access database. [`notebooks/fetch_unsoda.py`](notebooks/fetch_unsoda.py) downloads it,
exports the `soil_properties` and `general` tables, and writes the tidy extract used by the
notebook:

```bash
pixi run python notebooks/fetch_unsoda.py    # -> notebooks/data_temp/unsoda_bd_om.csv
```

- `notebooks/data_temp/` is **git-ignored** scratch (the UNSODA `.mdb` and derived CSV live there).
- Requires **`mdbtools`** on PATH to read the `.mdb` (not on conda-forge for osx-arm64; macOS:
  `brew install mdbtools`).
- If the extract is absent, Notebook 2 still runs — §1 skips its scatter plot and consistency
  check, and the blend falls back to the `OC_BD_SLOPE`/`OC_BD_INTERCEPT` constants in `_helpers.py`.

---

## Notebook 3 — Hydraulic conductivity & infiltration

Reads `rosetta_porosity_by_texture.csv` from Notebook 1 (which carries Rosetta's Ksat and the
Mualem–van Genuchten parameters) and visualizes conductivity and infiltration in stormwater units,
with bulk-density sliders and the same `implausible_bd` grey-dashed extrapolation flagging:

- **§1 Saturated K (Ksat)** — vs. bulk density on a log axis (in/hr).
- **§2 Unsaturated K(h)** — Rosetta's full **Mualem–van Genuchten** curve,
  K(h) = K0·Se^L·[1 − (1 − Se^(1/m))^m]² with Se(h) = [1 + (αh)^n]^(−m); uses `K0`/`L`
  (columns 5–6), *not* Ksat. Log–log, with field-capacity / wilting-point markers.
- **§3 Green–Ampt infiltration** — f(t) and cumulative F(t). Uses Rosetta's Ksat and moisture
  deficit (θₛ − θ_wp) with the **Rawls, Brakensiek & Miller (1983)** textural wetting-front
  suction ψ_f; deriving ψ_f from Rosetta's K(h) integral proved unreliable (its fitted `L`
  mis-orders the capillary drive — sand above clay), so the published textural values are used.

Caveats: Ksat/K0 are Rosetta's least-certain outputs, and the model is **matrix-only** (no
macropores/structure), so field infiltration is often 1–2 orders of magnitude higher.

---

## References

- Zhang & Schaap (2017), [rosetta-soil](https://github.com/usda-ars-ussl/rosetta-soil)
- van Genuchten (1980), *SSSAJ* 44:892–898
- Levi (2017), [Modified Centroid for Estimating Sand, Silt, and Clay from Soil Texture Class](https://acsess.onlinelibrary.wiley.com/doi/abs/10.2136/sssaj2016.09.0301), *SSSAJ*
- USDA-NRCS, [National Engineering Handbook, Part 630, Chapter 7 — Hydrologic Soil Groups](https://directives.nrcs.usda.gov/sites/default/files2/1712930597/11905.pdf)
- Saxton & Rawls (2006), [Soil Water Characteristic Estimates by Texture and Organic Matter for Hydrologic Solutions](https://www.researchgate.net/publication/43257423_Soil_Water_Characteristic_Estimates_by_Texture_and_Organic_Matter_for_Hydrologic_Solutions), *SSSAJ* 70:1569–1578
- Minasny & McBratney (2018), [Limited effect of organic matter on soil available water capacity](https://bsssjournals.onlinelibrary.wiley.com/doi/abs/10.1111/ejss.12475), *EJSS* 69:39–47
- Schaap & Leij (2000), Improved prediction of unsaturated hydraulic conductivity with the Mualem–van Genuchten model, *SSSAJ* 64:843–851
- Rawls, Brakensiek & Miller (1983), Green-Ampt infiltration parameters from soils data, *J. Hydraulic Engineering* 109:62–70
- Nemes, Schaap, Leij & Wösten (2001), [Description of the unsaturated soil hydraulic database UNSODA version 2.0](https://www.sciencedirect.com/science/article/abs/pii/S0022169401004656), *J. Hydrology* 251:151–162 — data: [UNSODA 2.0 on Ag Data Commons](https://agdatacommons.nal.usda.gov/articles/dataset/24851832)

---

## Developers & Contributors

### Environment (pixi)

Built with [pixi](https://pixi.sh) — see [pixi.toml](pixi.toml).

```bash
pixi install
pixi run jupyter lab    # open notebooks/ and run Notebook 1 first, then 2 and 3
```

Each notebook is paired with a [jupytext](https://jupytext.readthedocs.io/) `py:percent`
script (`*.py` alongside the `*.ipynb`) — the diff-friendly version to review and commit. After
editing either file, run e.g. `pixi run jupytext --sync notebooks/<name>.py` to keep the pair in
sync (the `.py` holds code/markdown; the `.ipynb` holds outputs). Plots are embedded as
self-contained HoloViews/Bokeh output, so the saved notebooks are interactive without a live
kernel.

Functions shared by more than one notebook live in
[`notebooks/_helpers.py`](notebooks/_helpers.py) — the scrollable-table display helper, the
van Genuchten–Mualem retention/conductivity functions, the extrapolation-aware bulk-density line
plot, and the ROSETTA + Minasny–McBratney organic-matter blend (`soil_water_bd_om_blend_table`),
its FAO band-diagram builder (`soil_water_texture_band_diagram`), and the storage-table HTML
(`soil_water_table_html`) used by the home page and Notebook 2. It is a plain importable module
(not a paired notebook); the notebooks `import` it at the top and run in `notebooks/`, so it
resolves both in Jupyter and during the Quarto build.

### Interactive website (GitHub Pages)

The notebooks are published as a static, **fully interactive** website — readers can move the
bulk-density / organic-matter sliders, hover, pan, and zoom right in the browser, with no install:

**<https://limnotech.github.io/soil-health-hydraulics/>**

#### How it works

The site is built with **[Quarto](https://quarto.org)** ([`_quarto.yml`](_quarto.yml),
[`index.qmd`](index.qmd)), with code folded by default. The home page (`index.qmd`) leads with a
**"Soil health at a glance"** summary — the headline two-slider blend figure and a slider-linked
storage table, composed with [Panel](https://panel.holoviz.org) (`embed=True`) so every slider
state works in static HTML. `_quarto.yml` lists the render inputs explicitly so the internal
`docs/` planning folder stays out of the published site.

- Build locally with **`pixi run render`** (or `pixi run preview` for a live server); output goes
  to `_site/`.
- Quarto **executes** the notebooks (run Notebook 1 first so its CSV exists) and **freezes** the
  results into [`_freeze/`](_freeze) — which **is committed**. When the cache matches, CI renders
  *from it* without re-executing (fast); otherwise CI re-executes (the committed
  `notebooks/data_temp/unsoda_bd_om.csv` lets Notebook 2 run fully). Either way it **deploys** — a
  stale cache is a non-fatal warning, never a blocked deploy.
- **Re-render and commit `_freeze/` after editing a notebook** to keep CI fast.
- [`.github/workflows/publish.yml`](.github/workflows/publish.yml) renders and deploys on every
  push to `main` using the same pixi environment.

> [!IMPORTANT]
> Quarto must **execute or freeze** the notebooks — that is what preserves the interactive
> HoloViews/Bokeh embeds (slider/hover/zoom). Rendering *pre-executed* `.ipynb` with execution
> disabled silently strips those embeds (dead plots). `execute-dir: file` is set so each notebook
> runs in `notebooks/` and its relative CSV read resolves.
