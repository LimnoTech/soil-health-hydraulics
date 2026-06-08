# soil_modeling

Estimating soil water-holding properties for every USDA texture class from pedotransfer
functions, with interactive HoloViews/Bokeh visualizations. All notebooks live in
[notebooks/](notebooks/).

## Notebooks

| Notebook | Purpose |
| --- | --- |
| [notebooks/1_rosetta_porosity_by_texture.ipynb](notebooks/1_rosetta_porosity_by_texture.ipynb) | **ROSETTA texture × bulk-density baseline** — total porosity, field capacity, wilting point, available & drainable water for each texture class across bulk densities 0.8–2.0 g/cm³. Writes `rosetta_porosity_by_texture.csv`. |
| [notebooks/2_organic_matter_water_holding.ipynb](notebooks/2_organic_matter_water_holding.ipynb) | **Organic-matter effects** — layers SOM/SOC sensitivity on top, two ways: a ROSETTA + Minasny–McBratney (2018) blend (anchored to a UNSODA-derived mineral-baseline OC), and Saxton–Rawls (2006). **Reads the CSV from the first notebook.** |

**Run order:** execute `1_rosetta_porosity_by_texture.ipynb` first (it generates
`rosetta_porosity_by_texture.csv`), then `2_organic_matter_water_holding.ipynb`, which loads that
CSV for the mineral baseline.

### Environment (pixi)

Built with [pixi](https://pixi.sh) — see [pixi.toml](pixi.toml).

```bash
pixi install
pixi run jupyter lab    # open notebooks/ and run the two notebooks in order
```

Each notebook is paired with a [jupytext](https://jupytext.readthedocs.io/) `py:percent`
script (`*.py` alongside the `*.ipynb`) — the diff-friendly version to review and commit. After
editing either file, run e.g. `pixi run jupytext --sync notebooks/<name>.py` to keep the pair in
sync (the `.py` holds code/markdown; the `.ipynb` holds outputs). Plots are embedded as
self-contained HoloViews/Bokeh output, so the saved notebooks are interactive without a live
kernel.

---

## Notebook 1 — ROSETTA porosity by texture and bulk density

Uses the [Rosetta](https://github.com/usda-ars-ussl/rosetta-soil) pedotransfer functions to
estimate, for every USDA soil texture class (using representative *median* sand/silt/clay values)
and for bulk densities from 0.8 to 2.0 g/cm³ in 0.1 steps:

- **Total porosity** — saturated volumetric water content θₛ
- **Field-capacity porosity** — volumetric water content at 33 kPa (330 cm suction)
- **Permanent-wilting-point porosity** — volumetric water content at 1500 kPa (15000 cm suction)
- **Available water capacity** — field capacity minus wilting point
- **Drainable water** — saturation minus field capacity (shown in the band plots)

Each output table also carries an **NRCS hydrologic soil group** (A–D) column, inferred from
texture class (see caveat below). Output is a 156-row DataFrame (12 texture classes × 13 bulk
densities), written to `notebooks/rosetta_porosity_by_texture.csv`. Sections 6–8 add interactive
plots (porosity/FC/WP vs. BD; a stacked-bar partition with a BD slider; and a transposed
FAO-style line/area diagram in inches of water per foot of soil).

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
  so the extremes (e.g. clay at 0.8 or sand at 2.0 g/cm³) are extrapolations: plausible-looking
  but outside the well-constrained range.
- **Hydrologic soil group (HSG)** — each row is tagged with an NRCS HSG (A–D) immediately to
  the right of `texture_class`, via the `HYDROLOGIC_SOIL_GROUP` mapping:

  | HSG | Texture classes | Runoff / infiltration |
  | --- | --- | --- |
  | A | sand, loamy sand, sandy loam | low runoff / high infiltration |
  | B | loam, silt loam, silt | moderate |
  | C | sandy clay loam | slow infiltration |
  | D | clay loam, silty clay loam, sandy clay, silty clay, clay | high runoff / very slow infiltration |

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

- **Section 1 — ROSETTA + Minasny & McBratney (2018) blend (recommended).** Keeps ROSETTA's
  texture + bulk-density skill for the mineral baseline, then adds M&M's empirical organic-carbon
  increments (Table 2 ΔWP / ΔAWC / ΔSAT slopes, by coarse/medium/fine texture group). We apply the
  WP, AWC and SAT slopes and derive FC = WP + AWC, so the blend reproduces M&M's headline AWC
  sensitivity exactly. **§1.1** estimates the mineral-baseline organic carbon from UNSODA 2.0
  (≈ 1 % OC; OC declines with BD, r ≈ −0.6) and the increments are applied *relative to that*
  baseline rather than from 0 % OC. The FAO-style diagram has **two sliders — mineral bulk density
  and organic matter (0–8 %)** — plus dedicated AWC- and drainable-vs-OC line plots.
- **Section 2 — Saxton & Rawls (2006).** An independent, self-contained PTF taking sand/clay/OM
  directly (calibrated for OM ≤ 8 %; higher OM flagged as extrapolation), with AWC and drainable
  line plots and an OM-slider diagram. **§2.1 validates** its OC sensitivity against M&M and shows
  it runs systematically lower — and negative for clays (a known Rawls-lineage behaviour) — which
  is why Section 1 builds the blend on the M&M increments rather than Saxton–Rawls.

### Key methodologies (organic matter)

- **M&M additive increments** (per +1 % organic carbon = +10 g C kg⁻¹, mm 100 mm⁻¹ ≡ vol %):
  coarse ΔWP/ΔAWC/ΔSAT = 0.86/1.94/4.59; medium 0.68/1.79/3.59; fine 0.54/1.41/3.23. Note
  1.16 mm 100 mm⁻¹ = 0.0116 cm³/cm³, so the per-1 %-OC AWC effect is ~0.01–0.02 cm³/cm³ — modest,
  largest in sands, smallest/negative in clays, while saturation (and drainable water) rises more.
- **Mineral-baseline OC (UNSODA)** — §1.1 reads `data_temp/unsoda_bd_om.csv` (367 UNSODA 2.0
  samples reporting both BD and OM) to justify the anchor: mineral all-horizon mean ≈ 0.9 % OC,
  topsoil median ≈ 1.1 %, with OC declining as BD rises (r ≈ −0.6). Default `OC_BASELINE_PCT = 1.0`.
- **Baseline & sliders** — ROSETTA's prediction is a *nominal* baseline at OC ≈ `OC_BASELINE_PCT`
  (not OC = 0); M&M increments are applied relative to it (values floored at ≥ 0), so the slider's
  0 % end is a truly organic-free mineral soil (drier than ROSETTA) and OC ≈ 1 % reproduces it. The
  BD and OM sliders are *independent what-if axes*; since organic matter physically lowers bulk
  density, the realistic region runs low-BD ↔ high-OM, and the extreme corner double-counts
  porosity. The modifier is linear (M&M found diminishing returns), so it may overstate gains at
  high OM; OM ≈ OC / 0.58 (van Bemmelen).

### Data & reproducibility (UNSODA)

The §1.1 mineral-baseline estimate uses **UNSODA 2.0** (Nemes et al., 2001), distributed as a
Microsoft Access database. [`notebooks/fetch_unsoda.py`](notebooks/fetch_unsoda.py) downloads it,
exports the `soil_properties` and `general` tables, and writes the tidy extract used by the
notebook:

```bash
pixi run python notebooks/fetch_unsoda.py    # -> notebooks/data_temp/unsoda_bd_om.csv
```

- `notebooks/data_temp/` is **git-ignored** scratch (the UNSODA `.mdb` and derived CSV live there).
- Requires **`mdbtools`** on PATH to read the `.mdb` (not on conda-forge for osx-arm64; macOS:
  `brew install mdbtools`).
- If the extract is absent, Notebook 2 still runs — §1.1 skips its scatter plot and falls back to
  the documented default `OC_BASELINE_PCT`.

---

## References

- Zhang & Schaap (2017), [rosetta-soil](https://github.com/usda-ars-ussl/rosetta-soil)
- van Genuchten (1980), *SSSAJ* 44:892–898
- Levi (2017), [Modified Centroid for Estimating Sand, Silt, and Clay from Soil Texture Class](https://acsess.onlinelibrary.wiley.com/doi/abs/10.2136/sssaj2016.09.0301), *SSSAJ*
- USDA-NRCS, [National Engineering Handbook, Part 630, Chapter 7 — Hydrologic Soil Groups](https://directives.nrcs.usda.gov/sites/default/files2/1712930597/11905.pdf)
- Saxton & Rawls (2006), [Soil Water Characteristic Estimates by Texture and Organic Matter for Hydrologic Solutions](https://www.researchgate.net/publication/43257423_Soil_Water_Characteristic_Estimates_by_Texture_and_Organic_Matter_for_Hydrologic_Solutions), *SSSAJ* 70:1569–1578
- Minasny & McBratney (2018), [Limited effect of organic matter on soil available water capacity](https://bsssjournals.onlinelibrary.wiley.com/doi/abs/10.1111/ejss.12475), *EJSS* 69:39–47
- Nemes, Schaap, Leij & Wösten (2001), [Description of the unsaturated soil hydraulic database UNSODA version 2.0](https://www.sciencedirect.com/science/article/abs/pii/S0022169401004656), *J. Hydrology* 251:151–162 — data: [UNSODA 2.0 on Ag Data Commons](https://agdatacommons.nal.usda.gov/articles/dataset/24851832)
