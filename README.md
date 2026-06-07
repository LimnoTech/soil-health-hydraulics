# soil_modeling

## Soil porosity by USDA texture class and bulk density (Rosetta)

[rosetta_porosity_by_texture.ipynb](rosetta_porosity_by_texture.ipynb) uses the
[Rosetta](https://github.com/usda-ars-ussl/rosetta-soil) pedotransfer functions to
estimate, for every USDA soil texture class (using representative *median* sand/silt/clay
values) and for bulk densities from 0.8 to 2.0 g/cm³ in 0.1 steps:

- **Total porosity** — saturated volumetric water content θₛ
- **Field-capacity porosity** — volumetric water content at 33 kPa (330 cm suction)
- **Permanent-wilting-point porosity** — volumetric water content at 1500 kPa (15000 cm suction)
- **Available water capacity** — field capacity minus wilting point

Each output table also carries an **NRCS hydrologic soil group** (A–D) column, inferred from
texture class (see caveat below).

Output is a 156-row pandas DataFrame (12 texture classes × 13 bulk densities), also written
to [rosetta_porosity_by_texture.csv](rosetta_porosity_by_texture.csv).

### Environment (pixi)

Built with [pixi](https://pixi.sh) — see [pixi.toml](pixi.toml).

```bash
pixi install
pixi run jupyter lab    # open and run rosetta_porosity_by_texture.ipynb
```

The notebook is paired with a [jupytext](https://jupytext.readthedocs.io/) `py:percent`
script, [rosetta_porosity_by_texture.py](rosetta_porosity_by_texture.py), which is the
diff-friendly version to review and commit. Editing either file and running
`pixi run jupytext --sync rosetta_porosity_by_texture.py` keeps the two in sync (the `.py`
holds the code/markdown; the `.ipynb` holds outputs).

### Key methodologies

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

### References

- Zhang & Schaap (2017), [rosetta-soil](https://github.com/usda-ars-ussl/rosetta-soil)
- van Genuchten (1980), *SSSAJ* 44:892–898
- Levi (2017), [Modified Centroid for Estimating Sand, Silt, and Clay from Soil Texture Class](https://acsess.onlinelibrary.wiley.com/doi/abs/10.2136/sssaj2016.09.0301), *SSSAJ*
- USDA-NRCS, [National Engineering Handbook, Part 630, Chapter 7 — Hydrologic Soil Groups](https://directives.nrcs.usda.gov/sites/default/files2/1712930597/11905.pdf)
