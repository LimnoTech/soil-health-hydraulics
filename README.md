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

Output is a 156-row pandas DataFrame (12 texture classes × 13 bulk densities), also written
to [rosetta_porosity_by_texture.csv](rosetta_porosity_by_texture.csv).

### Environment (pixi)

Built with [pixi](https://pixi.sh) — see [pixi.toml](pixi.toml).

```bash
pixi install
pixi run jupyter lab    # open and run rosetta_porosity_by_texture.ipynb
```

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

### Sanity checks

θₛ ≥ field capacity ≥ wilting point holds for every row, and per-class values track the
published [ROSETTA centroid reference](https://ncss-tech.github.io/AQP/aqp/water-retention-curves.html).

### References

- Zhang & Schaap (2017), [rosetta-soil](https://github.com/usda-ars-ussl/rosetta-soil)
- van Genuchten (1980), *SSSAJ* 44:892–898
- Levi (2017), [Modified Centroid for Estimating Sand, Silt, and Clay from Soil Texture Class](https://acsess.onlinelibrary.wiley.com/doi/abs/10.2136/sssaj2016.09.0301), *SSSAJ*
