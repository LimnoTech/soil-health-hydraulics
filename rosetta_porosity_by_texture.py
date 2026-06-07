# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Soil porosity by USDA texture class and bulk density (Rosetta)
#
# This notebook uses the **Rosetta** pedotransfer functions ([usda-ars-ussl/rosetta-soil](https://github.com/usda-ars-ussl/rosetta-soil)) to estimate, for **every USDA soil texture class** (using representative *median* sand/silt/clay values) and for **bulk densities from 0.8 to 2.0 g/cm³ in 0.1 steps**:
#
# | Output | Definition |
# |---|---|
# | **Total porosity** | Saturated volumetric water content θₛ (van Genuchten) |
# | **Field-capacity porosity** | Volumetric water content at 33 kPa (= 330 cm suction) |
# | **Permanent-wilting-point porosity** | Volumetric water content at 1500 kPa (= 15000 cm suction) |
#
# Rosetta predicts the **van Genuchten** water-retention parameters (θᵣ, θₛ, α, n) from texture + bulk density. We then evaluate the retention curve at the field-capacity and wilting-point suctions.
#
# **van Genuchten (1980) retention model:**
# $$\theta(h) = \theta_r + \frac{\theta_s - \theta_r}{\left[1 + (\alpha\,h)^{n}\right]^{m}}, \qquad m = 1 - \tfrac{1}{n}$$
#
# where $h$ is the suction head [cm] and $\alpha$ [1/cm], $n$ [-], $\theta_r$, $\theta_s$ are Rosetta outputs.
#
# **Environment:** built with [pixi](https://pixi.sh) (`pixi.toml`). Run `pixi install`, then open this notebook with the `soil_modeling` kernel (`pixi run jupyter lab`).

# %%
import numpy as np
import pandas as pd
from rosetta import rosetta, SoilData

pd.set_option("display.float_format", lambda v: f"{v:0.3f}")

# %% [markdown]
# ## 1. Representative (median) texture values for each USDA class
#
# Sand / silt / clay percentages below are the widely used central (representative) values for each of the 12 USDA texture classes. Each triplet sums to 100% and plots inside the correct region of the USDA texture triangle. Adjust these if you prefer a different convention (e.g. the geometric or modified centroids of Levi, 2017, SSSAJ).

# %%
# texture_class: (sand %, silt %, clay %)
TEXTURE_CLASSES = {
    "sand":            (92,  5,  3),
    "loamy sand":      (82, 12,  6),
    "sandy loam":      (65, 25, 10),
    "loam":            (40, 40, 20),
    "silt loam":       (20, 65, 15),
    "silt":            ( 7, 88,  5),
    "sandy clay loam": (60, 13, 27),
    "clay loam":       (30, 35, 35),
    "silty clay loam": (10, 56, 34),
    "sandy clay":      (50,  7, 43),
    "silty clay":      ( 7, 47, 46),
    "clay":            (20, 20, 60),
}

# NRCS hydrologic soil group (HSG) inferred from texture class. This is the standard
# texture-based approximation used when saturated hydraulic conductivity is unknown
# (e.g. HYSOGs250m, SWAT). Actual HSG also depends on Ksat, depth to a restrictive
# layer, and water-table depth. A = low runoff / high infiltration ... D = high runoff.
HYDROLOGIC_SOIL_GROUP = {
    "sand":            "A",
    "loamy sand":      "A",
    "sandy loam":      "A",
    "loam":            "B",
    "silt loam":       "B",
    "silt":            "B",
    "sandy clay loam": "C",
    "clay loam":       "D",
    "silty clay loam": "D",
    "sandy clay":      "D",
    "silty clay":      "D",
    "clay":            "D",
}

texture_df = (
    pd.DataFrame.from_dict(
        TEXTURE_CLASSES, orient="index", columns=["sand_pct", "silt_pct", "clay_pct"]
    )
    .rename_axis("texture_class")
    .reset_index()
)
# hydrologic soil group immediately to the right of texture_class
texture_df.insert(1, "hydrologic_soil_group", texture_df["texture_class"].map(HYDROLOGIC_SOIL_GROUP))

# sanity check: each class sums to 100%
assert (texture_df[["sand_pct", "silt_pct", "clay_pct"]].sum(axis=1) == 100).all()
texture_df

# %% [markdown]
# ## 2. Bulk-density range and suction set-points

# %%
# Bulk density 0.8 -> 2.0 g/cm3 in 0.1 steps (rounded to avoid float drift)
bulk_densities = np.round(np.arange(0.8, 2.0 + 1e-9, 0.1), 1)

# Suction heads (cm) at which to evaluate the retention curve
H_FIELD_CAPACITY = 330.0    # 33 kPa
H_WILTING_POINT = 15000.0   # 1500 kPa

print(f"{len(bulk_densities)} bulk densities: {bulk_densities}")


# %% [markdown]
# ## 3. van Genuchten retention helper

# %%
def vg_theta(h, theta_r, theta_s, alpha, n):
    """Volumetric water content at suction head h [cm] (van Genuchten, 1980).

    Parameters may be scalars or numpy arrays (broadcast together).
    alpha [1/cm], n [-], theta_r/theta_s [cm3/cm3]. m = 1 - 1/n.
    """
    theta_r = np.asarray(theta_r, dtype=float)
    theta_s = np.asarray(theta_s, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    n = np.asarray(n, dtype=float)
    m = 1.0 - 1.0 / n
    return theta_r + (theta_s - theta_r) / (1.0 + (alpha * h) ** n) ** m


# %% [markdown]
# ## 4. Run Rosetta and build the results DataFrame
#
# Rosetta input columns are `[sand%, silt%, clay%, bulk_density]`; model version **3** (2017 recalibration) is used. The mean output columns are `[θᵣ, θₛ, α, n, Ksat, K0, L]` with α in 1/cm and n dimensionless (linear, not log₁₀).

# %%
# Build every (texture class) x (bulk density) combination
rows = [
    {
        "texture_class": cls,
        "hydrologic_soil_group": HYDROLOGIC_SOIL_GROUP[cls],
        "sand_pct": sand,
        "silt_pct": silt,
        "clay_pct": clay,
        "bulk_density_g_cm3": bd,
    }
    for cls, (sand, silt, clay) in TEXTURE_CLASSES.items()
    for bd in bulk_densities
]
df = pd.DataFrame(rows)

# Rosetta expects [sand, silt, clay, bulk_density]
rosetta_input = df[["sand_pct", "silt_pct", "clay_pct", "bulk_density_g_cm3"]].to_numpy()
mean, stdev, codes = rosetta(3, SoilData.from_iter(rosetta_input.tolist()))
mean = np.asarray(mean)

theta_r, theta_s, alpha, n = mean[:, 0], mean[:, 1], mean[:, 2], mean[:, 3]

df["total_porosity"] = theta_s
df["field_capacity_porosity"] = vg_theta(H_FIELD_CAPACITY, theta_r, theta_s, alpha, n)
df["wilting_point_porosity"] = vg_theta(H_WILTING_POINT, theta_r, theta_s, alpha, n)
df["available_water_capacity"] = (
    df["field_capacity_porosity"] - df["wilting_point_porosity"]
)
df["rosetta_code"] = np.asarray(codes)  # 3 = texture + bulk density model

print(f"{len(df)} rows  ({len(TEXTURE_CLASSES)} classes x {len(bulk_densities)} bulk densities)")
df.head(15)

# %% [markdown]
# ## 5. Results

# %%
result = df[[
    "texture_class",
    "hydrologic_soil_group",
    "sand_pct",
    "silt_pct",
    "clay_pct",
    "bulk_density_g_cm3",
    "total_porosity",
    "field_capacity_porosity",
    "wilting_point_porosity",
    "available_water_capacity",
]].copy()

result

# %%
# Save to CSV
result.to_csv("rosetta_porosity_by_texture.csv", index=False)
print("Wrote rosetta_porosity_by_texture.csv")

# %% [markdown]
# ## 6. Quick visualization (optional)
#
# Interactive [hvPlot](https://hvplot.holoviz.org/) / [HoloViews](https://holoviews.org/) (Bokeh) line plots of total porosity, field-capacity porosity, and wilting-point porosity vs. bulk density, one line per texture class. Hover for values; use the toolbar to pan/zoom.

# %%
import hvplot.pandas  # noqa: F401  (registers the .hvplot accessor)

plot_opts = dict(
    x="bulk_density_g_cm3",
    by="texture_class",
    xlabel="Bulk density (g/cm³)",
    width=750,
    height=500,
    legend="right",
    grid=True,
)

result.hvplot.line(
    y="total_porosity",
    ylabel="Total porosity θₛ (cm³/cm³)",
    title="Rosetta total porosity vs. bulk density by USDA texture class",
    **plot_opts,
)

# %%
result.hvplot.line(
    y="field_capacity_porosity",
    ylabel="Field-capacity porosity at 33 kPa (cm³/cm³)",
    title="Rosetta field-capacity porosity vs. bulk density by USDA texture class",
    **plot_opts,
)

# %%
result.hvplot.line(
    y="wilting_point_porosity",
    ylabel="Wilting-point porosity at 1500 kPa (cm³/cm³)",
    title="Rosetta wilting-point porosity vs. bulk density by USDA texture class",
    **plot_opts,
)
