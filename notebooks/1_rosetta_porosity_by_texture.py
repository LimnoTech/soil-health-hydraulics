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
#     display_name: default
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Soil porosity by USDA texture class and bulk density (Rosetta)
#
# This notebook produces **interactive charts** showing how soil texture, compaction (bulk density), and organic matter change a soil's capacity to store and release water. Charts cover all 12 USDA texture classes across bulk densities from 0.8 to 1.9 g/cm³ — useful for stormwater design, soil-health assessment, and land-management decisions.
#
# Using the **Rosetta** pedotransfer functions ([usda-ars-ussl/rosetta-soil](https://github.com/usda-ars-ussl/rosetta-soil)), we estimate for each texture class and bulk density:
#
# | Output | Definition |
# |---|---|
# | **Total porosity** | Saturated volumetric water content θₛ (van Genuchten) |
# | **Field-capacity porosity** | Volumetric water content at 33 kPa (= 330 cm suction) |
# | **Permanent-wilting-point porosity** | Volumetric water content at 1500 kPa (= 15000 cm suction) |
# | **Saturated hydraulic conductivity** | Ksat (cm/day) — Rosetta's conductivity output |
#
# ::: {.callout-note collapse="true"}
# ## For researchers: van Genuchten–Mualem retention model
#
# Rosetta predicts the **van Genuchten** water-retention parameters (θᵣ, θₛ, α, n) from texture + bulk density. We then evaluate the retention curve at the field-capacity and wilting-point suctions.
#
# **van Genuchten (1980) retention model:**
# $$\theta(h) = \theta_r + \frac{\theta_s - \theta_r}{\left[1 + (\alpha\,h)^{n}\right]^{m}}, \qquad m = 1 - \tfrac{1}{n}$$
#
# where $h$ is the suction head [cm] and $\alpha$ [1/cm], $n$ [-], $\theta_r$, $\theta_s$ are Rosetta outputs.
# :::
#
# **Environment:** built with [pixi](https://pixi.sh) (`pixi.toml`). Run `pixi install`, then open this notebook with the `soil_modeling` kernel (`pixi run jupyter lab`).

# %%
import numpy as np
import pandas as pd
from rosetta import rosetta, SoilData

# Shared helpers (display table, van Genuchten–Mualem functions, extrapolation-aware line plot);
# importing _helpers also sets the shared pandas float_format. See notebooks/_helpers.py.
from _helpers import show, vg_theta, mualem_k, line_with_extrapolation, soil_water_texture_band_diagram

# %% [markdown]
# ## 1. Representative (median) texture values for each USDA class
#
# Each texture class is represented by a single sand/silt/clay triplet — the widely used central (representative) values. Each triplet sums to 100% and falls inside the correct region of the USDA texture triangle. Adjust these if you prefer a different convention (e.g. the geometric or modified centroids of Levi, 2017, SSSAJ).

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
show(texture_df)

# %% [markdown]
# ## 2. Bulk-density range and suction set-points
#
# Bulk density (g/cm³) is a direct measure of soil compaction: lower values indicate more pore space; higher values reflect a more compacted, denser soil. Field capacity (33 kPa) and permanent wilting point (1500 kPa) are the standard agronomic thresholds that bound plant-available water.

# %%
# Bulk density 0.8 -> 1.9 g/cm3 in 0.1 steps (rounded to avoid float drift). Capped at 1.9:
# BD 2.0 is physically implausible (θₛ > pore space) for most textures — see implausible_bd.
bulk_densities = np.round(np.arange(0.8, 1.9 + 1e-9, 0.1), 1)

# Suction heads (cm) at which to evaluate the retention curve
H_FIELD_CAPACITY = 330.0    # 33 kPa
H_WILTING_POINT = 15000.0   # 1500 kPa

print(f"{len(bulk_densities)} bulk densities: {bulk_densities}")


# %% [markdown]
# ## 3. van Genuchten–Mualem helpers
#
# ::: {.callout-note collapse="true"}
# ## For researchers: van Genuchten–Mualem retention model
#
# `vg_theta` gives the retention curve θ(h); `mualem_k` gives the **unsaturated** hydraulic
# conductivity K(h) from Rosetta's Mualem–van Genuchten parameters (matching-point `K0` and
# pore-connectivity `L`, columns 5–6 of the Rosetta output). Used for K at field capacity /
# wilting point (§4) and carried in the CSV for the hydraulic-conductivity notebook (Notebook 3).
# Both functions are shared with Notebook 3 and live in **`notebooks/_helpers.py`** (imported above).
# :::

# %% [markdown]
# ## 4. Run Rosetta and build the results DataFrame
#
# ::: {.callout-note collapse="true"}
# ## For researchers: how Rosetta is run
#
# Rosetta input columns are `[sand%, silt%, clay%, bulk_density]`; model version **3** (2017 recalibration) is used. The mean output columns are `[θᵣ, θₛ, α, n, Ksat, K0, L]` with α in 1/cm and n dimensionless (linear, not log₁₀). We keep **Ksat** (saturated hydraulic conductivity, column 4) as `ksat_cm_day` (cm/day) and `ksat_in_hr` (in/hr = cm/day ÷ 24 ÷ 2.54, the usual stormwater unit). We also retain the **Mualem–van Genuchten** parameters `theta_r, vg_alpha_1cm, vg_n, k0_cm_day, mualem_L` (used by the hydraulic-conductivity notebook, **Notebook 3**) and the unsaturated K at field capacity / wilting point (`k_fc_cm_day`, `k_wp_cm_day`).
# :::

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
ksat = mean[:, 4]  # saturated hydraulic conductivity, cm/day (linear, not log10)
k0 = mean[:, 5]    # Mualem matching-point conductivity, cm/day
Lpar = mean[:, 6]  # Mualem pore-connectivity exponent, dimensionless (often negative)

df["total_porosity"] = theta_s
df["field_capacity_porosity"] = vg_theta(H_FIELD_CAPACITY, theta_r, theta_s, alpha, n)
df["wilting_point_porosity"] = vg_theta(H_WILTING_POINT, theta_r, theta_s, alpha, n)
df["available_water_capacity"] = (
    df["field_capacity_porosity"] - df["wilting_point_porosity"]
)
df["ksat_cm_day"] = ksat  # Rosetta saturated hydraulic conductivity (cm/day)
df["ksat_in_hr"] = ksat / (24.0 * 2.54)  # same Ksat in inches/hour (stormwater units)
# van Genuchten–Mualem parameters kept for the unsaturated-K and infiltration sections (§10–§11)
df["theta_r"] = theta_r
df["vg_alpha_1cm"] = alpha
df["vg_n"] = n
df["k0_cm_day"] = k0
df["mualem_L"] = Lpar
# unsaturated K at field capacity (33 kPa) and wilting point (1500 kPa), cm/day
df["k_fc_cm_day"] = mualem_k(H_FIELD_CAPACITY, alpha, n, k0, Lpar)
df["k_wp_cm_day"] = mualem_k(H_WILTING_POINT, alpha, n, k0, Lpar)
df["rosetta_code"] = np.asarray(codes)  # 3 = texture + bulk density model

# Physical-plausibility flag. Saturated water content θₛ cannot exceed the pore space
# implied by the bulk density (1 - BD/ρ_particle, ρ_particle ≈ 2.65 g/cm³). Where Rosetta's
# θₛ exceeds it, the texture × BD combination is an extrapolation outside Rosetta's training
# domain (e.g. high-BD silt, which drives the spurious Ksat upturn). True = implausible.
PARTICLE_DENSITY = 2.65  # g/cm³, quartz-dominated mineral soil
df["implausible_bd"] = df["total_porosity"] > (1.0 - df["bulk_density_g_cm3"] / PARTICLE_DENSITY)

print(f"{len(df)} rows  ({len(TEXTURE_CLASSES)} classes x {len(bulk_densities)} bulk densities)")
show(df)

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
    "ksat_cm_day",
    "ksat_in_hr",
    "k_fc_cm_day",
    "k_wp_cm_day",
    "theta_r",
    "vg_alpha_1cm",
    "vg_n",
    "k0_cm_day",
    "mualem_L",
    "implausible_bd",
]].copy()

show(result)

# %%
# Save to CSV
result.to_csv("rosetta_porosity_by_texture.csv", index=False)
print("Wrote rosetta_porosity_by_texture.csv")

# %% [markdown]
# ## 6. Quick visualization of Water Storage
#
# Interactive [hvPlot](https://hvplot.holoviz.org/) / [HoloViews](https://holoviews.org/) (Bokeh)
# line plots of total porosity, field-capacity porosity, and wilting-point porosity vs. bulk
# density, one line per texture class. Hover for values; use the toolbar to pan/zoom.
# (Hydraulic conductivity and infiltration are in Notebook 3.)
#
# **Extrapolation greyed out.** Each line is solid only over physically plausible bulk densities;
# where Rosetta's θₛ exceeds the BD-implied pore space (the `implausible_bd` flag from §5) the curve
# continues as a faint **grey dashed** tail.

# %%
# line_with_extrapolation (shared with Notebook 3) is imported from _helpers: each texture's
# curve is solid over plausible bulk densities and grey-dashed where implausible_bd is True.
line_with_extrapolation(
    result,
    "total_porosity",
    "Total porosity θₛ (cm³/cm³)",
    "Rosetta total porosity vs. bulk density by USDA texture class",
)

# %%
line_with_extrapolation(
    result,
    "field_capacity_porosity",
    "Field-capacity porosity at 33 kPa (cm³/cm³)",
    "Rosetta field-capacity porosity vs. bulk density by USDA texture class",
)

# %%
line_with_extrapolation(
    result,
    "wilting_point_porosity",
    "Wilting-point porosity at 1500 kPa (cm³/cm³)",
    "Rosetta wilting-point porosity vs. bulk density by USDA texture class",
)

# %% [markdown]
# ::: {.callout-tip appearance="simple"}
# **Takeaway:** As bulk density increases — reflecting more compacted soil — total porosity, field capacity, and available water all decline across every texture class, so compaction directly reduces the water a soil can hold and release to plants and infiltration.
# :::

# %% [markdown]
# ## 7. Soil water partitioning by texture
#
# Horizontal stacked bars (one per texture class, in the order defined above) partitioning the pore space into three volumetric water states, in the style of the [AQP water-retention-curve plots](https://ncss-tech.github.io/AQP/aqp/water-retention-curves.html):
#
# | Band | Range | Color |
# | --- | --- | --- |
# | **Unavailable water** | 0 → wilting point (1500 kPa) | orange |
# | **Available water** | wilting point → field capacity (33 kPa) | green |
# | **Drainable water** | field capacity → total porosity (saturation) | blue |
#
# The three bands sum to total porosity θₛ. Use the **bulk-density slider** at the bottom to vary BD across the 0.8–1.9 g/cm³ range; the selected value is shown in the plot title. A red **⚠** marks bars at bulk densities that are physically implausible for that texture (`implausible_bd`: θₛ > 1 − BD/2.65 — extrapolation).

# %%
import holoviews as hv

# Place HoloViews widgets (the bulk-density slider) along the bottom of the plot.
hv.output(widget_location="bottom")

# Partition pore space into three water states for every bulk density, then plot as
# horizontal stacked bars (AQP water-retention-curve style) with a bulk-density slider.
part = result.copy()
part["unavailable"] = part["wilting_point_porosity"]                          # 0 -> PWP
part["available"] = part["available_water_capacity"]                          # PWP -> FC
part["drainable"] = part["total_porosity"] - part["field_capacity_porosity"]  # FC -> saturation

# y-axis label: texture class with its hydrologic soil group in parentheses, e.g. "sand (A)"
part["texture_label"] = part["texture_class"] + " (" + part["hydrologic_soil_group"] + ")"

water_states = ["unavailable", "available", "drainable"]
state_colors = ["orange", "green", "blue"]
label_order = [f"{cls} ({HYDROLOGIC_SOIL_GROUP[cls]})" for cls in TEXTURE_CLASSES]  # sand ... clay

part_long = part.melt(
    id_vars=["texture_label", "bulk_density_g_cm3"],
    value_vars=water_states,
    var_name="water_state",
    value_name="volumetric_fraction",
)
# reverse texture order so sand plots at the TOP of the horizontal (inverted) axis
part_long["texture_label"] = pd.Categorical(
    part_long["texture_label"], categories=label_order[::-1], ordered=True
)
part_long["water_state"] = pd.Categorical(
    part_long["water_state"], categories=water_states, ordered=True
)
part_long = part_long.sort_values(["bulk_density_g_cm3", "texture_label", "water_state"])

# groupby -> bulk-density slider; dynamic=False embeds every frame so the slider works
# in the saved/exported notebook without a live kernel. {dimensions} interpolates the
# selected slider value into the title, so it updates as the slider moves.
hmap = part_long.hvplot.barh(
    x="texture_label",
    y="volumetric_fraction",
    by="water_state",
    groupby="bulk_density_g_cm3",
    dynamic=False,
    stacked=True,
    color=state_colors,
    ylim=(0, 0.7),  # fixed value axis so bars don't rescale while sliding
    xlabel="texture class (hydrologic soil group)",
    ylabel="volumetric water content (cm³/cm³)",
    title="Soil water partitioning by texture (Rosetta) — {dimensions}",
    width=750,
    height=550,
    legend="top_right",
)
# Per-BD "⚠" markers at the end of bars for textures that are physically implausible at that
# BD (extrapolation, implausible_bd). Built as a HoloMap keyed like the bars so it tracks the slider.
def _implausible_marks(bd):
    sub = result[result["bulk_density_g_cm3"] == bd]
    # barh inverts axes, so in data space the categorical texture is the FIRST coordinate and the
    # numeric bar value is the SECOND; passing them the other way injects the numeric values as
    # spurious extra y-axis categories (the bug this fixes).
    marks = [
        hv.Text(f"{r.texture_class} ({r.hydrologic_soil_group})", r.total_porosity + 0.01, "⚠",
                halign="left").opts(text_color="red", text_font_size="9pt")
        for r in sub.itertuples() if r.implausible_bd
    ]
    # every frame must be the same type (Overlay); use an off-canvas blank when none are flagged
    if not marks:
        first = next(iter(TEXTURE_CLASSES))
        marks = [hv.Text(f"{first} ({HYDROLOGIC_SOIL_GROUP[first]})", -1.0, " ").opts(text_alpha=0)]
    return hv.Overlay(marks)

marks_hmap = hv.HoloMap({bd: _implausible_marks(bd) for bd in bulk_densities}, kdims=["bulk_density_g_cm3"])

# readable slider/title label + consistent one-decimal bulk density
(hmap * marks_hmap).redim(
    bulk_density_g_cm3=hv.Dimension("Bulk density, g/cm³ (higher is more compacted)", default=1.4, value_format=lambda v: f"{v:.1f}")
)

# %% [markdown]
# ::: {.callout-tip appearance="simple"}
# **Takeaway:** Fine-textured soils (clay, silty clay loam) hold more total and unavailable water, while medium-textured soils (loam, silt loam) typically offer the highest plant-available water — the green band — making them most productive for crops and vegetation.
# :::

# %% [markdown]
# ## 8. Soil water vs. texture (transposed line view)
#
# The same three water states as Section 7, transposed and drawn as a line/area diagram in the style of the classic FAO available-water figure: **texture on the x-axis** (coarse → fine), **water volume on the y-axis** expressed as **inches of water per foot of soil depth** (volumetric water content × 12), with the *permanent wilting point*, *field capacity*, and *total porosity* curves bounding the filled bands:
#
# - **Unavailable water** (orange): 0 → permanent wilting point
# - **Available water** (green): wilting point → field capacity
# - **Drainable water** (blue): field capacity → total porosity
#
# Use the **bulk-density slider** at the bottom to vary BD; the selected value appears in the title.
# Texture columns are **greyed out** at bulk densities where the combination is physically
# implausible (`implausible_bd`: Rosetta θₛ > pore space 1 − BD/2.65) — i.e. extrapolation.

# %%
# Transposed line/area view: texture on x (coarse -> fine), water content on y, with a
# bulk-density slider. Built as a HoloMap of overlays (one frame per bulk density) so it
# embeds every frame and works in the saved notebook without a live kernel. Columns flagged
# implausible_bd (θₛ > 1 − BD/2.65) are greyed as an extrapolation warning.
hv.output(widget_location="bottom")

texture_x = {cls: i for i, cls in enumerate(TEXTURE_CLASSES)}  # sand=0 ... clay=11
# x tick labels: texture class with its hydrologic soil group in parentheses, e.g. "sand (A)"
texture_ticks = [(i, f"{cls} ({HYDROLOGIC_SOIL_GROUP[cls]})") for cls, i in texture_x.items()]

# volumetric water content (cm³/cm³) -> inches of water per foot of soil depth (FAO style)
INCHES_PER_FOOT = 12.0


def _water_profile(bd):
    d = (
        result[result["bulk_density_g_cm3"] == bd]
        .set_index("texture_class")
        .reindex(list(TEXTURE_CLASSES))
        .reset_index()
    )
    x = d["texture_class"].map(texture_x).to_numpy()
    pwp = d["wilting_point_porosity"].to_numpy() * INCHES_PER_FOOT
    fc = d["field_capacity_porosity"].to_numpy() * INCHES_PER_FOOT
    por = d["total_porosity"].to_numpy() * INCHES_PER_FOOT
    return soil_water_texture_band_diagram(
        x, pwp, fc, por,
        texture_labels=d["texture_class"].astype(str).tolist(),
        implausible=d["implausible_bd"].to_numpy(),
    )


profiles = hv.HoloMap(
    {bd: _water_profile(bd) for bd in bulk_densities},
    kdims=[hv.Dimension("Bulk density, g/cm³ (higher is more compacted)", default=1.4, value_format=lambda v: f"{v:.1f}")],
)

profiles.opts(
    hv.opts.Overlay(
        width=820,
        height=520,
        legend_position="top_left",
        xlabel="texture class (hydrologic soil group); coarse → fine",
        ylabel="Water Storage Capacity (inches per foot of soil depth)",
        title="Soil water vs. texture (Rosetta) — {dimensions}",
    ),
    hv.opts.Curve(xticks=texture_ticks, xrotation=45, ylim=(0, 0.7 * INCHES_PER_FOOT)),
    hv.opts.Area(xticks=texture_ticks, xrotation=45, ylim=(0, 0.7 * INCHES_PER_FOOT)),
)

# %% [markdown]
# ::: {.callout-tip appearance="simple"}
# **Takeaway:** Loam and silt loam soils store the most plant-available water (roughly 2–2.5 inches per foot of depth), while sandy soils drain quickly and clay soils lock much of their water below the wilting point — a pattern that shifts noticeably as you move the bulk-density slider toward more compacted values.
# :::
