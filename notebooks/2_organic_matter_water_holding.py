# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.2
#   kernelspec:
#     display_name: default
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Organic-matter effects on soil water holding
#
# Companion to **`rosetta_porosity_by_texture.ipynb`**. That notebook estimates the ROSETTA
# texture × bulk-density baseline (total porosity, field capacity, wilting point) for every
# USDA texture class and writes **`rosetta_porosity_by_texture.csv`** — **run it first**.
#
# Here we layer **organic-matter** effects on top, with two approaches:
#
# - **Section 1 — ROSETTA + Minasny & McBratney (2018):** ROSETTA's mineral baseline plus the
#   empirical organic-carbon increments from M&M Table 2 (the OC sensitivity trusted over
#   Saxton–Rawls). **Two sliders: mineral bulk density and organic carbon.**
# - **Section 2 — Saxton & Rawls (2006):** an independent, self-contained PTF taking
#   sand/clay/OM directly, with a validation of its OC sensitivity against M&M (§2.1).
#
# Two outputs get equal billing for our **stormwater** audience: **available water**
# (FC − WP, plant-available) and **drainable water** (SAT − FC), the fast-draining pore
# space relevant to infiltration / detention storage.

# %%
import numpy as np
import pandas as pd
import hvplot.pandas  # noqa: F401  (registers the .hvplot accessor)
import holoviews as hv

pd.set_option("display.float_format", lambda v: f"{v:0.3f}")

# ROSETTA texture × bulk-density baseline produced by rosetta_porosity_by_texture.ipynb
result = pd.read_csv("rosetta_porosity_by_texture.csv")

# Reconstruct the shared constants from the baseline table (canonical sand -> clay order)
TEXTURE_CLASSES = {
    row.texture_class: (row.sand_pct, row.silt_pct, row.clay_pct)
    for row in result[["texture_class", "sand_pct", "silt_pct", "clay_pct"]]
    .drop_duplicates()
    .itertuples(index=False)
}
HYDROLOGIC_SOIL_GROUP = dict(
    result[["texture_class", "hydrologic_soil_group"]].drop_duplicates().itertuples(index=False, name=None)
)
bulk_densities = np.sort(result["bulk_density_g_cm3"].unique())

texture_x = {cls: i for i, cls in enumerate(TEXTURE_CLASSES)}  # sand=0 ... clay=11
texture_ticks = [(i, f"{cls} ({HYDROLOGIC_SOIL_GROUP[cls]})") for cls, i in texture_x.items()]
INCHES_PER_FOOT = 12.0  # volumetric water content (cm3/cm3) -> inches per foot of soil depth

print(f"loaded {len(result)} ROSETTA rows: {len(TEXTURE_CLASSES)} textures x {len(bulk_densities)} bulk densities")
result.head()

# %% [markdown]
# ## 1. ROSETTA + organic-matter modifier (Minasny & McBratney 2018)
#
# Keep **ROSETTA's** texture + bulk-density skill for the *mineral* soil baseline, then add the
# **empirical organic-carbon increments** from **Minasny & McBratney (2018) Table 2** — the OC
# sensitivity that Minnesota soil scientists trust over Saxton–Rawls. Per **+1 % organic carbon**
# (= +10 g C kg⁻¹), by USDA texture group:
#
# | M&M group | ΔWP | ΔAWC | ΔSAT  (mm 100 mm⁻¹ per 1 % OC) |
# | --- | --- | --- | --- |
# | Coarse | 0.86 | 1.94 | 4.59 |
# | Medium | 0.68 | 1.79 | 3.59 |
# | Fine | 0.54 | 1.41 | 3.23 |
#
# **Blend** (volumetric, cm³/cm³); ROSETTA gives the OC = 0 baseline at the **mineral bulk
# density set by the BD slider**:
#
# - WP(OC)  = WP\_ROSETTA + (ΔWP/100)·OC
# - AWC(OC) = AWC\_ROSETTA + (ΔAWC/100)·OC   (AWC\_ROSETTA = FC − WP from ROSETTA)
# - SAT(OC) = θₛ\_ROSETTA + (ΔSAT/100)·OC
# - FC(OC)  = WP + AWC  ;  **drainable = SAT − FC**
#
# We apply M&M's **WP, AWC and SAT** slopes — the three quantities that define the unavailable /
# available / drainable bands — and *derive* FC = WP + AWC. This reproduces M&M's headline AWC
# sensitivity exactly; the drainable response follows from ΔSAT − ΔFC. (Because M&M regressed each
# property independently, ΔAWC ≠ ΔFC − ΔWP; anchoring on ΔFC instead would understate the AWC
# response by ~25–50 %, especially in fine soils.)
#
# **Caveats.** (1) The OC = 0 curve is ROSETTA's prediction at the selected `BD` — a *nominal*
# mineral baseline, since ROSETTA's training data already carry some OC. (2) The BD and OC sliders
# are **independent "what-if" axes**; in reality organic matter *lowers* bulk density (the
# low-BD ↔ high-OC diagonal is the realistic region), and a low-BD + high-OC corner double-counts
# porosity, so don't read the extreme corners as coupled predictions. (3) The modifier is
# **linear**, whereas M&M found diminishing returns (largest gains 0→1 % OC), so it may overstate
# gains at high OC; their data span OC < 10 %. (4) OM ≈ OC / 0.58 (van Bemmelen), so the 0–8 % OC
# axis ≈ 0–14 % OM.

# %%
# ROSETTA mineral baseline (per bulk density) + additive Minasny & McBratney (2018) OC
# increments, computed across the full bulk-density range so BD is an interactive axis.
VB = 0.58  # van Bemmelen factor: OM ≈ OC / VB

# Minasny & McBratney (2018) Table 2 slopes, mm H2O 100 mm-1 per +1% OC (= +10 g C/kg)
MM_SLOPES = {
    "coarse": {"WP": 0.86, "AWC": 1.94, "SAT": 4.59},
    "medium": {"WP": 0.68, "AWC": 1.79, "SAT": 3.59},
    "fine":   {"WP": 0.54, "AWC": 1.41, "SAT": 3.23},
}
MM_GROUP = {
    "sand": "coarse", "loamy sand": "coarse", "sandy loam": "coarse", "sandy clay loam": "coarse",
    "loam": "medium", "silt loam": "medium", "silt": "medium", "clay loam": "medium", "silty clay loam": "medium",
    "sandy clay": "fine", "silty clay": "fine", "clay": "fine",
}

oc_values = np.round(np.arange(0.0, 8.0 + 1e-9, 0.5), 2)  # organic carbon %, 0–8 (≈ 0–14% OM)

base = result.set_index(["texture_class", "bulk_density_g_cm3"])
blend_rows = []
for bd in bulk_densities:
    for oc in oc_values:
        for cls in TEXTURE_CLASSES:
            s = MM_SLOPES[MM_GROUP[cls]]
            sat0 = base.loc[(cls, bd), "total_porosity"]
            fc0 = base.loc[(cls, bd), "field_capacity_porosity"]
            wp0 = base.loc[(cls, bd), "wilting_point_porosity"]
            awc0 = fc0 - wp0
            wp = wp0 + s["WP"] / 100 * oc
            awc = awc0 + s["AWC"] / 100 * oc   # M&M AWC slope applied directly
            sat = sat0 + s["SAT"] / 100 * oc
            fc = wp + awc                       # derive FC so AWC matches M&M exactly
            blend_rows.append(
                {
                    "texture_class": cls,
                    "hydrologic_soil_group": HYDROLOGIC_SOIL_GROUP[cls],
                    "mm_group": MM_GROUP[cls],
                    "bulk_density_g_cm3": bd,
                    "oc_pct": oc,
                    "om_pct_approx": round(oc / VB, 2),
                    "wilting_point_porosity": wp,
                    "field_capacity_porosity": fc,
                    "total_porosity": sat,
                    "available_water_capacity": awc,
                    "drainable_water": sat - fc,
                }
            )
blend_df = pd.DataFrame(blend_rows)
blend_df["texture_class"] = pd.Categorical(blend_df["texture_class"], categories=list(TEXTURE_CLASSES), ordered=True)

print(f"{len(blend_df)} rows  ({len(TEXTURE_CLASSES)} textures x {len(bulk_densities)} BD x {len(oc_values)} OC)")
blend_df[(blend_df["bulk_density_g_cm3"] == 1.5) & (blend_df["oc_pct"].isin([0.0, 2.0, 4.0]))].round(3).head(24)

# %%
# AVAILABLE water capacity vs. organic carbon, one line per texture class, with a BD slider.
hv.output(widget_location="bottom")

om8_oc = 8.0 * VB  # 8% OM ≈ 4.64% OC marks the top of the primary range

awc_lines = blend_df.hvplot.line(
    x="oc_pct",
    y="available_water_capacity",
    by="texture_class",
    groupby="bulk_density_g_cm3",
    dynamic=False,
    xlabel="soil organic carbon (% by weight)   [OM ≈ OC / 0.58]",
    ylabel="available water capacity (cm³/cm³)",
    title="ROSETTA + M&M blend: AVAILABLE water vs. organic carbon — {dimensions}",
    width=820,
    height=500,
    legend="right",
    grid=True,
    ylim=(0, 0.42),
)
(awc_lines * hv.VLine(om8_oc).opts(color="black", line_dash="dotted", line_width=1)).redim(
    bulk_density_g_cm3=hv.Dimension("Bulk density (g/cm³)", value_format=lambda v: f"{v:.1f}")
)

# %%
# DRAINABLE water (saturation − field capacity) vs. organic carbon — the rapidly draining pore
# space that matters for stormwater storage / infiltration. BD slider.
hv.output(widget_location="bottom")

drain_lines = blend_df.hvplot.line(
    x="oc_pct",
    y="drainable_water",
    by="texture_class",
    groupby="bulk_density_g_cm3",
    dynamic=False,
    xlabel="soil organic carbon (% by weight)   [OM ≈ OC / 0.58]",
    ylabel="drainable water  SAT − FC  (cm³/cm³)",
    title="ROSETTA + M&M blend: DRAINABLE water vs. organic carbon — {dimensions}",
    width=820,
    height=500,
    legend="right",
    grid=True,
    ylim=(0, 0.60),
)
(drain_lines * hv.VLine(om8_oc).opts(color="black", line_dash="dotted", line_width=1)).redim(
    bulk_density_g_cm3=hv.Dimension("Bulk density (g/cm³)", value_format=lambda v: f"{v:.1f}")
)

# %%
# FAO-style transposed diagram for the ROSETTA + M&M blend, with TWO sliders: mineral bulk
# density and organic carbon. dynamic=False embeds every (BD, OC) frame so it works in the
# saved notebook without a live kernel. (OC on a coarser 1% grid to keep the frame count sane.)
hv.output(widget_location="bottom")

oc_grid = np.round(np.arange(0.0, 8.0 + 1e-9, 1.0), 1)


def _blend_profile(bd, oc):
    d = (
        blend_df[(blend_df["bulk_density_g_cm3"] == bd) & (blend_df["oc_pct"] == oc)]
        .set_index("texture_class")
        .reindex(list(TEXTURE_CLASSES))
        .reset_index()
    )
    x = d["texture_class"].map(texture_x).to_numpy()
    pwp = d["wilting_point_porosity"].to_numpy() * INCHES_PER_FOOT
    fc = d["field_capacity_porosity"].to_numpy() * INCHES_PER_FOOT
    por = d["total_porosity"].to_numpy() * INCHES_PER_FOOT

    bands = (
        hv.Area((x, pwp, pwp * 0), vdims=["y", "y2"]).opts(color="orange", alpha=0.45, line_alpha=0)
        * hv.Area((x, fc, pwp), vdims=["y", "y2"]).opts(color="green", alpha=0.45, line_alpha=0)
        * hv.Area((x, por, fc), vdims=["y", "y2"]).opts(color="blue", alpha=0.40, line_alpha=0)
    )
    lines = (
        hv.Curve((x, pwp), label="Permanent wilting point").opts(color="black", line_width=2)
        * hv.Curve((x, fc), label="Field capacity").opts(color="black", line_width=2)
        * hv.Curve((x, por), label="Total porosity").opts(color="gray", line_width=1.5, line_dash="dashed")
    )
    labels = (
        hv.Text(8, pwp[8] * 0.5, "Unavailable\nwater").opts(text_color="saddlebrown", text_font_size="9pt")
        * hv.Text(5, (pwp[5] + fc[5]) / 2, "Available water").opts(text_color="darkgreen", text_font_size="10pt")
        * hv.Text(2.2, (fc[2] + por[2]) / 2, "Drainable\nwater").opts(text_color="navy", text_font_size="10pt")
    )
    return bands * lines * labels


blend_profiles = hv.HoloMap(
    {(bd, oc): _blend_profile(bd, oc) for bd in bulk_densities for oc in oc_grid},
    kdims=[
        hv.Dimension("Bulk density (g/cm³)", value_format=lambda v: f"{v:.1f}"),
        hv.Dimension("Organic carbon (%)", value_format=lambda v: f"{v:.1f}"),
    ],
)

blend_profiles.opts(
    hv.opts.Overlay(
        width=820,
        height=520,
        legend_position="top_left",
        xlabel="texture class (hydrologic soil group); coarse → fine",
        ylabel="Water volume (inches per foot of soil depth)",
        title="Soil water vs. texture — ROSETTA + M&M blend — {dimensions}",
    ),
    hv.opts.Curve(xticks=texture_ticks, xrotation=45, ylim=(0, 10)),
    hv.opts.Area(xticks=texture_ticks, xrotation=45, ylim=(0, 10)),
)

# %% [markdown]
# ## 2. Organic-matter sensitivity (Saxton–Rawls 2006)
#
# An independent alternative to Section 1's blend: the **Saxton & Rawls (2006)** pedotransfer
# functions, which take **sand, clay, and organic-matter %** directly and were developed from
# USDA/NRCS data for the continental USA. Self-contained (no ROSETTA baseline) — though (see §2.1)
# it gives a smaller, and for clays negative, OC effect than Minasny & McBratney.
#
# Calibrated for **OM ≤ 8 % by weight** (≈ 4.6 % organic carbon); values above that are shown but
# **flagged as extrapolation**.

# %%
def saxton_rawls(sand_frac, clay_frac, om_pct):
    """Saxton & Rawls (2006), SSSAJ 70:1569-1578 — soil-water characteristics from texture and
    organic matter. Inputs: sand & clay as decimal mass fractions (0-1), organic matter in % by
    weight. Returns volumetric water contents (cm³/cm³): theta_1500 (permanent wilting point),
    theta_33 (field capacity), theta_S (total porosity / saturation). Calibrated for OM up to
    ~8 %; higher OM is extrapolation. Scalars or numpy arrays (broadcast together).
    """
    S, C, OM = sand_frac, clay_frac, om_pct
    t1500t = -0.024 * S + 0.487 * C + 0.006 * OM + 0.005 * (S * OM) - 0.013 * (C * OM) + 0.068 * (S * C) + 0.031
    t1500 = t1500t + (0.14 * t1500t - 0.02)
    t33t = -0.251 * S + 0.195 * C + 0.011 * OM + 0.006 * (S * OM) - 0.027 * (C * OM) + 0.452 * (S * C) + 0.299
    t33 = t33t + (1.283 * t33t**2 - 0.374 * t33t - 0.015)
    tS33t = 0.278 * S + 0.034 * C + 0.022 * OM - 0.018 * (S * OM) - 0.027 * (C * OM) - 0.584 * (S * C) + 0.078
    tS33 = tS33t + (0.636 * tS33t - 0.107)
    tS = t33 + tS33 - 0.097 * S + 0.043
    return t1500, t33, tS


OM_VALID_MAX = 8.0  # Saxton–Rawls organic-matter calibration limit (% by weight)
om_values = np.round(np.arange(0.0, 15.0 + 1e-9, 1.0), 1)  # 0–15 %, 1 % steps

sr_rows = []
for om in om_values:
    for cls, (sand, silt, clay) in TEXTURE_CLASSES.items():
        pwp, fc, por = saxton_rawls(sand / 100, clay / 100, om)
        sr_rows.append(
            {
                "texture_class": cls,
                "hydrologic_soil_group": HYDROLOGIC_SOIL_GROUP[cls],
                "om_pct": om,
                "wilting_point_porosity": pwp,
                "field_capacity_porosity": fc,
                "total_porosity": por,
                "available_water_capacity": fc - pwp,
                "drainable_water": por - fc,
                "within_validity": om <= OM_VALID_MAX,
            }
        )
sr_df = pd.DataFrame(sr_rows)
sr_df["texture_class"] = pd.Categorical(sr_df["texture_class"], categories=list(TEXTURE_CLASSES), ordered=True)

print(f"{len(sr_df)} rows  ({len(TEXTURE_CLASSES)} textures x {len(om_values)} OM levels)")
sr_df[sr_df["om_pct"].isin([0.0, 2.0, 4.0, 8.0])].head(20)

# %%
# Saxton–Rawls AVAILABLE and DRAINABLE water vs. organic matter, one line per texture class.
# Shaded band / dashed line mark the 8% OM calibration limit (extrapolation beyond).
hv.output(widget_location="right")

sr_flag = (
    hv.VSpan(OM_VALID_MAX, float(sr_df["om_pct"].max())).opts(color="gray", alpha=0.12)
    * hv.VLine(OM_VALID_MAX).opts(color="black", line_dash="dashed", line_width=1)
)

sr_awc = sr_df.hvplot.line(
    x="om_pct", y="available_water_capacity", by="texture_class",
    xlabel="soil organic matter (% by weight)",
    ylabel="available water capacity (cm³/cm³)",
    title="Saxton–Rawls: AVAILABLE water vs. organic matter",
    width=820, height=460, legend="right", grid=True, ylim=(0, 0.22),
)
sr_drain = sr_df.hvplot.line(
    x="om_pct", y="drainable_water", by="texture_class",
    xlabel="soil organic matter (% by weight)",
    ylabel="drainable water  SAT − FC  (cm³/cm³)",
    title="Saxton–Rawls: DRAINABLE water vs. organic matter",
    width=820, height=460, legend="right", grid=True, ylim=(0, 0.45),
)
((sr_awc * sr_flag) + (sr_drain * sr_flag)).cols(1)

# %%
# FAO-style transposed diagram, Saxton–Rawls, with an organic-matter slider.
# Frames above 8% OM are annotated as extrapolation.
hv.output(widget_location="bottom")


def _sr_profile(om):
    d = (
        sr_df[sr_df["om_pct"] == om]
        .set_index("texture_class")
        .reindex(list(TEXTURE_CLASSES))
        .reset_index()
    )
    x = d["texture_class"].map(texture_x).to_numpy()
    pwp = d["wilting_point_porosity"].to_numpy() * INCHES_PER_FOOT
    fc = d["field_capacity_porosity"].to_numpy() * INCHES_PER_FOOT
    por = d["total_porosity"].to_numpy() * INCHES_PER_FOOT

    bands = (
        hv.Area((x, pwp, pwp * 0), vdims=["y", "y2"]).opts(color="orange", alpha=0.45, line_alpha=0)
        * hv.Area((x, fc, pwp), vdims=["y", "y2"]).opts(color="green", alpha=0.45, line_alpha=0)
        * hv.Area((x, por, fc), vdims=["y", "y2"]).opts(color="blue", alpha=0.40, line_alpha=0)
    )
    lines = (
        hv.Curve((x, pwp), label="Permanent wilting point").opts(color="black", line_width=2)
        * hv.Curve((x, fc), label="Field capacity").opts(color="black", line_width=2)
        * hv.Curve((x, por), label="Total porosity").opts(color="gray", line_width=1.5, line_dash="dashed")
    )
    labels = (
        hv.Text(8, pwp[8] * 0.5, "Unavailable\nwater").opts(text_color="saddlebrown", text_font_size="9pt")
        * hv.Text(5, (pwp[5] + fc[5]) / 2, "Available water").opts(text_color="darkgreen", text_font_size="10pt")
        * hv.Text(2.2, (fc[2] + por[2]) / 2, "Drainable\nwater").opts(text_color="navy", text_font_size="10pt")
    )
    overlay = bands * lines * labels
    if om > OM_VALID_MAX:
        overlay = overlay * hv.Text(6, 9.5, "⚠ extrapolated beyond 8% OM calibration").opts(
            text_color="red", text_font_size="9pt"
        )
    return overlay


sr_profiles = hv.HoloMap(
    {om: _sr_profile(om) for om in om_values},
    kdims=[hv.Dimension("Soil organic matter (% by weight)", value_format=lambda v: f"{v:.1f}")],
)

sr_profiles.opts(
    hv.opts.Overlay(
        width=820,
        height=520,
        legend_position="top_left",
        xlabel="texture class (hydrologic soil group); coarse → fine",
        ylabel="Water volume (inches per foot of soil depth)",
        title="Soil water vs. texture (Saxton–Rawls) — {dimensions}",
    ),
    hv.opts.Curve(xticks=texture_ticks, xrotation=45, ylim=(0, 10)),
    hv.opts.Area(xticks=texture_ticks, xrotation=45, ylim=(0, 10)),
)

# %% [markdown]
# ### 2.1 Validation: ΔAWC/ΔOC vs. Minasny & McBratney (2018)
#
# Do the two PTF families agree on how much available water organic carbon adds? We compute the
# Saxton–Rawls **ΔAWC per +1 % organic carbon** for each texture class — over the same
# **OC 0.5 % → 1.5 %** interval M&M used for PTF-derived slopes (OC = 0.58·OM) — then average by
# their coarse/medium/fine groups and compare against Table 2.
#
# Expect only *order-of-magnitude* agreement: both say the effect is small and decreases from
# coarse to fine textures, but **Saxton–Rawls is systematically lower** and turns **negative for
# clays** — a known feature of the Rawls/Saxton–Rawls lineage (M&M note their neural net "did not
# show a negative effect with an increase in OC for clay content larger than 60 %"). This is why
# Section 1 builds the blend on the M&M increments rather than Saxton–Rawls.

# %%
# Saxton–Rawls ΔAWC/ΔOC vs. Minasny & McBratney (2018) Table 2.
MM_AWC_SLOPE = {"general": 1.16, "coarse": 1.94, "medium": 1.79, "fine": 1.41}


def sr_awc_slope_per_pct_oc(sand_frac, clay_frac):
    """Saxton–Rawls ΔAWC over OC 0.5%->1.5% (M&M's interval), in mm H2O 100 mm-1 (= vol%)."""
    om_lo, om_hi = 0.5 / VB, 1.5 / VB  # OC% -> OM%
    p_lo, f_lo, _ = saxton_rawls(sand_frac, clay_frac, om_lo)
    p_hi, f_hi, _ = saxton_rawls(sand_frac, clay_frac, om_hi)
    return ((f_hi - p_hi) - (f_lo - p_lo)) * 100.0  # cm3/cm3 over 1% OC -> mm/100mm


val_df = pd.DataFrame(
    [
        {
            "texture_class": cls,
            "mm_group": MM_GROUP[cls],
            "saxton_rawls_dAWC": sr_awc_slope_per_pct_oc(sand / 100, clay / 100),
        }
        for cls, (sand, silt, clay) in TEXTURE_CLASSES.items()
    ]
)

grp_cmp = val_df.groupby("mm_group", sort=False)["saxton_rawls_dAWC"].mean().reset_index()
grp_cmp["mm_group"] = pd.Categorical(grp_cmp["mm_group"], categories=["coarse", "medium", "fine"], ordered=True)
grp_cmp = grp_cmp.sort_values("mm_group")
grp_cmp["minasny_mcbratney"] = grp_cmp["mm_group"].map(MM_AWC_SLOPE)

print("Mean ΔAWC/ΔOC by texture group (mm H₂O 100 mm⁻¹ per +1% OC):")
print(grp_cmp.round(2).to_string(index=False))
val_df.round(2)

# %%
# Grouped bars: Saxton–Rawls vs Minasny & McBratney mean ΔAWC/ΔOC, by texture group.
hv.output(widget_location="right")

grp_long = grp_cmp.melt(
    id_vars="mm_group",
    value_vars=["saxton_rawls_dAWC", "minasny_mcbratney"],
    var_name="method",
    value_name="dAWC_mm",
)
grp_long["method"] = grp_long["method"].map(
    {"saxton_rawls_dAWC": "Saxton–Rawls (2006)", "minasny_mcbratney": "Minasny & McBratney (2018)"}
).astype(str)

bars = grp_long.hvplot.bar(
    x="mm_group",
    y="dAWC_mm",
    by="method",
    color=["#4c78a8", "#f58518"],
    xlabel="texture group (Minasny & McBratney classes)",
    ylabel="ΔAWC / ΔOC  (mm H₂O 100 mm⁻¹ per +1% OC)",
    title="AWC sensitivity to organic carbon: Saxton–Rawls vs. Minasny & McBratney (2018)",
    width=760,
    height=470,
    legend="top_right",
    ylim=(-1.3, 2.3),
)
refs = (
    hv.HLine(0).opts(color="black", line_width=1)
    * hv.HLine(MM_AWC_SLOPE["general"]).opts(color="gray", line_dash="dotted", line_width=1)
)
bars * refs
