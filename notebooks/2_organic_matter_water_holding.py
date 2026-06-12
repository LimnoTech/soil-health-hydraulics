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
# # Organic-matter effects on soil water holding
#
# Companion to **`rosetta_porosity_by_texture.ipynb`**. That notebook estimates the ROSETTA
# texture × bulk-density baseline (total porosity, field capacity, wilting point) for every
# USDA texture class and writes **`rosetta_porosity_by_texture.csv`** — **run it first**.
#
# We first establish ROSETTA's **mineral-baseline organic carbon** (Section 1), then layer
# **organic-matter** effects on top, two ways:
#
# - **Section 1 — Mineral-baseline OC (UNSODA 2.0):** estimate the organic carbon implicit in
#   ROSETTA's mineral baseline, used as the anchor for the blend.
# - **Section 2 — ROSETTA + Minasny & McBratney (2018):** ROSETTA's mineral baseline plus M&M's
#   empirical organic-carbon increments (the OC sensitivity preferred over Saxton–Rawls).
#   **Two sliders: mineral bulk density and organic matter (0–8 %).**
# - **Section 3 — Saxton & Rawls (2006):** an independent, self-contained PTF taking
#   sand/clay/OM directly, with a validation of its OC sensitivity against M&M (§3.1).
#
# Two outputs get equal billing for our **stormwater** audience: **available water**
# (FC − WP, plant-available) and **drainable water** (SAT − FC), the fast-draining pore
# space relevant to infiltration / detention storage.

# %%
import numpy as np
import pandas as pd
import hvplot.pandas  # noqa: F401  (registers the .hvplot accessor)
import holoviews as hv
from scipy.stats import linregress

# Shared display helper (and the shared pandas float_format, set on import). See notebooks/_helpers.py.
from _helpers import show, soil_water_texture_band_diagram, soil_water_bd_om_blend_table, VB, OC_BASELINE_PCT, MM_SLOPES, MM_GROUP

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
show(result)

# %% [markdown]
# ## 1. Mineral-baseline organic carbon (from UNSODA 2.0)
#
# ROSETTA has no OC input, but its training samples (UNSODA + others) are mineral-dominated
# soils that *do* carry organic carbon — so the ROSETTA prediction is a **nominal baseline at
# OC > 0**, not a true organic-free (OC = 0) soil. Before blending in any organic-matter effect
# (Section 2), we estimate that baseline OC from the 367 UNSODA 2.0 samples that report both bulk
# density and organic-matter content (OC = 0.58·OM, van Bemmelen). OC falls clearly with bulk
# density (Pearson r ≈ −0.6); the scatter below overlays **OLS regressions of OC on bulk density**
# fit separately for topsoil, subsoil, and all mineral data, with the slopes/intercepts/R²/p
# reported in the accompanying table. For the mineral subset (OM ≤ 20 %): all-horizon mean ≈ 0.9 %
# OC; **mineral topsoil (≤15 cm) median ≈ 1.1 % OC**. We anchor ROSETTA at **`OC_BASELINE_PCT`
# (default 1.0 % OC)**; Section 2 then applies the M&M increments *relative to that* baseline.
#
# `data_temp/` is git-ignored; run `pixi run python notebooks/fetch_unsoda.py` to (re)create
# the UNSODA extract read below.

# %%
import os

# Organic-carbon anchor for the ROSETTA mineral baseline (estimated below from UNSODA).
# UNSODA mineral subset: all-horizon mean ~0.9 %, topsoil median ~1.1 % — round 1.0 % used here.
OC_BASELINE_PCT = 1.0  # % organic carbon; adjust to taste

_unsoda_path = "data_temp/unsoda_bd_om.csv"
if os.path.exists(_unsoda_path):
    unsoda = pd.read_csv(_unsoda_path)
    mineral = unsoda[unsoda["is_mineral"]].copy()
    mineral["horizon_group"] = np.where(
        mineral["depth_upper"] <= 15, "topsoil (≤15 cm)", "subsoil (>15 cm)"
    )
    print("UNSODA mineral subset (OM ≤ 20%) — organic carbon %, OC = 0.58·OM:")
    print(mineral.groupby("horizon_group")["OC_pct"].agg(["count", "mean", "median"]).round(2))
    print(f"OC~BD Pearson r = {mineral['OC_pct'].corr(mineral['bulk_density_g_cm3']):.2f}")

    # OLS regressions of organic carbon on bulk density: topsoil, subsoil, and all mineral data.
    _reg_groups = {
        "topsoil (≤15 cm)": mineral[mineral["horizon_group"] == "topsoil (≤15 cm)"],
        "subsoil (>15 cm)": mineral[mineral["horizon_group"] == "subsoil (>15 cm)"],
        "all mineral": mineral,
    }
    _reg_colors = {"topsoil (≤15 cm)": "#1f77b4", "subsoil (>15 cm)": "#ff7f0e", "all mineral": "black"}
    _bd_line = np.array([mineral["bulk_density_g_cm3"].min(), mineral["bulk_density_g_cm3"].max()])
    _reg_rows, _reg_lines = [], []
    for _name, _g in _reg_groups.items():
        _lr = linregress(_g["bulk_density_g_cm3"], _g["OC_pct"])
        _reg_rows.append({
            "regression": _name,
            "n": len(_g),
            "slope (%OC per g/cm³)": _lr.slope,
            "intercept (%OC)": _lr.intercept,
            "r": _lr.rvalue,
            "R²": _lr.rvalue ** 2,
            "p_value": _lr.pvalue,
        })
        _reg_lines.append(
            hv.Curve((_bd_line, _lr.intercept + _lr.slope * _bd_line), label=f"{_name} fit").opts(
                color=_reg_colors[_name], line_width=2,
                line_dash="solid" if _name == "all mineral" else "dashed",
            )
        )
    reg_df = pd.DataFrame(_reg_rows)

    _scatter = mineral.hvplot.scatter(
        x="bulk_density_g_cm3", y="OC_pct", by="horizon_group",
        xlabel="bulk density (g/cm³)", ylabel="organic carbon  OC = 0.58·OM  (%)",
        title="UNSODA 2.0: organic carbon vs. bulk density (mineral soils)",
        width=760, height=460, legend="top_right", alpha=0.6, size=25, ylim=(0, 10), grid=True,
    )
    _baseline = hv.HLine(OC_BASELINE_PCT).opts(color="firebrick", line_dash="dashed", line_width=2)
    display(_scatter * _baseline * hv.Overlay(_reg_lines))

    # Regression parameters table (p-values formatted in scientific notation to keep precision).
    _reg_show = reg_df.copy()
    _reg_show["p_value"] = _reg_show["p_value"].map(lambda p: f"{p:.2e}")
    print("OC = slope · bulk_density + intercept  (OLS):")
    display(show(_reg_show, height=160))
else:
    print(f"{_unsoda_path} not found — run `pixi run python notebooks/fetch_unsoda.py` to regenerate it.")
    print(f"Proceeding with the documented default OC_BASELINE_PCT = {OC_BASELINE_PCT} % OC.")

# %% [markdown]
# ## 2. ROSETTA + organic-matter modifier (Minasny & McBratney 2018)
#
# Keep **ROSETTA's** texture + bulk-density skill for the *mineral* soil baseline (Section 1),
# then add the **empirical organic-carbon increments** from **Minasny & McBratney (2018) Table 2**
# — an OC sensitivity derived from >50,000 measurements and preferred here over Saxton–Rawls
# (see §3.1). Per **+1 % organic carbon** (= +10 g C kg⁻¹), by USDA texture group:
#
# | M&M group | ΔWP | ΔAWC | ΔSAT  (mm 100 mm⁻¹ per 1 % OC) |
# | --- | --- | --- | --- |
# | Coarse | 0.86 | 1.94 | 4.59 |
# | Medium | 0.68 | 1.79 | 3.59 |
# | Fine | 0.54 | 1.41 | 3.23 |
#
# **Blend** (volumetric, cm³/cm³); ROSETTA gives the baseline at the **mineral bulk density set
# by the BD slider**, anchored at **OC = OC_base** (≈ 1 %, from UNSODA — see Section 1):
#
# - WP(OC)  = WP\_ROSETTA + (ΔWP/100)·(OC − OC_base)
# - AWC(OC) = AWC\_ROSETTA + (ΔAWC/100)·(OC − OC_base)   (AWC\_ROSETTA = FC − WP from ROSETTA)
# - SAT(OC) = θₛ\_ROSETTA + (ΔSAT/100)·(OC − OC_base)
# - FC(OC)  = WP + AWC  ;  **drainable = SAT − FC**
#
# We apply M&M's **WP, AWC and SAT** slopes — the three quantities that define the unavailable /
# available / drainable bands — and *derive* FC = WP + AWC. This reproduces M&M's headline AWC
# sensitivity exactly; the drainable response follows from ΔSAT − ΔFC. (Because M&M regressed each
# property independently, ΔAWC ≠ ΔFC − ΔWP; anchoring on ΔFC instead would understate the AWC
# response by ~25–50 %, especially in fine soils.)
#
# **Caveats.** (1) ROSETTA's prediction is a *nominal* baseline at **OC ≈ `OC_BASELINE_PCT`**
# (≈ 1 %, from UNSODA — Section 1), not OC = 0; the M&M increments are applied relative to it, so
# the slider's 0 % end is a truly organic-free mineral soil (drier than ROSETTA), clamped at ≥ 0.
# (2) The BD and OC sliders
# are **independent "what-if" axes**; in reality organic matter *lowers* bulk density (the
# low-BD ↔ high-OC diagonal is the realistic region), and a low-BD + high-OC corner double-counts
# porosity, so don't read the extreme corners as coupled predictions. (3) The modifier is
# **linear**, whereas M&M found diminishing returns (largest gains 0→1 % OC), so it may overstate
# gains at high OC; their data span OC < 10 %. (4) OM ≈ OC / 0.58 (van Bemmelen); the line-plot OC
# axis is capped at 5 % (≈ 8.6 % OM) and the diagram's OM slider spans 0–8 %. (5) In the diagram,
# texture columns are **greyed** where the blended
# saturation exceeds the BD-implied pore space (1 − BD/2.65) — physically impossible, i.e.
# extrapolation at that BD × OM (mirrors Notebook 1's `implausible_bd` flag).

# %%
# ROSETTA mineral baseline (per bulk density) + additive Minasny & McBratney (2018) OC
# increments, applied RELATIVE to OC_BASELINE_PCT, across the full bulk-density range.
# VB, MM_SLOPES, MM_GROUP are imported from _helpers; OC_BASELINE_PCT set in §1 above.

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
            d_oc = oc - OC_BASELINE_PCT          # increments relative to the mineral-baseline OC
            wp = max(wp0 + s["WP"] / 100 * d_oc, 0.0)
            awc = max(awc0 + s["AWC"] / 100 * d_oc, 0.0)   # M&M AWC slope; floor at 0
            sat = max(sat0 + s["SAT"] / 100 * d_oc, wp + awc)  # keep SAT >= FC
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
# order each (BD, texture) series ascending in OC, then flag extrapolation: a blended saturation
# exceeding the BD-implied pore space (1 - BD/2.65) is physically impossible (greyed in the plots).
blend_df = blend_df.sort_values(["bulk_density_g_cm3", "texture_class", "oc_pct"]).reset_index(drop=True)
blend_df["implausible"] = blend_df["total_porosity"] > (1.0 - blend_df["bulk_density_g_cm3"] / 2.65)
_blend_grp = blend_df.groupby(["bulk_density_g_cm3", "texture_class"], observed=True, sort=False)["implausible"]
_blend_extrap_mask = blend_df["implausible"] | _blend_grp.shift(-1).fillna(False)  # keep boundary so segments join

print(f"{len(blend_df)} rows  ({len(TEXTURE_CLASSES)} textures x {len(bulk_densities)} BD x {len(oc_values)} OC)")
show(blend_df[(blend_df["bulk_density_g_cm3"] == 1.5) & (blend_df["oc_pct"].isin([0.0, 2.0, 4.0]))].round(3))

# %%
# AVAILABLE water capacity vs. organic carbon, one line per texture class, with a BD slider.
# Each line is solid where plausible and grey-dashed where the (BD, OC) state is an extrapolation
# (blended saturation > 1 - BD/2.65); the greyed region grows as the BD slider increases.
hv.output(widget_location="bottom")

om8_oc = 8.0 * VB  # 8% OM ≈ 4.64% OC marks the top of the primary range


def blend_line(ycol, ylabel, title, ylim):
    common = dict(
        x="oc_pct", y=ycol, by="texture_class", groupby="bulk_density_g_cm3",
        dynamic=False, width=820, height=500, xlim=(0, 5), ylim=ylim,
    )
    solid = blend_df.assign(**{ycol: blend_df[ycol].where(~blend_df["implausible"])}).hvplot.line(
        xlabel="soil organic carbon (% by weight)   [OM ≈ OC / 0.58]",
        ylabel=ylabel, title=title, legend="right", grid=True, **common,
    )
    dashed = (
        blend_df.assign(**{ycol: blend_df[ycol].where(_blend_extrap_mask)})
        .hvplot.line(**common)
        .opts(hv.opts.Curve(color="lightgray", line_dash="dashed", alpha=0.9))
        .opts(show_legend=False)
    )
    y_lab = ylim[1]  # labels just below the top of the plot
    refs = (
        hv.VLine(OC_BASELINE_PCT).opts(color="firebrick", line_dash="dashed", line_width=1)  # baseline OC
        * hv.VLine(om8_oc).opts(color="black", line_dash="dotted", line_width=1)              # 8% OM
        * hv.Text(OC_BASELINE_PCT, y_lab, " ROSETTA baseline OC ≈ 1%", halign="left", valign="top").opts(
            text_color="firebrick", text_font_size="8pt")
        * hv.Text(om8_oc, y_lab, "8% OM ", halign="right", valign="top").opts(
            text_color="black", text_font_size="8pt")
    )
    return (solid * dashed * refs).redim(
        bulk_density_g_cm3=hv.Dimension("Bulk density, g/cm³ (higher is more compacted)", default=1.4, value_format=lambda v: f"{v:.1f}")
    )


blend_line(
    "available_water_capacity",
    "available water capacity (cm³/cm³)",
    "ROSETTA + M&M blend: AVAILABLE water vs. organic carbon — {dimensions}",
    (0, 0.42),
)

# %%
# DRAINABLE water (saturation − field capacity) vs. organic carbon — the rapidly draining pore
# space that matters for stormwater storage / infiltration. BD slider; extrapolated (BD, OC)
# states are grey-dashed (same flag as the AVAILABLE-water plot above).
hv.output(widget_location="bottom")

blend_line(
    "drainable_water",
    "drainable water  SAT − FC  (cm³/cm³)",
    "ROSETTA + M&M blend: DRAINABLE water vs. organic carbon — {dimensions}",
    (0, 0.60),
)

# %%
# FAO-style transposed diagram for the ROSETTA + M&M blend, with TWO sliders: mineral bulk
# density and organic MATTER (the way the audience thinks about it; OM ≈ OC / 0.58). Computed
# directly from the ROSETTA baseline + M&M increments so the slider can carry round OM values.
# dynamic=False embeds every (BD, OM) frame so it works without a live kernel. (OM capped at
# 8% on a 1% grid to bound the frame count / latency.)
hv.output(widget_location="bottom")

om_grid = np.round(np.arange(0.0, 8.0 + 1e-9, 1.0), 1)  # organic matter %, 0–8 (≈ 0–4.6% OC)
x_pos = np.array([texture_x[cls] for cls in TEXTURE_CLASSES])


def _blend_profile(bd, om):
    t = soil_water_bd_om_blend_table(result, bd, om)
    x = np.array([texture_x[cls] for cls in TEXTURE_CLASSES])
    pwp = t["wilting_point_porosity"].to_numpy() * INCHES_PER_FOOT
    fc = t["field_capacity_porosity"].to_numpy() * INCHES_PER_FOOT
    por = t["total_porosity"].to_numpy() * INCHES_PER_FOOT
    return soil_water_texture_band_diagram(
        x, pwp, fc, por,
        texture_labels=t["texture_class"].astype(str).tolist(),
        implausible=t["implausible"].to_numpy(),
    )


blend_profiles = hv.HoloMap(
    {(bd, om): _blend_profile(bd, om) for bd in bulk_densities for om in om_grid},
    kdims=[
        hv.Dimension("Bulk density, g/cm³ (higher is more compacted)", default=1.4, value_format=lambda v: f"{v:.1f}"),
        hv.Dimension("Organic matter (% by weight)", default=2.0, value_format=lambda v: f"{v:.1f}"),
    ],
)

blend_profiles.opts(
    hv.opts.Overlay(
        width=820,
        height=520,
        legend_position="top_left",
        xlabel="texture class (hydrologic soil group); coarse → fine",
        ylabel="Water Storage Capacity (inches per foot of soil depth)",
        title="Soil water vs. texture — ROSETTA + Minasny & McBratney blend\n{dimensions}",
    ),
    hv.opts.Curve(xticks=texture_ticks, xrotation=45, ylim=(0, 10)),
    hv.opts.Area(xticks=texture_ticks, xrotation=45, ylim=(0, 10)),
)

# %% [markdown]
# ## 3. Organic-matter sensitivity (Saxton–Rawls 2006)
#
# An independent alternative to Section 2's blend: the **Saxton & Rawls (2006)** pedotransfer
# functions, which take **sand, clay, and organic-matter %** directly and were developed from
# USDA/NRCS data for the continental USA. Self-contained (no ROSETTA baseline) — though (see §3.1)
# it gives a smaller, and for clays negative, OC effect than Minasny & McBratney.
#
# Restricted to its calibrated range, **OM ≤ 8 % by weight** (≈ 4.6 % organic carbon).

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
om_values = np.round(np.arange(0.0, OM_VALID_MAX + 1e-9, 1.0), 1)  # 0–8 %, 1 % steps (calibrated range)

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
            }
        )
sr_df = pd.DataFrame(sr_rows)
sr_df["texture_class"] = pd.Categorical(sr_df["texture_class"], categories=list(TEXTURE_CLASSES), ordered=True)

print(f"{len(sr_df)} rows  ({len(TEXTURE_CLASSES)} textures x {len(om_values)} OM levels)")
show(sr_df[sr_df["om_pct"].isin([0.0, 2.0, 4.0, 8.0])])

# %%
# Saxton–Rawls AVAILABLE and DRAINABLE water vs. organic matter, one line per texture class,
# over the calibrated 0–8 % OM range.
hv.output(widget_location="right")

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
(sr_awc + sr_drain).cols(1)

# %%
# FAO-style transposed diagram, Saxton–Rawls, with an organic-matter slider (calibrated 0–8 % OM).
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
    return soil_water_texture_band_diagram(
        x, pwp, fc, por, texture_labels=d["texture_class"].astype(str).tolist()
    )


sr_profiles = hv.HoloMap(
    {om: _sr_profile(om) for om in om_values},
    kdims=[hv.Dimension("Soil organic matter (% by weight)", default=2.0, value_format=lambda v: f"{v:.1f}")],
)

sr_profiles.opts(
    hv.opts.Overlay(
        width=820,
        height=520,
        legend_position="top_left",
        xlabel="texture class (hydrologic soil group); coarse → fine",
        ylabel="Water Storage Capacity (inches per foot of soil depth)",
        title="Soil water vs. texture (Saxton–Rawls) — {dimensions}",
    ),
    hv.opts.Curve(xticks=texture_ticks, xrotation=45, ylim=(0, 10)),
    hv.opts.Area(xticks=texture_ticks, xrotation=45, ylim=(0, 10)),
)

# %% [markdown]
# ### 3.1 Validation: ΔAWC/ΔOC vs. Minasny & McBratney (2018)
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
# Section 2 builds the blend on the M&M increments rather than Saxton–Rawls.

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
