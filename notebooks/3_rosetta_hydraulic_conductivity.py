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
# # Rosetta hydraulic conductivity & infiltration
#
# This notebook shows **how fast water moves through soil and how compaction slows it** — two
# questions central to stormwater management and soil health. Three interactive charts let you
# explore by texture class and bulk density:
#
# - **§1 Saturated hydraulic conductivity (Ksat)** — how fast a wet soil transmits water.
# - **§2 Unsaturated hydraulic conductivity K(h)** — how conductivity collapses as a soil dries.
# - **§3 Green–Ampt infiltration** f(t) / F(t) — how quickly a surface absorbs ponded water over time.
#
# All charts use stormwater units (in/hr) with bulk-density sliders. Physically implausible
# BD × texture combinations (`implausible_bd`, θₛ > 1 − BD/2.65) are greyed as extrapolation,
# as in Notebook 1.
#
# ::: {.callout-note collapse="true"}
# ## For researchers: data source
# Charts read `rosetta_porosity_by_texture.csv` (run Notebook 1 first), which carries Rosetta's
# saturated **Ksat** and the **Mualem–van Genuchten** unsaturated-conductivity parameters
# (`k0_cm_day`, `mualem_L`, `vg_alpha_1cm`, `vg_n`).
# :::

# %%
import numpy as np
import pandas as pd
import hvplot.pandas  # noqa: F401  (registers the .hvplot accessor)
import holoviews as hv

# Shared helpers (display table, Mualem K(h), extrapolation-aware line plot); importing _helpers
# also sets the shared pandas float_format. See notebooks/_helpers.py.
from _helpers import show, mualem_k, line_with_extrapolation

# Rosetta baseline + Mualem–van Genuchten parameters from Notebook 1
result = pd.read_csv("rosetta_porosity_by_texture.csv")

# Reconstruct shared constants from the table (canonical sand -> clay order)
TEXTURE_CLASSES = list(result["texture_class"].drop_duplicates())
HYDROLOGIC_SOIL_GROUP = dict(
    result[["texture_class", "hydrologic_soil_group"]].drop_duplicates().itertuples(index=False, name=None)
)
bulk_densities = np.sort(result["bulk_density_g_cm3"].unique())
texture_x = {cls: i for i, cls in enumerate(TEXTURE_CLASSES)}
texture_ticks = [(i, f"{cls} ({HYDROLOGIC_SOIL_GROUP[cls]})") for cls, i in texture_x.items()]


# mualem_k and line_with_extrapolation are shared with Notebook 1 and imported from _helpers.
# bd_slider_line (the slider-driven K(h) / Green–Ampt line plots) is specific to this notebook.
def bd_slider_line(dfL, x, y, xlabel, ylabel, title, **kw):
    """Line plot (one per texture) with a bulk-density slider; rows flagged implausible_bd are
    drawn grey-dashed. Returns a HoloMap keyed by bulk_density_g_cm3 (caller overlays/redims)."""
    common = dict(
        x=x, y=y, by="texture_class", groupby="bulk_density_g_cm3", dynamic=False,
        width=820, height=520, **kw,
    )
    solid = dfL.assign(**{y: dfL[y].where(~dfL["implausible_bd"])}).hvplot.line(
        xlabel=xlabel, ylabel=ylabel, title=title, legend="right", grid=True, **common
    )
    dashed = (
        dfL.assign(**{y: dfL[y].where(dfL["implausible_bd"])}).hvplot.line(**common)
        .opts(hv.opts.Curve(color="lightgray", line_dash="dashed", alpha=0.9))
        .opts(show_legend=False)
    )
    return solid * dashed


print(f"loaded {len(result)} rows: {len(TEXTURE_CLASSES)} textures x {len(bulk_densities)} bulk densities")
show(result)

# %% [markdown]
# ## 1. Saturated hydraulic conductivity (Ksat)
#
# Rosetta's saturated hydraulic conductivity vs. bulk density, one line per texture class, on a
# **log axis** in stormwater units (in/hr); grey-dashed where extrapolated (`implausible_bd`).
#
# Note the spurious **upturn at high BD for silt and other fine textures** — a neural-network
# *extrapolation artifact* (those dense fine-soil states are absent from Rosetta's training data),
# not a real rise in conductivity; it falls entirely inside the greyed region.

# %%
line_with_extrapolation(
    result,
    "ksat_in_hr",
    "Saturated hydraulic conductivity Ksat (in/hr, log scale)",
    "Rosetta saturated hydraulic conductivity vs. bulk density by USDA texture class",
    logy=True,
)

# %% [markdown]
# ::: {.callout-tip appearance="simple"}
# **Takeaway:** Sandy soils transmit water orders of magnitude faster than clays, and increasing bulk density (compaction) reduces Ksat sharply across all texture classes — a healthy, loose soil infiltrates far more water than a compacted one of the same texture.
# :::

# %% [markdown]
# ## 2. Unsaturated hydraulic conductivity K(h)
#
# As soil dries out, its ability to conduct water drops by many orders of magnitude. The
# log–log chart below shows how K falls with increasing suction for each texture class — use the
# bulk-density slider to see how compaction shifts every curve downward. Dotted verticals mark
# field capacity (FC, 330 cm suction) and wilting point (WP, 15 000 cm suction).
#
# ::: {.callout-note collapse="true"}
# ## For researchers: Mualem–van Genuchten K(h)
# Rosetta gives the full **Mualem–van Genuchten** parameter set, so conductivity is defined not
# just at saturation but at every suction: K(h) = K0·Se(h)^L·[1 − (1 − Se^(1/m))^m]² with
# Se(h) = [1 + (αh)^n]^(−m) (see `mualem_k`). Curves use Rosetta's **K0** (matching point) and **L**
# (columns 5–6), *not* Ksat. `k_fc_cm_day` / `k_wp_cm_day` in the table are K at those tensions.
# Log–log axes; bulk-density slider; implausible (`implausible_bd`) BD × texture combinations
# are grey-dashed.
# :::

# %%
hv.output(widget_location="bottom")

# K(h) over a log-spaced suction grid for every texture × bulk density
_suction = np.logspace(0.0, np.log10(15000.0), 36)  # 1 .. 15000 cm (log-spaced; 36 pts render identically on the log-log axis)
_kh_rows = []
for r in result.itertuples():
    K = mualem_k(_suction, r.vg_alpha_1cm, r.vg_n, r.k0_cm_day, r.mualem_L)  # cm/day
    for h, k in zip(_suction, K):
        _kh_rows.append((r.texture_class, r.bulk_density_g_cm3, h, k / (24.0 * 2.54), r.implausible_bd))
kh_df = pd.DataFrame(_kh_rows, columns=["texture_class", "bulk_density_g_cm3", "suction_cm", "k_in_hr", "implausible_bd"])
kh_df["texture_class"] = pd.Categorical(kh_df["texture_class"], categories=list(TEXTURE_CLASSES), ordered=True)

_kh_markers = (
    hv.VLine(330).opts(color="black", line_dash="dotted", line_width=1)
    * hv.VLine(15000).opts(color="gray", line_dash="dotted", line_width=1)
    * hv.Text(330, 30, " FC", halign="left", valign="top").opts(text_font_size="8pt")
    * hv.Text(15000, 30, "WP ", halign="right", valign="top").opts(text_font_size="8pt")
)
(
    bd_slider_line(
        kh_df, "suction_cm", "k_in_hr",
        "suction head h (cm, log)", "unsaturated K (in/hr, log)",
        "Rosetta unsaturated hydraulic conductivity K(h) — {dimensions}",
        logx=True, logy=True, ylim=(1e-7, 1e2),
    )
    * _kh_markers
).redim(bulk_density_g_cm3=hv.Dimension("Bulk density, g/cm³ (higher is more compacted)", default=1.4, value_format=lambda v: f"{v:.1f}"))

# %% [markdown]
# ::: {.callout-tip appearance="simple"}
# **Takeaway:** Between field capacity and wilting point — the range where plants can extract water — conductivity is already thousands of times lower than at saturation; compaction compresses this window further, making it harder for both water and roots to move through the soil profile.
# :::

# %% [markdown]
# ## 3. Green–Ampt infiltration
#
# Infiltration starts fast — especially in dry, coarse, or healthy soils — and slows toward a
# steady rate as the wetting front advances. The two charts below show the infiltration rate f(t)
# and cumulative depth F(t) over the first two hours of ponded conditions. Use the bulk-density
# slider to see how compaction cuts both curves dramatically, increasing runoff risk.
#
# ::: {.callout-note collapse="true"}
# ## For researchers: Green–Ampt parameterization
# The **Green–Ampt** model gives the infiltration rate f and cumulative depth F for ponded /
# intense-rain conditions:
#
# $$f = K_s\left(1 + \frac{\psi_f\,\Delta\theta}{F}\right), \qquad t = \frac{1}{K_s}\left[F - \psi_f\Delta\theta\,\ln\!\left(1 + \frac{F}{\psi_f\Delta\theta}\right)\right]$$
#
# **Parameterization:** Rosetta supplies the BD-sensitive **Ksat** and the **moisture deficit**
# Δθ = θₛ − θ_initial (initial = wilting point here, i.e. dry antecedent). The wetting-front suction
# **ψ_f** is taken from the **Rawls, Brakensiek & Miller (1983)** texture table (cm) — deriving ψ_f
# from Rosetta's Mualem K(h) integral proved unreliable (Rosetta's fitted L mis-orders the
# capillary drive, putting sand above clay), so the published textural values are more robust.
#
# Caveats: matrix-only (real infiltration is often higher via macropores/structure); Ksat is
# Rosetta's least-certain output; high-BD × texture extrapolations are grey-dashed.
# :::

# %%
hv.output(widget_location="bottom")

# Rawls, Brakensiek & Miller (1983) Green–Ampt wetting-front suction head, cm (silt ≈ silt loam)
RAWLS_PSI_F_CM = {
    "sand": 4.95, "loamy sand": 6.13, "sandy loam": 11.01, "loam": 8.89,
    "silt loam": 16.68, "silt": 16.68, "sandy clay loam": 21.85, "clay loam": 20.88,
    "silty clay loam": 27.30, "sandy clay": 23.90, "silty clay": 29.22, "clay": 31.63,
}
GA_T_MAX_HR = 2.0  # plot the first 2 hours
# Geometric (log) spacing keeps the early-time / small-F detail dense — where slow soils
# spend the whole 2-hour window — while using far fewer points than a uniform grid.
_F_grid = np.geomspace(0.2, 80.0, 120)  # cumulative infiltration, cm

_ga_rows = []
for r in result.itertuples():
    Ks = r.ksat_cm_day  # cm/day
    dtheta = max(r.total_porosity - r.wilting_point_porosity, 1e-6)
    Sf = RAWLS_PSI_F_CM[r.texture_class] * dtheta  # cm
    t_hr = (_F_grid - Sf * np.log1p(_F_grid / Sf)) / Ks * 24.0
    f_in_hr = Ks * (1.0 + Sf / _F_grid) / (24.0 * 2.54)
    F_in = _F_grid / 2.54
    keep = t_hr <= GA_T_MAX_HR
    for th, fi, Fi in zip(t_hr[keep], f_in_hr[keep], F_in[keep]):
        _ga_rows.append((r.texture_class, r.bulk_density_g_cm3, th, fi, Fi, r.implausible_bd))
ga_df = pd.DataFrame(_ga_rows, columns=["texture_class", "bulk_density_g_cm3", "t_hr", "f_in_hr", "F_in", "implausible_bd"])
ga_df["texture_class"] = pd.Categorical(ga_df["texture_class"], categories=list(TEXTURE_CLASSES), ordered=True)

bd_slider_line(
    ga_df, "t_hr", "f_in_hr",
    "time (hours)", "infiltration rate f (in/hr, log)",
    "Green–Ampt infiltration rate vs. time — {dimensions}",
    logy=True, ylim=(0.05, 300),
).redim(bulk_density_g_cm3=hv.Dimension("Bulk density, g/cm³ (higher is more compacted)", default=1.4, value_format=lambda v: f"{v:.1f}"))

# %%
# Cumulative infiltration depth F(t).
hv.output(widget_location="bottom")
bd_slider_line(
    ga_df, "t_hr", "F_in",
    "time (hours)", "cumulative infiltration F (inches)",
    "Green–Ampt cumulative infiltration vs. time — {dimensions}",
    ylim=(0, 25),
).redim(bulk_density_g_cm3=hv.Dimension("Bulk density, g/cm³ (higher is more compacted)", default=1.4, value_format=lambda v: f"{v:.1f}"))

# %% [markdown]
# ::: {.callout-tip appearance="simple"}
# **Takeaway:** Infiltration starts fast and settles toward the saturated rate — but compaction lowers the ceiling, so a compacted soil ponds and runs off far sooner than a healthy one of the same texture, even under the same rainfall intensity.
# :::
