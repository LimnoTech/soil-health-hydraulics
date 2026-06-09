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
# # Rosetta hydraulic conductivity & infiltration
#
# Companion to **`1_rosetta_porosity_by_texture.ipynb`** — reads `rosetta_porosity_by_texture.csv`
# (**run Notebook 1 first**), which carries Rosetta's saturated **Ksat** and the **Mualem–van
# Genuchten** unsaturated-conductivity parameters (`k0_cm_day`, `mualem_L`, `vg_alpha_1cm`, `vg_n`).
#
# - **§1 Saturated hydraulic conductivity (Ksat)** vs. bulk density.
# - **§2 Unsaturated hydraulic conductivity K(h)** from the Mualem–van Genuchten parameters.
# - **§3 Green–Ampt infiltration** f(t) / F(t).
#
# Stormwater units (in/hr); bulk-density sliders; physically implausible BD × texture combinations
# (`implausible_bd`, θₛ > 1 − BD/2.65) are greyed as extrapolation, as in Notebook 1.

# %%
import numpy as np
import pandas as pd
import hvplot.pandas  # noqa: F401  (registers the .hvplot accessor)
import holoviews as hv

pd.set_option("display.float_format", lambda v: f"{v:0.3f}")

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


def mualem_k(h, alpha, n, k0, L):
    """Unsaturated hydraulic conductivity K at suction head h [cm], Mualem–van Genuchten
    (Schaap & Leij, 2000), in the same units as k0 (Rosetta: cm/day):

        Se = [1 + (alpha h)^n]^(-m),  m = 1 - 1/n
        K(Se) = k0 * Se^L * [1 - (1 - Se^(1/m))^m]^2

    Use Rosetta's matching-point k0 (col 5) and exponent L (col 6, often negative) — not Ksat.
    """
    h = np.abs(np.asarray(h, dtype=float))
    alpha = np.asarray(alpha, dtype=float)
    n = np.asarray(n, dtype=float)
    m = 1.0 - 1.0 / n
    Se = (1.0 + (alpha * h) ** n) ** (-m)
    return k0 * Se**L * (1.0 - (1.0 - Se ** (1.0 / m)) ** m) ** 2


# --- plotting helpers (mirrors Notebook 1) ---
plot_opts = dict(
    x="bulk_density_g_cm3", by="texture_class", xlabel="Bulk density (g/cm³)",
    width=750, height=500, legend="right", grid=True,
)
# boundary row kept in both series so solid and grey-dashed (extrapolation) segments join
_next_implausible = result.groupby("texture_class", sort=False)["implausible_bd"].shift(-1).fillna(False)
_extrap_mask = result["implausible_bd"] | _next_implausible


def line_with_extrapolation(ycol, ylabel, title, **extra):
    """x = bulk density; solid over plausible BD, grey-dashed where implausible_bd."""
    solid = result.assign(**{ycol: result[ycol].where(~result["implausible_bd"])}).hvplot.line(
        y=ycol, ylabel=ylabel, title=title, **plot_opts, **extra
    )
    dashed = (
        result.assign(**{ycol: result[ycol].where(_extrap_mask)})
        .hvplot.line(y=ycol, **plot_opts, **extra)
        .opts(hv.opts.Curve(color="lightgray", line_dash="dashed", alpha=0.9))
        .opts(show_legend=False)
    )
    return solid * dashed


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
result.head()

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
    "ksat_in_hr",
    "Saturated hydraulic conductivity Ksat (in/hr, log scale)",
    "Rosetta saturated hydraulic conductivity vs. bulk density by USDA texture class",
    logy=True,
)

# %% [markdown]
# ## 2. Unsaturated hydraulic conductivity K(h)
#
# Rosetta gives the full **Mualem–van Genuchten** parameter set, so conductivity is defined not
# just at saturation but at every suction: K(h) = K0·Se(h)^L·[1 − (1 − Se^(1/m))^m]² with
# Se(h) = [1 + (αh)^n]^(−m) (see `mualem_k`). This is what governs **unsaturated** flow — as a soil
# dries, K drops by many orders of magnitude. Curves use Rosetta's **K0** (matching point) and **L**
# (columns 5–6), *not* Ksat. Dotted verticals mark field capacity (33 kPa = 330 cm) and wilting
# point (1500 kPa = 15000 cm); `k_fc_cm_day` / `k_wp_cm_day` in the table are K there. Log–log axes;
# bulk-density slider; implausible (`implausible_bd`) BD × texture combinations are grey-dashed.

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
).redim(bulk_density_g_cm3=hv.Dimension("Bulk density (g/cm³)", value_format=lambda v: f"{v:.1f}"))

# %% [markdown]
# ## 3. Green–Ampt infiltration
#
# Infiltration is a *process*, not just a conductivity. The **Green–Ampt** model gives the
# infiltration rate f and cumulative depth F for ponded / intense-rain conditions:
#
# $$f = K_s\left(1 + \frac{\psi_f\,\Delta\theta}{F}\right), \qquad t = \frac{1}{K_s}\left[F - \psi_f\Delta\theta\,\ln\!\left(1 + \frac{F}{\psi_f\Delta\theta}\right)\right]$$
#
# **Parameterization:** Rosetta supplies the BD-sensitive **Ksat** and the **moisture deficit**
# Δθ = θₛ − θ_initial (initial = wilting point here, i.e. dry antecedent). The wetting-front suction
# **ψ_f** is taken from the **Rawls, Brakensiek & Miller (1983)** texture table (cm) — deriving ψ_f
# from Rosetta's Mualem K(h) integral proved unreliable (Rosetta's fitted L mis-orders the
# capillary drive, putting sand above clay), so the published textural values are more robust.
# The rate starts high and decays toward Ksat as the wetting front advances.
#
# Caveats: matrix-only (real infiltration is often higher via macropores/structure); Ksat is
# Rosetta's least-certain output; high-BD × texture extrapolations are grey-dashed.

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
).redim(bulk_density_g_cm3=hv.Dimension("Bulk density (g/cm³)", value_format=lambda v: f"{v:.1f}"))

# %%
# Cumulative infiltration depth F(t).
hv.output(widget_location="bottom")
bd_slider_line(
    ga_df, "t_hr", "F_in",
    "time (hours)", "cumulative infiltration F (inches)",
    "Green–Ampt cumulative infiltration vs. time — {dimensions}",
    ylim=(0, 25),
).redim(bulk_density_g_cm3=hv.Dimension("Bulk density (g/cm³)", value_format=lambda v: f"{v:.1f}"))
