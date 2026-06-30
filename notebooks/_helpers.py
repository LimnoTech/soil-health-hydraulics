"""Shared helpers for the soil-health notebooks (1, 2, 3).

Imported as ``from _helpers import ...``. The notebooks execute with their working directory set
to ``notebooks/`` (Quarto ``execute-dir: file``), so this sibling module is importable. The leading
underscore also keeps Quarto from treating it as a page to render.

Importing this module has three global side-effects all three notebooks rely on: it sets the
pandas ``display.float_format`` (3 decimals, used by ``show``), registers the HoloViews **Bokeh**
backend, and sets a wheel-zoom-off ``hv.opts.defaults(... active_tools=[])`` so scroll doesn't
zoom a figure until the user toggles it.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import hvplot.pandas  # noqa: F401  (registers the .hvplot accessor used by line_with_extrapolation)
import holoviews as hv
from bokeh.models import HoverTool
from IPython.display import HTML

# All three notebooks display volumetric water contents to 3 decimals.
pd.set_option("display.float_format", lambda v: f"{v:0.3f}")

# Wheel zoom stays in the toolbar but is INACTIVE on load on every figure (users toggle it on).
hv.extension("bokeh")  # ensure the Bokeh backend is registered before setting opts defaults
hv.opts.defaults(
    hv.opts.Curve(active_tools=[]),
    hv.opts.Area(active_tools=[]),
    hv.opts.Scatter(active_tools=[]),
    hv.opts.Bars(active_tools=[]),
    hv.opts.Overlay(active_tools=[]),
    hv.opts.Rectangles(active_tools=[]),
)


def show(df, height=360):
    """Display the *full* DataFrame in a fixed-height, scrollable box — renders the same in
    JupyterLab and in the exported HTML site (`to_html` emits every row and respects the
    float_format set above)."""
    return HTML(
        "<style>.scroll-df thead th{position:sticky;top:0;background:#fff;"
        "box-shadow:inset 0 -1px 0 #ccc;}</style>"
        f'<div class="scroll-df" style="max-height:{height}px;overflow:auto;'
        'border:1px solid #ddd;border-radius:4px;">'
        f"{df.to_html()}</div>"
    )


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


def mualem_k(h, alpha, n, k0, L):
    """Unsaturated hydraulic conductivity K at suction head h [cm], Mualem–van Genuchten
    (Schaap & Leij, 2000), in the same units as k0 (Rosetta: cm/day):

        Se = [1 + (alpha h)^n]^(-m),  m = 1 - 1/n
        K(Se) = k0 * Se^L * [1 - (1 - Se^(1/m))^m]^2

    k0 = matching-point conductivity (Rosetta col 5), L = pore-connectivity exponent
    (col 6, often negative). At h = 0 (Se = 1) this returns k0 — note Rosetta's k0 is the
    fitted matching point and is typically < Ksat (col 4). Scalars or numpy arrays.
    """
    h = np.abs(np.asarray(h, dtype=float))
    alpha = np.asarray(alpha, dtype=float)
    n = np.asarray(n, dtype=float)
    m = 1.0 - 1.0 / n
    Se = (1.0 + (alpha * h) ** n) ** (-m)
    return k0 * Se**L * (1.0 - (1.0 - Se ** (1.0 / m)) ** m) ** 2


# Default options for line_with_extrapolation: x = bulk density, one line per texture class.
LINE_PLOT_OPTS = dict(
    x="bulk_density_g_cm3",
    by="texture_class",
    xlabel="Bulk density (g/cm³)",
    width=750,
    height=500,
    legend="right",
    grid=True,
)


def line_with_extrapolation(result, ycol, ylabel, title, plot_opts=None, **extra):
    """Line plot of `ycol` vs. bulk density, one line per texture class, drawn solid over
    physically plausible bulk densities and **grey-dashed** where ``implausible_bd`` is True
    (Rosetta θₛ > BD-implied pore space — extrapolation). The boundary row is kept in both
    series so the solid and dashed segments join.

    `result` must carry ``texture_class``, ``bulk_density_g_cm3``, ``implausible_bd`` and `ycol`.
    `plot_opts` overrides the LINE_PLOT_OPTS defaults; `extra` (e.g. ``logy=True``) is forwarded
    to ``hvplot.line``.
    """
    if plot_opts is None:
        plot_opts = LINE_PLOT_OPTS
    next_implausible = result.groupby("texture_class", sort=False)["implausible_bd"].shift(-1).fillna(False)
    extrap_mask = result["implausible_bd"] | next_implausible
    solid = result.assign(**{ycol: result[ycol].where(~result["implausible_bd"])}).hvplot.line(
        y=ycol, ylabel=ylabel, title=title, **plot_opts, **extra
    )
    dashed = (
        result.assign(**{ycol: result[ycol].where(extrap_mask)})
        .hvplot.line(y=ycol, **plot_opts, **extra)
        .opts(hv.opts.Curve(color="lightgray", line_dash="dashed", alpha=0.9))
        .opts(show_legend=False)
    )
    return solid * dashed


# --- Minasny & McBratney (2018) organic-matter blend (shared by NB2 and the home page) ---
VB = 0.58  # van Bemmelen factor: OM ≈ OC / VB
OC_BASELINE_PCT = 1.0  # reference mineral-baseline OC (%) at the mean mineral BD ≈ 1.4 (NB2 §1)

# ROSETTA mineral-baseline organic carbon as a FUNCTION of bulk density, from the UNSODA 2.0
# all-mineral OC~BD regression (NB2 §1): OC% = OC_BD_SLOPE·BD + OC_BD_INTERCEPT (r ≈ −0.63).
# The M&M increments are applied relative to this BD-dependent baseline, so the ROSETTA prediction
# carries the OC mineral soils typically have at that BD (≈ 2.9% at BD 0.8, 1.0% at BD 1.4, 0 by
# BD ≈ 1.72). NB2 §1 recomputes these from UNSODA and asserts they still match these constants.
OC_BD_SLOPE = -3.1602      # %OC per g/cm³
OC_BD_INTERCEPT = 5.4266   # %OC at BD = 0 (extrapolated intercept)


def oc_baseline_for_bd(bd):
    """ROSETTA mineral-baseline organic carbon (%) at bulk density `bd` (g/cm³), from the UNSODA
    all-mineral OC~BD regression (NB2 §1), floored at 0 — negative OC is nonsensical and the fit
    crosses zero near BD ≈ 1.72. At the mean mineral BD ≈ 1.4 this ≈ ``OC_BASELINE_PCT``."""
    return np.maximum(OC_BD_SLOPE * np.asarray(bd, dtype=float) + OC_BD_INTERCEPT, 0.0)

# Table 2 slopes, mm H2O per 100 mm per +1% OC (= vol %), by USDA texture group
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


def soil_water_bd_om_blend_table(result, bd, om, oc_baseline=None):
    """ROSETTA mineral baseline + Minasny & McBratney (2018) organic-matter increments for a
    single (bulk density `bd`, organic matter `om` %) state — one row per texture class in the
    canonical sand→clay order of `result`. Volumetric water contents (cm³/cm³).

    Applies M&M's WP/AWC/SAT slopes relative to `oc_baseline` and derives FC = WP + AWC so AWC
    matches M&M exactly (mirrors NB2's _blend_profile). `oc_baseline` defaults to the BD-dependent
    mineral baseline ``oc_baseline_for_bd(bd)`` (the UNSODA all-mineral OC~BD fit); pass a scalar to
    override. `implausible` flags blended SAT exceeding the BD-implied pore space (1 − BD/2.65). The
    single source feeding both the FAO figure and the home-page / NB2 linked table.
    """
    if oc_baseline is None:
        oc_baseline = float(oc_baseline_for_bd(bd))
    base = result.set_index(["texture_class", "bulk_density_g_cm3"])
    texture_classes = list(result["texture_class"].drop_duplicates())
    hsg = dict(
        result[["texture_class", "hydrologic_soil_group"]]
        .drop_duplicates().itertuples(index=False, name=None)
    )
    d_oc = om * VB - oc_baseline  # OM% -> OC%, relative to the BD-dependent mineral-baseline OC
    rows = []
    for cls in texture_classes:
        s = MM_SLOPES[MM_GROUP[cls]]
        sat0 = base.loc[(cls, bd), "total_porosity"]
        fc0 = base.loc[(cls, bd), "field_capacity_porosity"]
        wp0 = base.loc[(cls, bd), "wilting_point_porosity"]
        wp = max(wp0 + s["WP"] / 100 * d_oc, 0.0)
        awc = max((fc0 - wp0) + s["AWC"] / 100 * d_oc, 0.0)
        fc = wp + awc
        sat = max(sat0 + s["SAT"] / 100 * d_oc, fc)
        rows.append({
            "texture_class": cls,
            "hydrologic_soil_group": hsg[cls],
            "wilting_point_porosity": wp,
            "field_capacity_porosity": fc,
            "total_porosity": sat,
            "available_water_capacity": awc,
            "drainable_water": sat - fc,
            "implausible": sat > (1.0 - bd / 2.65),
        })
    df = pd.DataFrame(rows)
    df["texture_class"] = pd.Categorical(df["texture_class"], categories=texture_classes, ordered=True)
    return df


def soil_water_texture_band_diagram(x, pwp, fc, por, *, texture_labels=None, implausible=None, hover=True):
    """FAO-style soil-water band diagram for ONE profile (one slider frame).

    `x` = texture x-positions (0..11); `pwp`/`fc`/`por` = wilting point / field capacity / total
    porosity arrays **already in the plot's y-units** (inches per foot here). `texture_labels` =
    per-column texture names shown first in the hover tooltip (defaults to "texture <i>" if omitted).
    Returns an hv.Overlay of: filled bands (orange unavailable / green available / blue drainable),
    the three boundary curves, the three text labels, optional grey extrapolation spans
    (`implausible` boolean array), and (if `hover`) an invisible per-texture hover layer reporting
    the texture plus available / drainable / total stormwater capacity. Visible geometry is identical
    to the previous per-notebook _*_profile code.
    """
    x = np.asarray(x); pwp = np.asarray(pwp); fc = np.asarray(fc); por = np.asarray(por)
    bands = (
        hv.Area((x, pwp, pwp * 0), vdims=["y", "y2"]).opts(color="orange", alpha=0.45, line_alpha=0)
        * hv.Area((x, fc, pwp), vdims=["y", "y2"]).opts(color="green", alpha=0.45, line_alpha=0)
        * hv.Area((x, por, fc), vdims=["y", "y2"]).opts(color="blue", alpha=0.40, line_alpha=0)
    )
    lines = (
        hv.Curve((x, por), label="Saturation").opts(color="gray", line_width=2, line_dash="solid")
        * hv.Curve((x, fc), label="Field Capacity").opts(color="black", line_width=2, line_dash="dashed")
        * hv.Curve((x, pwp), label="Permanent Wilting Point").opts(color="black", line_width=2, line_dash="dotted")
    )
    labels = (
        hv.Text(8, pwp[8] * 0.5, "Unavailable\nwater").opts(text_color="saddlebrown", text_font_size="12pt")
        * hv.Text(5, (pwp[5] + fc[5]) / 2, "Available water").opts(text_color="darkgreen", text_font_size="12pt")
        * hv.Text(2.2, (fc[2] + por[2]) / 2, "Drainable\nwater").opts(text_color="navy", text_font_size="12pt")
    )
    overlay = bands * lines * labels
    if implausible is not None:
        flagged = x[np.asarray(implausible, dtype=bool)]
        if len(flagged):
            overlay = overlay * hv.Overlay(
                [hv.VSpan(xi - 0.5, xi + 0.5).opts(color="gray", alpha=0.2) for xi in flagged]
            )
    if hover:
        available = fc - pwp
        drainable = por - fc
        total = por - pwp
        if texture_labels is None:
            texture_labels = [f"texture {int(xi)}" for xi in x]
        rects = [
            (xi - 0.5, 0.0, xi + 0.5, max(por_i, 1e-6), tx, av, dr, tot)
            for xi, por_i, tx, av, dr, tot in zip(x, por, texture_labels, available, drainable, total)
        ]
        hover_tool = HoverTool(tooltips=[
            ("Texture", "@texture"),
            ("Available water", "@available{0.00} in/ft"),
            ("Drainable water", "@drainable{0.00} in/ft"),
            ("Total stormwater capacity", "@total{0.00} in/ft"),
        ])
        hover_layer = hv.Rectangles(rects, vdims=["texture", "available", "drainable", "total"]).opts(
            fill_alpha=0, line_alpha=0, tools=[hover_tool]
        )
        overlay = overlay * hover_layer
    # Pin the x-range to the texture extent (sand at the left edge, clay at the right edge),
    # removing Bokeh's default auto-range padding. The ±0.5 hover rects / VSpans clip at the edges.
    # Also bump the title / axis-label / tick / legend fonts up from Bokeh's small defaults.
    overlay = overlay.opts(
        hv.opts.Overlay(
            xlim=(float(np.min(x)), float(np.max(x))),
            fontsize={"title": 14, "labels": 12, "xticks": 11, "yticks": 11, "legend": 12},
        )
    )
    return overlay


def soil_water_table_html(df):
    """Styled HTML for the slider-linked storage table, for use inside a Panel ``pn.pane.HTML``
    (Panel embeds it per widget state, so it stays slider-linked — unlike a Bokeh DataTable, this
    gives full CSS control). Formatting: **bold, word-wrapped** column headers; data cells and the
    first (texture) column do NOT wrap, so each texture label stays on one line; the row index is
    hidden; all rows render with no scrolling."""
    table_html = df.to_html(index=False, border=0, classes="swt")
    css = (
        "<style>"
        "table.swt{border-collapse:collapse;font-size:0.95rem;font-family:inherit;}"
        "table.swt th{font-weight:700;white-space:normal;vertical-align:bottom;text-align:right;"
        "padding:4px 9px;border-bottom:1px solid #bbb;max-width:96px;}"
        "table.swt th:first-child{text-align:left;white-space:nowrap;max-width:none;}"
        "table.swt td{white-space:nowrap;text-align:right;padding:3px 9px;}"
        "table.swt td:first-child{text-align:left;}"
        "</style>"
    )
    return css + table_html


# Canonical water-storage columns shared by the home-page table and the per-BD CSV exports:
# source column -> exported (display) name. drainable_water is derived (SAT − FC) if absent.
WATER_STORAGE_COLUMNS = {
    "wilting_point_porosity": "wilting point",
    "field_capacity_porosity": "field capacity",
    "total_porosity": "saturation",
    "available_water_capacity": "available water",
    "drainable_water": "drainable water",
}


def export_water_storage_tables(df, outdir, *, prefix="rosetta", round_to=3):
    """Write one CSV snapshot per bulk-density step — rows = ``texture (HSG)``, columns =
    wilting point / field capacity / saturation / available water / drainable water (cm³/cm³,
    rounded to `round_to`). Files are named ``{prefix}_bd_{bd:.1f}.csv`` in `outdir` (created if
    needed). Reusable across notebooks (ROSETTA baseline in NB1, blended outputs in NB2).

    `df` must carry ``texture_class``, ``hydrologic_soil_group``, ``bulk_density_g_cm3`` and the
    porosity columns in ``WATER_STORAGE_COLUMNS``; ``drainable_water`` is derived as
    total_porosity − field_capacity_porosity when not already present. Texture row order follows
    `df` (canonical sand → clay). Returns the list of written ``Path``s.
    """
    df = df.copy()
    if "drainable_water" not in df.columns:
        df["drainable_water"] = df["total_porosity"] - df["field_capacity_porosity"]
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    texture_order = list(df["texture_class"].drop_duplicates())
    written = []
    for bd, g in df.groupby("bulk_density_g_cm3", sort=True):
        g = g.set_index("texture_class").loc[texture_order].reset_index()  # canonical sand → clay
        table = pd.DataFrame(
            {"texture (HSG)": g["texture_class"].astype(str) + " (" + g["hydrologic_soil_group"].astype(str) + ")"}
        )
        for src, name in WATER_STORAGE_COLUMNS.items():
            table[name] = g[src].to_numpy().round(round_to)
        path = outdir / f"{prefix}_bd_{bd:.1f}.csv"
        table.to_csv(path, index=False)
        written.append(path)
    return written
