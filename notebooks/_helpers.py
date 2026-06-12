"""Shared helpers for the soil-health notebooks (1, 2, 3).

Imported as ``from _helpers import ...``. The notebooks execute with their working directory set
to ``notebooks/`` (Quarto ``execute-dir: file``), so this sibling module is importable. The leading
underscore also keeps Quarto from treating it as a page to render.

Importing this module sets the pandas ``display.float_format`` that all three notebooks share
(3 decimals) — the same setting ``show`` relies on for its HTML table rendering.
"""

import numpy as np
import pandas as pd
import hvplot.pandas  # noqa: F401  (registers the .hvplot accessor used by line_with_extrapolation)
import holoviews as hv
from IPython.display import HTML

# All three notebooks display volumetric water contents to 3 decimals.
pd.set_option("display.float_format", lambda v: f"{v:0.3f}")


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
