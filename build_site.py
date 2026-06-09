#!/usr/bin/env python3
"""Build the static website published to GitHub Pages.

Each notebook is rendered with ``jupyter nbconvert --to html --no-input`` (code
hidden). nbconvert preserves the notebooks' *self-contained* interactive
Bokeh/HoloViews output verbatim, so the bulk-density / organic-matter sliders,
hover, pan, and zoom all work client-side with no kernel or server. (Quarto's
pandoc pipeline strips those embeds for pre-executed notebooks, which is why we
render with nbconvert directly.)

A shared top navigation bar and matching styles are injected into every page,
and a landing page (``index.html``) is generated, so the pages read as one site.

Output: ``_site/`` (index.html + notebooks/*.html) -> uploaded by the GitHub
Action to GitHub Pages.

Run locally:  ``pixi run build-site``
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.resolve()
SITE = ROOT / "_site"
NB_SRC = ROOT / "notebooks"
NB_OUT = SITE / "notebooks"

SITE_TITLE = "Soil Health & Water Storage"
REPO_URL = "https://github.com/LimnoTech/soil-health-hydraulics"

# (filename stem, short nav label, landing-page heading, landing-page blurb)
NOTEBOOKS = [
    (
        "1_rosetta_porosity_by_texture",
        "1 · Porosity by texture",
        "1 · Porosity by texture &amp; bulk density",
        "Total porosity, <strong>field capacity</strong>, and <strong>wilting point</strong> "
        "for every USDA texture class across bulk densities 0.8–1.9 g/cm³ — with "
        "plant-available and drainable water in inches per foot of soil, plus an "
        "FAO-style storage diagram. The foundation for the other two notebooks.",
    ),
    (
        "2_organic_matter_water_holding",
        "2 · Organic matter",
        "2 · Organic-matter effects",
        "How adding <strong>organic matter</strong> (0–8%, up to 15%) shifts available and "
        "drainable water, via two independent methods — Rosetta + Minasny &amp; McBratney "
        "(2018), and Saxton–Rawls (2006) — with a mineral-baseline calibration from the "
        "UNSODA database.",
    ),
    (
        "3_rosetta_hydraulic_conductivity",
        "3 · Conductivity & infiltration",
        "3 · Conductivity &amp; infiltration",
        "Saturated <strong>Ksat</strong>, unsaturated <strong>K(h)</strong> from the "
        "Mualem–van Genuchten curve, and <strong>Green–Ampt infiltration</strong> rates — "
        "all in stormwater units (in/hr), with bulk-density sliders.",
    ),
]

# Shared look-and-feel injected into every page (notebook pages + landing page).
SITE_CSS = """
<style>
  :root { --shw-bar: #2c6e8f; --shw-bar-text: #ffffff; --shw-accent: #1b4a61; }
  body { margin: 0; }
  .shw-nav {
    position: sticky; top: 0; z-index: 1000;
    display: flex; flex-wrap: wrap; align-items: center; gap: .25rem 1.25rem;
    background: var(--shw-bar); color: var(--shw-bar-text);
    padding: .55rem 1.25rem;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    font-size: 15px; box-shadow: 0 1px 4px rgba(0,0,0,.18);
  }
  .shw-nav .shw-brand { font-weight: 700; margin-right: .5rem; }
  .shw-nav a { color: var(--shw-bar-text); text-decoration: none; opacity: .9; }
  .shw-nav a:hover { opacity: 1; text-decoration: underline; }
  .shw-nav a.shw-active { text-decoration: underline; font-weight: 600; opacity: 1; }
  .shw-nav .shw-spacer { flex: 1 1 auto; }
</style>
"""


def nav_html(prefix: str, active: str) -> str:
    """Top navigation bar. ``prefix`` makes links relative to the page location
    ('' for the root index, '../' for pages under notebooks/). ``active`` is the
    href of the current page so it can be highlighted."""
    links = [("index.html", "Home")] + [
        (f"notebooks/{stem}.html", label) for stem, label, _, _ in NOTEBOOKS
    ]
    items = [f'<span class="shw-brand">{SITE_TITLE}</span>']
    for href, label in links:
        cls = " class=\"shw-active\"" if href == active else ""
        items.append(f'<a href="{prefix}{href}"{cls}>{label}</a>')
    items.append('<span class="shw-spacer"></span>')
    items.append(f'<a href="{REPO_URL}">GitHub ↗</a>')
    return '<nav class="shw-nav">' + "\n  ".join(items) + "</nav>"


def inject(html: str, prefix: str, active: str) -> str:
    """Insert the shared CSS into <head> and the nav bar just after <body>."""
    if "</head>" in html:
        html = html.replace("</head>", SITE_CSS + "</head>", 1)
    # Insert nav after the opening <body ...> tag.
    body_start = html.find("<body")
    if body_start != -1:
        tag_end = html.find(">", body_start)
        if tag_end != -1:
            insert_at = tag_end + 1
            html = html[:insert_at] + "\n" + nav_html(prefix, active) + html[insert_at:]
    return html


def render_notebooks() -> None:
    NB_OUT.mkdir(parents=True, exist_ok=True)
    for stem, _, _, _ in NOTEBOOKS:
        src = NB_SRC / f"{stem}.ipynb"
        if not src.exists():
            sys.exit(f"ERROR: missing notebook {src}")
        print(f"  rendering {src.name} …")
        subprocess.run(
            [
                "jupyter", "nbconvert", "--to", "html", "--no-input",
                "--output-dir", str(NB_OUT),
                "--output", f"{stem}.html",
                str(src),
            ],
            check=True,
        )
        out = NB_OUT / f"{stem}.html"
        html = out.read_text(encoding="utf-8")
        active = f"notebooks/{stem}.html"
        out.write_text(inject(html, prefix="../", active=active), encoding="utf-8")


def cards_html() -> str:
    cards = []
    for stem, _, heading, blurb in NOTEBOOKS:
        cards.append(
            f'<a class="shw-card" href="notebooks/{stem}.html">'
            f"<h3>{heading}</h3><p>{blurb}</p></a>"
        )
    return "\n".join(cards)


def write_index() -> None:
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{SITE_TITLE}</title>
{SITE_CSS}
<style>
  .shw-main {{
    max-width: 960px; margin: 0 auto; padding: 1.5rem 1.25rem 3rem;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    line-height: 1.55; color: #1c2b33;
  }}
  .shw-main h1 {{ font-size: 1.9rem; line-height: 1.2; }}
  .shw-main h1 .shw-sub {{ display:block; font-size: 1.05rem; font-weight: 400; color:#557; margin-top:.4rem; }}
  .shw-tip {{ background:#eef6fb; border-left:4px solid var(--shw-bar); padding:.8rem 1rem; border-radius:4px; margin:1.5rem 0; }}
  .shw-cards {{ display:grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap:1rem; margin:1.5rem 0; }}
  .shw-card {{ display:block; text-decoration:none; color:inherit; border:1px solid #d6e1e8; border-radius:8px; padding:1rem 1.1rem; transition:box-shadow .15s, border-color .15s; }}
  .shw-card:hover {{ box-shadow:0 3px 12px rgba(0,0,0,.1); border-color: var(--shw-bar); }}
  .shw-card h3 {{ margin:.1rem 0 .5rem; color: var(--shw-accent); font-size:1.08rem; }}
  .shw-card p {{ margin:0; font-size:.93rem; color:#3a4a52; }}
  .shw-note {{ font-size:.9rem; color:#557; border-top:1px solid #e3e9ec; margin-top:2rem; padding-top:1rem; }}
  .shw-main a {{ color: var(--shw-accent); }}
</style>
</head>
<body>
{nav_html(prefix="", active="index.html")}
<main class="shw-main">
  <h1>Soil Health Impacts on Water Storage &amp; Hydraulic Conductivity
    <span class="shw-sub">Interactive results from the Rosetta pedotransfer functions</span>
  </h1>

  <p>How do <strong>soil texture</strong>, <strong>bulk density</strong> (compaction), and
  <strong>organic matter</strong> control how much water a soil can store and how fast water
  moves through it? These pages turn the
  <a href="https://github.com/usda-ars-ussl/rosetta-soil">Rosetta</a> pedotransfer functions
  into <strong>interactive charts</strong> you can explore right in your browser.</p>

  <div class="shw-tip">
    <strong>How to explore:</strong> every chart is live — <strong>drag the sliders</strong>
    (bulk density, organic matter) and <strong>hover</strong> for exact values; use the toolbar
    to <strong>zoom and pan</strong>. Nothing to install. Grey, dashed regions flag physically
    implausible texture × bulk-density combinations (extrapolation beyond Rosetta's training range).
  </div>

  <h2>The notebooks</h2>
  <div class="shw-cards">
{cards_html()}
  </div>

  <p class="shw-note">All results come from the <strong>Rosetta v3</strong> pedotransfer
  functions, which predict van Genuchten–Mualem retention and conductivity parameters from
  sand/silt/clay percentages and bulk density (texture classes use representative <em>median</em>
  values). These are <strong>model estimates</strong> for matrix (non-macropore) flow, intended
  for relative comparison and planning — not a substitute for site-specific measurements. See
  each notebook for full methodology and caveats, or the
  <a href="{REPO_URL}#readme">project README</a>.</p>
</main>
</body>
</html>
"""
    (SITE / "index.html").write_text(html, encoding="utf-8")


def main() -> None:
    if SITE.exists():
        shutil.rmtree(SITE)
    SITE.mkdir(parents=True)
    print(f"Building site into {SITE} …")
    render_notebooks()
    write_index()
    # Disable Jekyll so GitHub Pages serves files (e.g. site_libs) as-is.
    (SITE / ".nojekyll").write_text("", encoding="utf-8")
    print("Done.")


if __name__ == "__main__":
    main()
