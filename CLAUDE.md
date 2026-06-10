# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A small scientific-computing project: three Jupyter notebooks that turn the **Rosetta v3**
pedotransfer functions into interactive soil water-storage / hydraulic-conductivity charts for a
stormwater & soil-health audience, published as a static interactive website on GitHub Pages.
There is no application code or test suite — the "product" is the executed notebooks and the site
built from them.

## Environment & commands

Environment is managed by **pixi** (`pixi.toml`); never use bare `pip`/`conda`.

```bash
pixi install                      # create/refresh the env from pixi.lock
pixi run jupyter lab              # work on notebooks interactively
pixi run render                   # Quarto: execute + freeze + build _site/ (the published website)
pixi run preview                  # Quarto live-reload preview server
```

Notebook edit/run cycle (the `.py` is the source of truth to edit; the `.ipynb` is the generated pair):

```bash
cd notebooks
pixi run jupytext --sync <name>.py     # propagate .py edits into the paired .ipynb
```

Then `pixi run render` (from repo root) — Quarto executes the changed notebook and refreshes the
`_freeze/` cache. **Commit `_freeze/`**; CI renders from it without re-executing.

One-off data regeneration (rarely needed — see UNSODA note below):

```bash
pixi run python notebooks/fetch_unsoda.py   # needs `mdbtools` on PATH (brew install mdbtools)
```

## Architecture: a 3-notebook pipeline joined by one CSV

Run order is **1 → (2 and 3)**. Notebook 1 is the source of truth; 2 and 3 are independent consumers.

- **`notebooks/1_rosetta_porosity_by_texture`** runs Rosetta over 12 USDA texture classes ×
  bulk densities 0.8–1.9 g/cm³ and **writes `notebooks/rosetta_porosity_by_texture.csv`** — the
  data contract for the whole project (porosity/FC/WP, available & drainable water, Ksat, and the
  Mualem–van Genuchten params `k0_cm_day`, `mualem_L`, `vg_alpha_1cm`, `vg_n`).
- **`notebooks/2_organic_matter_water_holding`** reads that CSV and layers organic-matter effects
  (Rosetta + Minasny & McBratney 2018 blend; Saxton–Rawls 2006).
- **`notebooks/3_rosetta_hydraulic_conductivity`** reads that CSV for saturated/unsaturated K and
  Green–Ampt infiltration.

If you change Notebook 1's outputs, re-run Notebook 1 (to rewrite the CSV) **then** re-run 2 and 3,
then `pixi run render` to refresh `_freeze/`. The CSV is committed (so notebooks 2/3 open standalone),
and the site builds from the committed `_freeze/` cache — see the website-build section.

## Critical, non-obvious constraints

**Interactive sliders need their widget states *baked in* to work without a kernel.** Plots use
`hvplot`/HoloViews `HoloMap`s with `dynamic=False`, which embeds every slider frame as self-contained
Bokeh (`holoviews_exec`). Three things preserve this in static HTML; one silently destroys it:

- ✅ Quarto **executing** the doc, or rendering from its **`_freeze/`** cache — preserves embeds.
- ✅ Panel composition (`pn.Column(md, plot)`) saved with **`embed=True`** (e.g. `.save(..., embed=True)`)
  or marked `.servable()` and rendered by Quarto — preserves embeds. (So you *can* add captions/controls
  around a slider this way — verified.)
- ❌ Feeding Quarto (or pandoc/nbconvert) a **pre-executed `.ipynb` with execution disabled**, or
  displaying a Panel layout in *live* mode (no `embed`) — strips the frames → dead plots. This is the
  trap; both were verified empirically.

**Slider defaults / labels** are set via `hv.Dimension(name, default=<value>, value_format=...)` on
the HoloMap key dimension; the default must be an actual key value (it snaps otherwise). Current
defaults: bulk density **1.4**, organic matter **2%**. The bulk-density dimension label doubles as
on-slider help text (`"Bulk density, g/cm³ (higher is more compacted)"`). That label is rendered as
**HTML** at the slider, so **never put raw `<` or `>` in a dimension label** — use words/arrows.
The label also appears in plot titles via `{dimensions}`.

**`implausible_bd`** marks physically impossible texture×BD combinations (θₛ > 1 − BD/2.65,
i.e. Rosetta extrapolation). It is a CSV column and all plots grey-dash those regions; preserve this
flag/treatment when adding plots.

## Website build (Quarto → GitHub Pages)

- Built with **Quarto** (`_quarto.yml`, `index.qmd`), code folded by default. `pixi run render`
  → `_site/`. Config: `execute: { enabled: true, freeze: auto }` and `execute-dir: file`.
- **`execute-dir: file`** runs each notebook in `notebooks/` so its relative CSV read resolves.
- **Freeze workflow:** `pixi run render` executes notebooks locally and writes `_freeze/`, which
  **is committed**. CI runs `pixi run render`: if the committed cache matches, it renders without
  executing (~15 s). **Caveat: the freeze does not reliably "hit" across macOS→Linux** (and any
  change to a notebook — including a `jupytext_version` header bump — invalidates it), so CI often
  re-executes (~4 min). That is fine: re-execution succeeds (rosetta installs on Linux; the
  committed `notebooks/data_temp/unsoda_bd_om.csv` lets NB2 run fully) and **always deploys**. The
  freeze-staleness step is a **non-fatal warning** — never let it `exit 1` and block the deploy.
- `_site/` and `.quarto/` are git-ignored; **`_freeze/` is committed** (render cache / local speed-up).
  (`linux-64` is in `pixi.toml` platforms so the Ubuntu runner can solve the env.)
- Single-file `quarto render <file>` bypasses the project (no `_site`/freeze) — always run the
  **project** render (`pixi run render`).
- One-time per repo: **Settings → Pages → Source → GitHub Actions**, or the deploy job fails.

## Gotchas

- **UNSODA data.** The derived `notebooks/data_temp/unsoda_bd_om.csv` **is committed** (an
  exception to the otherwise git-ignored `data_temp/`) so CI can fully execute Notebook 2 — incl.
  its UNSODA scatter — on a freeze miss. NB2 also falls back to `OC_BASELINE_PCT = 1.0` if the file
  is ever absent. Only regenerating that CSV via `fetch_unsoda.py` needs Homebrew `mdbtools` (not on
  conda-forge for osx-arm64); reading it does not.
- **Renaming the repo directory breaks the pixi env** (conda console-script shebangs hardcode the
  old absolute path — `jupytext`/`jupyter` fail though `python` works). Fix: `rm -rf .pixi && pixi install`.
- **NB3 page weight**: the K(h) suction grid (36 pts) and Green–Ampt `geomspace(…, 120)` were tuned
  down to keep the embedded page small with no visual change; large grids bloat the static HTML.

## Workflow conventions

- **Commits are staged manually by the user in GitHub Desktop — do not run `git commit`.** Make and
  verify the file changes; leave staging/committing to the user.
- Both the `.py` and `.ipynb` of each notebook are committed and must stay in sync (`jupytext --sync`).
