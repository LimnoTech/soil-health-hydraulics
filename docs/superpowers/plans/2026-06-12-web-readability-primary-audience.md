# Web Readability & Usability for the Primary Audience — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Elevate primary-audience content (interactive charts + plain-language takeaways) across the site, headlined by the two-slider FAO blend figure + a slider-linked storage table on the home page, while keeping researcher/dev detail available but demoted.

**Architecture:** Extract the M&M blend math and a shared FAO band-diagram builder (with category hover) into `notebooks/_helpers.py`; reuse them on `index.qmd` and in NB1/NB2; add a global wheel-zoom-off default; then rewrite notebook prose and restructure the README. Verification is by rendering the Quarto site and inspecting the emitted HTML — there is no unit-test suite.

**Tech Stack:** Python, pandas, hvPlot/HoloViews/Bokeh, Quarto, jupytext, pixi.

**Spec:** `docs/superpowers/specs/2026-06-12-web-readability-primary-audience-design.md`

**Scope:** This plan covers **Phases 1–3** (the committed work). **Phase 4** (per-notebook section reordering) is *optional and contingent* on the Phase 3 review gate; it will get its own plan only if that review calls for it.

---

## Conventions for EVERY task (read once, apply throughout)

These replace the generic "write test / commit" steps in the writing-plans template, adapting to this project:

1. **Edit the `.py` source of truth** (never the `.ipynb` directly) for notebooks. `index.qmd`, `_helpers.py`, `README.md` are edited directly.
2. **Sync** after editing a notebook `.py`:
   `cd notebooks && pixi run jupytext --sync <name>.py` (run from repo root with the `cd` inline).
3. **Render** to execute + refresh the freeze: `pixi run render` (from repo root).
4. **Verify** with the concrete command(s) in the task, AND for any interactive change, **open the HTML** (`_site/...html`) in a browser and confirm sliders move, plots are live, hover shows the right fields, and wheel-zoom is OFF until toggled. A clean render log is necessary but NOT sufficient — embeds can render "successfully" yet be dead.
5. **Do NOT run `git commit`.** Staging/committing is done manually by the user in GitHub Desktop (CLAUDE.md hard rule). Each task ends at a **Checkpoint**: report what changed and what to stage (including updated `notebooks/*.ipynb` and `_freeze/...`), then pause.
6. **Pure-Python functions can be tested directly** without rendering: `pixi run python -c "..."`. Use this for `soil_water_bd_om_blend_table`.
7. Preserve CLAUDE.md constraints: `dynamic=False` embeds, **no raw `<`/`>` in slider dimension labels**, `implausible_bd` grey treatment, keep `.py`/`.ipynb` synced, commit `_freeze/`.

---

# PHASE 1 — Shared code + cross-cutting interactions

### Task 1: Add M&M blend constants + `soil_water_bd_om_blend_table()` to `_helpers.py`

**Files:**
- Modify: `notebooks/_helpers.py` (append after `line_with_extrapolation`)

- [ ] **Step 1: Add the blend constants and function.** Append to `notebooks/_helpers.py`:

```python
# --- Minasny & McBratney (2018) organic-matter blend (shared by NB2 and the home page) ---
VB = 0.58  # van Bemmelen factor: OM ≈ OC / VB
OC_BASELINE_PCT = 1.0  # ROSETTA mineral-baseline organic carbon anchor (%), from UNSODA (NB2 §1)

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


def soil_water_bd_om_blend_table(result, bd, om, oc_baseline=OC_BASELINE_PCT):
    """ROSETTA mineral baseline + Minasny & McBratney (2018) organic-matter increments for a
    single (bulk density `bd`, organic matter `om` %) state — one row per texture class in the
    canonical sand→clay order of `result`. Volumetric water contents (cm³/cm³).

    Applies M&M's WP/AWC/SAT slopes relative to `oc_baseline` and derives FC = WP + AWC so AWC
    matches M&M exactly (mirrors NB2's _blend_profile). `implausible` flags blended SAT exceeding
    the BD-implied pore space (1 − BD/2.65). The single source feeding both the FAO figure and the
    home-page / NB2 linked table.
    """
    base = result.set_index(["texture_class", "bulk_density_g_cm3"])
    texture_classes = list(result["texture_class"].drop_duplicates())
    hsg = dict(
        result[["texture_class", "hydrologic_soil_group"]]
        .drop_duplicates().itertuples(index=False, name=None)
    )
    d_oc = om * VB - oc_baseline  # OM% -> OC%, relative to the mineral-baseline OC
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
```

- [ ] **Step 2: Verify it imports and matches NB2's current math.** The values must equal NB2's `_blend_profile` for a sample state. Run:

```bash
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics/notebooks && pixi run python -c "
import pandas as pd
from _helpers import soil_water_bd_om_blend_table
r = pd.read_csv('rosetta_porosity_by_texture.csv')
t = soil_water_bd_om_blend_table(r, bd=1.4, om=2.0)
row = t[t.texture_class=='loam'].iloc[0]
print('loam @ BD1.4 OM2%:', round(row.wilting_point_porosity,4), round(row.field_capacity_porosity,4), round(row.total_porosity,4), 'impl=', row.implausible)
assert abs((row.field_capacity_porosity - row.wilting_point_porosity) - row.available_water_capacity) < 1e-9
assert len(t) == 12
print('OK')
"
```

Expected: prints four numbers + `impl= False` + `OK` (no assertion error). Sanity: FC ≥ WP, SAT ≥ FC.

- [ ] **Step 3: Checkpoint.** Report the printed values; stage `notebooks/_helpers.py`. Do NOT git commit.

---

### Task 2: Add the shared `soil_water_texture_band_diagram()` builder (bands + lines + labels + grey + hover)

**Files:**
- Modify: `notebooks/_helpers.py` (add a `from bokeh.models import HoverTool` near the imports; append `soil_water_texture_band_diagram`)

- [ ] **Step 1: Add the Bokeh import.** In `notebooks/_helpers.py`, after `import holoviews as hv` add:

```python
from bokeh.models import HoverTool
```

- [ ] **Step 2: Append the builder:**

```python
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
    return overlay
```

- [ ] **Step 3: Smoke-test the import + return type:**

```bash
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics/notebooks && pixi run python -c "
import numpy as np, holoviews as hv
from _helpers import soil_water_texture_band_diagram
x=np.arange(12); pwp=np.linspace(1,3,12); fc=pwp+2; por=fc+2
ov = soil_water_texture_band_diagram(x,pwp,fc,por, implausible=(x>9))
print(type(ov), 'OK')
"
```

Expected: prints an `Overlay` type + `OK`, no exception.

- [ ] **Step 4: Checkpoint.** Note: hover hit-testing on `fill_alpha=0` rectangles is verified later when wired into a rendered plot (Task 4). If hover does not trigger there, change `fill_alpha=0` to `fill_alpha=0.01` in `soil_water_texture_band_diagram`. Stage `notebooks/_helpers.py`; do NOT git commit.

---

### Task 3: Add the global wheel-zoom-off default to `_helpers.py`

**Files:**
- Modify: `notebooks/_helpers.py` (after the `pd.set_option(...)` line near the top)

- [ ] **Step 1: Add the default.** After the existing `pd.set_option("display.float_format", ...)` line, add:

```python
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
```

- [ ] **Step 2: Verify import is clean:**

```bash
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics/notebooks && pixi run python -c "import _helpers; print('import OK')"
```

Expected: `import OK` with no traceback. (A HoloViews "extension already loaded" notice is harmless.)

- [ ] **Step 3: Checkpoint.** The *effect* (scroll doesn't zoom until toggled) is confirmed visually in Task 4's render. Stage `notebooks/_helpers.py`; do NOT git commit.

---

### Task 4: Re-point NB1 §8 to `soil_water_texture_band_diagram` + new y-axis label; verify hover & wheel-zoom

**Files:**
- Modify: `notebooks/1_rosetta_porosity_by_texture.py` (`_water_profile`, the `profiles.opts` ylabel)

- [ ] **Step 1: Add the import.** In NB1's top cell, change the helper import to include `soil_water_texture_band_diagram`:

```python
from _helpers import show, vg_theta, mualem_k, line_with_extrapolation, soil_water_texture_band_diagram
```

- [ ] **Step 2: Replace the body of `_water_profile`** so it delegates to the shared builder. Replace the existing function (from `def _water_profile(bd):` through its `return overlay * lines * labels`) with:

```python
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
```

- [ ] **Step 3: Update the y-axis label.** In the `profiles.opts(...)` block, change:

```python
        ylabel="Water volume (inches per foot of soil depth)",
```
to
```python
        ylabel="Water Storage Capacity (inches per foot of soil depth)",
```

- [ ] **Step 4: Sync + render:**

```bash
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics/notebooks && pixi run jupytext --sync 1_rosetta_porosity_by_texture.py
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics && pixi run render 2>&1 | tail -15
```

Expected: render completes, "Output created: _site/index.html", NB1 cells all `Done`.

- [ ] **Step 5: Verify content + interactivity.**

```bash
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics
grep -c "Water Storage Capacity (inches per foot" _site/notebooks/1_rosetta_porosity_by_texture.html
grep -c "Total stormwater capacity" _site/notebooks/1_rosetta_porosity_by_texture.html
```

Expected: first grep ≥ 1 (new label present); second grep ≥ 1 (hover tooltip text embedded). **Then open `_site/notebooks/1_rosetta_porosity_by_texture.html`** and confirm: §8 diagram renders identically to before; the BD slider works; **hovering a texture column shows Available / Drainable / Total stormwater capacity**; **scrolling does NOT zoom** until you click the wheel-zoom toolbar button. If hover doesn't fire, apply the `fill_alpha=0.01` fix from Task 2 and re-render.

- [ ] **Step 6: Checkpoint.** Stage `notebooks/1_rosetta_porosity_by_texture.py` + `.ipynb` + changed `_freeze/...`. Do NOT git commit.

---

### Task 5: Re-point NB2 §2 and §3 FAO diagrams to `soil_water_texture_band_diagram` + blend table; new y-labels

**Files:**
- Modify: `notebooks/2_organic_matter_water_holding.py` (`_blend_profile`, `_sr_profile`, the two `*.opts` ylabels, and the blend constants now sourced from `_helpers`)

- [ ] **Step 1: Import shared blend + builder.** In NB2's top cell, change the helper import to:

```python
from _helpers import show, soil_water_texture_band_diagram, soil_water_bd_om_blend_table, VB, OC_BASELINE_PCT, MM_SLOPES, MM_GROUP
```

- [ ] **Step 2: Remove the now-duplicated constants.** Delete NB2's local definitions of `VB`, `MM_SLOPES`, and `MM_GROUP` (they now come from `_helpers`). **Keep** `OC_BASELINE_PCT`'s narrative in §1, but delete the bare re-assignment line `OC_BASELINE_PCT = 1.0` only if §1 no longer needs to override it — otherwise leave §1 as-is (it may set its own). Keep `base = result.set_index([...])` (still used by `blend_line` and `_sr_profile` paths) and `om_grid`.

- [ ] **Step 3: Replace `_blend_profile` body** to delegate to the shared builder via the shared blend table. Replace the whole function with:

```python
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
```

(Note: this requires `texture_x` and `INCHES_PER_FOOT` to exist in NB2. NB2 already defines `texture_x` and `INCHES_PER_FOOT` in its setup cell — confirm with `grep -n "texture_x\|INCHES_PER_FOOT" notebooks/2_organic_matter_water_holding.py`; if absent, add `INCHES_PER_FOOT = 12.0` and `texture_x = {cls: i for i, cls in enumerate(TEXTURE_CLASSES)}` to the setup cell.)

- [ ] **Step 4: Replace `_sr_profile` body** to delegate to the shared builder (no greying — Saxton–Rawls has no implausible flag here):

```python
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
    return soil_water_texture_band_diagram(x, pwp, fc, por)
```

- [ ] **Step 5: Update both y-axis labels.** In `blend_profiles.opts(...)` and `sr_profiles.opts(...)`, change `ylabel="Water volume (inches per foot of soil depth)"` to `ylabel="Water Storage Capacity (inches per foot of soil depth)"`.

- [ ] **Step 6: Sync + render:**

```bash
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics/notebooks && pixi run jupytext --sync 2_organic_matter_water_holding.py
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics && pixi run render 2>&1 | tail -15
```

Expected: NB2 executes all cells `Done`; site builds.

- [ ] **Step 7: Verify.**

```bash
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics
grep -c "Water Storage Capacity (inches per foot" _site/notebooks/2_organic_matter_water_holding.html
grep -c "Total stormwater capacity" _site/notebooks/2_organic_matter_water_holding.html
```

Expected: both ≥ 1. **Open the HTML**: the §2 two-slider diagram and §3 OM-slider diagram render as before; hover shows the three categories; grey extrapolation columns still appear in §2; sliders work; wheel-zoom off by default. Confirm the §2 diagram is **visually unchanged** from the pre-refactor version (the blend table feeds identical numbers).

- [ ] **Step 8: Checkpoint.** Stage NB2 `.py` + `.ipynb` + `_freeze/...`. Do NOT git commit.

---

### Task 6: Phase 1 gate — full render + cross-page interactive verification

- [ ] **Step 1:** `cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics && pixi run render 2>&1 | tail -20` → clean build, all notebooks `Done`.
- [ ] **Step 2:** Open `_site/notebooks/1_...html` and `_site/notebooks/2_...html`. Confirm, on every FAO diagram: (a) category hover works, (b) y-axis reads "Water Storage Capacity …", (c) wheel-zoom is off until toggled, (d) `dynamic=False` sliders still move and the plot is live without a kernel.
- [ ] **Step 3: Checkpoint.** Phase 1 is the foundation for Phases 2–3. Report results; the user stages the Phase 1 changes. STOP for review before Phase 2.

---

# PHASE 2 — Home-page summary section (+ NB2 linked table)

> **AS-BUILT (2026-06-13):** Tasks 7–8 below describe the original `HoloMap` + `hv.Table` approach.
> It was **superseded during implementation by a Panel `embed=True` layout** — `pn.Column(figure,
> pn.Row(bd_slider, om_slider), pn.pane.HTML(soil_water_table_html(...)))` with sliders **between**
> the figure and table — because `hv.Table` can't format the table as required and `hv.Div` doesn't
> embed inline in Quarto. Same on `index.qmd` and NB2 §2. See the spec's "As-built notes" for the
> full rationale and parameters (figure width 720, toolbar right, fonts, page weight).

### Task 7: Build the home-page summary (headline figure + slider-linked table)

**Files:**
- Modify: `index.qmd` (insert a new section + an executed Python cell above "## The notebooks")

- [ ] **Step 1: Insert the summary section.** In `index.qmd`, immediately **above** the `## The notebooks` line, insert:

````markdown
## Soil health at a glance: how compaction & organic matter change water storage

Drag the two sliders — **bulk density** (compaction) and **organic matter** — and watch the
soil-water bands respond. Easing compaction and adding organic matter shift water out of the
unavailable band into **plant-available** and **drainable** storage. The table below the chart
lists the same values per soil texture class (volume fraction, cm³/cm³).

```{python}
#| echo: false
import sys
sys.path.insert(0, "notebooks")  # index.qmd executes at the repo root (execute-dir: file)

import numpy as np
import pandas as pd
import hvplot.pandas  # noqa: F401
import holoviews as hv
from _helpers import soil_water_texture_band_diagram, soil_water_bd_om_blend_table, show  # also sets float_format + wheel-zoom default

hv.output(widget_location="bottom")

result = pd.read_csv("notebooks/rosetta_porosity_by_texture.csv")
TEXTURE_CLASSES = list(result["texture_class"].drop_duplicates())
texture_x = {cls: i for i, cls in enumerate(TEXTURE_CLASSES)}
HSG = dict(result[["texture_class", "hydrologic_soil_group"]].drop_duplicates().itertuples(index=False, name=None))
texture_ticks = [(i, f"{cls} ({HSG[cls]})") for cls, i in texture_x.items()]
INCHES_PER_FOOT = 12.0

bulk_densities = np.sort(result["bulk_density_g_cm3"].unique())
om_grid = np.round(np.arange(0.0, 8.0 + 1e-9, 1.0), 1)

BD_DIM = hv.Dimension("Bulk density, g/cm³ (higher is more compacted)", default=1.4, value_format=lambda v: f"{v:.1f}")
OM_DIM = hv.Dimension("Organic matter (% by weight)", default=2.0, value_format=lambda v: f"{v:.1f}")

def _fig(bd, om):
    t = soil_water_bd_om_blend_table(result, bd, om)
    x = np.array([texture_x[c] for c in TEXTURE_CLASSES])
    pwp = t["wilting_point_porosity"].to_numpy() * INCHES_PER_FOOT
    fc = t["field_capacity_porosity"].to_numpy() * INCHES_PER_FOOT
    por = t["total_porosity"].to_numpy() * INCHES_PER_FOOT
    return soil_water_texture_band_diagram(
        x, pwp, fc, por,
        texture_labels=t["texture_class"].astype(str).tolist(),
        implausible=t["implausible"].to_numpy(),
    )

def _tbl(bd, om):
    t = soil_water_bd_om_blend_table(result, bd, om)
    out = pd.DataFrame({
        "texture (HSG)": [f"{c} ({HSG[c]})" for c in t["texture_class"]],
        "wilting point": t["wilting_point_porosity"].to_numpy(),
        "field capacity": t["field_capacity_porosity"].to_numpy(),
        "saturation": t["total_porosity"].to_numpy(),
        "available water": t["available_water_capacity"].to_numpy(),
        "drainable water": t["drainable_water"].to_numpy(),
    }).round(3)
    return hv.Table(out)

fig_hmap = hv.HoloMap({(bd, om): _fig(bd, om) for bd in bulk_densities for om in om_grid}, kdims=[BD_DIM, OM_DIM])
tbl_hmap = hv.HoloMap({(bd, om): _tbl(bd, om) for bd in bulk_densities for om in om_grid}, kdims=[BD_DIM, OM_DIM])

fig_hmap = fig_hmap.opts(
    hv.opts.Overlay(width=820, height=520, legend_position="top_left",
                    xlabel="texture class (hydrologic soil group); coarse → fine",
                    ylabel="Water Storage Capacity (inches per foot of soil depth)",
                    title="Soil water vs. texture — ROSETTA + Minasny & McBratney blend\n{dimensions}"),
    hv.opts.Curve(xticks=texture_ticks, xrotation=45, ylim=(0, 10)),
    hv.opts.Area(xticks=texture_ticks, xrotation=45, ylim=(0, 10)),
)

(fig_hmap + tbl_hmap).cols(1)
```

::: {.callout-note appearance="simple"}
Figure axis is in **inches per foot of soil depth**; the table is the same quantities as
**volume fraction (cm³/cm³)**.
:::
````

- [ ] **Step 2: Render:** `cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics && pixi run render 2>&1 | tail -20`. Expected: `index.qmd` now shows a kernel start + cell executions; "Output created: _site/index.html".

- [ ] **Step 3: Verify shared widget + content:**

```bash
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics
grep -c "Soil health at a glance" _site/index.html
grep -c "Total stormwater capacity" _site/index.html
ls -lh _site/index.html
```

Expected: heading present (≥1); hover text present (≥1); note the file size. **Open `_site/index.html`** and confirm: ONE pair of sliders (bulk density + organic matter) controls **both** the figure and the table; both update together; figure hover works; wheel-zoom off. If the figure and table render with **separate** slider sets, switch the layout to a single `hv.HoloMap` of a `Layout` value instead of `fig_hmap + tbl_hmap` (see fallback note below).

- [ ] **Step 4: Page-weight check.** If `_site/index.html` is heavy (rule of thumb: > ~8–10 MB), reduce frames on the **home page only**: coarsen `om_grid` to `np.arange(0,8.1,2.0)` and/or restrict `bulk_densities` to a representative subset (e.g. `[1.0,1.2,1.4,1.6,1.8]`). Re-render and re-check. Record the chosen grid.

- [ ] **Step 5: Checkpoint.** Stage `index.qmd` + `_freeze/index/...` + `_site` is git-ignored (not staged). Do NOT git commit. **Fallback note (shared widget):** if `+` doesn't share one widget, replace the last line with a HoloMap of Layouts:
  `hv.HoloMap({(bd,om): (_fig(bd,om).opts(...) + _tbl(bd,om)).cols(1) for ...}, kdims=[BD_DIM, OM_DIM])` and apply figure opts inside `_fig`.

---

### Task 8: Add the slider-linked table to NB2 §2

**Files:**
- Modify: `notebooks/2_organic_matter_water_holding.py` (the §2 `blend_profiles` cell)

- [ ] **Step 1: Add a table HoloMap alongside `blend_profiles`.** In the §2 cell, after `blend_profiles` is built and `.opts(...)` applied, replace the final expression so the figure and a linked table render together. Add:

```python
def _blend_tbl(bd, om):
    t = soil_water_bd_om_blend_table(result, bd, om)
    out = pd.DataFrame({
        "texture (HSG)": [f"{c} ({HYDROLOGIC_SOIL_GROUP[c]})" for c in t["texture_class"]],
        "wilting point": t["wilting_point_porosity"].to_numpy(),
        "field capacity": t["field_capacity_porosity"].to_numpy(),
        "saturation": t["total_porosity"].to_numpy(),
        "available water": t["available_water_capacity"].to_numpy(),
        "drainable water": t["drainable_water"].to_numpy(),
    }).round(3)
    return hv.Table(out)

_blend_tbl_hmap = hv.HoloMap(
    {(bd, om): _blend_tbl(bd, om) for bd in bulk_densities for om in om_grid},
    kdims=[
        hv.Dimension("Bulk density, g/cm³ (higher is more compacted)", default=1.4, value_format=lambda v: f"{v:.1f}"),
        hv.Dimension("Organic matter (% by weight)", default=2.0, value_format=lambda v: f"{v:.1f}"),
    ],
)

(blend_profiles + _blend_tbl_hmap).cols(1)
```

(The two HoloMaps share kdims with `blend_profiles`, so they render under one shared slider pair. `bulk_densities`, `om_grid`, `HYDROLOGIC_SOIL_GROUP` already exist in NB2's namespace.)

- [ ] **Step 2: Sync + render:**

```bash
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics/notebooks && pixi run jupytext --sync 2_organic_matter_water_holding.py
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics && pixi run render 2>&1 | tail -15
```

- [ ] **Step 3: Verify.** Open `_site/notebooks/2_organic_matter_water_holding.html`; confirm §2 now shows the figure **and** a table beneath it driven by the same BD+OM sliders. Same shared-widget fallback as Task 7 applies if needed.

- [ ] **Step 4: Checkpoint.** Stage NB2 `.py` + `.ipynb` + `_freeze/...`. Do NOT git commit.

---

### Task 9: Phase 2 gate

- [ ] **Step 1:** Full `pixi run render`; open `_site/index.html` and NB2 HTML; confirm both linked figure+table work with single shared sliders, hover, wheel-zoom-off, and acceptable page weight.
- [ ] **Step 2: Checkpoint.** Report; user stages. STOP for review before Phase 3.

---

# PHASE 3 — Plain language + in-place collapse + README (NO reordering)

> Phase 3 does **not** reorder notebook sections. It rewrites prose for the primary audience, adds a one-line takeaway under each chart, and collapses heavy method **in place** using Quarto collapsible callouts. Code stays folded (already configured).

**Collapsible callout syntax** (use throughout Phase 3 to demote method in place):

```markdown
::: {.callout-note collapse="true"}
## For researchers: <topic>
<the existing technical paragraph(s), moved verbatim into here>
:::
```

**Takeaway styling** (one line under a chart):

```markdown
::: {.callout-tip appearance="simple"}
**Takeaway:** <one plain-language sentence about what the chart shows for soil health / stormwater>.
:::
```

### Task 10: Notebook 1 — plain language + takeaways + in-place collapse

**Files:**
- Modify: `notebooks/1_rosetta_porosity_by_texture.py` (markdown cells only; no code reordering)

- [ ] **Step 1: Add chart takeaways.** Under each of §6, §7, §8's chart cells, add a `# %% [markdown]` cell with a `.callout-tip` takeaway. Example for §8 (write equivalents for §6 and §7):

```python
# %% [markdown]
# ::: {.callout-tip appearance="simple"}
# **Takeaway:** Compaction (higher bulk density) squeezes the blue *drainable* band first, then
# the green *plant-available* band — so a compacted soil both stores less plant water and drains
# storm flows more slowly.
# :::
```

- [ ] **Step 2: Collapse heavy method in place.** Wrap the dense passages in collapsible callouts **without moving them**: the van Genuchten equation block in the §3 markdown, and the §4 Rosetta-mechanics paragraph. Convert each `# ## N. …` method intro's technical body into:

```python
# %% [markdown]
# ::: {.callout-note collapse="true"}
# ## For researchers: van Genuchten–Mualem retention model
# <move the existing θ(h) equation + parameter explanation here, unchanged>
# :::
```

- [ ] **Step 3: Plain-language pass.** Lightly rewrite the §1, §2 intro sentences and the top-of-notebook header prose to lead with what the reader gets (charts of how compaction changes water storage), deferring "Rosetta pedotransfer functions predict van Genuchten parameters…" to a collapsed "How this works" note. Keep all numbers/claims identical.

- [ ] **Step 4: Sync + render + verify:**

```bash
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics/notebooks && pixi run jupytext --sync 1_rosetta_porosity_by_texture.py
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics && pixi run render 2>&1 | tail -10
grep -c "Takeaway" _site/notebooks/1_rosetta_porosity_by_texture.html
grep -c 'callout-note' _site/notebooks/1_rosetta_porosity_by_texture.html
```

Expected: takeaways present (≥3), collapsible callouts present. Open the HTML: confirm callouts collapse/expand, charts unchanged, no broken markdown.

- [ ] **Step 5: Checkpoint.** Stage NB1 `.py` + `.ipynb` + `_freeze/...`. Do NOT git commit.

---

### Task 11: Notebook 2 — plain language + takeaways + in-place collapse

**Files:**
- Modify: `notebooks/2_organic_matter_water_holding.py` (markdown cells only)

- [ ] **Step 1: Add takeaways** under the §2 blend figure/table, the §2 AWC and drainable line plots. Example for the headline §2:

```python
# %% [markdown]
# ::: {.callout-tip appearance="simple"}
# **Takeaway:** Adding organic matter lifts both the plant-available (green) and drainable (blue)
# bands — most strongly in coarse soils — so soil-health practices buy measurable storm and drought
# resilience. Slide organic matter up and watch the bands grow.
# :::
```

- [ ] **Step 2: Collapse the researcher detail in place.** Wrap §1 (UNSODA mineral-baseline OC + the OLS regressions narrative) and §3 + §3.1 (Saxton–Rawls and its M&M validation) in `::: {.callout-note collapse="true"}` blocks titled "For researchers: …", leaving their code cells where they are (code stays folded). Do not move sections.

- [ ] **Step 3: Plain-language pass** on the notebook header + §2 intro: lead with the practical result; defer "Minasny & McBratney increments / van Bemmelen / mineral baseline" mechanics into the collapsed §1 note.

- [ ] **Step 4: Sync + render + verify** (same commands as Task 10, with NB2 paths). Confirm takeaways + collapsible callouts present; figure/table from Phase 2 still work.

- [ ] **Step 5: Checkpoint.** Stage NB2 `.py` + `.ipynb` + `_freeze/...`. Do NOT git commit.

---

### Task 12: Notebook 3 — plain language + takeaways + in-place collapse

**Files:**
- Modify: `notebooks/3_rosetta_hydraulic_conductivity.py` (markdown cells only)

- [ ] **Step 1: Add takeaways** under §1 (Ksat), §2 (K(h)), §3 (Green–Ampt). Example for §3:

```python
# %% [markdown]
# ::: {.callout-tip appearance="simple"}
# **Takeaway:** Infiltration starts fast and settles toward the saturated rate (Ksat); compaction
# (higher bulk density) lowers both, so compacted soils pond and run off sooner in a storm.
# :::
```

- [ ] **Step 2: Collapse the dense method paragraphs in place** — the Mualem–van Genuchten K(h) explanation in §2 and the Rawls/Brakensiek/Miller ψ_f rationale in §3 — into `::: {.callout-note collapse="true"}` "For researchers" blocks.

- [ ] **Step 3: Sync + render + verify** (NB3 paths). Confirm takeaways + callouts; sliders/plots still live; wheel-zoom-off inherited from Phase 1.

- [ ] **Step 4: Checkpoint.** Stage NB3 `.py` + `.ipynb` + `_freeze/...`. Do NOT git commit.

---

### Task 13: README — Developers & Contributors section

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Move the dev sections to the bottom.** Cut the `### Environment (pixi)` subsection and the `## Interactive website (GitHub Pages)` section from their current mid-document positions and paste them, unchanged, **after** the `## References` section, under a new top-level heading:

```markdown
## Developers & Contributors

### Environment (pixi)
<existing Environment content>

### Interactive website (GitHub Pages)
<existing Interactive website content>
```

Resulting order: intro → Notebooks table → Notebook 1/2/3 detail → References → **Developers & Contributors**. Leave the per-notebook "Key methodologies" blocks where they are (out of scope).

- [ ] **Step 2: Verify structure:**

```bash
cd /Users/aaufdenkampe/Documents/git_Limno/soil-health-hydraulics
grep -nE "^#{1,3} (Environment|Interactive website|References|Developers & Contributors)" README.md
```

Expected: `References` appears **before** `Developers & Contributors`, and `Environment` + `Interactive website` appear **after** it. No duplicate headings.

- [ ] **Step 3: Checkpoint.** Stage `README.md`. Do NOT git commit.

---

### Task 14: Phase 3 review gate (decides whether Phase 4 happens)

- [ ] **Step 1:** Full `pixi run render`; open `_site/index.html` and all three notebook pages.
- [ ] **Step 2: Review with the user against the goal:** does the home-page summary + clearer prose + collapsed method read well for the primary audience *without* reordering notebook sections? 
- [ ] **Step 3: Decision.** If yes → the project is done (Phases 1–3 stand on their own). If no → spin up a **separate Phase 4 plan** for per-notebook section reordering (charts lead, method pooled into a bottom section for NB1 and NB2). Record the decision.

---

## Self-review (completed by plan author)

- **Spec coverage:** headline figure + linked table on home page (Task 7) and NB2 (Task 8) ✓; shared FAO builder + blend (Tasks 1–2) ✓; category hover with Available/Drainable/Total stormwater capacity, no unavailable (Task 2) ✓; y-axis "Water Storage Capacity" on all FAO plots (Tasks 4,5,7) ✓; wheel-zoom off globally (Task 3, verified 4/6) ✓; plain language + in-place collapse, no reorder (Tasks 10–12) ✓; README Developers & Contributors (Task 13) ✓; Phase 4 deferred behind a review gate (Task 14) ✓; `index.qmd` execution path handled (Task 7 sys.path + `notebooks/` CSV) ✓; page-weight tuning (Task 7 Step 4) ✓; embed/CLAUDE.md constraints in Conventions ✓.
- **Placeholders:** none — `<topic>`/`<existing … content>` markers in Phase 3 denote *verbatim existing prose to relocate*, not unwritten code; all code steps contain complete code.
- **Type/name consistency:** `soil_water_bd_om_blend_table` columns (`wilting_point_porosity`, `field_capacity_porosity`, `total_porosity`, `available_water_capacity`, `drainable_water`, `implausible`) are used identically in Tasks 1, 5, 7, 8; `soil_water_texture_band_diagram(x, pwp, fc, por, *, implausible=, hover=)` signature is called consistently in Tasks 4, 5, 7.

## Key risks (carried from the spec)

- **Embeds can render "successfully" yet be dead** — always open the HTML, never trust the render log alone.
- **Shared-widget layout** (`fig_hmap + tbl_hmap`) — verify one slider pair drives both; fallback = HoloMap-of-Layout (Tasks 7/8 notes).
- **Hover hit-testing on transparent rectangles** — fallback `fill_alpha=0.01` (Task 2/4).
- **Wheel-zoom global default** may be overridden by per-figure opts — verify visually (Task 4).
- **Home-page weight** from 108 (BD×OM) frames carrying plot + table — measure and coarsen the grid if needed (Task 7 Step 4).
