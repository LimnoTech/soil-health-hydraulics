# Design: Improve web readability & usability for the primary audience

**Date:** 2026-06-12
**Branch:** `improve_docs` (do not merge to `main` until the broader doc work is complete)
**Status:** Approved design — ready for implementation plan

## Problem & goal

The published site (https://limnotech.github.io/soil-health-hydraulics/) and the repo README
currently present content in **computational-pipeline order**, not **reader-interest order**. The
primary audience hits methodology, math, and code before reaching the interactive charts and the
"what this means for soil health / stormwater" payoff.

**Primary audience:** people with a reasonable soil-science and stormwater background and an
interest in soil health — *not* graduate soil physicists and *not* Python programmers.
**Secondary audiences (keep, but demote below the primary content):** (1) soil-science
researchers, (2) contributors to the Python repo.

**Goal:** elevate primary-audience content (interactive charts + plain-language takeaways) and
demote secondary content (derivations, validation, environment/build mechanics) without deleting
it.

## Guiding pattern (applied to every page)

1. Lead with the **interactive chart** + a **one-line plain-language takeaway**.
2. Keep **short caveats** as inline **collapsed callouts** near their chart.
3. Pool **heavy methodology / derivations / validation** into a bottom
   **"Methods, assumptions & references"** section (collapsible).
4. **Code stays folded** (Quarto `code-fold: true`, already configured).

**For the notebooks, this pattern is split across two phases.** Items 2 and 4 (in-place collapse,
folded code) plus a plain-language prose rewrite happen in **Phase 3** *without reordering sections*.
Items 1 and 3 — actually **reordering** so charts lead and method pools at the bottom — are deferred
to an **optional Phase 4**, to be decided only after reviewing whether Phase 3 alone is sufficient.
The home page (Phase 2) realizes the full pattern from the start.

## Decisions (locked during brainstorming)

- Barrier to fix: **both** ordering (charts buried) **and** prose density/jargon.
- Demotion mechanism: **hybrid** — charts + takeaways lead; short caveats inline-collapsed; heavy
  method pooled at the bottom per page.
- Headline figure: **NB2 §2 two-slider FAO blend plot** (bulk density + organic matter), shown on
  the **home page** in a new summary section above "The Notebooks", **and kept in NB2 §2**.
- A **slider-linked table** sits below the headline figure (home page) and is **also added to
  NB2 §2**: rows = texture class (HSG) in canonical sand→clay order; columns = wilting point,
  field capacity, saturated/total, available water, drainable water; all **volume fraction
  (cm³/cm³)**; driven by the **same** BD + OM sliders as the figure.
- Code sharing: **extract a shared FAO-diagram builder + the M&M blend math** into shared code
  (revises the earlier "shared functions only" call now that the diagram is reused and gains hover).
  Keep it in **`notebooks/_helpers.py`** (one shared module) unless implementation shows it should
  split.
- **Hover by category** added to all "Soil water vs. texture" FAO-style plots (additive overlay;
  visuals unchanged). Hover reports **available**, **drainable**, and **total stormwater capacity**
  (= available + drainable); **unavailable water is omitted** (not relevant to the audience).
- **Y-axis label** on all FAO-style plots changes from "Water volume …" to **"Water Storage
  Capacity (inches per foot of soil depth)"**.
- **Wheel zoom off by default** on **every** figure, toggleable from the toolbar (global default).
- **README.md:** move the **Environment (pixi)** and **Interactive website (GitHub Pages)**
  subsections to the **very bottom**, below **References**, into a new **`## Developers &
  Contributors`** section.
- Summary-section heading: **"Soil health at a glance: how compaction & organic matter change
  water storage."**
- Visual-companion tooling: declined (token budget) — text-only brainstorming.

## Architecture

### Shared code (extend `notebooks/_helpers.py`)

New shared API (names indicative; finalize in the plan):

- `mm_blend_table(result, bd, om, oc_baseline=OC_BASELINE_PCT) -> DataFrame`
  Per-texture wilting point / field capacity / saturated / available / drainable (volume
  fraction), in canonical sand→clay order, plus the `implausible` flag (blended SAT > 1 − BD/2.65).
  **Single source** feeding both the FAO figure and the linked table so they cannot drift.
  Subsumes the current NB2 `_blend_profile` math and the `MM_SLOPES` / `MM_GROUP` / `VB` /
  `OC_BASELINE_PCT` constants.
- `fao_diagram(profile, *, units="in_per_ft", implausible_mask=None, hover=True, labels=True) -> hv.Overlay`
  The `bands * lines * labels (* grey extrapolation spans) (* hover overlay)` diagram. Consumed by
  NB1 §8 (BD only), NB2 §2 (BD + OM, the headline), NB2 §3 (Saxton–Rawls, OM only), and the home
  page. Inputs are the per-texture pwp/fc/por arrays so each caller supplies its own model values.
- **Global plot default:** wheel zoom present-but-inactive on load, set once on import (e.g. via a
  small `hv.opts.defaults(...)` / per-builder `active_tools=[]`). Verify it does not clobber other
  per-figure opts.

NB2-only code that **stays** in NB2: `saxton_rawls`, `blend_line`, the §1 UNSODA/regression logic.
NB3-only `bd_slider_line` stays in NB3.

### Hover by category (FAO plots)

Direct hover on `hv.Area` only exposes a single band, so add **one invisible hover marker per
texture column** (transparent element at each texture x-position) carrying the category volumes as
hover fields, with a custom Bokeh `HoverTool` tooltip showing, for the hovered texture:
**Available water**, **Drainable water**, and **"Total stormwater capacity"**
(= available + drainable), all in **inches/foot** (matching the axis). **Unavailable water is
omitted** — it is
not relevant to the stormwater/soil-health audience. The visible bands/lines/labels are untouched
(the unavailable band still renders; it just isn't in the tooltip). Implemented inside
`fao_diagram` so all four call sites get it for free. Tooltip formatting to be verified against the
rendered HTML.

The FAO-plot **y-axis label** also changes from "Water volume (inches per foot of soil depth)" to
**"Water Storage Capacity (inches per foot of soil depth)"** across all four call sites.

## Page-by-page changes

### Home page (`index.qmd`) — Phase 2

- New section **"Soil health at a glance: how compaction & organic matter change water storage"**
  inserted **above** "## The notebooks".
- Becomes an **executed** document: a code cell reads `notebooks/rosetta_porosity_by_texture.csv`
  and imports the shared module. Because `execute-dir: file` runs `index.qmd` at the **repo root**,
  the plan must handle the path (add `notebooks/` to `sys.path` and read `notebooks/…csv`).
- Contents: heading → 1–2 sentence plain framing (e.g. "Drag the sliders: adding organic matter and
  easing compaction shift water from unavailable toward plant-available and drainable storage.") →
  **headline two-slider FAO figure** → **slider-linked table** beneath it.
- **Figure + table share one BD+OM slider pair**: a single `HoloMap` keyed by (bulk density,
  organic matter) whose value is a `Layout(fao_diagram(...), hv.Table(...))` with `dynamic=False`,
  so every frame embeds both and one widget pair drives both (preserves the CLAUDE.md embed rule).
- **Units:** figure y-axis is inches/foot; table is volume fraction (cm³/cm³) per the spec — add a
  one-line caption noting the difference rather than duplicating columns.
- The existing "How to explore" callout already documents sliders/hover/zoom and covers this figure.

> The FAO-diagram swap to the shared `fao_diagram` (hover + new y-axis label) happens in **Phase 1**
> for NB1 §8 and NB2 §2/§3; it is not a reorder. The items below are the prose/structure work.

### Notebook 1

- **Phase 3 (no reordering):** plain-language rewrite of the §1–§8 prose; add a one-line takeaway
  under each chart (§6–§8); collapse heavy method **in place** — VG equation and Rosetta-run detail
  into collapsed callouts; code stays folded.
- **Phase 4 (optional, pending review):** move the visualizations (§6 line plots, §7 stacked-bar
  partition, §8 FAO band diagram) up to **lead** the page, and relocate §1 (texture values), §2 (BD
  range / suction set-points), §3 (VG–Mualem helpers), §4 (Rosetta run), §5 (results table) into a
  bottom **"How these numbers are computed (methods)"** section.

### Notebook 2

- **Phase 2:** add the **slider-linked table** to §2 (same shared builder/mechanism as the home page).
- **Phase 3 (no reordering):** plain-language rewrite + a plain takeaway on the §2 blend figure/line
  plots; collapse **§1 (UNSODA mineral-baseline OC + OLS regressions)** and **§3 + §3.1
  (Saxton–Rawls and its validation vs. M&M)** **in place** into collapsed callouts.
- **Phase 4 (optional, pending review):** reorder so **§2 blend leads**, and relocate §1 and §3/§3.1
  into a bottom **"Calibration & alternative method (for researchers)"** section.

### Notebook 3

- **Phase 3 only:** already chart-led (§1 Ksat, §2 K(h), §3 Green–Ampt), so **light touch** — add
  plain takeaways and fold dense method paragraphs into collapsed callouts. The wheel-zoom default
  is inherited from Phase 1. No FAO plots → **no hover work**. No reordering expected (no Phase 4).

### README.md

- Move **Environment (pixi)** and **Interactive website (GitHub Pages)** subsections to the bottom,
  below **References**, under a new **`## Developers & Contributors`** section.
- Resulting order: intro → Notebooks table → Notebook 1/2/3 detail → References → **Developers &
  Contributors** (Environment + Interactive website).
- **Out of scope (unless requested):** the per-notebook "Key methodologies" blocks stay where they
  are.

## Phasing

1. **Phase 1 — Shared code + cross-cutting interactions.** Extract `mm_blend_table` + `fao_diagram`
   into `_helpers.py`; add category hover; set the global wheel-zoom-off default. Re-point NB1 §8,
   NB2 §2/§3 at the shared builder. **Gate:** `pixi run render`, open the HTML, confirm sliders,
   embeds, hover, and wheel-zoom-default all work before proceeding.
2. **Phase 2 — Home-page summary section (+ NB2 linked table).** Headline figure + slider-linked
   table + framing + heading in `index.qmd`; wire `index.qmd` to the shared module/CSV; add the
   same linked table to NB2 §2. **Gate:** render + open; confirm the shared widget drives both
   figure and table; measure page weight.
3. **Phase 3 — Plain language + in-place collapse + README (no notebook reordering).** Rewrite
   notebook prose for the primary audience, add chart takeaways, collapse heavy method **in place**
   (folded code + collapsed callouts), and restructure the README (Developers & Contributors
   section). **Review gate:** with the home-page summary live and the prose clearer, evaluate
   whether the notebooks read well enough **before** deciding to do Phase 4.
4. **Phase 4 — Per-notebook section reordering (OPTIONAL, pending the Phase 3 review).** Reorder
   NB1 and NB2 so charts lead and method pools into a bottom section. Only undertaken if the Phase 3
   review concludes the lighter touch is insufficient. NB3 is already chart-led and likely needs no
   reordering.

## Risks & verification

- **Embed preservation is the top risk.** The linked table (HoloMap-of-Layout), the hover overlay,
  and the shared builder must all survive as static HTML. Verify by rendering and **opening the
  HTML** at each phase gate — not just by a clean render log.
- **Home-page weight.** The headline embeds every (BD × OM) frame — currently 12 × 9 = 108 frames,
  now each carrying a plot **and** a table. Measure the emitted HTML size; if heavy, trim the OM
  grid or BD range **on the home-page copy only** (NB2 may keep the finer grid). Mirrors the
  existing NB3 page-weight tuning.
- **`index.qmd` execution path.** `execute-dir: file` runs it at the repo root; the import + CSV
  read must be made robust to that.
- **CLAUDE.md constraints preserved throughout:** `dynamic=False` embeds, no raw `<`/`>` in slider
  dimension labels, `implausible_bd` grey-dash/grey-column treatment, commit `_freeze/` after
  notebook changes, keep `.py`/`.ipynb` in sync, **do not `git commit`** (user stages manually).

## Out of scope

- New figures/tables beyond the headline FAO figure + linked table (other ideas deferred — "focus
  on NB2/the headline for now").
- Hover on non-FAO plots (NB3, NB1 line/bar plots).
- Demoting README per-notebook "Key methodologies" blocks.
- Any change to the underlying Rosetta / M&M / Saxton–Rawls science.
