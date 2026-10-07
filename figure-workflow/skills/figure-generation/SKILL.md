---
name: figure-generation
description: "Workflow for making a deliverable figure that follows the class data-visualization rules: write a figure spec, pick the chart from the question, build from saved analysis results with the helper functions, self-lint, save script + vector + png, then hand to the figure-review agent and revise (max 3 rounds). Use when asked to make, plot, draw or redo any figure that will be shown, reported or kept. Follows the checkpoint style of data-analysis-workflow (figure checkpoint after group comparison and interpretation)."
---

# Figure generation (checkpoint F)

Load `figure-rubric` first (rules, IDs, helpers); load `figure-style` for layout mechanics and `figure-composer` for multi-panel figures. Do not duplicate them. This skill is the order of work.

## Step 0: raw data first becomes saved results
Raw tables are rarely plottable as they are. If the input is raw (counts, intensities, readings, a sample-by-feature table) or no results file exists, do the analysis with the week-1 checkpoint skills BEFORE any figure, following `data-analysis-workflow` (single-cell data: `scrna-seq-workflow`). Load each checkpoint skill when you reach it, save its restart file and results table, and stop on FAIL.

| Figure needs | Run (skill) | Saved result to plot |
|---|---|---|
| Any figure | `data-intake-integrity`, `data-qc-filtering` | clean table, n per group, exclusions |
| Comparable values across samples | `data-normalization-features` | normalized matrix, chosen features, log base |
| Clusters, gradients, batch effects (scatter/PCA/UMAP, heatmap) | `data-structure-discovery` (single-cell: `scrna-dimensionality-reduction`, `scrna-clustering`) | PCA/embedding coordinates, cluster labels, stability |
| Group differences (volcano, MA, effect-size or paired plots) | `data-group-comparison` (single-cell: `scrna-condition-comparison`) | effect sizes, intervals, adjusted p-values |
| What features mean (enrichment, annotation) | `omics-biological-interpretation`, `scrna-cell-type-annotation` | enrichment table, labels |
| Trust in the claim | `data-interpretation-robustness` | robust/sensitive/unstable labels |

Rules: (1) plot only what the analysis actually supports; carry its WARNs into the spec and caption (E7). (2) If a summary would hide structure (E9), explore first with PCA or clustering, check cluster stability, and show the raw points or embedding rather than only means. (3) Choose the analysis from the question, not the other way around: a clustering question needs structure discovery; a difference question needs a replicate-aware comparison, not a mean-and-SD bar. (4) The data summary sent to the reviewer lists the checkpoints run, their verdicts and warnings, and what was not run.

## Steps
1. **Inspect** the saved results file (units, ranges, missingness, sample composition, what one n is). Never plot from in-memory variables (P1, P2).
2. **Write the spec** with `figure_spec_template()`: one question, chart type and why (C1), n and denominators (E1), error-bar meaning (E3), reference/control (E6), baseline and range reasons (H3), missing-data handling (H8), color plan (K1-K7), purpose (exploration or explanation). Save it as `figN.spec.md`. If the question is unclear, ask the user; do not guess.
3. **Choose the chart from the question** using rubric C1. No 3D, no explode/shadow pies (H1, H2).
4. **Decide scale deliberately.**
   - Bars start at zero. Points/lines may use a focused range, with a stated reason (H3).
   - Linear vs log: if the question is not clear, make BOTH, show the user, and record which was chosen for later figures (H6).
   - Binned data: try several bin widths and keep the one that preserves the structure (C3).
   - Panels to be compared share axes (H9).
5. **Build** with `apply_figure_rules()`, `categorical_colors`, `colormap_for`, `bubble_area`, `break_at_gaps`. Show raw points when n is small (E4); define error bars and n (E1, E3); add control/reference and a specific annotation (E6, E8). More than about 6 series -> split or bin (S1).
6. **Self-check:** `audit_figure(fig)` before saving, `lint_script(script)` after. Fix every FAIL. Fix or justify every WARN in the spec.
7. **Save** `figN.py`, `figN.svg` (or pdf), `figN.png`, `figN.spec.md` via the artifact tool.
8. **Automatic review (mandatory).** As soon as the figure is saved, delegate it to the `FIGURE_REVIEWER` agent from the `repl` tool with `host.delegate({"name": "Review", "profile": "FIGURE_REVIEWER", "task": ...})`. The task must contain the image, script and spec as literal artifact markers plus a data summary (n, ranges, filters). The reviewer runs `lint_script` itself and returns a card. Revise on every FAIL and fix or justify every WARN, re-save, and review again. Stop when there is no FAIL, or after 3 rounds and report what remains. The maker never grades its own figure.
9. **Report** in the week-1 style: what each rule checked, what was found, PASS/WARN/FAIL with reason, a rating, what was not checked, and the final `figN.review.md`.

## Rules the generator must never break
No 3D. Bars from zero. No rainbow or neon colors. Size by area. No line across a gap or between unordered points. No means without n and defined error bars. No claim stronger than the design supports.
