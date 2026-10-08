---
name: figure-generation
description: "Workflow for making a deliverable figure that follows the class data-visualization rules: write a figure spec, pick the chart from the question, build from saved analysis results with the helper functions, self-lint, save script + vector + png, then hand to the figure-review agent and revise (max 3 rounds). Use when asked to make, plot, draw or redo any figure that will be shown, reported or kept. Follows the checkpoint style of data-analysis-workflow (figure checkpoint after group comparison and interpretation)."
---

# Figure generation: the decision walkthrough

Load `figure-rubric` first (rules, IDs, helpers). Load `figure-style` for layout mechanics and `figure-composer` for multi-panel figures. This skill is the order of decisions. Do not skip a stage; when a stage's question cannot be answered from the request or the data, ask the user.

## Stage 1: What is the question?
1. State the ONE question or claim the figure answers, and whether it is for exploration or explanation (P1, S1). Ask the user if unclear.
2. Name the audience and the final size (single column, full width, slide) so fonts can be judged at that size.

## Stage 2: Is the data ready to plot?
3. Inspect the input: raw or already analyzed? units, ranges, missingness, what one n is (`data-intake-integrity`).
4. If it is raw, or no saved results exist, run the week-1 analysis checkpoints BEFORE any figure (`data-analysis-workflow`; single-cell: `scrna-seq-workflow`). Stop on FAIL. Save every restart file and results table.

| The figure needs | Run (skill) | Saved result to plot |
|---|---|---|
| Any figure | `data-intake-integrity`, `data-qc-filtering` | clean table, n per group, exclusions |
| Comparable values across samples | `data-normalization-features` | normalized matrix, chosen features, log base |
| Clusters, gradients, batch effects (scatter, PCA, UMAP, heatmap) | `data-structure-discovery` (single-cell: `scrna-dimensionality-reduction`, `scrna-clustering`) | coordinates, cluster labels, stability |
| Group differences (volcano, MA, effect-size, paired plots) | `data-group-comparison` (single-cell: `scrna-condition-comparison`) | effect sizes, intervals, adjusted p-values |
| What features mean (enrichment, annotation) | `omics-biological-interpretation`, `scrna-cell-type-annotation` | enrichment table, labels |
| Trust in the claim | `data-interpretation-robustness` | robust, sensitive or unstable labels |

5. Choose the analysis from the question, not the reverse. A summary that would hide structure (E9) means explore first with PCA or clustering, check cluster stability, and show the raw points or embedding. A difference question needs a replicate-aware comparison, not mean and SD bars. Carry every checkpoint warning and every not-run check into the spec and the caption. Plot only what the analysis supports.

## Stage 3: Choose the style
6. List the variable types (categorical, ordered, continuous, time, paired) and map the question to a chart with rubric C1. Write down the chosen style AND the best alternative you rejected and why (C6).
7. For complicated data keep a clear overview (heatmap, small multiples, overview plus detail) and split into several figures rather than overload one (S1). Never simplify by removing evidence.

## Stage 4: Decide the scale and encoding
8. Baseline and range: bars start at zero; points and lines may use a focused range with a stated reason (H3). Area, not radius, for size (H4). Shading only with meaning (H5, K4).
9. Linear or log: if the question does not settle it, render BOTH and let the user choose, then record the choice (H6). Binned data: try several bin widths and keep the one that preserves the structure (C3). Compared panels share axes, bins and order (H9, U3).
10. Lines only for meaningful order; break at gaps; never fill missing values (H7, H8).

## Stage 5: Decide color
11. Categories get distinct readable hues (`categorical_colors`, at most 6; never neon or cyan or yellow); ordered magnitude gets a sequential map; change around a meaningful zero gets a diverging map centered there (`colormap_for`, `symmetric_limits`); never rainbow. Never color alone; keep the same color for the same group across figures (K1-K7).

## Stage 6: Decide typography, title and legend (T1-T4)
12. Call `apply_figure_rules()` first. Size ladder at the base of 10 pt: title 13 bold > axis labels 11 regular > ticks 10 = legend 10. Nothing below 8 pt at final size; check at the final print size, not on screen.
13. Title: exactly one, at the TOP (default) or the BOTTOM (`title_bottom`), bold and larger than the axis labels, never inside the plot area. It describes what is plotted (outcome, units, comparison) and does not persuade (U2). Panel titles bold at least axis-label size.
14. Axis labels regular weight, larger than ticks, with units. Legend: at least tick size, outside the axes or in empty space, not covering data, entries in visual order, at most 6 (else label directly), no frame.

## Stage 7: Evidence and unbiased presentation (E1-E9, U1-U5)
15. Show n and denominators, define error bars (SD, SE or 95% CI) and the independent n, show raw points when n is small, show pairing, add the control or reference, annotate the specific comparison.
16. Same filters, scales, bins and ordering rule for every group; exclusions stated with counts; uncertainty and contrary or null results shown; limits and analysis warnings in the caption with a one-line takeaway the data support. As clear as the data honestly allow, never clearer.

## Stage 8: Build, check, save
17. Write the figure spec first with `figure_spec_template()` (question, data source, variable types, style and rejected alternative, n, error bars, reference, baseline reason, missing-data handling, color plan, title and legend placement, analysis checkpoints run). Save it as `figN.spec.md`. Also plan an on-figure caption (see 20).
18. Build from the saved results file, never in-memory variables (P2), with fixed seeds (P3). Make the script self-contained: do not `exec` helper files. If helpers must be shared across figures, list every dependency file for the reviewer.
19. Run `audit_figure(fig)` before saving and `lint_script(script)` after. Fix every FAIL; fix or justify every WARN in the spec. Save the audit and lint output as `figN.audit.txt` so the reviewer can see it.
20. Put a short caption ON the figure (`fig.text`, at least 8 pt): n and exclusions, data origin (for example simulated), the analysis warnings that limit the claim, and a one-line takeaway the data support (U1, U5). A reader of the image alone must see the limits. Then save `figN.py`, `figN.svg` (or pdf), `figN.png`, `figN.spec.md` and `figN.audit.txt` as artifacts.

## Stage 9: Automatic independent review (mandatory)
21. As soon as the figure is saved, delegate it to the `FIGURE_REVIEWER` agent from the `repl` tool with `host.delegate({"name": "Review", "profile": "FIGURE_REVIEWER", "task": ...})`. The task holds the image, script and spec as literal artifact markers plus a data summary: n, ranges, filters, and the analysis checkpoints run with their verdicts, warnings and not-run list. Include the audit/lint output and every dependency file. If you are running as a sub-agent and cannot dispatch the reviewer, say so, label the figures UNREVIEWED DRAFTS and list exactly what the parent must send; never present them as reviewed.
22. The reviewer returns a card with Style fit, Bias check and Typography lines. Revise on every FAIL. For every WARN either fix it or write a one-line justification in the spec, then re-save and review again after any change. Stop when there is no FAIL and every WARN is fixed or justified, or after 3 rounds, and report what remains. The maker never grades its own figure.

## Stage 10: Report
23. Report in the week-1 style: what each rule checked, what was found, PASS/WARN/FAIL with reason, the analysis verdicts that fed the figure, what was not checked, and the final `figN.review.md`.

## Never break
No 3D. Bars from zero. No rainbow or neon. Size by area. No line across a gap or between unordered points. No means without n and defined error bars. No text under 8 pt. Title bold and larger than axis labels. Legend legible and off the data. A style that fits the question. Nothing hidden, nothing favored, uncertainty shown. No claim stronger than the design supports.
