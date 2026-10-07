---
name: figure-rubric
description: "Shared, numbered rubric of figure rules from the Lecture 2 data-visualization class (readability, not skewing results, readable color, chart choice, context, simplicity, process), each with its source slide, PASS/WARN/FAIL severity and how it is checked. Load when generating a figure (via figure-generation) or reviewing one (figure-review agent). Ships helpers: lint_script, audit_figure, review_card, apply_figure_rules, categorical_colors, colormap_for, symmetric_limits, bubble_area, break_at_gaps, figure_spec_template."
---

# Figure rubric (single source of truth)

Both the generator and the reviewer read this file, so they cannot drift apart. Principle behind every rule: **a figure is a contract**. Visual magnitude represents numerical magnitude, every visual dimension encodes something stated, and the reader can see n, uncertainty and context.

**Origin tags:** `NOTES` = verdict from the student's own slide notes; `SLIDE` = stated on the slide itself; `ADDED` = added by the assistant from standard practice, not from class (student should accept or change). Severity: FAIL blocks the figure, WARN must be fixed or justified in the spec, N/A when not applicable.
**Check:** CODE = `lint_script`/`audit_figure`; IMG = look at the rendered image; SPEC = compare to the figure spec; DATA = needs the data summary.

## A. Do not distort magnitude
| ID | Rule | Slide | Origin | Sev | Check |
|---|---|---|---|---|---|
| H1 | Pie charts are flat 2D only: no depth, explode or shadow. Use a pie only for a few parts of a whole with a stated total; otherwise bars. | 2,3 | NOTES (pie choice: ADDED) | FAIL | CODE, IMG |
| H2 | No 3D of any kind (3D bars, cylinders, 3D lines/ribbons). Every dimension must encode data. | 4,6 | NOTES | FAIL | CODE, IMG |
| H3 | Bars and filled areas start at zero (length = magnitude). Points and lines may use a focused range (all values far from zero) but the range must still contain the relevant context; if zero lies inside or near the range, include it. 49 vs 51 must look close. | 15-18 | NOTES + SLIDE | FAIL for bars; WARN otherwise | CODE, IMG |
| H4 | Size encodes AREA, not radius (30 vs 10 is 3x, not 9x). Shaded area is proportional to value. Use `bubble_area`. | 13,21 | NOTES | FAIL | CODE, IMG |
| H5 | Shading and fill carry meaning only (density, magnitude); never decorative. | 13,39 | NOTES | WARN | IMG |
| H6 | Linear vs log: choose by the question (absolute differences vs ratios/fold change). If unsure, render BOTH and let the user choose; label the scale and transformation. Ratio plots have a reference at 1 (log2 ratio 0), centered. Equal tick spacing must mean equal quantitative steps. | 14,19,20 | NOTES + SLIDE | WARN | CODE, IMG |
| H7 | Connect with a line only when order is meaningful (time, dose, paired visits). Do not connect independent observations or rows in file order. | 23,46 | NOTES | FAIL if line implies false order | IMG, SPEC |
| H8 | Missing is not measured: break lines at gaps (`break_at_gaps`), show gaps explicitly, do not fill or interpolate unless stated and marked as estimated. | 7,46 | NOTES | WARN | CODE, IMG |
| H9 | Things to be compared share axes (same x and y limits, same bins, same order) across panels and across figures of the same quantity. | 5,28,49,54 | NOTES | WARN | CODE, IMG |

## B. Evidence and context
| ID | Rule | Slide | Origin | Sev | Check |
|---|---|---|---|---|---|
| E1 | Show the total n and denominators (e.g. 57/90 recovered), exclusions and transformations. Tiny n must be visible so readers can judge representativeness. | 12,44,55 | NOTES + SLIDE | FAIL if n absent | IMG, SPEC |
| E2 | Check hidden subgroup/case mix before trusting a pooled comparison; show stratified view when a confounder exists (pooled A<B but A>B within each severity). | 45 | SLIDE | WARN | DATA |
| E3 | Error bars/intervals are defined (SD, SE or 95% CI) and n is the independent unit (cultures/subjects, not cells/rows). | 33,34 | SLIDE | FAIL if undefined | CODE, IMG |
| E4 | Show raw observations when n is small (points over bars; means alone are insufficient). | 29,10 | SLIDE + NOTES | WARN | CODE, IMG |
| E5 | Paired data show the pairing (lines between before/after) and the within-pair difference. | 32,36 | SLIDE | WARN | IMG |
| E6 | Include reference points: control, normal response, baseline, no-change line, comparison group. | 35 | NOTES | WARN | IMG |
| E7 | Title/claim says only what the data and design support (no causal claim from association). | 35,55 | SLIDE | WARN | SPEC |
| E8 | Annotate the specific comparison with outcome, units and direction, not generic labels. | 52 | SLIDE | WARN | IMG |
| E9 | Summary statistics do not replace the plot: mean/SD/r can hide structure (Datasaurus). Plot the data first; use a test or summary that fits the structure (nonlinear, clustered, paired). A correlation matrix still needs scatterplots. | 10,43 | NOTES + SLIDE | WARN | DATA, SPEC |

## C. Right chart for the question
| ID | Rule | Slide | Origin | Sev | Check |
|---|---|---|---|---|---|
| C1 | Chart follows the question: how many/what fraction -> bars/proportion/mosaic; how it varies -> points, histogram, box, violin, ECDF; relationships/clusters -> scatter, binned counts; change -> line/paired lines; many features -> heatmap/small multiples. | 11,56 | NOTES + SLIDE | FAIL if mismatched | SPEC, IMG |
| C2 | Box plots hide multiple modes; use violin (or points) when the shape of the distribution matters, and say smoothing may affect shape. | 31 | NOTES | WARN | IMG |
| C3 | Binning is acceptable but must reveal the true dynamics: try several bin widths and keep the one that survives; use common bins when comparing groups. | 30,49 | NOTES | WARN | IMG |
| C4 | Overplotting: use transparency or binned counts with a density scale (grays or blues); state the scale. | 37 | NOTES | WARN | IMG |
| C5 | Scale axes so clusters and structure are visible (not squashed into a corner). | 34 | NOTES | WARN | IMG |

## D. Color
| ID | Rule | Slide | Origin | Sev | Check |
|---|---|---|---|---|---|
| K1 | Unordered groups get clearly distinct hues, never similar shades of one color. Four groups need four distinguishable colors. | 6,39 | NOTES | FAIL | IMG, CODE |
| K2 | Readable colors only: no neon or pale colors such as cyan or yellow on white. | 6 | NOTES | FAIL | CODE, IMG |
| K3 | Match color to variable type: categorical -> distinct hues; ordered magnitude -> sequential (blues/grays); change around a meaningful zero -> diverging centered on it. No rainbow: it invents boundaries. | 30,38,40 | NOTES + SLIDE | FAIL | CODE, IMG |
| K4 | Shade changes only when they mean something (density, magnitude). Otherwise use different colors. | 39 | NOTES | WARN | IMG |
| K5 | Never color alone: add direct labels, shapes, line styles or panels; check color-vision-deficiency safety. | 41 | SLIDE | WARN | IMG |
| K6 | A color scale is defined (units, limits, center). State any row/column scaling (z-score hides level differences). | 33,42 | SLIDE | WARN | IMG |
| K7 | The same group keeps the same color in every figure of the analysis. | - | ADDED | WARN | SPEC |

## E. Simplicity and legibility
| ID | Rule | Slide | Origin | Sev | Check |
|---|---|---|---|---|---|
| S1 | One question per figure. If lines cannot be followed or the legend needs effort, simplify: bin, color, or split into several figures or small multiples. Exploration may be dense; explanation shows one comparison at a time. Heuristic: more than about 6 colored series is too many. | 24,26,53 | NOTES + SLIDE (limit 6: ADDED) | FAIL if unreadable | IMG, CODE |
| S2 | Readable at final size: tick text at least 8 pt, labeled axes with units, no overlapping text, upright 2D layout. | 9 | NOTES (8 pt: ADDED) | WARN | CODE, IMG |
| S3 | Consistent order, labels and layout across aligned panels. Layout changes must preserve the evidence. | 54 | SLIDE | WARN | IMG |

## F. Process and reproducibility
| ID | Rule | Slide | Origin | Sev | Check |
|---|---|---|---|---|---|
| P1 | Order: inspect data (units, ranges, missingness, composition) -> explore broadly -> state ONE question -> refine (labels, common scales, uncertainty, context). Write the figure spec before plotting. | 47 | SLIDE | WARN | SPEC |
| P2 | Figure reads saved results, not in-memory variables; saved as vector (svg/pdf) plus png with the script. | - | ADDED (week-1 convention) | WARN | CODE |
| P3 | Fixed random seeds; data path recorded in the spec. | - | ADDED | WARN | CODE |

## Helper use
- `lint_script(path_or_text)` -> findings (static, code-level; cannot see the image).
- `audit_figure(fig)` -> findings on a live matplotlib figure before saving.
- `review_card(findings)` -> markdown table.
- Generation helpers: `apply_figure_rules`, `categorical_colors(n)`, `colormap_for(kind)`, `symmetric_limits`, `bubble_area`, `break_at_gaps`, `figure_spec_template`.

## Open items for the student (not decided here)
- Is a pie ever acceptable beyond being flat 2D (current rule: only a few parts of a stated whole)?
- Confirm or change every rule tagged ADDED.
