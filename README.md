# Reproducible-Science

Reusable, tested workflows for reproducible analysis. Every workflow runs as ordered checkpoints that save a restart file, stop on a failed check, and report in plain English what was checked, what was found, and why it passed or failed, with a quality rating.

| Folder | What it is |
|---|---|
| `data-analysis-workflow/` | The general rules and skills for any quantitative dataset: intake, QC, normalization, structure, group comparison, biological meaning, robustness. Start here for proteomics, CRISPR screens, bulk omics or assays. |
| `scrna-workflow/` | The single-cell RNA-seq workflow: a restartable notebook, engine, tests and eight single-cell skills. It reuses the general rules. |
| `figure-workflow/` | Rules and agents for making and reviewing figures: a shared rubric, a generation workflow, a figure-maker agent and a read-only figure-reviewer agent. Use it after the analysis checkpoints. |

The skills in both folders are also published to Claude as skills, so a new conversation can load them by name (`data-analysis-workflow`, `scrna-seq-workflow`, ...). The folder copies and the published skills are kept in step by hand.

## Figure workflow

Built from Lecture 2 (data visualization: "a figure is a contract"). It has four parts:

- **`figure-rubric` (skill):** the single list of 35 numbered rules, each with its source slide, a PASS/WARN/FAIL severity and an origin tag (`NOTES` from class notes, `SLIDE` stated on a slide, `ADDED` added by the assistant and awaiting confirmation). The rules cover not distorting magnitude (no 3D, bars from zero, area-proportional sizes), evidence (show n, define error bars, show pairing), chart choice, readable color, simplicity and reproducibility. It also ships code helpers: `lint_script` (checks a plotting script) and `audit_figure` (checks a live figure).
- **`figure-generation` (skill):** the order of work. Raw data first goes through the analysis checkpoints (for example clustering or a group comparison) and is saved as results. Then a figure spec is written, the figure is built and self-checked, and it is sent to review.
- **`FIGURE_MAKER` (agent):** builds the figure by following the two skills above and sends every saved figure to the reviewer automatically. It revises until there is no FAIL, for at most 3 rounds, and it never grades its own work.
- **`FIGURE_REVIEWER` (agent):** a read-only critic. It runs the code checks first, then judges the image, script, spec and data summary, and returns a PASS/WARN/FAIL/N/A card per rule.

```
raw data -> analysis checkpoints -> saved results -> figure spec -> build -> self-lint -> FIGURE_REVIEWER -> revise (max 3 rounds) -> final figure + review card
```

Details are in `figure-workflow/README.md`. As with the other folders, the live versions are the published skills and agents in Claude, and the copies here are kept in step by hand.
