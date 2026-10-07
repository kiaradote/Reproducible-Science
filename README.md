# Reproducible-Science

Reusable, tested workflows for reproducible analysis. Every workflow runs as ordered checkpoints that save a restart file, stop on a failed check, and report in plain English what was checked, what was found, and why it passed or failed, with a quality rating.

| Folder | What it is |
|---|---|
| `data-analysis-workflow/` | The general rules and skills for any quantitative dataset: intake, QC, normalization, structure, group comparison, biological meaning, robustness. Start here for proteomics, CRISPR screens, bulk omics or assays. |
| `scrna-workflow/` | The single-cell RNA-seq workflow: a restartable notebook, engine, tests and eight single-cell skills. It reuses the general rules. |

The skills in both folders are also published to Claude as skills, so a new conversation can load them by name (`data-analysis-workflow`, `scrna-seq-workflow`, ...). The folder copies and the published skills are kept in step by hand.
