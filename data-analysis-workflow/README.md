# General data-analysis workflow

A reusable set of rules and skills for reproducible analysis of any quantitative dataset (omics, screens, assays, tables). The single-cell workflow in `../scrna-workflow` specialises it and reuses these rules instead of repeating them.

## Checkpoints, in order
| # | Checkpoint | Skill |
|---|---|---|
| 1 | Intake and integrity | `skills/data-intake-integrity` |
| 2 | QC and filtering | `skills/data-qc-filtering` |
| 3 | Normalization and features | `skills/data-normalization-features` |
| 4 | Structure discovery | `skills/data-structure-discovery` |
| 5 | Group comparison (if there is a group question) | `skills/data-group-comparison` |
| 6 | Biological meaning (enrichment, expectations, other studies) | `skills/omics-biological-interpretation` |
| 7 | Interpretation and robustness | `skills/data-interpretation-robustness` |

The master rule is `skills/data-analysis-workflow`: run in order, one restart file per checkpoint, stop on FAIL, a plain-English card after every checkpoint (what it checked, what we found, result and why), a quality rating, a running `RUN_SUMMARY.md`, and a brief finished summary.
Helpers: `skills/data-analysis-workflow/kernel.py` (standard library only: `da_record`, `da_card`, `da_rating`, `da_write_report`, `da_run_summary`, `da_is_current`) and `skills/omics-biological-interpretation/kernel.py` (numpy, scipy, pandas: `bio_ora`, `bio_gsea`, `bio_read_gmt`, `bio_expectations`, `bio_concordance`).

## Quality scale
Good / Usable, with N caveats / Accepted despite a failed check / Not reliable until fixed / Not applicable / Not done. The run is rated by its worst checkpoint.

## Worked example
`examples/proteomics_worked_example.md` shows the cards produced for a synthetic proteomics table with planted effects (33 of 40 planted proteins found, 0 of 20 null proteins called); `examples/proteomics_example.py` generates it and doubles as a test.

## Tests
`pytest` from this folder (the helpers need numpy, scipy, pandas and pytest). `tests/da_kernel.py` and `tests/bio_kernel.py` are copies of the two skill kernels; if you change a kernel, change both copies or republish the skill.
Only the proteomics example has been run on data; thresholds for CRISPR, plate-assay and survey data are untested defaults.
