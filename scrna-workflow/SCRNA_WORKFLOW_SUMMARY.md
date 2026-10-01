# scRNA-seq workflow: checkpoints at a glance

Run `scrna_checkpoints.ipynb` top to bottom (set `NB_DATASET`). Every stage saves a restartable file, writes a report, and updates `runs/<dataset>/RUN_SUMMARY.md`.

| # | Stage | Question | Skill | Restart file |
|---|---|---|---|---|
| 0-1 | Fetch, assemble | Right files, complete, aligned? | scrna-data-integrity | 00_fetch.json, 01_assemble.h5ad |
| 2 | QC | Which cells are damaged or empty? | scrna-qc-filtering | 02_qc.h5ad |
| 3 | Features | Normalised, informative genes chosen reproducibly? | scrna-normalization-features | 03_features.h5ad |
| 4 | Dimensionality reduction | Structure captured, no technical driver? | scrna-dimensionality-reduction | 04_dimred.h5ad |
| 5 | Clustering | Groups stable and distinct? | scrna-clustering | 05_cluster.h5ad |
| 6 | Annotation | Which cell types, how sure? | scrna-cell-type-annotation | 06_annotate.h5ad |
| 7 | Condition comparison | What differs between conditions (if any)? | scrna-condition-comparison | 07_compare.h5ad |
| 8 | Interpretation, validation | What does it mean, what survives changed choices? | scrna-interpretation-validation | 08_interpret.h5ad |

Master rule: skill `scrna-seq-workflow` (order, stop on FAIL, checkpoint card, quality scale, finished summary).

Every checkpoint reports: verdict and quality rating, tests run (what each tests, result), conclusion, Not run list, caveats carried, how to go back.
Quality scale: Good (all checks passed) / Usable, with caveats (a check warned) / Accepted despite a failed check / Not reliable until fixed (a check failed) / Not done. The run is rated by its worst stage.

Go back: change a parameter and rerun the notebook (unchanged earlier stages load from disk), or `wf.invalidate("<stage>")`.
Tested on: pancreas donor 1 (GSE84133, counts), PBMC control vs interferon-beta (GSE96583 batch 2, 8 paired donors). Large GSE120575 not run beyond loading (compute).
Not implemented: transcription-factor activity, cell-cell communication, spatial association; ribosomal percentage and doublet detection in QC.
