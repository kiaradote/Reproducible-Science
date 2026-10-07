---
name: scrna-seq-workflow
description: "Master rule for any single-cell RNA-seq (scRNA-seq) analysis request. It specialises data-analysis-workflow: load that skill first for the shared rules (ordered checkpoints, a restart file per step, stop on FAIL, plain-English test cards with what was checked, what was found and the result and why, quality ratings, running summary, brief finished summary), then follow the single-cell checkpoints in notebook order: fetch and integrity, QC, normalization and variable genes, dimensionality reduction, clustering, cell type annotation, optional condition comparison, interpretation and validation. Use whenever the user asks to analyse, process, cluster, annotate or compare scRNA-seq or single-cell data, or wants a reproducible, restartable single-cell workflow."
---

# scRNA-seq workflow rule

**Load `data-analysis-workflow` first and follow its rules for everything below.** This skill adds only what is specific to single-cell data: which skill handles each checkpoint, the notebook and engine, and a few single-cell rules. If a rule here differs from the general one, this one wins.
Fallback if the general skill cannot be loaded: run the checkpoints in order, save a restart file after each, stop on FAIL, report each checkpoint in the chat as a table of *what it checked / what we found / result and why* with a conclusion, the checks that did not run and how to go back, and rate it Good, Usable with caveats, Accepted despite a failed check, Not reliable, Not applicable or Not done.

## The checkpoints, in order

| Notebook stage | General skill it extends | Single-cell skill (adds what is specific) |
|---|---|---|
| 0 fetch, 1 assemble | `data-intake-integrity` | `scrna-data-integrity` |
| 2 qc | `data-qc-filtering` | `scrna-qc-filtering` (mitochondrial genes, nuclear pseudogenes, doublets) |
| 3 features | `data-normalization-features` | `scrna-normalization-features` |
| 4 dimred | `data-structure-discovery` | `scrna-dimensionality-reduction` |
| 5 cluster | `data-structure-discovery` | `scrna-clustering` |
| 6 annotate | none (single-cell only) | `scrna-cell-type-annotation` |
| 7 compare (only with a condition and replicate donors) | `data-group-comparison` | `scrna-condition-comparison` (pseudobulk, abundance) |
| 8 interpret | `omics-biological-interpretation`, `data-interpretation-robustness` | `scrna-interpretation-validation` |

Stage 7 is recorded as not applicable when there is no condition. Stage 8 always comes last and uses stage 7's differential expression when it exists.

## Single-cell rules

1. **Notebook and engine.** If `scrna_checkpoints.ipynb` with the `scrna_workflow/` folder is available, run it. It runs the stages in this order, saves `runs/<dataset>/checkpoints/NN_<stage>.h5ad`, writes `reports/<stage>.md`, rebuilds `RUN_SUMMARY.md` with the quality table after every stage, skips stages whose parameters and input are unchanged, and goes back with `wf.invalidate("<stage>")`. Otherwise follow the single-cell skills by hand and write the same files.
2. **The unit is the sample, not the cell.** Donors or patients are the replicates; thousands of cells from one donor are not. Significance in stage 7 always uses sample-level (pseudobulk) totals.
3. **Ask only what the data cannot tell.** The species, which column is the condition, which is the sample (donor) and which is the batch or library. Choose other defaults and say what was chosen.
4. **Published labels are for evaluation only.** Never use a dataset's own author labels as the marker panel or the reference for that dataset.
5. **Large data.** Load sparse and in chunks, record peak memory, run sensitivity analyses on a capped subsample and say so, and pause a heavy run when the user's compute is limited instead of shrinking the analysis silently.
6. **Single-cell wording.** A cluster name means "matches this marker panel", not "proven cell type". When each condition is its own library, a treatment effect is "treatment or library". A cell-type name that follows the condition may be a change of name, not of cells.
7. **Finished summary.** Use the general template with these rows: 0-1 Fetch/assemble, 2 QC, 3 Features, 4 Dimensionality reduction, 5 Clustering, 6 Annotation, 7 Comparison (or "not requested"), 8 Interpretation; add the dataset, number of cells and species to the first line, and point to `runs/<dataset>/RUN_SUMMARY.md` and the checkpoint files.
