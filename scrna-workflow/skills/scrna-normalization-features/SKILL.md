---
name: scrna-normalization-features
description: "Checkpoint 3 of a reproducible scRNA-seq workflow: normalization and feature selection, run after scrna-qc-filtering. Use when making single-cell expression profiles comparable between cells and choosing informative genes: picks the normalization from the detected data scale (counts, TPM, log2 TPM, log-normalized), records the log base, selects 2000 highly variable genes (seurat_v3 on counts, seurat on log data, optionally per batch), and checks per-cell totals, depth spread, dominant-gene cells, technical genes among the selection, and selection stability. Produces a PASS/WARN/FAIL report."
---

# scRNA-seq normalization and feature selection (checkpoint 3)

Question answered: can the expression profiles be compared between cells, and which genes carry the informative variation?
Scope stops at variable-gene selection. Scaling, PCA, batch correction and clustering belong to the next checkpoint.
Start only after checkpoints 1 and 2 have no FAIL. Input: the QC checkpoint (cells x genes, sparse) and the scale detected at
checkpoint 1 (`counts`, `linear_tpm`, `log2tpm`, `log1p_cp10k`; stored as `adata.uns["scale"]`).
Records use the same `{check, status, detail}` format as the other checkpoints. Call `write_feat_report(records, "features_report")`,
show the overall status, and do not continue on FAIL without asking the user. Carry every WARN into the final report.

Helpers (loaded from `kernel.py`; need scanpy, numpy, scipy, pandas; scikit-misc for `seurat_v3`):
`feat_run`, `feat_normalize`, `feat_check_normalization`, `feat_select`, `feat_check_selection`, `feat_plot`, `write_feat_report`,
`feat_catalog`, `linear_stats`, `feat_record`.

## Workflow

1. **One call:** `adata, records = feat_run(adata, scale, species="human", n_top_genes=2000, batch_key=None, sample_col=None, plot_path=...)`.
   It normalizes, checks, selects genes, checks again. Then `write_feat_report(records, prefix)` and save `adata` as the checkpoint.
2. **Method by scale (the standard choices):**
   - `counts`: normalize to 10,000 per cell, natural log1p; raw counts kept in `layers["counts"]`; `seurat_v3` on the raw counts.
   - `linear_tpm`: natural log1p; input kept in `layers["input"]`; `seurat` on the logged values.
   - `log2tpm` and `log1p_cp10k`: already normalized, so `X` is left unchanged; the log base is recorded in `uns["log1p"]`
     (base 2 for log2 TPM) and `seurat` is used.
   - scaled or unrecognized values: refused (FAIL). Go back to the unscaled matrix.
3. **Number of genes:** 2000 is the common default. Use `batch_key` (donor, sample or batch column) whenever there is more than one
   batch: genes are then ranked by how many batches call them variable, which keeps single-batch artefacts out. Pass `sample_col`
   to get the depth-spread check.
4. **Read the Not run list.** A check that did not run has not passed (for example batch robustness needs `batch_key`).

## Checks (thresholds are defaults, not literature values; change them with the function arguments)

- Normalization: method fits the scale; per-cell linear totals have CV below 0.1; no NaN, infinite or negative values; sparsity unchanged;
  input values kept; median depth within 5-fold across samples; no more than 1% of cells with a single gene above 50% of their total.
- Selection: requested number of unique genes; mitochondrial, ribosomal and haemoglobin genes at most 10% of the selection; at most 10%
  of selected genes detected in under 1% of cells; selected genes not concentrated in the lowest-expressed 10% of genes (limit 30%);
  Jaccard overlap of selections on two random halves of the cells at least 0.5 (WARN only); with `batch_key`, at least half of the
  selected genes variable in at least half the batches.

## Pitfalls

- **Log base.** Scanpy's `seurat` flavor assumes natural-log data unless `adata.uns["log1p"]["base"]` is set (checked in scanpy 1.10.4).
  Log2(TPM+1) data without the base set is treated as natural log. In a test on log2(TPM+1) values derived from a real count matrix
  (GSE84133 donor 1, 388 cells) the wrong base changed 259 of 2000 selected genes. The Sade-Feldman example object set only
  `uns["data_type"]`, so its selection probably had this problem; its saved object was not available to confirm.
- `seurat_v3` needs raw counts; never run it on TPM or log data.
- Select genes on normalized, unscaled values; scale only afterwards, on the selected genes.
- Library-size normalization cannot remove real composition or depth differences between samples; the depth and dominant-gene checks flag them.
- Plate-based TPM data has no recoverable read depth; the depth check then uses genes detected per cell.

## Tests

`test_scrna_features.py`: counts path on synthetic and real data, log2 TPM and linear TPM paths, base handling against explicit natural
log, scaled data refused, and the warnings (depth spread, dominant-gene cells, technical genes, low-detection genes, unstable selection),
batch check and report text.
