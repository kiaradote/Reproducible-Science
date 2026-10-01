---
name: scrna-qc-filtering
description: "Checkpoint 2 of a reproducible scRNA-seq workflow: quality control and filtering, run after scrna-data-integrity has passed. Use when computing QC metrics and filtering cells and genes in a single-cell dataset: validating the mitochondrial gene set (MT- prefix, nuclear MTRNR2L pseudogenes excluded, datasets with no mitochondrial genes) and reporting cells removed per filter, in total and per sample. Produces a PASS/WARN/FAIL report."
---

# scRNA-seq QC and filtering (checkpoint 2)

Scope: decide which cells and genes to remove, and show exactly what each filter did. Start only after checkpoint 1
(`scrna-data-integrity`) has no FAIL. Records use the same `{check, status, detail}` format as checkpoint 1, so both
reports can be concatenated. Call `write_qc_report(records, "qc_report")` and show the overall status;
do not continue to normalization or clustering on FAIL without asking the user. Carry every WARN into the final report.

Helpers (loaded from `kernel.py`; need pandas): `check_mito`, `check_mito_signal`, `qc_filter_report`, `write_qc_report`, `qc_record`, `qc_catalog`.

## Workflow

1. **Mitochondrial genes.** `records, mito = check_mito(var_names, species)`. Uses the `MT-` (human) or `mt-` (mouse)
   prefix only; `MTRNR2L*` nuclear pseudogenes are excluded and reported. FAIL if no mitochondrial genes exist
   (Ensembl IDs, wrong species, or genes removed by the authors). In that case a mito-percentage filter is meaningless
   and must be dropped or justified. Compute `pct_mt` only from the returned `mito` list, then call `check_mito_signal(pct_mt, max_pct_mt)`: WARN when almost no cell has any
   mitochondrial counts although the genes are listed (the matrix was filtered upstream; found in the Kang 2018 PBMC release, where all 13 MT- genes are present with zero counts). The mito filter then removes nothing and cannot find stressed cells.
2. **Choose thresholds before filtering** from the QC plots (genes per cell, total counts or TPM, pct_mt), per sample as well
   as pooled, and write them into the parameter file. Do not tune thresholds to obtain a desired cluster result.
3. **Filter report.** `qc_filter_report(obs, {name: removed_mask}, sample_col=...)` gives cells removed per filter,
   total loss (WARN above 20%) and per-sample loss (WARN above 50% in any sample). The 20% and 50% limits are defaults,
   not literature values; change them with `max_total_frac` and `max_sample_frac`.
4. **Report.** `write_qc_report(records, prefix)` writes `<prefix>.json` and `<prefix>.md`. The report opens with the
   verdict and a one-sentence **Why**, then a table of every check that ran with what it tests and its result, then a
   **Not run** list of catalogued checks that were never called (for example per-sample loss when no sample column was
   given). A check that did not run has not passed. `qc_catalog()` lists the checks and their purposes.

## Reference: what the Sade-Feldman example did (GSE120575, Smart-seq2, log2(TPM+1))

- Authors had already kept cells with at least 1,093 genes (0 cells below 250 genes), so no gene-count filter was applied.
- Mito set: names starting `MT-` or `MTRNR`, 50 genes. Likely 37 mitochondrial plus 13 nuclear MTRNR2L genes, which this skill excludes.
- pct_mt (computed on linear TPM): median 12.74, 95th percentile 24.75, maximum 90.44. Filter: pct_mt below 40 kept,
  removing 231 cells (16,291 to 16,060). Genes seen in fewer than 3 cells removed (55,737 to 45,692).
- QC plots pooled all cells; no per-patient or per-plate view.

## Pitfalls

- A "MT-" plus "MTRNR" prefix rule picks up nuclear MTRNR2L genes (37 + 13 = 50 genes in the example above).
- A naive "MT" prefix also matches metallothioneins and unrelated genes (MT1A, MTHFD1, MTCH1); keep the hyphen.
- Some deposited matrices contain no mitochondrial genes at all (for example GSE84133 human donor 1).
- Thresholds that look fine pooled can remove most of one sample; check per-sample loss.

## Tests

`test_scrna_qc.py`: mito prefix rules (pseudogenes, metallothioneins, Ensembl IDs, mouse), a real matrix with no mito genes,
filter reports with a sample losing most of its cells, and report writing.
