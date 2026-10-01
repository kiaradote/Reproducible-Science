---
name: scrna-data-integrity
description: "Checkpoint 1 of a reproducible scRNA-seq workflow: verify data integrity before any QC or analysis. Use right after downloading or assembling a single-cell expression matrix (GEO, 10x, h5ad, TPM or counts): file checksums and truncated downloads, header/column alignment against metadata, matrix scale (counts vs TPM vs log), NaN and duplicate IDs, and published shape and label counts. Produces a PASS/WARN/FAIL report. QC and filtering belong to checkpoint 2 (scrna-qc-filtering)."
---

# scRNA-seq data integrity (checkpoint 1)

Scope: is the data what it claims to be? No filtering or biological judgment happens here; QC is checkpoint 2
(`scrna-qc-filtering`). Run this BEFORE QC, normalization, clustering or any statistics. Each check returns a record
`{check, status, detail}` with status PASS, WARN, FAIL or INFO. Collect the records into one list,
call `write_integrity_report(records, "integrity_report")`, and show the user the overall status.
**Do not continue to the next checkpoint while the overall status is FAIL**; tell the user what failed and
let them decide. A WARN must be mentioned in the final report.

Helpers (loaded from `kernel.py` when this skill is loaded; need numpy, pandas, scipy):
`check_file`, `check_alignment`, `check_values`, `check_expected`, `write_integrity_report`, `integrity_record`, `integrity_catalog`.

## Workflow

1. **Files.** `check_file(path, expected_sha256=None, expected_bytes=None)` for every raw download.
   It streams the whole gzip, which catches truncated downloads. If no reference hash exists, record the
   printed sha256 in `DATA_SOURCES.md`.
2. **Alignment.** `check_alignment(matrix_ids, meta_ids, matrix_labels, meta_labels)`.
   Always pass a label row/column that comes from the matrix file itself (for example a sample-label row)
   and the same label from the metadata. Without it the check can only compare ID sets, and a header shifted
   by one position still has identical sets. Any disagreement is a FAIL, not a percentage: a one-position
   shift only disagrees at sample boundaries.
3. **Values.** `records, detected = check_values(X, declared=...)` on the matrix with ALL genes (the sum tests are
   invalid after HVG subsetting or scaling). Detects `counts`, `log2tpm`, `log1p_cp10k`, `linear_tpm`, or
   `scaled_or_unknown`. Declare what the source says; a mismatch is a FAIL. The detected scale decides which methods
   are valid: count-based methods (`seurat_v3` HVGs, scVI) need `counts`.
4. **Expected numbers.** `check_expected(name, observed, expected)` for shape, number of samples/patients and label
   counts, with the expected values taken from the source paper or GEO page (retrieve them; do not recall from memory).
5. **Report.** `write_integrity_report(records, prefix)` writes `<prefix>.json` and `<prefix>.md`; save the .md as an
   artifact and record the overall status in the run log. The report opens with the verdict and a one-sentence
   **Why** (which checks failed or warned, with their details), then a table of every check that ran with what it
   tests and its result, then a **Not run** list of catalogued checks that were never called. Read the Not run list:
   a check that did not run has not passed. `integrity_catalog()` lists the checks and their purposes. When a stage runs only some of the checks (for example a fetch stage runs only the file checks), pass `scope=[...]` with the catalogue names that stage should run so the Not run list is not padded with checks that belong to a later stage. `check_values` works on sparse matrices and inspects only a sample of cells, so it is safe on large datasets.

## Pitfalls seen in practice

- Row/column offsets: the Sade-Feldman GSE120575 TPM matrix has a trailing empty column and a header offset by one
  position, so every column name was wrong until relabeled by position. Verified only by set equality, which cannot
  detect misordering.
- Raw-file parsing can need unusual options (comment lines before the header, latin-1 encoding, misspelled column names).
  Select columns by name and assert they exist.
- Reading a large dense matrix can need many times the file size in RAM (the example peaked near 25 GB for a 127 MB gz).

## Tests

`test_scrna_integrity.py` exercises every check on a real count matrix (GSE84133, GSM2230757) plus injected faults:
shifted header, duplicate IDs, NaN, wrong declared scale, truncated gzip, wrong checksum.
