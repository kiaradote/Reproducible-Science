---
name: data-intake-integrity
description: "Checkpoint 1 of the general data-analysis workflow (data-analysis-workflow): verify the data before any analysis. Use right after downloading or receiving any dataset (omics counts or intensities, assay or plate readings, screens, surveys, measurement tables). Checks that files are complete and readable, records checksums, that rows and columns match what was expected, that sample and feature IDs are unique and agree with the sample sheet in content and order, how much is missing, what scale the numbers are on (counts, intensities, logged, normalised, scaled), that nothing impossible is present, and that groups, replicates and batches are described. Every test is reported in plain English with its reason and a PASS/WARN/FAIL verdict."
---

# Data intake and integrity (checkpoint 1)

Question answered: are these the intended files, complete, readable and consistent with what was expected, and what kind of numbers are in them?
Follow the rules of `data-analysis-workflow` (stop on FAIL, plain-English card with *what it checked / what we found / result and why*, quality rating, restart file). Nothing is filtered or changed here: this checkpoint only looks and records.

## Tests, plain definitions and defaults

| Test | What it checks, in plain words | Default rule |
|---|---|---|
| File complete | The file opens to the end (compressed files are not cut short) and its size matches the expected size, if one is known | FAIL if it cannot be read to the end or the size differs from the reference |
| Checksum recorded | A fingerprint (SHA-256) of every input file is stored so a later reader can tell if the data changed | INFO when no reference exists (store it); FAIL if it differs from a given reference |
| Shape as expected | The number of samples and features equals what the source says | FAIL on any difference when an expected number is given; INFO otherwise |
| Unique IDs | No sample ID and no feature ID appears twice | FAIL if any duplicate |
| Sample sheet agrees with the data | The samples in the table and in the sample sheet are the same set, and in the same order | FAIL if the sets differ; FAIL if the order differs (labels would be attached to the wrong samples) |
| Missing values | How much of the table is empty, per sample and per feature | WARN if over 30% of all values are missing or any sample is over 50% missing; FAIL if a whole sample or the whole table is empty |
| Impossible values | No infinite values, no negative values where only zero or more can occur (counts, intensities before logging), no non-integer counts if the data are declared as counts | FAIL on any |
| Scale detected | What the numbers look like: non-negative integers (counts), positive continuous with a long right tail (intensities), roughly symmetric around small values (logged), centred with negative values (scaled) | FAIL if the declared scale disagrees with the detected one; INFO records the detection |
| Duplicated samples | No two samples have identical values | FAIL if two samples are exact copies; WARN if correlation is above 0.9999 |
| Design described | The sample sheet names the group, the replicate or subject, and the batch (or says there is none) | WARN if group, replicate or batch is missing; INFO shows group sizes |

## Choosing values for other data types

- **Count data (RNA-seq, CRISPR guide counts, ATAC peaks):** values are whole numbers. A non-integer or negative value is a FAIL. Feature IDs should be genes, guides or peaks from a known annotation.
- **Intensity data (proteomics, metabolomics, imaging features):** values are positive and skewed; missing values are expected and not random (low-abundance features go missing). The 30% missing default is a start; for label-free proteomics 40 to 60% can be normal, so raise the limit and say why.
- **Assay or plate readings:** the plate layout (row, column, plate) belongs in the sample sheet; check that controls (blanks, positives, negatives) are present and identified.
- **CRISPR screens:** check that every guide maps to one target, that non-targeting and essential-gene control guides exist, and that the library design file lists each guide once.
- **Survey or clinical tables:** categorical columns have the allowed levels only; dates and units are consistent.

## Pitfalls

- A matrix whose header is shifted by one column still has the right shape and the right IDs; only comparing the order with the sample sheet catches it.
- Declared and detected scale can disagree when someone saved a log-transformed table as "counts". Normalising twice gives wrong results without any error.
- Do not fix a problem here silently. Report it, and let the user decide how to fix it.

## Not run is a result

If an expected number, a reference checksum or a sample sheet is not available, list the test under *Did not run* with that reason.
