---
name: data-qc-filtering
description: "Checkpoint 2 of the general data-analysis workflow (data-analysis-workflow), run after data-intake-integrity passes. Use to find and, where justified, remove damaged, empty, outlying or duplicated samples and uninformative features in any quantitative dataset (omics, screens, assays, imaging, tables). Checks total signal and features detected per sample, outliers by robust distance, duplicates, replicate agreement, technical feature share (for example mitochondrial reads, contaminants, non-targeting controls), balance of groups across batches, how much each filter removes overall and per group, and whether filtering is biased against one group. Plain-English tests with reasons and a PASS/WARN/FAIL verdict."
---

# QC and filtering (checkpoint 2)

Question answered: which samples or measurements are damaged, empty, outlying or duplicated, and are the filters sensible and fair to every group?
Follow the rules of `data-analysis-workflow`. Compute the metrics first, show them, then filter; report what every filter removed.

## Tests, plain definitions and defaults

| Test | What it checks, in plain words | Default rule |
|---|---|---|
| Total signal per sample | Whether any sample has far less (or far more) total signal than typical | WARN for samples under 0.2 or over 5 times the median; FAIL if more than 30% of samples are flagged |
| Features detected per sample | Whether any sample detects far fewer features than typical | WARN for samples under 0.4 times the median |
| Outliers | Whether any sample sits far from the rest on the QC metrics or in the first principal components | Flag robust z above 4 (median and MAD); WARN if more than 5% of samples are flagged |
| Duplicates | Whether any two samples are nearly identical | WARN above correlation 0.999 between different samples |
| Replicate agreement | Whether replicates of the same group resemble each other more than other groups | WARN if a sample's average correlation with its group is below 0.8, or lower than its correlation with another group |
| Technical feature share | Share of signal from features that mark a technical problem (mitochondrial or ribosomal reads in RNA, contaminants or reverse hits in proteomics, non-targeting guides in screens, blanks in assays) | Check the feature set is actually present and has signal before filtering on it; WARN if no sample has any |
| Batch and group balance | Whether every batch contains every group | WARN if a batch lacks a group; FAIL if every batch holds one group only (the group effect cannot be separated from batch) |
| Filter impact | How many samples and features each filter removes, overall and per group | WARN if a filter removes over 10% of samples; FAIL if over 50% |
| Fair filtering | Whether one group loses much more than another | WARN if group loss differs by more than 20 percentage points |
| Feature filter | Features seen in too few samples are dropped (at least as many samples as the smallest group) | INFO with count removed; WARN if over 80% of features are removed |
| Controls behave | Positive, negative and blank controls give the expected signal | FAIL if negative controls show strong signal or positive controls none |

## Choosing values for other data types

- **RNA-seq (bulk):** total counts, genes detected, mitochondrial and ribosomal share; samples under 1 million counts or 10,000 genes detected are worth a look.
- **Proteomics / metabolomics:** number of identified features, share missing per sample, and total ion current; remove contaminants and reverse-database hits as a feature filter.
- **CRISPR screens:** reads per guide and the number of guides with zero reads; the library should be covered (typically at least 90% of guides above 30 reads in the starting sample); plasmid or time-zero distribution should be even.
- **Plate assays:** edge effects (compare outer and inner wells), Z-prime of controls (0.5 or more is good), and drift across plates.

## Pitfalls

- A "technical feature" set can be empty. For example, a dataset can list mitochondrial genes that have zero counts everywhere, and then a "mitochondrial genes found" test passes while the filter removes nothing. Always check that the set has signal.
- Filtering can remove a whole real condition (damaged by the treatment itself). Always compare what each group loses.
- Every flagged sample needs a decision on the card (removed or kept, and why); a sample can be flagged only because of its batch and then be kept and modelled.
- Do not remove samples to improve a result. Remove them for a reason stated before looking at the comparison.

## Not run is a result

Metrics that cannot be computed (no replicates, no batch column, no controls) go under *Did not run* with the reason.
