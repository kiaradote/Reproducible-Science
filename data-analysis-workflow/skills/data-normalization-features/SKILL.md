---
name: data-normalization-features
description: "Checkpoint 3 of the general data-analysis workflow (data-analysis-workflow), run after data-qc-filtering. Use to make samples comparable and choose informative features in any quantitative dataset: picks the normalisation from the detected data scale (counts, intensities, already logged), records log base and pseudocount, handles missing values on purpose and records it, then selects features by variability or other stated rule. Checks that scale factors are not extreme, that sample totals or medians become comparable, that no single feature dominates, that selected features are not technical or rarely observed, and that the selection is stable across halves of the samples and across batches. Plain-English tests with reasons and a PASS/WARN/FAIL verdict."
---

# Normalization and feature selection (checkpoint 3)

Question answered: are values comparable between samples, and are the informative features chosen in a reproducible way?
Follow the rules of `data-analysis-workflow`. The method is chosen from the scale found at intake, never assumed, and everything chosen is written down.

## Method by data scale (defaults)

| Detected scale | Default normalisation | Then |
|---|---|---|
| Counts (RNA-seq, guide counts) | Size factors (median of ratios) or counts per million, then log2(x + 1) for exploration; keep the raw counts for count-based tests | Select the top 2000 variable features, or all features with enough counts for small panels |
| Intensities (proteomics, metabolomics) | log2 transform, then median centring per sample (or quantile normalisation if distributions should be identical) | Handle missing values on purpose (below) |
| Already logged or normalised | Do not log again; check the distribution is comparable between samples | Record that it was supplied as normalised |
| Plate readings | Normalise to plate controls (percent of control, or B-score for plate position effects) | Check edge effects after |

## Tests, plain definitions and defaults

| Test | What it checks, in plain words | Default rule |
|---|---|---|
| Method matches scale | The normalisation fits the detected scale and the log base and pseudocount are recorded | FAIL if a transform is applied twice or a count method is used on non-counts |
| Scale factors | Whether any sample needed an extreme correction | WARN if the largest factor is over 10 times the smallest |
| Samples comparable after | Whether sample medians (or totals) are now similar | WARN if their spread (coefficient of variation) is over 10% after normalisation |
| No dominant feature | Whether one feature takes a large share of a sample's signal | WARN if a single feature holds over 10% of any sample's total |
| Missing values | How missing values were handled and how many were filled in | INFO records the method; WARN if over 20% of values were imputed; imputation is never applied to a column that is mostly missing |
| Selection size | How many features were selected, and by what rule | INFO |
| Technical features among the selection | Whether features that mark a technical problem were selected | WARN if over 10% of the selection |
| Rarely observed features | Whether the selected features are seen in very few samples | WARN if over 10% are detected in under 5% of samples |
| Selection stability | Whether two random halves of the samples choose the same features | WARN if their overlap (Jaccard index) is below 0.5 |
| Batch robustness | Whether the selection agrees when chosen within each batch | WARN if the average overlap between batches is below 0.5 |

## Choosing values for other data types

- **Proteomics:** missing values are often not random; prefer filtering to features seen in at least 70% of samples in one group, and impute (for example from a low-abundance distribution) only when a method needs a complete matrix. Report the share imputed per sample.
- **CRISPR screens:** normalise guide counts to total reads, then compute log fold changes against the starting sample; do not select features by variability, because all guides are kept.
- **Small feature panels (under a few hundred):** skip variable-feature selection and keep all features; say so.

## Pitfalls

- Choose features without using group labels at this checkpoint. Selecting the features that differ between groups here is circular and makes the comparison look better than it is. If no stable set of variable features exists (most features vary by similar amounts), say so and keep every feature that passes the missing-value rule.
- Normalising total signal removes real differences when most features change in one direction (for example a global transcription shift). If the design could cause that, use spike-ins or control features and say so.
- Variable-feature selection run on all samples can pick up batch differences. Compare with the per-batch result.
- A stable selection can still be dominated by technical features. Check both.

## Not run is a result

Stability and batch tests need enough samples (at least 6) and a batch column; if either is missing, list them under *Did not run*.
