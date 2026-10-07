---
name: data-structure-discovery
description: "Checkpoint 4 of the general data-analysis workflow (data-analysis-workflow), run after data-normalization-features. Use to find the main patterns in any normalised quantitative dataset: dimension reduction (PCA), how many components carry real structure compared with shuffled data, whether components track technical factors (depth, missing values, batch, plate, run order) instead of biology, optional batch correction, and optional grouping (clustering) with stability checks. Plain-English tests with reasons, figures with one-line captions, and a PASS/WARN/FAIL verdict. Use before comparing groups so that a technical pattern is not mistaken for a group effect."
---

# Structure discovery (checkpoint 4)

Question answered: what are the main patterns (groups, gradients, batch effects), and are they real or technical?
Follow the rules of `data-analysis-workflow`. This checkpoint describes patterns; it does not test hypotheses.

## Tests, plain definitions and defaults

| Test | What it checks, in plain words | Default rule |
|---|---|---|
| Components above noise | How many principal components carry more variation than the same data with each feature shuffled across samples | Keep the smallest number that captures 90% of the structure above the shuffled level; WARN if fewer than 2 are above noise (no clear structure) |
| Technical drivers | Whether a leading component tracks a technical factor (total signal, features detected, share missing, plate, run date) | WARN if a component holding over 30% of the retained variation correlates above 0.5 with a technical factor |
| Batch structure | How much of the structure above noise is explained by batch, and how much by group (variance explained by each, weighted by the share each component carries) | WARN if batch explains more than 20% of it, because batch then has to be handled in the comparison; FAIL if batch and group are the same variable |
| Batch correction effect | If a correction (Harmony, ComBat, removal of batch means) was applied: batch mixing improves and group differences are not erased | WARN if group separation drops by more than half; never correct when batch and group are confounded |
| Neighbour graph | The nearest-neighbour graph is connected and each sample has enough neighbours | WARN if samples are disconnected |
| Grouping stability (if clustering) | Whether groups are found again under resampling, with different settings | WARN if the median agreement (adjusted Rand index) is below 0.8 or a group is recovered in under 60% of resamples |
| Group distinctness | Whether neighbouring groups differ in enough features to be called separate | WARN if two groups differ in fewer than 10 features |

## Choosing values for other data types

- **Few samples (under 20):** PCA and a heat map of sample correlations are enough; skip clustering or report it as exploratory.
- **Proteomics:** run PCA after imputation, and also colour by share of missing values; missingness often forms its own component.
- **CRISPR screens:** the structure should separate by time point or treatment, not by guide library batch; check replicates of the same time point cluster together.
- **Plate assays:** colour PCA by plate, row and column to reveal edge or drift effects.

## Pitfalls

- Measure batch (a category) by the variance it explains, not by correlating components with batch numbers like 0, 1, 2: the correlation depends on the arbitrary order of the labels.
- A component can be perfectly "real" and still technical. Only the correlation with technical factors tells them apart.
- Batch correction can create structure or hide the group effect when batch and group overlap. Check before and after, and record the choice.
- A 2D plot (UMAP, t-SNE) shows neighbourhoods, not distances. Do not read the distance between groups from it.

## Not run is a result

With no batch column or fewer than 6 samples, batch and stability tests go under *Did not run* with the reason.
