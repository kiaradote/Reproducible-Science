---
name: data-group-comparison
description: "Checkpoint 5 of the general data-analysis workflow (data-analysis-workflow), run after data-structure-discovery, only when a group or condition question exists. Use to find what differs between two or more groups in any quantitative dataset. Starts with a design gate (independent replicates per group, pairing, group confounded with batch, groups that follow a technical factor), then fits replicate-aware models (negative binomial for counts, linear models for logged intensities or readings, paired designs with the subject as a blocking factor), corrects for multiple testing, checks a label-swap null and direction consistency across replicates, and reports effect sizes with the question first. Plain-English tests with reasons and a PASS/WARN/FAIL verdict."
---

# Group comparison (checkpoint 5)

Question answered: what differs between the groups, by how much, and can this dataset support that conclusion at all?
Follow the rules of `data-analysis-workflow`. Write the question in one sentence first (which groups, which direction, for which features). Skip this checkpoint, and record it as not applicable, when there is no group question or no replication.

## Design gate (first, before any test)

| Test | What it checks, in plain words | Default rule |
|---|---|---|
| Independent replicates | Each group has several independent samples (subjects, animals, plates); repeated measurements of one sample are not replicates | FAIL under 2 per group; WARN under 3 |
| Pairing | Whether the same subjects appear in every group | PASS when every subject is in all groups (use the subject as a blocking factor); INFO when unpaired |
| Group vs batch | Whether the group is perfectly mixed up with a batch | WARN when every batch holds one group only: the effect is "group or batch" and must be reported that way |
| Group follows a technical factor | Whether the group tracks depth, plate or date | WARN when the group is almost the same as a technical factor |
| Power note | How large an effect could be found with this many replicates | INFO |

## Tests after the gate

| Test | What it checks, in plain words | Default rule |
|---|---|---|
| Model fits the data | Counts use a negative-binomial model; logged intensities or readings use a linear model; paired designs include the subject | WARN if a simpler fallback model was used and why |
| Multiple testing | p-values are corrected across all features tested (Benjamini-Hochberg) | A feature counts when its adjusted p is below 0.05 and its change is at least 1.4-fold (or the stated effect threshold) |
| Label-swap null | Swapping the group labels (within subjects when paired) finds almost nothing | WARN if the swapped data give more than 5 features, or more than 5% of the real count |
| Direction consistency | Features called different move the same way in most individual replicates | WARN if the median share of replicates agreeing is below 75% |
| Effect sizes reported | The size of each change and its uncertainty are given, not just p-values | INFO |
| Biological meaning | What the changed features mean (pathways, categories, other studies) | Done in the next checkpoint with `omics-biological-interpretation`; this checkpoint only produces the ranked table |

## Choosing values for other data types

- **Proteomics:** limma-style linear models on log2 intensities; features with many missing values need a presence/absence test or are listed as not testable.
- **CRISPR screens:** compare guide log fold changes against non-targeting controls; summarise guides per gene (MAGeCK or a robust rank) and check that essential-gene controls drop out.
- **Plate assays:** model plate as a blocking factor; compare to controls on the same plate.
- **Survey or clinical tables:** adjust for stated confounders, report effect sizes with intervals, and keep the same multiple-testing rule.

## Pitfalls

- Testing single cells, wells or pixels as if they were independent replicates inflates significance enormously; the sample is the unit.
- A group difference confounded with batch is not a treatment effect. Say "treatment or batch".
- Do not choose the threshold after seeing the list. State it first.

## Not run is a result

Direction consistency needs paired or replicated samples; the null check needs at least 6 samples. Otherwise list them under *Did not run*.
