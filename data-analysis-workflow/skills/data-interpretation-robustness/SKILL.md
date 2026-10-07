---
name: data-interpretation-robustness
description: "Checkpoint 6 (final) of the general data-analysis workflow (data-analysis-workflow). Use to say what an analysis result means and how far to trust it, for any quantitative dataset. Tests expectations written before looking (known biology or physics, positive and negative controls), validates findings on held-out samples or batches, re-runs the analysis under reasonable changes (subsets, seeds, thresholds, normalisation, feature counts) and labels each conclusion robust, sensitive or unstable, and carries every warning and unrun check from earlier checkpoints into the conclusions. Produces a findings table with evidence, robustness and caveats, and a PASS/WARN/FAIL verdict with plain-English reasons."
---

# Interpretation and robustness (checkpoint 6)

Question answered: what do the results mean, and which conclusions survive changed analysis choices?
Follow the rules of `data-analysis-workflow`. This checkpoint always comes last and reads every earlier report.

## Tests, plain definitions and defaults

| Test | What it checks, in plain words | Default rule |
|---|---|---|
| Expectations met | Things that must be true (known markers, positive controls) are found, and things that must not be true (negative controls) are not | Written before the analysis; WARN if any is not met; an expectation that cannot be tested is listed, not counted as met |
| Held-out validation | A finding or a signature holds in samples (or batches) the analysis did not use | Leave one sample or batch out; WARN if the held-out performance (area under the curve) is below 0.8; samples with too few measurements are dropped and reported |
| Sensitivity to choices | The result is found again after changing the seed, using 80% of samples, a different number of features, another threshold, another normalisation | Label a conclusion robust (found in 90% of changes with the worst case 75%), sensitive (70% or more) or unstable (below) |
| Controls in the final result | Negative controls give no effect and positive controls give the expected one | FAIL if a negative control shows a clear effect |
| Caveats carried | Every warning and every unrun check from earlier checkpoints is listed next to the conclusions it affects | INFO; the findings table shows the caveats |
| Names are hints | Category, pathway or cluster names are described as shared members, not shown activity | Wording check on the report |

## Findings table (always produced)

One row per conclusion: the finding in one sentence, the evidence (numbers), the robustness label, and the caveats that apply. A conclusion labelled sensitive or unstable depends on a choice that could reasonably have been made differently, and the report says so.

## Choosing values for other data types

- **Proteomics:** expect abundant housekeeping proteins to be unchanged and known pathway members to move together; validate on a second batch if one exists.
- **CRISPR screens:** essential genes (ribosomal, proteasome) should drop out and non-targeting guides should not; the strongest hits should replicate across guides of the same gene.
- **Plate assays:** controls on every plate; reproduce a hit on a second plate or day.
- **Any data:** if no held-out sample or batch exists, say that validation did not run. Do not report an in-sample score as validation.

## Pitfalls

- Expectations chosen after seeing the result test nothing. Write them first or label them as post hoc.
- A perfect held-out score with very few samples proves little. Report how many samples were held out.
- Robustness checks on a subsample show how much the answer depends on choices, not what the full-data answer is. State the subsample.

## Not run is a result

If there is no second batch, no controls and no expectations file, validation, controls and expectations go under *Did not run*, and the overall rating says so.
