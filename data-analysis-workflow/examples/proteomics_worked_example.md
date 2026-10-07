# Worked example: the general workflow on a synthetic proteomics table

Synthetic label-free proteomics data: 24 samples (12 control, 12 treated) in 3 batches, 800 proteins. Planted on purpose: 40 proteins raised 2-fold by treatment, batch offsets, one weak sample, one duplicated sample, intensity-dependent missing values, and 20 proteins with no effect used as negative controls. Produced by `proteomics_example.py` using the helpers of `data-analysis-workflow`.

The six cards below are what a conversation following the rule should print in the chat, one after each checkpoint.

## Card 1

**Checkpoint 1: intake: WARN, usable, with 2 caveat(s)**
*Asks:* are these the intended data, complete and consistent?

| What it checked | What we found | Result and why |
|---|---|---|
| how much of the table is empty | 19.1% of all values missing; the emptiest sample is 73% missing | WARN: we warn above 30% overall or 50% in one sample; one sample is above 50% |
| whether two samples are copies of each other | 1 pair(s) with correlation above 0.9999: S10, S11 | WARN: an exact copy is a FAIL, a near copy is a WARN |
| no sample or protein ID appears twice | 24 samples and 800 proteins, all unique | PASS: any duplicate would be a FAIL; there are none |
| the samples in the table and the sample sheet are the same set, in the same order | 24 of 24 match in the same order | PASS: a different set or order would be a FAIL; both agree |
| no negative or infinite values | 0 non-positive, 0 infinite | PASS: any would be a FAIL |
| what the numbers look like (counts, intensities, logged, scaled) | all positive, not whole numbers, largest 1.4e+09 | INFO: this pattern means linear-scale intensities, so a log2 transform is needed; recorded, not scored |
| the sample sheet names group, subject and batch | groups {'control': 12, 'treated': 12}; batches {'b0': 8, 'b1': 8, 'b2': 8} | INFO: recorded, not scored |

*Conclusion:* readable and consistent; one near-duplicate pair and sparse missing values need handling at QC
*Did not run:* file checksum (no reference given)
*Go back:* outputs/01_intake.csv

## Card 2

**Checkpoint 2: QC and filtering: WARN, usable, with 2 caveat(s)**
*Asks:* which samples are damaged or duplicated, and are the filters fair?

| What it checked | What we found | Result and why |
|---|---|---|
| whether any sample has far less signal than typical | 1 sample(s) under 0.2 of the median signal or 0.4 of the median proteins: S07 | WARN: we warn for any such sample |
| whether any sample sits far from the others on median intensity | 2 sample(s) beyond robust z of 4: S07, S17 | WARN: we warn for any outlier |
| every batch contains every group | {'control': {'b0': 4, 'b1': 4, 'b2': 4}, 'treated': {'b0': 4, 'b1': 4, 'b2': 4}} | PASS: FAIL would mean a group effect could not be separated from batch |
| how many samples the filters remove and whether one group loses more | removing S07, S11 = 8% of samples; share lost: control 0%, treated 17% | PASS: we warn above 10% removed or a 20-point gap between groups |

*Conclusion:* removed S07, S11 (weak signal; near-copy of S10); kept S17, which is flagged only by its batch offset and is retained with the batch modelled later; groups and batches stay balanced
*Did not run:* replicate agreement (no technical replicates); control behaviour (no control samples)
*Carried caveats:* missing values are intensity-dependent
*Go back:* outputs/02_qc.csv

## Card 3

**Checkpoint 3: normalization and features: WARN, usable, with 1 caveat(s)**
*Asks:* are samples comparable and are features chosen reproducibly?

| What it checked | What we found | Result and why |
|---|---|---|
| whether two halves of the samples choose the same 100 most variable proteins (no group labels used) | overlap between halves is 0.14 | WARN: we warn below 0.5 |
| the normalisation fits the detected scale | linear intensities: log2 then median-centred per sample | PASS: intensities need a log transform and centring; counts methods were not used |
| whether sample medians are now similar | spread of sample medians is 0.00% | PASS: we warn above 10% |
| how many values were filled in | kept 630 of 800 proteins seen in 70% or more of samples; 4.9% of kept values were filled in | PASS: we warn above 20% imputed |

*Conclusion:* samples are comparable; a stable set of most-variable proteins does not exist here because most proteins vary by similar amounts, so all proteins that passed the missing-value rule are kept
*Did not run:* technical features among the selection (no contaminant list given); batch robustness of the selection
*Go back:* outputs/03_normalized.csv

## Card 4

**Checkpoint 4: structure discovery: WARN, usable, with 1 caveat(s)**
*Asks:* what are the main patterns and are they technical?

| What it checked | What we found | Result and why |
|---|---|---|
| how much of the structure above noise is explained by batch, and by group | batch explains 75% of it and group 20% | WARN: we warn when batch explains more than 20%, because it then has to be handled in the comparison |
| how many components carry more variation than shuffled data | 3 components exceed the shuffled-data level (together 32% of all variation) | PASS: we warn if fewer than 2 are above noise |

*Conclusion:* there is clear structure; batch accounts for a large share of it, so batch has to be a covariate in the comparison (it is)
*Did not run:* batch correction effect (none applied); grouping stability (no clustering run)
*Go back:* outputs/04_structure.csv

## Card 5

**Checkpoint 5: group comparison: PASS, good**
*Asks:* which proteins differ between control and treated?

| What it checked | What we found | Result and why |
|---|---|---|
| each group has several independent samples | 12 control and 10 treated samples after QC | PASS: we fail under 2 per group and warn under 3 |
| the group is not perfectly mixed up with a batch | every batch holds both groups | PASS: a batch holding one group only would be a WARN |
| proteins that differ between groups, replicate-aware, corrected for multiple testing | 33 proteins at adjusted p below 0.05 and at least 1.4-fold; 33 of 40 planted found, 0 of 20 null proteins called | PASS: linear model on log2 intensities with batch as covariate; Benjamini-Hochberg across all proteins |
| shuffled group labels should find almost nothing | median 0 proteins with shuffled labels against 34 with real labels | PASS: we warn above 5 proteins or 5% of the real count |
| called proteins move the same way in each batch | median 100% of batches agree | PASS: we warn below 75% |
| whether subjects appear in every group | each sample is a different subject: unpaired | INFO: recorded, not scored; batch is included as a covariate |

*Conclusion:* a clear and consistent treatment effect on the planted proteins; batch was modelled
*Did not run:* pathway enrichment (no protein sets given)
*Go back:* outputs/05_comparison.csv

## Card 6

**Checkpoint 6: interpretation and robustness: PASS, good**
*Asks:* what do the results mean and how far to trust them?

| What it checked | What we found | Result and why |
|---|---|---|
| written before the analysis: the 40 spiked proteins are higher in treated | 94% of them were called | PASS: met when at least 80% are found |
| 20 proteins with no spiked effect should not be called | 0 of 20 called | PASS: any called null protein is a FAIL |
| a 40-protein signature trained on two batches separates the groups in the third | area under the curve per held-out batch: 1.00, 1.00, 1.00 | PASS: we warn below 0.8 |
| the called proteins are found again using 80% of samples | median overlap with the full result 0.82 over 5 subsamples | PASS: robust when the overlap is 0.6 or more |

*Conclusion:* the treatment effect on the planted proteins is validated across batches and stable to subsampling
*Carried caveats:* intensity-dependent missing values; one weak sample and one duplicate removed
*Go back:* outputs/06_findings.csv

## Run summary

## Quality check

**Overall: Usable, with caveats: read the caveats before drawing conclusions.**

| checkpoint | quality | why | to go back to this point |
|---|---|---|---|
| 1. Intake and integrity | Usable, with 2 caveat(s) | 2 check(s) warned | outputs/intake.csv |
| 2. QC and filtering | Usable, with 2 caveat(s) | 2 check(s) warned | outputs/qc.csv |
| 3. Normalization and features | Usable, with 1 caveat(s) | 1 check(s) warned | outputs/normalize.csv |
| 4. Structure discovery | Usable, with 1 caveat(s) | 1 check(s) warned | outputs/structure.csv |
| 5. Group comparison | Good | every scored check passed | outputs/compare.csv |
| 6. Interpretation and robustness | Good | every scored check passed | outputs/interpret.csv |
