# Statistical review: Bacsik et al. 2023, eLife 12:RP86852 (DOI 10.7554/eLife.86852)
Reviewer: STAT_REVIEWER | Review date: 2026-10-13 | Data level: analysis-ready per-cell tables from the authors' GitHub repo (jbloomlab/barcoded_flu_pdmH1N1, `main` branch, files all_samples.csv sha256 prefix 6c81c7ac973fd236 and all_samples_complete_measurements_cells_data.csv prefix f32596cab16a51a7; downloaded 2026-10-13). Raw FASTQ (GEO GSE214938), the Snakemake pipeline and replicate-level barcode tables were NOT re-processed.

Study type: single-cell RNA-seq plus barcoded-virus progeny sequencing, in vitro (MDCK-SIAT1-TMPRSS2 cells); one experiment at low MOI (~0.15) and one at higher MOI (~0.6). Closest guideline: MDAR (the paper includes an MDAR checklist); no omics-specific statistical guideline applies.

## 1. Verdict
The central finding, that per-cell viral transcription is a poor predictor of per-cell progeny, is supported by the data as far as I could recompute it. Every number I tested reproduced from the public table: Pearson r of -0.14 and 0.39 (my values -0.139 and 0.387), the Gini coefficients, the 6-cell and 23-cell concentration figures, and both NS-gene permutation tests. The main reservations are about inference, not arithmetic. Each MOI is a single experiment, so the effective replication for any claim of generality is n = 1. The analysed cells (91 of 254 and 290 of 357) are a selected subset. Pearson r is a fragile summary for these heavy-tailed, zero-inflated data. A few secondary statements are slightly stronger than the data allow.

## 2. Claims table
| id | claim (paper) | estimand | true replication unit | rating | evidence |
|---|---|---|---|---|---|
| C1 | Transcription and progeny production are poorly correlated in single cells | cell-level correlation, viral mRNA fraction vs progeny barcode frequency | cells within 1 experiment per MOI | Supported (association weak/absent) with caveats M1, M2, M3 | Verified |
| C2 | Progeny production is more heterogeneous than transcription | Gini coefficients | cells | Supported | Verified |
| C3 | Highest-transcribing cells are often NS-missing aberrant infections | difference in mean viral fraction, NS-missing vs all-genes | 13 (low) and 15 (high) NS-missing cells | Supported but small n; direction firm at low MOI, weaker at high MOI | Verified |
| C4 | No trend of progeny with transcription even at high MOI | correlation | 290 cells | Partly: weak positive association present (Pearson 0.165, p=0.005) | Verified |
| C5 | Physical and infectious progeny are much more correlated than transcription and progeny | Pearson 0.39 vs -0.14 | 91 cells | Supported in direction, understated: Spearman is 0.87 | Verified |
| C6 | Wide range of progeny among cells infected by complete unmutated virions | descriptive | cells | Supported descriptively; not separated from measurement noise (see M3) | Inferred from text + partial recompute |

## 3. Major concerns
### M1. Replication unit: one experiment per MOI (Inferred from text)
All inference is across cells inside a single infection and single 10x run per MOI. Cell-level CIs and p-values do not capture run-to-run (batch, ambient RNA, virus-stock) variation. The paper itself notes the two MOIs give consistent broad trends, which is useful but is two non-identical conditions, not replicates (infectious progeny exist only at low MOI; thresholds for calling infection differ by MOI). Resolution: a biological replicate of at least the low-MOI experiment, or explicit limitation language.

### M2. Selection of the analysed cells (Verified)
Progeny analysis is restricted to cells expressing both barcoded genes (HA and NA) and, at low MOI, with complete sequence. That is 91/254 (36%) low MOI and 290/357 (81%) high MOI. The retained cells are slightly more transcriptionally active (median viral UMI fraction 0.061 vs 0.040 low MOI; 0.252 vs 0.227 high) and more often have all eight genes (73.6% vs 57.1%; 84.5% vs 70.9%). Cells with absent HA or NA, the aberrant infections the paper highlights, are therefore under-represented, and the NS-missing group has only 13 and 15 cells. This probably biases toward a stronger transcription-progeny link, so the headline is conservative, but the NS comparison and the "range among complete unmutated virions" claim rest on subsets. No sensitivity analysis to the exclusions is reported.

### M3. Pearson r on heavy-tailed, zero-inflated, floor-censored data (Verified)
Progeny frequency is highly skewed (Gini 0.78 to 0.88) with 23 to 55 percent of cells at or below the limit of detection depending on the readout (low-MOI physical 40%, low-MOI infectious 55%, high-MOI physical 23%). Pearson r is dominated by a handful of cells. Evidence from my recomputation: physical vs infectious progeny Pearson 0.39 but Spearman 0.87; transcription vs infectious progeny Pearson -0.14 (bootstrap 95% CI -0.22 to -0.07, which sits oddly against p=0.19, showing instability) but Spearman -0.28 (CI -0.47 to -0.06); at high MOI Pearson 0.165, Spearman 0.128 (CI -0.005 to 0.256); dropping the top 5 percent of progeny producers moves low-MOI physical r from -0.006 to -0.19. The qualitative conclusion (weak, not strongly positive) is robust across all variants, but the reported values depend on the statistic chosen, and a rank-based or log-scale/censored-aware analysis should be primary. The low-MOI relation is, if anything, modestly negative in rank terms.

### M4. Measurement noise and attenuation not quantified for the correlation (Verified, partial)
Progeny frequency is a bottlenecked sequencing measurement (the authors show replicate-frequency scatter and set a detection limit at 1e-5, Inferred from text). Using the two independent barcode readouts per cell as parallel measures, HA vs NA supernatant frequency has Spearman 0.60 (low MOI) and 0.43 (high MOI). Spearman-Brown reliability of the averaged measure is about 0.75 and 0.61. Dividing the observed rho by the square root of reliability gives an upper bound on the transcription-progeny |rho| of about 0.26 (low) and 0.16 (high) if transcription were error-free. So progeny noise does not explain the poor correlation, which strengthens C1. Caveat: HA/NA discordance may include real biology (segment reassortment, barcode swaps; authors estimate 4-8 percent of cells could be affected), so this reliability is a lower bound. Reliability of the transcription measure (UMI fraction per cell, 10x dropout) was not assessable, so full disattenuation is not possible.

## 4. Minor concerns
1. Statement "no trend at high MOI" (C4) is stronger than the data: Pearson 0.165 (p=0.005, 95% CI 0.05-0.27), though rank correlation is weaker (0.128) and ~0.10 after excluding top 5% producers.
2. Infected-cell calls use different virus-fraction thresholds per MOI (1% vs 2.5%) chosen to handle ambient RNA; no sensitivity sweep reported, and thresholds' prespecification is not stated. Cells carrying both virus tags (51 low MOI, 17 high MOI) are dropped from the 254/357 counts without sensitivity analysis.
3. NS-missing vs all-genes tests: the paper reports p=0.0002 and 0.004 (permutation, 5000 draws); mine are 0.0004 and 0.0036 (same test, different seed). The Mann-Whitney p at high MOI is 0.047, so the high-MOI result is sensitive to the test chosen. Effect sizes and CIs are not reported with these p-values; n=13 and 15.
4. The paper reports p-values as p<0.01 in places (Fig 5C legend, Methods) rather than exact values or effect sizes with CIs; the correlation coefficients carry no intervals.
5. Progeny frequencies are compositional (shares of one pooled library). No dependence adjustment is discussed, though with hundreds of contributing cells the induced dependence is small.
6. Total recovered cells: the table holds 3129 unique low-MOI cells and 1243 high-MOI cells, while the text gives ~1600 and ~800 (used to compute empirical MOI). Infected counts reconcile exactly (254 and 357 wt-tagged), so the discrepancy is probably cells carrying the second virus tag, but it is not stated in the paper and I did not resolve it.
7. Gini coefficients are reported without uncertainty (bootstrap would be easy), and the Gini of viral_UMIs (0.38, 0.40) differs from that of the viral fraction (0.46, 0.40) at low MOI; the paper's values match the latter.

## 5. Reproduction results (all Verified from the analysis-ready table)
| quantity | reported | recomputed | match |
|---|---|---|---|
| analysed cells low / high MOI | 91 / 290 | 91 / 290 | yes |
| wt-tagged infected cells low / high | 254 / 357 | 254 / 357 (from all_samples.csv) | yes |
| Pearson r transcription vs infectious progeny, low | -0.14 | -0.139 | yes |
| Pearson r physical vs infectious progeny, low | 0.39 | 0.387 | yes |
| Gini physical / infectious, low | 0.78 / 0.88 | 0.780 / 0.875 | yes |
| Gini physical, high | 0.66 | 0.661 | yes |
| Gini transcription low / high | 0.46 / 0.40 | 0.464 / 0.403 | yes |
| cells for half of physical progeny / transcripts, low | 6 / 23 | 6 / 23 | yes |
| NS-missing vs all-genes permutation p, low / high | 0.0002 / 0.004 | 0.0004 / 0.0036 | yes (seed-level) |
Full correlation table with bootstrap CIs and subsets: correlation_reproduction.csv.

## 6. Sensitivity analyses
See correlation_reproduction.csv. Subsets: all-8-genes cells only, missing-gene cells only, detected-only (above LOD), top 5 percent of progeny removed, Pearson vs Spearman. The conclusion "weak or absent positive association" holds in every variant. Not run: threshold sweep for infection calling, leave-one-batch-out (no batches), censored (Tobit-type) model, reprocessing of raw reads.

## 7. Reporting-guideline benchmark (MDAR-style, statistical items)
Present: data and code availability (GEO, GitHub with Snakemake pipeline, interactive plots), exact software versions, permutation test details, technical replicate plots. Missing or weak: biological replicates and n per claim; justification of sample size; effect sizes with CIs; exact p-values; prespecification status; sensitivity to exclusion criteria; uncertainty on Gini coefficients and correlations.

## 8. What was done well
Complete public data and a reproducible pipeline; every reported number I checked reproduced; the authors describe technical replicates, set an empirical detection limit, address barcode swaps and ambient RNA, use a permutation test rather than a parametric test, replicate the main finding in a second MOI condition, and report cells without progeny rather than hiding them.

## 9. Requests to authors
1. A biological replicate of at least the low-MOI experiment, or explicit limitation text.
2. Report Spearman (or censored-model) correlations with CIs alongside Pearson.
3. Sensitivity of main results to the infection threshold and to the HA/NA expression and sequencing requirements; describe the excluded cells.
4. Effect sizes with CIs for the NS-missing comparison and a note that n=13/15.
5. Soften "no trend" at high MOI.

## 10. Limitations of this review and reproducibility log
- Not assessed: raw reads (GEO GSE214938), replicate-level barcode frequencies, long-read genome calling, the single-cell clustering and ambient-RNA correction, per-cell transcription measurement reliability, the figure supplements and the claims about specific mutations.
- Data-level substitution: analysis-ready repo tables used, not raw data; repo `main` branch may differ from the archived copy cited by the paper (counts reconcile, so this seems unlikely to matter).
- Full text read from PMC (PMC10484525) text extraction; figures could not be downloaded, so figure-based claims were checked against data, not images.
- Environment: Python 3 with numpy, pandas, scipy; bootstrap B=4000, seed 0; permutation 5000, seed 1. Helpers from the statistical-review-workflow skill: permutation_diff.
- The authors' interactive plots and pipeline were not run.
