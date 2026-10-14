---
name: bias-missingness-alignment-screen
description: "Stage 4b of a statistical claim-verification review: screen for bias, missing values, selective reporting / p-hacking risk, and misalignment between conclusions and the data and tests actually run. Quantifies missingness and selection (analysed vs eligible), maps forking paths and runs specification curves, screens reported p-values, and compares each conclusion across abstract, results, discussion, figures, supplement, preprint/published versions, peer review and corrections. Findings are risk flags unless demonstrated; p-hacking cannot be demonstrated from published results alone. Ships missingness_report, funnel_compare, p_curve_summary, spec_curve, alignment_row."
---

# Bias, missingness and conclusion-alignment screen

Run after the design gate (statistical-methodology-audit). Record every finding with `concern_record` (falsification attempt required). Everything here is a **risk flag or a quantified sensitivity**, not an accusation. In particular, p-hacking, selective reporting and bias cannot be demonstrated from published numbers alone; the strongest honest statement is "the published record is consistent / not consistent with X, and here is what would settle it".

## A. Missing values and attrition
1. Find every place units or values drop out: eligible -> screened -> analysed (flow diagram, Methods, legends). Build the funnel with counts at each step.
2. Quantify missingness in the obtained data with `missingness_report(df, by=group_col)`. A group-dependent pattern means the data are not missing completely at random; a null test never proves randomness.
3. State the handling used by the authors (complete-case, imputation, censoring at a detection limit, exclusion) and whether the claim changes under another reasonable handling (sensitivity).
4. Values at a floor or ceiling (detection limit, saturation) are censored, not missing: check whether the statistic copes with them.

## B. Selection and other bias
1. `funnel_compare(eligible, analysed, cols)` on covariates and on the outcome-related variables: which direction would the selection move the headline estimate? Say whether the bias would strengthen or weaken the conclusion.
2. Consider by type and keep them separate: selection (who/what enters), measurement (instrument, proxy, differential error), confounding (design gate), survivorship, attrition, time-varying or batch effects, outcome-dependent filtering, circularity (selecting on the outcome then testing it).
3. Label each: demonstrated (shown in data), supported concern, or Inferred from text.

## C. Selective reporting and p-hacking risk (forking paths)
1. List the analytic choices the authors made: thresholds, exclusions, covariates, transformations, subgroups, model family, outcome definition, multiple outcomes. Mark which are justified in the paper, which look data-dependent, and which are prespecified (registration, protocol, preprint v1).
2. Count the tests that could have been run versus reported (supplements, other outcomes, other cutoffs); apply `multiplicity_table` to what is reported.
3. Run `spec_curve(df, estimate_fn, choices)` over the reasonable alternatives. Report the range of estimates and the fraction that keep the paper's sign and conclusion, next to the paper's own specification. A conclusion that survives most specifications is robust; one that holds only at the chosen cutoff is fragile (still not proof of intent).
4. `p_curve_summary(pvals)` only if there are about 10+ independent significant p-values on one question; otherwise write Not assessable. A pile-up just under 0.05 is a flag only.
5. Exact p-values versus "p<0.05", rounding, one-sided tests, and unreported non-significant results are reporting findings.
6. Never write "p-hacked". Write what was checked, what was found, and what would resolve it (registration, full set of analyses, raw data).

## D. Conclusion-to-data alignment across everything published
Build one `alignment_row` per central claim, comparing the wording and numbers in: abstract, title, results text, figure/table (as printed and as recomputed), discussion, supplement, and, where they exist, the preprint versus published version, peer-review reports and author responses, corrections or retractions, and the authors' code/README.
- Does the data (recomputed) support the strength of wording (descriptive, associative, causal, predictive, generalised)?
- Do the tests that were run actually address the stated question? A non-significant test is not evidence of no effect; a correlation is not a mechanism; a subset result is not a population result.
- Are numbers consistent across locations (n, effect, p, CI)? Name any disagreement.
- Does a later section walk back or strengthen a claim made earlier?
- Has anything changed between versions that bears on the claim?
Rate each as aligned, overstated, understated, inconsistent or not assessable; overstated and inconsistent must name the locations that disagree.

## Output
- Funnel and missingness tables, the specification-curve table, the alignment table (save as CSV artifacts).
- Concerns in the concern ledger with type (methodology / interpretation / reporting / generalization), confidence, impact, falsification attempt.
- A short paragraph for the verdict: where bias could move the conclusions, in which direction, and whether the conclusions align with the data and the tests that were run.
