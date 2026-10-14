---
name: statistical-methodology-audit
description: "Stage 4 of a statistical claim-verification review: assess whether the study design and statistical methods support the claims. Requires a design gate first (observational unit, independent sampling unit, replication unit, assumptions), then assesses independence, pseudoreplication, confounding, selection bias, multiplicity, missing data, model assumptions, circular analysis and interpretation. Every concern must be recorded with concern_record, which demands a falsification attempt and refuses to label a concern a demonstrated error without a performed check. Ships design_gate and concern_record."
---

# Statistical methodology audit

A plausible methodological concern is not a demonstrated statistical error. This stage turns candidate concerns into calibrated findings.

## 1. Design gate (before any criticism)
Call `design_gate(observational_unit, independent_sampling_unit, replication_unit, assumptions, source)`. State what the paper's design is, where you read it (section or figure), and what the method assumes. If the design cannot be established from the paper and data, write that, and treat every design-dependent concern as Not assessable or Inferred from text.

## 2. Candidate concerns, by type
Keep error types separate because their consequences differ:
- **arithmetic**: a reported number does not follow from the data or formula.
- **methodology**: the method is invalid or inefficient for the design (pseudoreplication, ignored clustering or pairing, wrong model family, circular analysis, leakage, multiplicity without control, unmodelled missing data, confounding with batch).
- **interpretation**: the conclusion goes beyond the estimate (causal wording, equating non-significance with no effect, practical significance).
- **reporting**: information needed to evaluate is missing (n, test, exact p, CI, thresholds, exclusions, software).
- **generalization**: claims about populations, conditions or datasets beyond what was sampled.

Checks to consider: independence and unit of replication, randomization/blinding, sample-size rationale, prespecification, selection and exclusion steps (compare analysed vs eligible units), batch confounding, assumptions and diagnostics, correlation-specific pitfalls (outliers, floor or censored values, heteroscedasticity, noise attenuation), multiplicity, researcher degrees of freedom, effect sizes with intervals, reporting guideline items (references/reporting_guidelines.md in statistical-review-workflow).

## 3. For each candidate concern
1. Name the claim(s) affected (claim_ids), and quote the source evidence.
2. Write the statistical reasoning: why the design plus method could produce a problem.
3. State a verification procedure and run it when the data are in the manifest (a recomputation, a sensitivity analysis, a documented simulation).
4. **Attempt to falsify it**: does a different analysis, subset, tolerance or reading remove the concern? Is there an alternative explanation that fits equally well? Does the concern change the paper's conclusion at all?
5. Record with `concern_record(...)`. Result is one of: demonstrated (requires Verified evidence), supported concern, withdrawn, not assessable. Withdrawn concerns stay in the ledger with the reason.
6. Rate confidence and impact separately (low/medium/high). Priority = confidence x impact. Order the report by priority.

## 4. Language rules
- Use "could", "may", "would need X to resolve" for supported concerns; use "error" only for demonstrated ones.
- Do not infer intent, misconduct or overall scientific integrity from statistical limitations.
- A concern that does not affect any central conclusion is minor even if it is real.
