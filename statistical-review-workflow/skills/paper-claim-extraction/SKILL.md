---
name: paper-claim-extraction
description: "Stage 1 of a statistical claim-verification review: extract every statistical claim from a paper (estimate, n, test, p-value, CI, figure, conclusion) with its exact location, and build a claim ledger with the six-part structure author claim, required evidence, independent calculation, methodological assessment, conclusion, remaining uncertainty. Ships claim_record and ledger_to_csv, which refuse a Verified or Contradicted verdict without a performed reproduction."
---

# Paper and claim extraction

Output: a **claim ledger**, one row per claim, built with `claim_record(...)`.

## Steps
1. Read the full text (fetch_article_fulltext, publisher HTML or PDF). Note which version you read (preprint, accepted, published) and its DOI.
2. List headline claims first (abstract, title, last paragraph of the introduction, discussion), then supporting claims in results, figure legends and tables.
3. For each claim record: a short verbatim quote (under 20 words) and its location (section, page, figure panel, table row). Quote numbers exactly as printed: estimate, n, test, p, CI, df.
4. Classify each claim: descriptive, associative, causal, predictive, or methodological ("method X performs well").
5. Fill **required_evidence** for each claim BEFORE looking at data: which data file, which sample labels or group assignments, which formula or model, which threshold or cutoff. If the paper does not state one of these, write "not stated in paper"; that gap is itself a reporting finding.
6. Leave reproduction, method assessment and conclusion at their defaults ("not attempted", "not assessed", "Not assessable") until later stages fill them.
7. Count tests reported anywhere (main text, legends, supplements) for later multiplicity assessment.

## The six-part ledger row
1. Author's claim: what exactly is reported and where.
2. Required evidence: data, labels, formula, threshold needed to verify.
3. Independent calculation: filled in stage 3 (code path, inputs, outputs).
4. Methodological assessment: filled in stage 4, after the design gate.
5. Conclusion: Verified, Contradicted, Supported concern, or Not assessable.
6. Remaining uncertainty: what missing evidence could change the conclusion.

## Verdict meanings (use only these)
- **Verified**: recomputed from obtained data and agrees at reported precision (or within a stated tolerance).
- **Contradicted**: recomputed from obtained data and disagrees, with the size of the discrepancy stated.
- **Supported concern**: a reasoned concern, partly tested, that does not demonstrate an error.
- **Not assessable**: the required evidence was not obtained or the method cannot be applied; say what is missing.

Save the ledger as `claims_table.csv` via `ledger_to_csv`.
