---
name: statistical-review-workflow
description: "Orchestrator for statistical claim verification of a paper and its public data. Runs ordered stages: inventory, claim extraction (six-part claim ledger), data retrieval and reproducibility gate (evidence manifest), independent statistical verification, methodology audit (design gate, falsification of each concern), adversarial review, and a calibrated report. Distinguishes demonstrated errors from plausible concerns; verdicts are Verified, Contradicted, Supported concern, Not assessable. Ships helpers: recompute_p, grim_check, cluster_check, multiplicity_table, smd_ci, permutation_diff, p_near_threshold, ci_p_consistent."
---

# Statistical claim-verification workflow (orchestrator)

Domain-agnostic. Never carry conclusions from one paper to another. Stage skills: `paper-claim-extraction`, `data-retrieval-reproducibility`, `independent-statistical-verification`, `statistical-methodology-audit`, `adversarial-review-synthesis`. Supporting: `data-intake-integrity`, `data-group-comparison`, `data-analysis-workflow`, `literature-review`, `mcp-pubmed`, `figure-rubric`. Load each stage skill when you reach it (loading also installs its helpers).

## Quality rules (apply at every stage)
1. A plausible methodological concern is not a demonstrated statistical error. Use "error" only for demonstrated, Verified findings.
2. Every major criticism names the claim, the source evidence, the statistical reasoning, and a verification procedure, and records an attempt to falsify it.
3. Reproduce key numerical claims whenever data and methods are available; preserve code, inputs and outputs so the analysis can be repeated. Never claim a reproduction using data that were not obtained.
4. Before criticising a method, establish the study design, observational unit, independent sampling unit and the method's assumptions (design gate).
5. Separate arithmetic, methodology, interpretation, reporting and generalization issues; state each one's consequence separately.
6. State missing data, inaccessible files and unverified analyses explicitly. Never describe an unperformed check as verified.
7. Consider alternative explanations; try to falsify each proposed criticism before accepting it.
8. Simulations and sensitivity analyses: document null, assumptions, randomization scheme, iterations, seed, uncertainty, limitations. A null-model result is not an estimate of the proportion of genuinely affected units without justification.
9. Prioritise concerns by confidence and by impact on the main conclusions.
10. Finish with a calibrated verdict per central claim (supported, weakened, contradicted, unresolved). No overall judgement of scientific integrity from statistical limitations alone.
11. Environment: run on the existing environment. Do not create environments or install packages; if a package is missing, use the standard library or ask the user first.

## Evidence tags and verdicts
Findings: **Verified** (recomputed, agrees), **Contradicted** (recomputed, disagrees; discrepancy stated), **Supported concern** (reasoned, partly tested, not a demonstrated error), **Not assessable** (evidence missing; say what). Text-only reading is "Inferred from text" and can never support Verified or Contradicted. Declare the data level reproduced from (raw, processed, summary) and any substitution.

## Stages
0. **Inventory.** Title, DOI, version read, study type, matching reporting guideline (references/reporting_guidelines.md), links to data, code and supplements; check corrections or retractions.
1. **Claim extraction** (`paper-claim-extraction`): claim ledger with location, required evidence (data, labels, formula, threshold).
2. **Data retrieval and reproducibility gate** (`data-retrieval-reproducibility`): data availability statement, publisher supplement links, accession inspection, file roles, label reconciliation, evidence manifest with source/size/checksum/parse status, alternative sources before declaring unavailable. Output: evidence-to-claim map. Stop and report if no claim has obtainable evidence.
3. **Independent verification** (`independent-statistical-verification`): recompute, compare, consistency checks, sensitivity, documented simulations.
4. **Methodology audit** (`statistical-methodology-audit`): design gate, then concerns with falsification attempts, typed by error kind, rated by confidence and impact.
4b. **Bias, missingness and alignment screen** (`bias-missingness-alignment-screen`): attrition and missing-value quantification, selection and bias direction, forking-paths and p-hacking risk (flags only), and a per-claim comparison of conclusions across abstract, results, discussion, figures, supplement, versions, peer review and corrections.
5. **Adversarial review and synthesis** (`adversarial-review-synthesis`): challenge pass (independent sub-agent when available), downgrade or withdraw weak concerns, calibrated verdicts.
6. **Report** (references/report_template.md): verdict summary; claim ledger (six-part); evidence manifest; reproduction table; concern ledger ordered by priority with withdrawn concerns listed; sensitivity and simulation records; reporting-guideline gaps; what was done well; requests to authors; a dedicated section 'Limitations of this review and items for the reader to oversee' (intent, unreported analyses, domain judgment, evidence not obtained, data level, cross-publication evidence, simulation meaning, self-pass vs independent pass, automation caveat) followed by the reproducibility log (versions, accessions, checksums, environment, anything not assessed). Fill every bullet for the paper at hand.

Save: report (.md), claims_table.csv, evidence_manifest.csv, concern ledger CSV, scripts and result tables as artifacts. Be even-handed and report what the authors did well.
