# Statistical review: <paper title> (<DOI>)
Reviewer: STAT_REVIEWER | Date: <date> | Version of paper read: <preprint/published> | Data level reproduced from: <raw/processed/summary>

## 1. Calibrated verdict
Per central claim: Supported / Weakened / Contradicted / Unresolved, with one sentence of reason and the evidence tag. No statement about overall scientific integrity.

## 2. Claim ledger (six-part)
| id | 1 author claim + location | 2 required evidence (data, labels, formula, threshold) | 3 independent calculation (script, inputs, outputs) | 4 methodological assessment | 5 conclusion (Verified / Contradicted / Supported concern / Not assessable) | 6 remaining uncertainty |

## 3. Evidence manifest and data retrieval
Table from evidence_manifest.csv: file, source, size, sha256, retrieved, parsed. Label reconciliation results. Sources tried before declaring anything unavailable. Evidence-to-claim map.

## 4. Reproduction results
| claim id | reported | recomputed | class (match/minor/major) | notes |

## 5. Design gate
Observational unit, independent sampling unit, replication unit, assumptions, where established.

## 6. Concern ledger (ordered by priority = confidence x impact)
### <id> <title> [error type: arithmetic | methodology | interpretation | reporting | generalization]
- Claims affected / Source evidence / Statistical reasoning / Verification procedure and result / Falsification attempt / Alternative explanations / Result (demonstrated | supported concern | not assessable) / Confidence / Impact / Evidence tag
### Withdrawn or downgraded concerns (with reasons)

## 6b. Bias, missingness and conclusion alignment
Funnel (eligible -> analysed) and missingness tables; direction in which selection or bias could move each headline estimate; forking-paths map and specification-curve result; selective-reporting risk flags (flags only, never 'p-hacked'); alignment table comparing each central conclusion across abstract, results, discussion, figures, supplement, versions, peer review and corrections.

## 7. Sensitivity and simulation records
For each: null hypothesis, assumptions, randomization, iterations, seed, statistic, result, interval, limitations, and what the result does NOT show.

## 8. Reporting-guideline gaps

## 9. What was done well

## 10. Requests to authors

## 11. Limitations of this review and items for the reader to oversee
State plainly that this is an automated statistical review and that the reader should independently judge the items below. Fill every bullet for the paper at hand (write "none identified" only after checking), and add paper-specific items.
- **Intent and integrity:** the review cannot establish why an analysis was chosen. Risk flags about selective reporting or p-hacking are not findings of misconduct and make no statement about the authors' integrity. Reader to oversee: whether any flag warrants asking the authors or the editor.
- **Unreported analyses:** only published, deposited or registered analyses were visible. Analyses that were run but not reported are invisible except through registrations, protocols or earlier versions. Reader to oversee: whether a registration or protocol exists and was consulted.
- **Domain judgment:** whether thresholds, exclusions, controls, assays or models are scientifically reasonable was assessed from the paper and general statistical principles, not by a subject expert. Reader to oversee: field-specific validity of each choice named in the concern ledger.
- **Evidence not obtained:** list each file, source or analysis that could not be retrieved, parsed or run, and which claims depend on it (copy from the evidence-to-claim map). Those claims are Not assessable, not confirmed.
- **Data level and substitutions:** state the level reproduced from (raw, processed, summary) and any substituted source; anything upstream of that level (raw reads, preprocessing, alignment, image analysis, instrument output) was not re-checked.
- **Cross-publication evidence:** related papers, replication attempts, corrections and peer-review files were checked only if listed in the manifest. Reader to oversee: anything else published on the same question.
- **Simulation and sensitivity results:** they describe what chance or alternative choices produce under stated assumptions, not the true share of affected units or the true effect. Reader to oversee: whether the assumptions fit the study.
- **Reviewer self-check:** state whether the adversarial pass was an independent sub-agent or a self-pass; a self-pass is weaker.
- **Automation caveat:** calculations were scripted and are in the saved artifacts; reasoning in the ledgers is judgment and should be spot-checked, starting with the highest-priority concerns and every Verified or Contradicted verdict.

### Reproducibility log
Versions of the paper and files, accessions, checksums, environment (no environments created, no packages installed), seeds, and anything not assessed.
