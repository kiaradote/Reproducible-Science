---
name: adversarial-review-synthesis
description: "Stage 5 of a statistical claim-verification review: challenge the draft findings before they are reported. Looks for unsupported accusations, wrong or mismatched calculations, alternative explanations, concerns that do not affect the central conclusions, and language that overstates certainty; then produces a calibrated verdict (supported, weakened, contradicted, unresolved) per central claim. Can be run by an independent sub-agent that receives only the ledgers, manifest and artifacts."
---

# Adversarial review and synthesis

Run after stages 1 to 4. The aim is to find what is wrong with the review, not with the paper.

## Procedure
1. Assemble the claim ledger, concern ledger, evidence manifest, scripts and outputs.
2. **Independent challenge pass.** If sub-agents are available, delegate this pass to a fresh agent that sees only those files (not your reasoning) and ask it to attack each concern. Otherwise do it yourself as a separate pass using the checklist below. Record which was done.
3. Checklist per concern:
   - Is every number in the concern traceable to a script and input in the manifest?
   - Does the recomputation use the paper's own subset, labels, thresholds and preprocessing?
   - Does the verification actually test the concern, or only a related quantity?
   - Is a plausible concern worded as a proven error? Downgrade the wording or the result field.
   - Is there an equally good alternative explanation? State it.
   - Does the concern, if true, change a central conclusion? If not, move it to minor.
   - Is a simulation or null-model result being read as a proportion of affected units without justification?
   - Is any "Verified" row built on data that were not obtained, or on a substituted source that is not declared?
   - Are arithmetic, methodology, interpretation, reporting and generalization issues kept apart?
4. Resolve each challenge: concern upheld, downgraded, reworded, or withdrawn (with the reason kept in the ledger).
5. **Calibrated verdict.** For each central claim: Supported, Weakened, Contradicted, or Unresolved, each with the evidence tags and remaining uncertainty. Do not issue an overall verdict on the paper's scientific integrity from statistical limitations alone.
6. Report counts: concerns raised, upheld, downgraded, withdrawn.
