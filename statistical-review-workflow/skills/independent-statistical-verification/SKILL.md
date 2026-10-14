---
name: independent-statistical-verification
description: "Stage 3 of a statistical claim-verification review: recompute reported statistics from the obtained data, check p-values, confidence intervals, degrees of freedom and effect sizes for internal consistency, run sensitivity analyses and documented simulations, and classify each reproduction as match, minor or major. Ships repro_compare and sim_record (documents null, assumptions, randomization, iterations, seed, uncertainty, limitations). Uses the helpers loaded by statistical-review-workflow."
---

# Independent statistical verification

Work only from files in the evidence manifest. For each claim with obtained evidence:

1. **Preserve inputs and code.** Write the recomputation as a script or cell that reads the manifest files, records versions and seeds, and saves its outputs (tables, figures). Another person must be able to rerun it.
2. **Recompute** the exact quantity the claim states, with the paper's own definition, threshold and subset. If the paper's definition is ambiguous, run each plausible reading and report all of them.
3. **Compare** with `repro_compare(reported, recomputed, decimals)`: match, minor, or major. Small differences from random seeds, software versions or ties are minor; state the likely cause only if tested.
4. **Internal consistency** from the text alone (Inferred from text, but computed): `recompute_p`, `ci_p_consistent`, `grim_check`, df against n.
5. **Sensitivity**: alternative model or statistic (rank-based, permutation), leave-one-unit-out, threshold sweep, outlier handling, analysis at the true replication unit, bootstrap resampling at that unit.
6. **Simulations and null models**: document each with `sim_record(...)`: null hypothesis, assumptions, randomization scheme, iterations, seed, statistic, result, interval, limitations. A null-model result says what chance alone produces under that null. Do not read it as the proportion of genuinely affected units without a separate justification; the record carries this guard.
7. **Conclusion per claim**: Verified (agrees), Contradicted (disagrees, discrepancy stated), Not assessable (inputs missing or method not applicable). A disagreement is Contradicted only after you have checked that you used the same subset, threshold, labels and preprocessing as the paper; otherwise report it as an unresolved discrepancy and list the readings tried.
8. Update the ledger rows (`claim_record`) with the reproduction description and the conclusion.

Report what could not be run and why. Never describe an unperformed check as verified.
