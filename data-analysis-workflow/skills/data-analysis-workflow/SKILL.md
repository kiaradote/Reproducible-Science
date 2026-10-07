---
name: data-analysis-workflow
description: "Master rule for any quantitative data analysis the user wants reproducible, restartable or presentable: omics (bulk or single-cell RNA-seq, proteomics, metabolomics, CRISPR screens), assay or plate readings, imaging measurements, surveys or any table of samples by measurements. Run the analysis as ordered checkpoints (intake and integrity, QC and filtering, normalization and features, structure discovery, group comparison, biological meaning, interpretation and robustness), each with its own skill. Every checkpoint saves a restart file, stops on FAIL, and reports in plain English: what each test checked, what was found, and the result and why, with a conclusion, a quality rating, what did not run and how to go back. Keeps a running summary and ends with a brief finished summary. Use a specialised master rule first when one exists (for example scrna-seq-workflow for single-cell RNA-seq)."
---

# General data-analysis workflow rule

Use this skill first for any analysis that should be reproducible and checkable, then load the skill of each checkpoint when you reach it. If a specialised master rule exists for the data type (single-cell RNA-seq: `scrna-seq-workflow`), follow that one; it extends this rule.
Do not improvise a different order and do not run a later checkpoint on the output of an earlier one that has not been reported.

## The checkpoints, in order

| # | Checkpoint | Question it answers | Skill |
|---|---|---|---|
| 1 | intake | Are these the intended files, complete, readable, consistent with what was expected, and what kind of numbers are in them? | `data-intake-integrity` |
| 2 | qc | Which samples or measurements are damaged, empty, outlying or duplicated, and are the filters sensible? | `data-qc-filtering` |
| 3 | normalize | Are values comparable between samples, and are the informative features chosen in a reproducible way? | `data-normalization-features` |
| 4 | structure | What are the main patterns (groups, gradients, batch effects), and are they real or technical? | `data-structure-discovery` |
| 5 | compare | What differs between groups or conditions (only when a design with replicates exists)? | `data-group-comparison` |
| 6 | biology | What do the changed or marker features mean biologically, and how well is that supported (enrichment, expectations, agreement with other studies)? | `omics-biological-interpretation` |
| 7 | interpret | Which conclusions survive changed analysis choices, held-out samples and controls? | `data-interpretation-robustness` |

Checkpoint 5 is skipped, and recorded as not applicable, when there is no group question or no replicates. Checkpoint 6 is skipped, and recorded as not applicable, when the data are not biological features with a usable annotation. Checkpoint 7 always comes last and carries every caveat forward.

## Rules

1. **Order.** Run the checkpoints in the order above. If a notebook or pipeline for this data type exists, run it; otherwise do the steps yourself and write the same files.
2. **One restart file per checkpoint.** After each checkpoint save its output in an open format (CSV, Parquet, h5ad, RDS), its report (`reports/<checkpoint>.md` and `.json`), its figures, and a manifest entry holding the parameters used, the input file hash and the output file hash (`da_hash_file`, `da_hash_params`, `da_is_current` in this skill's helpers). A checkpoint is rerun only when its parameters or its input changed; to go back, delete or invalidate it and everything after it.
3. **Stop on FAIL.** A FAIL stops the run. Show the reason and ask the user to fix the cause or to continue on purpose; continuing is recorded with a written reason. A WARN never stops the run but travels forward as a caveat. Never change a threshold just to turn a FAIL into a PASS; change it only when you can say why the data justify it, and record that.
4. **Thresholds are defaults, and they are written down.** Each checkpoint skill gives default thresholds. Use them unless the data type needs another value, say which value you used and why, and keep it in the manifest. When a threshold is a judgement call, say so.
5. **Report every checkpoint in plain English, in the chat, the moment it finishes.** No code and no function names. For each checkpoint write:
   - **Heading:** checkpoint number and name, verdict (PASS, WARN or FAIL) and quality rating.
   - **What this checkpoint asks:** one sentence.
   - **Tests** as a short table with three plain columns per test: *What it checked* (a definition in everyday words, written as "Whether ..." so every row reads the same way), *What we found* (the numbers), and *Result and why* (PASS, WARN, FAIL or INFO plus the reason in words: the rule or threshold and how the data compared with it). A result is never given without its reason. An INFO row says why it is not scored.
   - **Conclusion:** one to three sentences on how far the result can be trusted and what happens next, naming the biggest caveat.
   - **Did not run:** each check that did not run, with the reason. A check that did not run counts for nothing.
   - **Carried caveats** from earlier checkpoints that still apply.
   - **Go back:** the output file and how to rerun from here.
   In the chat show every WARN and FAIL and at least the five most informative PASS rows, and say how many other tests passed; the full list stays in the report file. `da_card` builds this text from the test records.
6. **Quality scale**, for every checkpoint and for the whole run:
   - *Good*: every scored check met its threshold.
   - *Usable, with N caveats*: at least one check warned; everything built on this checkpoint inherits the caveat.
   - *Accepted despite a failed check*: a FAIL was continued on purpose, with the reason on record.
   - *Not reliable until fixed*: a check failed; nothing downstream is trusted.
   - *Not applicable*: no scored check could run for this data (for example, no groups to compare).
   - *Not done*: the checkpoint has not run.
   The run is rated by its worst checkpoint. Thresholds are defaults for typical data, so Good means the checks passed, not that the conclusion is proven.
7. **Keep a running summary.** After every checkpoint rebuild `RUN_SUMMARY.md` with the quality table at the top (`da_run_summary`), then each checkpoint's report, parameters, software versions, input file hashes and how to resume.
8. **Be brief at the end.** When the run finishes, give a finished summary of at most one screen: the table of checkpoints with their quality rating and a one-line conclusion each, the overall rating, the caveats that matter most, and the paths of the summary and output files. Do not repeat the cards.
9. **Every flagged item gets a decision.** A sample, feature or group that a test flags is either removed or kept, and the card says which and why. A flag that is simply left unexplained is not allowed.
10. **Defaults and questions.** Choose sensible defaults and say what was chosen. Ask the user only for what they alone know: which columns identify the groups, replicates and batches, what biology is expected, and whether a heavy run should proceed when time or memory is limited. Never use the data's own published labels to build the thing they are later used to evaluate.
11. **Honest wording.** A group name or cluster is a label, not a proven category. A difference between groups is an effect of the treatment only if the design check says the groups are not mixed up with a batch. Pathway or category names describe shared members, not demonstrated activity. State every skipped sample, subsample or unrun check next to the number it affects.
12. **Large data.** Read in chunks, keep matrices sparse where most values are zero, record peak memory, run sensitivity analyses on a capped subsample and say so, and pause a heavy run when compute is limited instead of silently shrinking the analysis.

## Helpers (kernel.py)

`da_record` (one test with its plain-English definition, numbers and reason), `da_overall`, `da_rating`, `da_card` (the chat card), `da_write_report`, `da_run_summary` (quality table across checkpoints), `da_hash_file`, `da_hash_params`, `da_is_current`. They need only the Python standard library.

## Example checkpoint card

```
Checkpoint 2: QC and filtering. WARN, usable with 1 caveat
Asks: which measurements are damaged, empty or outlying, and are the filters sensible?

| What it checked | What we found | Result and why |
|---|---|---|
| Whether any sample has too few detected features to be trusted | 2 of 24 samples have under 40% of the median feature count | WARN: more than one sample is far below typical; we expect none |
| Whether any sample is a duplicate of another | no pair is more correlated than 0.999 | PASS: no near-identical samples |
| Whether the batch is balanced across groups | each batch holds both groups | PASS: a group effect cannot be mistaken for a batch effect |

Conclusion: usable, but the two weak samples were kept; later results should be checked with them removed.
Did not run: replicate-concordance (only one replicate per group).
Go back: outputs/02_qc.csv; rerun from qc.
```

## Finished summary template

```
Analysis: <dataset>, <n samples> samples, <n features> features; overall quality: <Good | Usable, with caveats | Not reliable | Incomplete>
1 Intake          <rating>  <one line>
2 QC              <rating>  <one line>
3 Normalization   <rating>  <one line>
4 Structure       <rating>  <one line>
5 Comparison      <rating or "not requested">  <one line>
6 Biology         <rating or "not applicable">  <one line>
7 Interpretation  <rating>  <one line>
Main caveats: <up to three>
Go back to any step: <output files>; details in RUN_SUMMARY.md
```
