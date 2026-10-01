---
name: scrna-seq-workflow
description: "Master rule for any single-cell RNA-seq (scRNA-seq) analysis request: run the analysis as an ordered series of checkpoints (fetch and integrity, QC, normalization and variable genes, dimensionality reduction, clustering, cell type annotation, optional condition comparison, interpretation and validation), each with its own skill, in the order of the reusable notebook. Use whenever the user asks to analyse, process, cluster, annotate or compare scRNA-seq / single-cell data, or asks for a reproducible, restartable or presentable single-cell workflow. Requires every checkpoint to save a restartable file, stop on FAIL, and report in a fixed plain-language card: tests run and what each tests, the conclusion, a quality rating, what did not run, and how to go back; and requires a running RUN_SUMMARY.md plus a brief finished summary."
---

# scRNA-seq workflow rule

Use this skill first for any scRNA-seq analysis, then load the skill of each checkpoint when you reach it. Do not improvise a different order, and do not run a later checkpoint on the output of an earlier one that has not been reported.

## The checkpoints, in order

| # | Stage (notebook name) | Question it answers | Skill to load |
|---|---|---|---|
| 0-1 | fetch, assemble | Are these the intended files, complete, and do matrix, genes and cell metadata line up? | `scrna-data-integrity` |
| 2 | qc | Which cells are damaged or empty, and are the filters sensible? | `scrna-qc-filtering` |
| 3 | features | Is expression normalised and are the informative genes chosen reproducibly? | `scrna-normalization-features` |
| 4 | dimred | Is the main structure captured (PCs against a shuffled null, batch effect checked) without technical drivers? | `scrna-dimensionality-reduction` |
| 5 | cluster | Are the groups of cells stable and distinct? | `scrna-clustering` |
| 6 | annotate | Which known cell types do the clusters resemble, and how sure are we? | `scrna-cell-type-annotation` |
| 7 | compare | What differs between two conditions (only when a condition and replicate samples exist)? | `scrna-condition-comparison` |
| 8 | interpret | What do the results mean and which conclusions survive changed analysis choices? | `scrna-interpretation-validation` |

Checkpoint 7 is skipped, and recorded as skipped, when there is no condition to compare. Checkpoint 8 always comes last and uses checkpoint 7's results when they exist.

## Rules

1. **Order.** Run the checkpoints in the order above. If the reusable notebook (`scrna_checkpoints.ipynb` with the `scrna_workflow/` folder) is available, run it: it already executes the stages in this order, skips stages whose parameters and inputs have not changed, and writes the files below. Otherwise follow the skills in order and produce the same files by hand.
2. **One restartable file per checkpoint.** After each checkpoint save its result (`runs/<dataset>/checkpoints/NN_<stage>.h5ad`, `.json` for fetch), its report (`reports/<stage>.md` and `.json`) and its figures (`figures/<stage>/`, with `figure_captions.json` of one-line captions and `figure_legends.md` of plain-language legends). A later stage reads only the previous stage's file, so any stage can be restarted by invalidating it and rerunning.
3. **Stop on FAIL.** A FAIL stops the run. Show the reason and ask the user to fix the cause or to continue on purpose; continuing is recorded with a written reason (`accept`). A WARN never stops the run but is carried forward as a caveat. Never lower a threshold to turn a FAIL into a PASS.
4. **Report every checkpoint as a card** the moment it finishes, in this order and in plain language (a student must be able to read it without opening code):
   - **Verdict and quality rating** (see the scale below), with the reason in one sentence.
   - **Tests run:** a table of status, test, what it tests (one phrase), and the result with its numbers.
   - **Conclusion:** one or two sentences saying what the checkpoint found and what it means for the next step.
   - **Not run:** the checks that did not run and why. A check that did not run counts for nothing.
   - **Caveats carried forward** from earlier checkpoints that still apply.
   - **How to go back:** the checkpoint file and the command that reruns from this stage.
5. **Quality scale** (shown for every checkpoint and for the whole run, in the quality check table at the top of `RUN_SUMMARY.md`):
   - *Good*: every scored check met its threshold.
   - *Usable, with N caveats*: at least one check warned; every result built on this checkpoint inherits the caveat.
   - *Accepted despite a failed check*: a FAIL was continued on purpose, with the reason on record.
   - *Not reliable until fixed*: a check failed; nothing downstream is trusted.
   - *Not done*: the stage has not run.
   The whole run is rated by its worst stage. Say plainly that thresholds are defaults for typical data and a Good rating means the checks passed, not that the biology is proven.
6. **Keep `RUN_SUMMARY.md` current.** Rebuild it after every checkpoint: quality check table, verdict and reason per stage, every check table with what it tests, the Not run list, parameters, software versions, timing and memory, input file hashes, figures with captions, folder map and resume instructions. This file is the record the user returns to.
7. **Be brief at the end.** When the run finishes, give a finished summary of at most one screen: the table of checkpoints with their quality rating, one line of conclusion each, the overall rating, the caveats that matter most, and the paths of `RUN_SUMMARY.md` and the checkpoint files. Do not repeat the cards.
8. **Defaults and questions.** Choose universal defaults and say what was chosen. Ask the user only for what they alone know: the species, which obs column is the condition and which is the sample (donor, patient), the expected biology for the expectations file, and whether a large run should proceed when memory or time is limited. Never use the dataset's own authors' labels as the marker panel or reference for that dataset; use them only to evaluate.
9. **Large data.** Load matrices sparse and in chunks, record peak memory per stage, work on a capped subsample for sensitivity runs and say so, and pause a heavy run when the user's compute is limited instead of shrinking the analysis silently.
10. **Honest wording.** A cluster name means "matches this marker panel", not "proven cell type". Pathway names describe shared genes, not activity. A difference between conditions is a treatment effect only if the design check says the condition is not confounded with a batch. State any subsampling, skipped sample or unrun check next to the number it affects.

## Finished summary template

```
scRNA-seq run: <dataset>, <n cells> cells, <species>; overall quality: <Good | Usable, with caveats | Not reliable | Incomplete>
0-1 Fetch/assemble  <rating>  <one line>
2 QC                <rating>  <one line>
3 Features          <rating>  <one line>
4 Dimension red.    <rating>  <one line>
5 Clustering        <rating>  <one line>
6 Annotation        <rating>  <one line>
7 Comparison        <rating or "not requested">  <one line>
8 Interpretation    <rating>  <one line>
Main caveats: <up to three>
Go back to any step: runs/<dataset>/checkpoints/NN_<stage>.h5ad; details in runs/<dataset>/RUN_SUMMARY.md
```
