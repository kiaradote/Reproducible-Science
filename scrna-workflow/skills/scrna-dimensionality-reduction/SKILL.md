---
name: scrna-dimensionality-reduction
description: "Checkpoint 4 of a reproducible scRNA-seq workflow: dimensionality reduction, run after scrna-normalization-features. Use when summarizing single-cell expression with PCA and UMAP: scales the selected genes, runs PCA, chooses the number of PCs from the data (smallest number capturing a set share of the real structure above shuffled-data noise, default 90%), tests whether PCs track depth, genes detected, mitochondrial share or batch, optionally corrects batch with Harmony, builds the neighbor graph and UMAP, and writes four presentation figures with plain-language legends for non-expert audiences. Produces a PASS/WARN/FAIL report."
---

# scRNA-seq dimensionality reduction (checkpoint 4)

Question answered: how many independent axes of variation does this dataset have, are they biology or measurement, and what does the cell landscape look like?
Start only after checkpoints 1 to 3 have no FAIL. Input: the checkpoint 3 object (normalized `X`, `var["highly_variable"]`, QC columns `n_genes`, `total`, `pct_mt` in `obs`).
Records use the same `{check, status, detail}` format as the other checkpoints. Call `write_dr_report(records, "dimred_report")`, show the overall status, and do not
continue on FAIL without asking the user. Carry every WARN into the final report. Scope stops at the neighbor graph and UMAP; clustering is the next checkpoint.

Helpers (loaded from `kernel.py`; need scanpy, scikit-learn, numpy, scipy, pandas, matplotlib; `harmonypy` for batch correction):
`dr_run`, `dr_choose_pcs`, `dr_pcs_for_threshold`, `dr_pc_corr`, `dr_batch_r2`, `dr_batch_share`, `dr_mixing`, `dr_harmony`, `dr_legends`, `dr_figures`, `write_dr_report`, `dr_catalog`.

## Workflow

1. **One call:** `adata, records = dr_run(adata, batch_key=None, correct_batch="auto", variance_threshold=0.9, variance_basis="signal", fig_dir=...)`.
   Then `write_dr_report(records, prefix)` and save `adata` as the checkpoint (`obsm["X_pca"]`, `obsm["X_umap"]`, optional `obsm["X_pca_harmony"]`, neighbor graph, `uns["dimred"]`).
2. **Number of PCs (chosen from the data, never fixed):** the smallest number of PCs that together capture `variance_threshold` of the *real structure*, meaning the
   variation that rises above what shuffled gene values would give (parallel analysis, 10 shuffles). Only the leading run of PCs that beat the shuffled data counts.
   `variance_basis="total"` uses all variation instead; with 2000 selected genes the first 50 PCs hold only about a third of total variation, so a high threshold on
   the total is usually unreachable. The shuffled-data count and the elbow are reported as cross-checks. `min_pcs` (5) and `max_pcs` (50) bound the result.
3. **Batch:** pass `batch_key` (donor, sample or batch column) to get the batch checks. `correct_batch="auto"` runs Harmony only when the batch explains more than 20% of the
   retained PC variation; `True` forces it, `False` forbids it. The report gives batch variation and neighbor mixing before and after. Correction cannot tell batch from
   real biology: if batches are confounded with the biology of interest (for example patient and response), check that known groups stay separate.
4. **Figures with legends:** `fig_dir` receives `pca_variance.png`, `pca_covariates.png`, `pca_scatter.png`, `umap.png` and `figure_legends.md`. Each legend is plain language,
   built from the computed numbers, and makes no biological claim. Show each figure with its legend when presenting to a general audience.
5. **Read the Not run list.** A check that did not run has not passed (batch checks need `batch_key`).

## Checks (thresholds are defaults, not literature values; change them with the function arguments)

PCA ran with finite results; PC count and the share captured; threshold reached within the PCs computed; every PC kept beats shuffled data; at most 30% of retained variation
in PCs that correlate above 0.5 with depth, genes detected or mitochondrial share; with `batch_key`, batch share of retained variation (limit 20%) and, if corrected, whether
batch variation fell and neighbor mixing rose; the neighbor graph has at most 1% of cells stranded in pieces under 10 cells (well separated cell groups can form their own large
pieces, which is not a fault); UMAP coordinates are finite.

## Pitfalls

- UMAP is for viewing. Distances between separated groups and group sizes are not quantitative; use the PCs or the neighbor graph for analysis.
- A raw "explain X% of total variance" rule does not transfer between datasets: noise level sets how many PCs are needed. Thresholding the structure above the shuffled
  noise floor does.
- The shuffled-gene test alone is generous on real single-cell data (all 50 PCs passed on a 1,937-cell pancreas dataset where the elbow was at 9); the threshold rule keeps it usable.
- Depth and genes detected are plotted and tested but not regressed out; a PC that tracks them is flagged, not removed.
- Dense scaled matrix of the selected genes is held in memory (cells x genes selected, float32); very large cell counts need more RAM.

## Tests

`test_scrna_dimred.py`: structure detection on synthetic data, threshold helper (monotone, exact, leading-run rule), noise fallback, threshold-not-reached warning, determinism with a seed,
depth-driven PCs flagged, batch effect detected and corrected (Harmony), correction declined or not needed, forced correction without a batch column, figures and legends with computed numbers,
report text, and a real count matrix.
