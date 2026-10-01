---
name: scrna-clustering
description: "Checkpoint 5 of a reproducible scRNA-seq workflow: clustering, run after scrna-dimensionality-reduction. Use when grouping single cells into clusters: clusters on the PCs chosen at checkpoint 4 (Harmony PCs if batch was corrected), scans Leiden resolutions and picks the finest one that is stable under cell resampling, has no tiny or dominant clusters and whose neighboring clusters still differ by enough genes, checks per-cluster stability and batch mixing, finds marker genes for every cluster, and writes four presentation figures with plain-language legends. Produces a PASS/WARN/FAIL report. Does not name cell types."
---

# scRNA-seq clustering (checkpoint 5)

Question answered: how many reproducible groups of cells does this dataset contain, and what genes define each?
Start only after checkpoints 1 to 4 have no FAIL. Input: the checkpoint 4 object (`obsm["X_pca"]` or `obsm["X_pca_harmony"]`, `uns["dimred"]` giving the representation and the PC count, normalized `X`, `var["highly_variable"]`, optional batch column in `obs`).
Records use the same `{check, status, detail}` format as the other checkpoints. Call `write_cl_report(records, "clustering_report")`, show the overall status, and do not continue on FAIL without asking the user. Carry every WARN into the final report. Scope stops at numbered clusters and their marker genes; naming cell types is a separate, later checkpoint.

Helpers (loaded from `kernel.py`; need scanpy, python-igraph, scikit-learn, numpy, scipy, pandas, matplotlib):
`cl_run`, `cl_graph`, `cl_leiden`, `cl_stability`, `cl_pick_resolution`, `cl_silhouette`, `cl_best_jaccard`, `cl_nearest_pairs`, `cl_pair_de`, `cl_markers`, `cl_batch_mixing`, `cl_legends`, `cl_figures`, `write_cl_report`, `cl_catalog`.

## Workflow

1. **One call:** `adata, records = cl_run(adata, batch_key=None, resolution=None, fig_dir=...)`. Then `write_cl_report(records, prefix)` and save `adata` as the checkpoint
   (`obs["cluster"]`, one `obs["leiden_r<res>"]` column per scanned resolution, `uns["clustering"]`, `uns["cluster_resolution_scan"]`, `uns["cluster_markers"]`).
2. **Representation:** the same PCs checkpoint 4 chose (`uns["dimred"]["rep"]`, first `n_pcs`); the neighbor graph is rebuilt on them (15 neighbors). Clusters are never computed on the UMAP.
3. **Resolution (chosen from the data unless `resolution=` is given):** Leiden (modularity) at 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0. For each, stability is the median adjusted Rand index between the full clustering and
   the clustering of 5 resamples that keep 80% of cells. Starting from the finest resolution, the first one that has stability of at least 0.8, at least 2 clusters, no cluster under max(10 cells, 0.5% of cells) and
   neighboring clusters that differ by at least 10 genes is used. If none qualifies, the most stable resolution is used and the report says WARN. Silhouette (up to 5000 cells) is reported as context, not used to choose.
4. **Checks on the chosen clustering:** stability, plausible cluster count (2 to min(60, cells/20)), no tiny clusters, no cluster holding more than 50% of cells, per-cluster recovery (best-match Jaccard of at least 0.6 under resampling),
   neighboring clusters distinguishable (each cluster against its nearest by centroid; Wilcoxon on the selected genes, at least 10 genes with |log2FC| of at least 0.5 and adjusted p under 0.01), markers (at least 5 per cluster,
   genes detected in at least 25% of the cluster, up to 500 cells per cluster), and with `batch_key`, no cluster of 10 or more cells that is at least 90% one batch.
5. **Figures with legends:** `fig_dir` receives `cluster_resolution.png`, `cluster_umap.png`, `cluster_sizes.png`, `cluster_markers.png`, `cluster_markers.csv`, `resolution_scan.csv` and `figure_legends.md`.
   Legends are plain language, built from the computed numbers, and name no cell types. Show each figure with its legend when presenting to a general audience.
6. **Read the Not run list.** A check that did not run has not passed (batch mixing needs `batch_key`).

## Pitfalls

- Thresholds are defaults, not literature values; change them with the function arguments. Record the values used when reporting.
- A real rare cell type can trip the tiny-cluster and per-cluster stability checks: on a 1,937-cell pancreas sample an 8-cell group with a clear mast-cell signature gave WARN on both although the clustering was otherwise stable
  (adjusted Rand index 0.84 against the authors' labels, which were not used). Inspect the cluster and its markers before dropping it; do not merge it just to turn the report green. With 80% resampling, a cluster of under 10 cells is
  hard to recover by construction.
- Splitting at the highest "stable" resolution can create halves that differ by almost no genes; that is why distinguishability is part of the resolution choice, not only a final check.
- Resolution is a modelling choice: clusters are numbered groups, not cell types. Marker genes found within the clusters are descriptive, since the same data defined the groups (double dipping), so p-values are optimistic.
- Windows: scanpy's Leiden wrapper routes igraph through numpy's random generator, which on Windows raises thousands of silent errors per call (about 20 s instead of 0.02 s on 2,000 cells). `cl_leiden` calls igraph directly with a seeded Python generator; keep it that way.
- Runs are reproducible for a fixed `seed`; changing igraph or scanpy versions can change cluster membership slightly.
- Memory: the graph is held in a light AnnData and marker tests use at most 500 cells per cluster, so cost grows with the number of clusters more than with the number of cells.

## Tests

`test_scrna_cluster.py`: synthetic groups recovered and stable, outputs and legends, determinism with a seed, tiny and dominant clusters flagged, pure noise not called stable, fixed resolution honoured,
resolution-choice rules, Jaccard and label ranking helpers, single-batch clusters flagged, Harmony representation used when present, report text, and a real pancreas count matrix.
