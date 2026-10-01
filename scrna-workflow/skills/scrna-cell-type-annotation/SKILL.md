---
name: scrna-cell-type-annotation
description: "Checkpoint 6 of a reproducible scRNA-seq workflow: marker-based cell type annotation of clusters, run after scrna-clustering. Use when naming clusters: scores each cluster against an editable marker panel with positive and negative markers (panel can be built from the CellGuide connector), optionally compares with a labelled reference dataset, and gives every cluster a confidence tier (confident, probable, ambiguous, unassigned) with the reason, so clusters that are mixed, transitional or unmatched are kept as ambiguous instead of being forced into a type. Writes three presentation figures with plain-language legends and a PASS/WARN/FAIL report."
---

# scRNA-seq cell type annotation (checkpoint 6)

Question answered: which known cell types do the clusters resemble, how sure are we, and which clusters do not fit a single type?
Start only after checkpoints 1 to 5 have no FAIL. Input: the checkpoint 5 object (`obs["cluster"]`, normalized log `X` with gene names, `obsm["X_umap"]`, optionally `uns["cluster_markers"]`) and a marker panel.
Records use the same `{check, status, detail}` format as the other checkpoints. Call `write_ann_report(records, "annotation_report")`, show the overall status, and do not continue on FAIL without asking the user.
A cluster is not a cell type: a label here means "matches this marker panel", not "is proven to be this cell type". Scope stops at cluster names; group comparisons are a later checkpoint.

Helpers (loaded from `kernel.py`; need scanpy-compatible AnnData, numpy, scipy, pandas, matplotlib):
`ann_run`, `ann_read_panel`, `ann_clean_cellguide`, `ann_suggest_negatives`, `ann_gene_stats`, `ann_score_table`, `ann_call`, `ann_reference_labels`, `ann_match_label`, `ann_legends`, `ann_figures`, `write_ann_report`, `ann_catalog`.

## Workflow

1. **Marker panel** (a table with columns `cell_type`, `positive`, `negative`, `source`; genes separated by `;`; saved with the run as `marker_panel_used.csv`). Three ways to get one:
   - From the CellGuide connector (calls must run in the `repl` tool): for each cell type name run `search_cell_types` (limit 1), then `get_marker_genes` with `marker_type="canonical"` and `"computational"`; collect
     `{name: {"id", "name", "canonical": [...], "computational": [...]}}` into a dict, write it to a JSON file under `handoff/`, and in the analysis kernel call `ann_clean_cellguide(raw)`. It keeps literature markers first, then data-derived markers for one species by score,
     removes duplicates, mitochondrial and ribosomal genes, and genes listed for more than two of the requested types.
   - Write the CSV by hand. This is the right route when the tissue is unusual.
   - Edit a CellGuide panel. **Always read it before use:** the report's panel strength and panel checks name the types that need work.
2. **Negative markers:** genes that should be off in a type (for example a lineage marker of an unrelated type). They are used only if you give them. `ann_suggest_negatives(panel)` fills empty ones from the other types' positives, but shared or leaky genes then trigger false conflicts
   (on the pancreas test it turned 5 correct names into ambiguous), so it is opt-in (`suggest_negatives=True`).
3. **One call:** `adata, records = ann_run(adata, panel, reference=None, ref_label_key="cell_type", label_map=None, fig_dir=...)`, then `write_ann_report(records, prefix)` and save `adata` as the checkpoint (`obs["cell_type"]`, `obs["annotation_tier"]`, `obs["cell_type_candidates"]`, `uns["annotation"]`, `uns["annotation_panel"]`).
4. **Scoring:** for every cluster and type, score = the mean over the type's positive genes of that gene's cluster mean rescaled across clusters (0 lowest cluster, 1 highest); positives present = share of the type's positive genes detected in the cluster, meaning in at least 25% of its cells or, when a gene is never that common (sparse droplet data), at least half as often as in the gene's best cluster, and never under 5%. Without the relative rule, B cells in a 10x PBMC set (markers detected in 10 to 20% of cells) stayed unassigned.
5. **Tiers (thresholds are defaults, not literature values):**
   - *confident:* at least half the positives present, score of 0.3 or more, leads the next type by 0.15 or more, none of its negatives on, and the reference (if given) agrees.
   - *probable:* lead of 0.05 to 0.15, or fewer than 10 cells, or the reference disagrees with an otherwise confident call.
   - *ambiguous:* top two types within 0.05, or negatives of the leading type are on, or the reference disagrees with a modest lead. The label stays "ambiguous" and `cell_type_candidates` holds the two types.
   - *unassigned:* no type has enough positive markers detected. Transition or mixed clusters, doublet-like clusters and types missing from the panel end up here or in ambiguous; that is intended.
6. **Reference comparison (optional):** pass a labelled AnnData (`reference`, labels in `ref_label_key`). Each cell is named by rank correlation with every reference type's average profile over the genes that differ most between reference types (a simplified SingleR), so the two datasets' normalisation does not need to match. The cluster takes the majority name and its agreement share.
   Reference names are matched to panel names by `label_map`, then exact, then substring match; unmatched names are not compared. Reference types with fewer than 10 cells are dropped.
7. **Figures with legends:** `fig_dir` receives `annotation_scores.png` (cluster by type score matrix with tier), `annotation_umap.png`, `annotation_markers.png` (panel genes by cluster), `annotation_table.csv` and `figure_legends.md`.
8. **Read the Not run list.** Pretrained-model tools (for example CellTypist) are never run by this checkpoint; the reference comparison is Not run without a reference.

## Pitfalls

- CellGuide computational markers mix species and repeat genes; endocrine cell types share most of their data-derived genes, and some types return nothing usable (on the pancreas test: acinar returned only mitochondrial rRNA genes, ductal returned one data-derived gene). The panel is a starting point; add markers for such types.
- A reference that lacks a cell type will name its cells after the nearest other type, sometimes with high agreement (the pancreas test reference had 3 acinar cells, so acinar cells were called ductal at 99%). Agreement shows consistency with the reference, not truth. Mast cells (9 reference cells) were called macrophages by the reference while their markers CPA3 and TPSAB1 were clear.
- Cluster means are rescaled across clusters, so a gene that is on in every cluster scores low everywhere. With very few clusters (two or three) the rescaling is coarse; treat tiers with caution.
- Tiny clusters are capped at probable. Marker and reference evidence for a cluster of fewer than 10 cells is weak by construction.
- Marker-based calls reflect expression at the cluster level; a cluster can hold two cell types if it was under-split. Compare clusters sharing a name before merging.
- Author or reference labels of the same dataset must never be passed as the panel or the reference for that dataset; use them only to evaluate.

## Tests

`test_scrna_annotate.py`: clear clusters named confident, mixed cluster ambiguous with both candidates, no-signal cluster unassigned, negative markers conflict, case-insensitive genes, missing panel genes (WARN and FAIL), single-gene panel flagged, missing cluster column fails, call rules (pure), negative suggestion only for empty types,
label matching, reference agreement across scales, reference disagreement downgrade, too few shared genes, CellGuide cleaning, figures and report text, and a real pancreas run in which every named cluster matches the authors' majority label.
