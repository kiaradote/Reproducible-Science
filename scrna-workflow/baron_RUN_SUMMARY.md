# Run summary

- **Run folder:** `C:\Users\Kiara\repos\Reproducible-Science\scrna-workflow\runs\baron`
- **Species:** human
- **Generated:** 2026-10-01 19:03:22 UTC
- **Software:** python 3.11.16, numpy 2.2.6, pandas 2.3.3, scipy 1.14.1, anndata 0.12.19, scanpy 1.10.4
- **Stages finished:** 9 of 9; **worst verdict so far:** FAIL

**Input files (SHA-256 recorded at download):**

- `matrix`: 5,670,823 bytes, `6fdee6c7f0c2c3098f851c761d7fafa4f9a359af200703258f7fda255ddd3e69`

## Quality check

**Overall: Usable, with caveats: read the caveats before drawing conclusions.** 3 stage(s) good, 4 usable with caveats, 1 accepted despite a failed check, 0 not reliable or needing attention, 0 not done, 1 not applicable.

How to read it: *Good* means every scored check met its default threshold; *Usable* means at least one check warned and the caveat applies to everything built on that stage; *Not reliable* means a check failed and nothing after it should be trusted until it is fixed. Thresholds are defaults chosen for typical data, not guarantees of a correct biological answer, and a check that did not run counts for nothing (see the last section).

| stage | what it makes sure of | quality | why | to go back to this point |
|---|---|---|---|---|
| 0. Fetch | the files are the intended ones and complete | Good | every scored check passed | `wf.invalidate("fetch")` then rerun; checkpoint `00_fetch.json` |
| 1. Assemble and data integrity | matrix, genes and cell metadata line up and match what was expected | Good | every scored check passed | `wf.invalidate("assemble")` then rerun; checkpoint `01_assemble.h5ad` |
| 2. QC and filtering | damaged or empty cells are removed and the filters are sensible | Accepted despite a failed check | accepted with the reason: GSE84133 has no mitochondrial genes, so the mitochondrial filter does not apply | `wf.invalidate("qc")` then rerun; checkpoint `02_qc.h5ad` |
| 3. Normalization and feature selection | counts are normalised and informative genes are chosen reproducibly | Usable, with 2 caveat(s) | 2 warning(s), no failures: low-detection HVGs (27.0% of selected genes are detected in under 1% of cells; limit 10%); HVG stability (two ... | `wf.invalidate("features")` then rerun; checkpoint `03_features.h5ad` |
| 4. Dimensionality reduction | the main structure is captured without being driven by technical effects | Good | every scored check passed | `wf.invalidate("dimred")` then rerun; checkpoint `04_dimred.h5ad` |
| 5. Clustering | groups of cells are stable under resampling and each is distinct | Usable, with 3 caveat(s) | 3 warning(s), no failures: resolution chosen from the data (resolution 0.1 (no resolution reached the stability target; the most stable ... | `wf.invalidate("cluster")` then rerun; checkpoint `05_cluster.h5ad` |
| 6. Cell type annotation | cell-type names are supported by marker genes and say when they are not | Usable, with 3 caveat(s) | 3 warning(s), no failures: marker panel (12 of 13 cell types scored; left out for missing genes: pancreatic acinar cell (0% present); ... | `wf.invalidate("annotate")` then rerun; checkpoint `06_annotate.h5ad` |
| 7. Condition comparison | the design can answer the condition question, and the differences are statistically valid | Not applicable | no scored check ran for this dataset: all 0 scored checks passed (1 informational). 17 catalogued check(s) did not run. | `wf.invalidate("compare")` then rerun; checkpoint `07_compare.h5ad` |
| 8. Biological interpretation and validation | the biological conclusions hold up when analysis choices change | Usable, with 1 caveat(s) | 1 warning(s), no failures: robustness: cell-type calls (5 robust, 1 sensitive, 0 unstable of 6 named cell types: mast cell (sensitive, ... | `wf.invalidate("interpret")` then rerun; checkpoint `08_interpret.h5ad` |

## Stages at a glance

| stage | state | verdict | why | time (s) | peak MB | checkpoint |
|---|---|---|---|---|---|---|
| 0. Fetch | current | PASS | all 1 scored checks passed, none failed or warned (2 informational) | 5.2 | 207 | `00_fetch.json` (0.0 MB) |
| 1. Assemble and data integrity | current | PASS | all 8 scored checks passed, none failed or warned (2 informational) | 17.7 | 694 | `01_assemble.h5ad` (29.5 MB) |
| 2. QC and filtering | current | FAIL (accepted: GSE84133 has no mitochondrial genes, so the mitochondrial filter does not apply) | 1 of 4 checks failed: mito genes found (no 'MT-' genes: Ensembl IDs, wrong species, or genes filtered?). 2 catalogued check(s) did not run. | 1.1 | 694 | `02_qc.h5ad` (29.3 MB) |
| 3. Normalization and feature selection | current | WARN | 2 warning(s), no failures: low-detection HVGs (27.0% of selected genes are detected in under 1% of cells; limit 10%); HVG stability (two halves) (overlap (Jaccard) 0.49 between selections on two random halves of 968 cells; limit 0.5). 2 cat | 7.2 | 694 | `03_features.h5ad` (58.2 MB) |
| 4. Dimensionality reduction | current | PASS | all 6 scored checks passed, none failed or warned (1 informational). 2 catalogued check(s) did not run. | 57.2 | 694 | `04_dimred.h5ad` (62.1 MB) |
| 5. Clustering | current | WARN | 3 warning(s), no failures: resolution chosen from the data (resolution 0.1 (no resolution reached the stability target; the most stable one): 9 clusters, stability 0.98, ); no tiny clusters (below 10 cells: 8 (8 cells)); per-cluster stabili | 7.5 | 738 | `05_cluster.h5ad` (62.2 MB) |
| 6. Cell type annotation | current | WARN | 3 warning(s), no failures: marker panel (12 of 13 cell types scored; left out for missing genes: pancreatic acinar cell (0% present); sources: CellGuid); panel strength (fewer than 2 positive genes: pancreatic acinar cell (0 of 0 positive g | 2.3 | 738 | `06_annotate.h5ad` (62.2 MB) |
| 7. Condition comparison | current | PASS | all 0 scored checks passed (1 informational). 17 catalogued check(s) did not run. | 0.3 | 747 | `07_compare.h5ad` (62.2 MB) |
| 8. Biological interpretation and validation | current | WARN | 1 warning(s), no failures: robustness: cell-type calls (5 robust, 1 sensitive, 0 unstable of 6 named cell types: mast cell (sensitive, median 100%)). 2 catalogued check(s) did not run. | 34.3 | 810 | `08_interpret.h5ad` (62.3 MB) |

**Checks across all stages:** 39 passed, 9 warned, 1 failed, 16 informational. A WARN is a point to review, not an automatic fault; a FAIL stops the run until fixed or accepted with a reason.

## 0. Fetch

Parameters used (stage-specific): 

```json
{
 "files": [
  {
   "name": "matrix",
   "filename": "GSM2230757_human1_umifm_counts.csv.gz",
   "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2230nnn/GSM2230757/suppl/GSM2230757_human1_umifm_counts.csv.gz"
  }
 ]
}
```

### Data integrity report: PASS

**Why:** all 1 scored checks passed, none failed or warned (2 informational)

**Checks run (3):** 1 pass, 0 warn, 0 fail, 2 info

| status | check | what it tests | result |
|---|---|---|---|
| INFO | file size | file size equals the expected size | matrix: 5670823 bytes; no reference given |
| INFO | file sha256 | file checksum equals the reference, or is recorded | matrix: 6fdee6c7f0c2c3098f851c761d7fafa4f9a359af200703258f7fda255ddd3e69 (5670823 bytes); no reference given, record this in DATA_SOURCES.md |
| PASS | gzip stream | gzip file decompresses fully (not truncated) | matrix: decompresses fully (78333463 bytes) |

## 1. Assemble and data integrity

Parameters used (stage-specific): 

```json
{
 "adapter": "baron_csv",
 "declared_scale": "counts",
 "expected": {
  "shape": [
   1937,
   20125
  ],
  "obs_counts": {
   "assigned_cluster": {
    "beta": 872,
    "alpha": 236,
    "delta": 214,
    "endothelial": 130,
    "ductal": 120,
    "acinar": 110,
    "quiescent_stellate": 92,
    "gamma": 70,
    "activated_stellate": 51,
    "macrophage": 14,
    "epsilon": 13,
    "mast": 8,
    "schwann": 5,
    "t_cell": 2
   }
  }
 }
}
```

### Data integrity report: PASS

**Why:** all 8 scored checks passed, none failed or warned (2 informational)

**Checks run (10):** 8 pass, 0 warn, 0 fail, 2 info

| status | check | what it tests | result |
|---|---|---|---|
| PASS | unique cell IDs | no duplicated cell IDs | 0 duplicated IDs |
| PASS | cell ID sets match metadata | matrix and metadata describe the same cells | 0 only in matrix, 0 only in metadata |
| PASS | column order vs metadata | matrix columns are in the metadata's order (label-row test) | 1937 columns agree in order |
| INFO | alignment note | additional check | labels and metadata come from the same table, so the order test here is trivially true |
| PASS | NaN/inf values | no missing or infinite values | 0 in 1937 sampled cells |
| PASS | negative values | no negative values (scaled data legitimately has some) | 0 (scaled data legitimately has negatives) |
| INFO | detected scale | what the matrix holds: counts, TPM, log-normalized or scaled | counts; max 4318.00, integer-valued True |
| PASS | scale matches declaration | detected scale equals the scale the source states | declared counts, detected counts |
| PASS | shape (cells x genes) | shape and label counts equal the published numbers | (1937, 20125) |
| PASS | counts of assigned_cluster | shape and label counts equal the published numbers | {'beta': 872, 'alpha': 236, 'delta': 214, 'endothelial': 130, 'ductal': 120, 'acinar': 110, 'quiescent_stellate': 92, 'gamma': 70, 'activated_stellate': 51, 'macrophage': 14, 'epsilon': 13, 'mast': 8, 'schwann': 5, 't_cell': 2} |

## 2. QC and filtering

Parameters used (stage-specific): 

```json
{
 "min_genes": 200,
 "max_genes": null,
 "max_pct_mt": 20,
 "min_cells": 3,
 "sample_col": null
}
```

### QC report: FAIL

**Why:** 1 of 4 checks failed: mito genes found (no 'MT-' genes: Ensembl IDs, wrong species, or genes filtered?). 2 catalogued check(s) did not run.

**Checks run (4):** 1 pass, 0 warn, 1 fail, 2 info

| status | check | what it tests | result |
|---|---|---|---|
| FAIL | mito genes found | a valid mitochondrial gene set exists in the matrix | no 'MT-' genes: Ensembl IDs, wrong species, or genes filtered? |
| INFO | filter: n_genes<200 | cells removed by each QC filter | removes 0 of 1937 cells (0.0%) |
| PASS | total cells removed | combined cell loss is within the limit | 0 of 1937 (0.0%); limit 20% |
| INFO | genes removed (<min_cells) | cells removed by each QC filter | 5387 of 20125 genes seen in fewer than 3 cells |

**Not run:**

- mito signal present: the mitochondrial genes carry counts, so the mitochondrial filter can actually detect stressed cells
- per-sample loss: no sample loses more than the limit

## 3. Normalization and feature selection

Parameters used (stage-specific): 

```json
{
 "n_top_genes": 2000,
 "batch_key": null,
 "sample_col": null
}
```

### Normalization and feature selection report: WARN

**Why:** 2 warning(s), no failures: low-detection HVGs (27.0% of selected genes are detected in under 1% of cells; limit 10%); HVG stability (two halves) (overlap (Jaccard) 0.49 between selections on two random halves of 968 cells; limit 0.5). 2 catalogued check(s) did not run.

**Checks run (11):** 9 pass, 2 warn, 0 fail, 0 info

| status | check | what it tests | result |
|---|---|---|---|
| PASS | method matches scale | the normalization fits what the matrix holds (scaled or unknown data is refused) | counts: normalize to 10000 per cell, then natural log1p; raw counts kept in layers['counts'] |
| PASS | per-cell totals comparable | after normalization every cell has about the same linear-scale total | linear-scale total median 10000, CV 0.000 (limit 0.1) |
| PASS | finite and non-negative values | no NaN, infinite or negative values after the transform | 0 NaN/inf, 0 negative |
| PASS | sparsity preserved | the matrix was not densified or padded | 3724704 stored values, 3724704 before |
| PASS | input values kept | the original values survive (counts layer, or X unchanged) | layers['counts'] holds 3724704 stored values, integer-valued |
| PASS | dominant-gene cells | few cells are dominated by a single gene | 7 of 1937 cells have one gene above 50% of their total (limit 1% of cells) |
| PASS | HVG count | the requested number of unique variable genes was selected | 2000 selected, 2000 requested, gene names unique: True |
| PASS | technical genes among HVGs | mitochondrial, ribosomal and haemoglobin genes do not dominate the selection | 2 of 2000 (0.1%): mitochondrial 0, ribosomal 2, haemoglobin 0; limit 10% |
| WARN | low-detection HVGs | selected genes are detected in enough cells | 27.0% of selected genes are detected in under 1% of cells; limit 10% |
| PASS | mean expression coverage | selected genes are not concentrated among the lowest-expressed genes | 14.8% of selected genes are in the lowest-expressed 10% of genes; limit 30% |
| WARN | HVG stability (two halves) | two random halves of the cells select similar genes | overlap (Jaccard) 0.49 between selections on two random halves of 968 cells; limit 0.5 |

**Not run:**

- depth spread across samples: median depth does not differ wildly between samples
- batch robustness: selected genes are variable in many batches, not one

## 4. Dimensionality reduction

Parameters used (stage-specific): 

```json
{
 "batch_key": null,
 "variance_threshold": 0.9
}
```

### Dimensionality reduction report: PASS

**Why:** all 6 scored checks passed, none failed or warned (1 informational). 2 catalogued check(s) did not run.

**Checks run (7):** 6 pass, 0 warn, 0 fail, 1 info

| status | check | what it tests | result |
|---|---|---|---|
| PASS | PCA on selected genes | PCA ran on the scaled selected genes and gave finite results | 2000 selected genes, 1937 cells, 50 components computed |
| INFO | number of PCs chosen from the data | the PC count is the smallest number capturing a set share of the real structure in this dataset | 12 PCs used: the smallest number capturing 90% of the real structure (variation above shuffled data); 50 PCs beat shuffled data (10 shuffles), elbow at 9 |
| PASS | variance threshold reached | the chosen share of variation is reached within the PCs computed | 91% of the real structure (variation above shuffled data) captured by 12 PCs (target 90%) |
| PASS | chosen PCs beat shuffled data | every PC kept rises above what shuffled gene values would give | 12 of 12 PCs rise above shuffled data |
| PASS | PCs driven by technical covariates | retained PCs are not mainly depth, genes-detected or mitochondrial effects | 0 of 12 PCs (0% of retained variation) correlate above 0.5 with a measurement; limit 30% |
| PASS | neighbor graph fragments | no cells are stranded in tiny disconnected pieces of the neighbor graph | 1 connected piece(s), largest holds 100.0% of cells; 0 cells (0.0%) sit in pieces under 10 cells. Well separated cell groups can form their own large pieces, which is not a fault |
| PASS | UMAP finite | UMAP coordinates are finite for every cell | all coordinates finite |

**Not run:**

- batch structure in PCs: how much of the retained PC variation is explained by the batch column
- batch correction effect: correction lowered batch structure and improved mixing

## 5. Clustering

Parameters used (stage-specific): 

```json
{
 "batch_key": null,
 "resolution": null
}
```

### Clustering report: WARN

**Why:** 3 warning(s), no failures: resolution chosen from the data (resolution 0.1 (no resolution reached the stability target; the most stable one): 9 clusters, stability 0.98, ); no tiny clusters (below 10 cells: 8 (8 cells)); per-cluster stability (recovered under 0.6: 8 (0.22)). 1 catalogued check(s) did not run.

**Checks run (9):** 5 pass, 3 warn, 0 fail, 1 info

| status | check | what it tests | result |
|---|---|---|---|
| INFO | representation used | clusters are found on the same (batch-corrected if needed) PCs chosen at the previous checkpoint | X_pca, first 12 components; neighbor graph rebuilt with 15 neighbors |
| WARN | resolution chosen from the data | the resolution is the finest one whose clusters stay stable when cells are resampled and whose neighboring clusters still differ by enough genes | resolution 0.1 (no resolution reached the stability target; the most stable one): 9 clusters, stability 0.98, silhouette 0.55; scanned 0.1 to 2.0 |
| PASS | chosen resolution is stable | clusters at the chosen resolution are reproduced when 20% of cells are left out | median adjusted Rand index 0.98 (lowest 0.95) over 5 resamples keeping 80% of cells; target 0.8 |
| PASS | cluster count | the number of clusters is plausible for the number of cells | 9 clusters for 1937 cells (plausible range 2 to 60) |
| WARN | no tiny clusters | no cluster is too small to trust | below 10 cells: 8 (8 cells) |
| PASS | no dominant cluster | no single cluster holds most of the cells | largest cluster 0 holds 40% of cells; limit 50% |
| WARN | per-cluster stability | every cluster is recovered when cells are resampled | recovered under 0.6: 8 (0.22) |
| PASS | neighboring clusters distinguishable | clusters next to each other differ by enough genes to be kept apart | each cluster's nearest neighbor differs by at least 10 genes (/log2FC/ at least 0.5, adjusted p under 0.01; selected genes) |
| PASS | cluster markers | every cluster has marker genes that define it | every cluster has at least 5 markers; median 591 per cluster. Tested on up to 500 cells per cluster and genes detected in at least 25% of a cluster's cells |

**Not run:**

- batch mixing within clusters: no cluster is made almost entirely of one batch

## 6. Cell type annotation

Parameters used (stage-specific): 

```json
{
 "panel": "C:\\Users\\Kiara\\repos\\Reproducible-Science\\scrna-workflow\\runs\\baron\\marker_panel.csv",
 "reference": null,
 "reference_label_key": "cell_type",
 "label_map": {
  "quiescent_stellate": "pancreatic stellate cell",
  "activated_stellate": "pancreatic stellate cell",
  "t_cell": "T cell",
  "schwann": "Schwann cell",
  "gamma": "pancreatic PP cell",
  "epsilon": "pancreatic epsilon cell",
  "delta": "pancreatic delta cell"
 },
 "panel_sha256": "64d4a21dfd009f71b3e34effa09eb0bb8c0f768153db3c54fabbec28943a2a28",
 "reference_sha256": null
}
```

### Cell type annotation report: WARN

**Why:** 3 warning(s), no failures: marker panel (12 of 13 cell types scored; left out for missing genes: pancreatic acinar cell (0% present); sources: CellGuid); panel strength (fewer than 2 positive genes: pancreatic acinar cell (0 of 0 positive genes in the data); add markers to the pa); enough usable markers (too few detectable positive markers to score these types, so no cluster can be named after them: Schwann cell ). 2 catalogued check(s) did not run.

**Checks run (11):** 3 pass, 3 warn, 0 fail, 5 info

| status | check | what it tests | result |
|---|---|---|---|
| WARN | marker panel | each cell type in the panel has enough of its positive marker genes present in the data to be scored | 12 of 13 cell types scored; left out for missing genes: pancreatic acinar cell (0% present); sources: CellGuide canonical + Homo sapiens computational markers; genes shared by more than 2 types dropped |
| WARN | panel strength | every cell type has at least 2 positive markers, so no call rests on a single gene | fewer than 2 positive genes: pancreatic acinar cell (0 of 0 positive genes in the data); add markers to the panel file |
| INFO | negative markers | genes that should be off in a cell type are defined (given, or suggested from the other types) | none given: only positive evidence is used (add negative markers to the panel to enable the conflict check) |
| INFO | sparse panel genes | panel genes too rarely detected in this dataset (dropout) to count as evidence are set aside and listed | set aside because they are detected in under 5% of cells in every cluster: pancreatic stellate cell: HOXC9; mast cell: IL5RA; T cell: CCL4; Schwann cell: NRXN3, PLP1, SOX10, MPZ, CADM2 |
| WARN | enough usable markers | every panel type keeps at least 3 detectable positive markers (or all it was given, if fewer), so a call never rests on one or two noisy genes | too few detectable positive markers to score these types, so no cluster can be named after them: Schwann cell (2 of 8 usable); add markers that are detected in this kind of data |
| INFO | annotation tiers | how many clusters could be named with confidence, and how many are ambiguous or unassigned | 4 confident, 3 probable, 1 ambiguous, 1 unassigned of 9 clusters |
| PASS | cells named | most cells sit in clusters that could be named (confident or probable) | 79% of cells are in named clusters (limit 70%); not named: 1 (ambiguous: pancreatic PP cell and pancreatic delta cell score almost the same (difference 0.02)), 6 (unassigned: no cell type has at least half of its positive markers detected with a clear signal) |
| PASS | negative marker conflicts | no named cluster switches on genes that should be off in its assigned type | no leading type has negative markers switched on |
| PASS | data-driven marker support | the cluster's own top marker genes (found without the panel) include genes from its assigned type's panel | every named cluster has at least one panel gene among its own top markers |
| INFO | same type in several clusters | clusters sharing a name are reported, since they may be subtypes, states or over-splitting | pancreatic beta cell: clusters 0, 7 (subtypes, states or over-splitting; compare their markers) |
| INFO | panel types not found | cell types in the panel that no cluster was assigned to | no cluster was named: pancreatic delta cell, pancreatic PP cell, pancreatic epsilon cell, macrophage, T cell, Schwann cell |

**Not run:**

- reference agreement: an independent labelled reference dataset gives the same name as the marker panel
- pretrained annotation tools: automated annotation with pretrained models (for example CellTypist) was compared

## 7. Condition comparison

Parameters used (stage-specific): 

```json
{
 "condition": null,
 "gene_sets_sha256": null
}
```

### Condition comparison report: PASS

**Why:** all 0 scored checks passed (1 informational). 17 catalogued check(s) did not run.

## Checks

| status | check | what it tests | result |
|---|---|---|---|
| INFO | groups | exactly two conditions are compared and both contain cells | no comparison requested for this dataset: set 'condition' and 'sample' in the stage parameters to compare two groups |

**Not run:**

- replicates: each condition has enough independent samples (donors, patients, animals) for a test; cells of one sample are not replicates
- pairing: samples measured in both conditions are compared with themselves (paired design)
- condition vs batch: the condition is not perfectly confounded with a batch or library, which would make a treatment effect indistinguishable from a technical one
- clusters follow condition: clusters are not almost entirely made of one condition, which would mean the treatment itself drives the clustering
- cell types follow condition: a cell-type name is not almost entirely used in one condition only
- testable cell types: each compared cell type has enough cells from enough samples in both conditions
- unlabelled cells: the share of cells without a cell-type name is small
- abundance: the share of each cell type per sample is compared between conditions by two different models that agree
- differential expression: genes are compared between conditions within each cell type on sample-level (pseudobulk) totals, not on single cells
- null check: swapping the condition labels within samples gives (almost) no discoveries
- sample consistency: differentially expressed genes move in the same direction in most individual samples
- state change: cells of one type can be told apart by condition in samples the classifier never saw
- pathway enrichment: ranked differential-expression results are tested against gene sets with a pre-ranked method
- transcription factor activity: regulator activity inferred from a transcription-factor network (not implemented in this version: needs the OmniPath network and the decoupler package)
- cell-cell communication: ligand-receptor signalling between cell types (not implemented in this version: needs a ligand-receptor resource)
- spatial association: cell types or genes that co-vary in tissue position (needs coordinates in obsm['spatial'])
- gene modules: co-expressed gene modules that change with the condition (use the gene programs of checkpoint 8 for the descriptive version)

## 8. Biological interpretation and validation

Parameters used (stage-specific): 

```json
{
 "gene_sets": "C:\\Users\\Kiara\\repos\\Reproducible-Science\\scrna-workflow\\runs\\baron\\gene_sets\\kegg_hsa.gmt",
 "gene_sets_kegg": "hsa",
 "expectations": "C:\\Users\\Kiara\\repos\\Reproducible-Science\\scrna-workflow\\runs\\baron\\expectations.csv",
 "signature": null,
 "perturb": true,
 "max_cells": 3000,
 "gene_sets_sha256": null,
 "expectations_sha256": "6d8d5259ea1e839e20bb4f26df33f39294a3bc932b0941d38fa47edbe72b0cfc"
}
```

### Biological interpretation and validation report: WARN

**Why:** 1 warning(s), no failures: robustness: cell-type calls (5 robust, 1 sensitive, 0 unstable of 6 named cell types: mast cell (sensitive, median 100%)). 2 catalogued check(s) did not run.

**Checks run (9):** 6 pass, 1 warn, 0 fail, 2 info

## Findings, their evidence and how far to trust them

28 findings; 22 carry a robust or validated label. A finding labelled sensitive or unstable depends on a choice that could reasonably have been made differently.

| finding | evidence | robustness | caveats |
|---|---|---|---|
| endothelial cell cells are present | 129 cells; clusters 4 (confident); same name in a median 100% of cells across settings (worst 99%) | robust | none recorded |
| mast cell cells are present | 8 cells; clusters 8 (probable); same name in a median 100% of cells across settings (worst 0%) | sensitive | cluster 8: only 8 cells |
| pancreatic alpha cell cells are present | 243 cells; clusters 2 (confident); same name in a median 100% of cells across settings (worst 99%) | robust | none recorded |
| pancreatic beta cell cells are present | 864 cells; clusters 0, 7 (confident, probable); same name in a median 100% of cells across settings (worst 100%) | robust | cluster 7: modest lead over pancreatic PP cell (0.08) |
| pancreatic ductal cell cells are present | 120 cells; clusters 5 (probable); same name in a median 100% of cells across settings (worst 99%) | robust | cluster 5: modest lead over T cell (0.14) |
| pancreatic stellate cell cells are present | 163 cells; clusters 3 (confident); same name in a median 100% of cells across settings (worst 90%) | robust | none recorded |
| pancreatic beta cell: marker genes are enriched for Maturity onset diabetes of the young [hsa04950] | 6 of 23 pathway genes among the markers, fold enrichment 38.4, adjusted p 2.0e-07 | robust (term set overlap 1.00 across settings) | enrichment describes shared genes, not pathway activity; KEGG disease or infection pathway: it contains generic immune genes, so the name does not mean the disease is present |
| pancreatic beta cell: marker genes are enriched for Insulin secretion [hsa04911] | 7 of 70 pathway genes among the markers, fold enrichment 14.7, adjusted p 5.9e-06 | robust (term set overlap 1.00 across settings) | enrichment describes shared genes, not pathway activity; KEGG disease or infection pathway: it contains generic immune genes, so the name does not mean the disease is present |
| pancreatic alpha cell: marker genes are enriched for Antigen processing and presentation [hsa04612] | 5 of 56 pathway genes among the markers, fold enrichment 13.2, adjusted p 1.3e-03 | robust (term set overlap 1.00 across settings) | enrichment describes shared genes, not pathway activity |
| pancreatic alpha cell: marker genes are enriched for Cell adhesion molecule (CAM) interaction [hsa04514] | 6 of 112 pathway genes among the markers, fold enrichment 7.9, adjusted p 2.0e-03 | robust (term set overlap 1.00 across settings) | enrichment describes shared genes, not pathway activity |
| pancreatic stellate cell: marker genes are enriched for Cytoskeleton in muscle cells [hsa04820] | 16 of 173 pathway genes among the markers, fold enrichment 13.6, adjusted p 2.1e-12 | robust (term set overlap 0.90 across settings) | enrichment describes shared genes, not pathway activity |
| pancreatic stellate cell: marker genes are enriched for Protein digestion and absorption [hsa04974] | 11 of 74 pathway genes among the markers, fold enrichment 21.9, adjusted p 6.7e-11 | robust (term set overlap 0.90 across settings) | enrichment describes shared genes, not pathway activity; KEGG disease or infection pathway: it contains generic immune genes, so the name does not mean the disease is present |
| endothelial cell: marker genes are enriched for Leukocyte transendothelial migration [hsa04670] | 7 of 97 pathway genes among the markers, fold enrichment 10.6, adjusted p 1.9e-04 | robust (term set overlap 1.00 across settings) | enrichment describes shared genes, not pathway activity |
| endothelial cell: marker genes are enriched for Cell adhesion molecule (CAM) interaction [hsa04514] | 7 of 112 pathway genes among the markers, fold enrichment 9.2, adjusted p 1.9e-04 | robust (term set overlap 1.00 across settings) | enrichment describes shared genes, not pathway activity |
| pancreatic ductal cell: marker genes are enriched for Cornified envelope formation [hsa04382] | 15 of 110 pathway genes among the markers, fold enrichment 20.1, adjusted p 3.7e-14 | robust (term set overlap 1.00 across settings) | enrichment describes shared genes, not pathway activity |
| pancreatic ductal cell: marker genes are enriched for Leukocyte transendothelial migration [hsa04670] | 10 of 97 pathway genes among the markers, fold enrichment 15.2, adjusted p 2.7e-08 | robust (term set overlap 1.00 across settings) | enrichment describes shared genes, not pathway activity |
| program_1 is used mostly by pancreatic ductal cell | 91% of its usage; top genes TMSB4X, ANXA2, KRT7, KRT19, S100A6 | reproducible across restarts | programs are descriptive and not tied to a pathway unless tested |
| program_3 is used mostly by pancreatic beta cell | 70% of its usage; top genes PCSK1, SCG2, HSP90B1, HSPA5, VGF | reproducible across restarts | programs are descriptive and not tied to a pathway unless tested |
| program_4 is used mostly by pancreatic stellate cell | 91% of its usage; top genes TIMP1, IGFBP4, FTH1, IFITM3, SPARC | reproducible across restarts | programs are descriptive and not tied to a pathway unless tested |
| program_6 is used mostly by pancreatic stellate cell | 61% of its usage; top genes FTH1, CTSD, TMSB4X, CTSB, CD68 | reproducible across restarts | programs are descriptive and not tied to a pathway unless tested |
| program_7 is used mostly by pancreatic beta cell | 74% of its usage; top genes PPP1R1A, FTH1, G6PC2, CDKN1C, TTR | reproducible across restarts | programs are descriptive and not tied to a pathway unless tested |
| program_8 is used mostly by endothelial cell | 94% of its usage; top genes PLVAP, HLA-B, TMSB4X, TMSB10, ENG | reproducible across restarts | programs are descriptive and not tied to a pathway unless tested |
| program_11 is used mostly by pancreatic alpha cell | 89% of its usage; top genes GCG, TTR, CLU, CHGA, TM4SF4 | reproducible across restarts | programs are descriptive and not tied to a pathway unless tested |
| Expectation: Beta cells express insulin and islet genes | 4 of 4 genes higher in pancreatic beta cell (mean difference; * = passes): INS +3.95*, IAPP +3.42*, MAFA +0.80*, PDX1 +0.54* | n/a | none |
| Expectation: Alpha cells express glucagon | 3 of 3 genes higher in pancreatic alpha cell (mean difference; * = passes): GCG +5.24*, ARX +0.66*, IRX2 +1.81* | n/a | none |
| Expectation: Ductal cells express ductal genes | 3 of 3 genes higher in pancreatic ductal cell (mean difference; * = passes): CFTR +1.30*, KRT19 +2.92*, ANXA4 +1.84* | n/a | none |
| Expectation: Negative control: glucagon is not higher in endothelial cells | 1 of 1 genes lower in endothelial cell (mean difference; * = passes): GCG -0.68* | n/a | none |
| Expectation: Beta cells are 30 to 60% of islet donor cells | 44.6% of cells (expected 30 to 60%) | n/a | none |

## Checks

| status | check | what it tests | result |
|---|---|---|---|
| PASS | gene sets | at least one gene-set collection is loaded and usable against the genes detected in this dataset | 352 of 372 sets keep 10 to 500 genes among the 14,738 detected genes |
| INFO | background | the genes tested for enrichment are compared with the genes actually detected, not with the whole genome | 14,738 genes detected in at least 3 cells are the background for every enrichment test |
| PASS | pathway analysis | cell-type marker genes are tested for over-representation of gene sets with multiple-testing correction | 83 non-redundant pathways below FDR 0.05 in 5 of 5 cell types; markers per type: median 100 |
| PASS | gene programs | data-driven expression programs are reproducible across random restarts | 12 programs; restart similarity 0.95 (target 0.8); 7 used mostly (50% or more) by one cell type |
| PASS | expectations | known biology listed as positive and negative controls is reproduced by the data | 5 met, 0 not met, 0 not testable |
| PASS | robustness: clusters | cluster membership is stable when seeds, cells, variable genes, PCs and resolution change | median adjusted Rand index 0.99 across 10 changed settings on 1,937 cells (robust) |
| WARN | robustness: cell-type calls | each cell-type name survives the same changes | 5 robust, 1 sensitive, 0 unstable of 6 named cell types: mast cell (sensitive, median 100%) |
| PASS | robustness: pathways | the enriched pathways are found again after the same changes | median overlap of significant pathway sets with the baseline 1.00 across settings (0.6 or more = robust) |
| INFO | caveats carried | warnings and unrun checks from earlier checkpoints are carried into the conclusions | 8 warning(s), 1 failure(s) and 26 unrun check(s) from earlier checkpoints |

**Not run:**

- pre-ranked enrichment: a ranked differential-expression result (checkpoint 7) is tested against the gene sets
- signature validation: a gene signature separates two groups in samples it was not trained on

## Caveats carried from earlier checkpoints

These limit every conclusion above.

- [FAIL] qc: mito genes found: no 'MT-' genes: Ensembl IDs, wrong species, or genes filtered?
- [WARN] annotate: marker panel: 12 of 13 cell types scored; left out for missing genes: pancreatic acinar cell (0% present); sources: CellGuide canonical + Homo sapiens computational markers; genes shared by more than 2 types droppe
- [WARN] annotate: panel strength: fewer than 2 positive genes: pancreatic acinar cell (0 of 0 positive genes in the data); add markers to the panel file
- [WARN] annotate: enough usable markers: too few detectable positive markers to score these types, so no cluster can be named after them: Schwann cell (2 of 8 usable); add markers that are detected in this kind of data
- [WARN] cluster: resolution chosen from the data: resolution 0.1 (no resolution reached the stability target; the most stable one): 9 clusters, stability 0.98, silhouette 0.55; scanned 0.1 to 2.0
- [WARN] cluster: no tiny clusters: below 10 cells: 8 (8 cells)
- [WARN] cluster: per-cluster stability: recovered under 0.6: 8 (0.22)
- [WARN] features: low-detection HVGs: 27.0% of selected genes are detected in under 1% of cells; limit 10%
- [WARN] features: HVG stability (two halves): overlap (Jaccard) 0.49 between selections on two random halves of 968 cells; limit 0.5
- [not run] annotate: reference agreement
- [not run] annotate: pretrained annotation tools
- [not run] cluster: batch mixing within clusters
- [not run] compare: replicates
- [not run] compare: pairing
- [not run] compare: condition vs batch
- [not run] compare: clusters follow condition
- [not run] compare: cell types follow condition
- [not run] compare: testable cell types
- [not run] compare: unlabelled cells
- [not run] compare: abundance
- [not run] compare: differential expression
- [not run] compare: null check
- [not run] compare: sample consistency
- [not run] compare: state change
- [not run] compare: pathway enrichment
- [not run] compare: transcription factor activity
- [not run] compare: cell-cell communication
- [not run] compare: spatial association
- [not run] compare: gene modules
- [not run] dimred: batch structure in PCs
- [not run] dimred: batch correction effect
- [not run] features: depth spread across samples
- [not run] features: batch robustness
- [not run] qc: mito signal present
- [not run] qc: per-sample loss

## Checks that did not run

A check that did not run has not passed. Each is listed with the stage that owns it.

- 2. QC and filtering: mito signal present
- 2. QC and filtering: per-sample loss
- 3. Normalization and feature selection: depth spread across samples
- 3. Normalization and feature selection: batch robustness
- 4. Dimensionality reduction: batch structure in PCs
- 4. Dimensionality reduction: batch correction effect
- 5. Clustering: batch mixing within clusters
- 6. Cell type annotation: reference agreement
- 6. Cell type annotation: pretrained annotation tools
- 7. Condition comparison: replicates
- 7. Condition comparison: pairing
- 7. Condition comparison: condition vs batch
- 7. Condition comparison: clusters follow condition
- 7. Condition comparison: cell types follow condition
- 7. Condition comparison: testable cell types
- 7. Condition comparison: unlabelled cells
- 7. Condition comparison: abundance
- 7. Condition comparison: differential expression
- 7. Condition comparison: null check
- 7. Condition comparison: sample consistency
- 7. Condition comparison: state change
- 7. Condition comparison: pathway enrichment
- 7. Condition comparison: transcription factor activity
- 7. Condition comparison: cell-cell communication
- 7. Condition comparison: spatial association
- 7. Condition comparison: gene modules
- 8. Biological interpretation and validation: pre-ranked enrichment
- 8. Biological interpretation and validation: signature validation

## Figures

**figures/annotate/**

Legends: see `figures/annotate/figure_legends.md`

![annotation_markers](figures/annotate/annotation_markers.png)

*annotation markers: Dot plot of panel marker genes (columns, grouped by cell type) across 9 clusters (rows): dot size = share of cells expressing the gene, colour = mean expression scaled per gene.*

![annotation_scores](figures/annotate/annotation_scores.png)

*annotation scores: Heat map of marker-panel scores for 9 clusters (rows) against 12 cell types (columns); black outlines mark the assigned type and the right-hand label gives the confidence tier.*

![annotation_umap](figures/annotate/annotation_umap.png)

*annotation umap: UMAP of 1,937 cells coloured by assigned cell type; grey = ambiguous or unassigned. 7 of 9 clusters are named.*

**figures/cluster/**

Legends: see `figures/cluster/figure_legends.md`

![cluster_markers](figures/cluster/cluster_markers.png)

*cluster markers: Dot plot of the top marker genes (columns) for each of 9 clusters (rows): dot size = share of the cluster's cells expressing the gene, colour = mean expression scaled per gene.*

![cluster_resolution](figures/cluster/cluster_resolution.png)

*cluster resolution: Number of clusters (top), stability under resampling (middle; line = median adjusted Rand index, band = range, dotted line = target) and silhouette (bottom) at each tested resolution; the red line marks the chosen resolution 0.1.*

![cluster_sizes](figures/cluster/cluster_sizes.png)

*cluster sizes: Left: number of cells in each of 9 clusters (red = under the minimum size). Right: how well each cluster is recovered when 20% of cells are left out (dotted line = threshold).*

![cluster_umap](figures/cluster/cluster_umap.png)

*cluster umap: UMAP of 1,937 cells coloured by 9 clusters, with cluster numbers at the cluster centres; clusters were computed on the first 12 principal components, not on the UMAP.*

**figures/dimred/**

Legends: see `figures/dimred/figure_legends.md`

![pca_covariates](figures/dimred/pca_covariates.png)

*pca covariates: Absolute correlation (colour and number) between each of the first 12 principal components (columns) and each per-cell measurement (rows).*

![pca_scatter](figures/dimred/pca_scatter.png)

*pca scatter: 1,937 cells on the first two principal components, coloured by total signal per cell on a log scale.*

![pca_variance](figures/dimred/pca_variance.png)

*pca variance: Left: variation captured by each of 50 principal components, with the level expected from shuffled data as a dashed line and the 12 components used marked by a red line. Right: cumulative share of real structure reached by the first components.*

![umap](figures/dimred/umap.png)

*umap: UMAP of 1,937 cells, coloured by total signal per cell on a log scale.*

**figures/interpret/**

Legends: see `figures/interpret/figure_legends.md`

![int_pathways](figures/interpret/int_pathways.png)

*int pathways: Dot plot of pathways (columns) enriched among the marker genes of each cell type (rows); dot size = marker genes in the pathway, colour = -log10 adjusted p-value.*

![int_programs](figures/interpret/int_programs.png)

*int programs: Heat map of the usage of 12 data-driven gene programs (rows) in each cell type (columns), with each program's top genes.*

![int_robustness](figures/interpret/int_robustness.png)

*int robustness: Left: adjusted Rand index of the clusters against the baseline for each changed setting. Right: share of each cell type's cells keeping the same name under each setting.*

![int_validation](figures/interpret/int_validation.png)

*int validation: Left: outcome of each listed expectation. Right: leave-one-sample-out area under the curve of the gene signature.*

**figures/**

![mean_variance](figures/mean_variance.png)

*mean variance: Mean expression (x, log scale) against normalised variability (y) for 14,738 genes; the 2,000 selected variable genes are red and the rest grey.*

![qc_before_filtering](figures/qc_before_filtering.png)

*qc before filtering: Histograms of genes detected per cell, total signal per cell and mitochondrial percentage for 1,937 cells before filtering; red lines mark the cut-offs (200 genes minimum). The mitochondrial panel is empty because no mitochondrial percentage could be computed.*

## Folder map

- `checkpoints/`: one restartable file per stage (`NN_stage.h5ad`, or `.json` for fetch)
- `reports/`: `stage.md` (read this) and `stage.json` (same content for programs) for every stage
- `figures/`: plots, with `figure_legends.md` beside each stage's figures
- `manifest.json`: parameters hash, input and output hashes, timing, memory, versions per stage
- `RUN_SUMMARY.md`: this file, rebuilt after every stage

## Going back and resuming

- Rerun the notebook: finished stages whose parameters and input checkpoint are unchanged load from disk; the first changed stage and everything after it recompute.
- Force one stage: `run_stage("<stage>", ..., force=True)`; forget a stage and all later ones: `wf.invalidate("<stage>")`.
- A stage that reported FAIL stops the run with its checkpoint kept; fix the cause, or continue on purpose with `wf.accept("<stage>", "reason")` (the reason is recorded above).
- Stage state tracks parameters and input files, not code: after editing code, `invalidate` the stage.
