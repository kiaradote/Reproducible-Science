# Restartable scRNA-seq workflow (checkpoints 0-8)

See `SCRNA_WORKFLOW_SUMMARY.md` for the one-page overview.

Open `scrna_checkpoints.ipynb` from this folder (it expects `scrna_workflow/` next to it). Set the environment variable `NB_DATASET` to
`baron` (GSE84133 donor 1, UMI counts, small), `kang` (GSE96583 batch 2, 10x PBMC, control vs interferon-beta, 8 donors) or `sade_feldman`
(GSE120575, Smart-seq2 TPM, large). Run all cells, or change one parameter and rerun: finished stages whose inputs and parameters are unchanged load from
their checkpoint file instead of recomputing; everything after a changed stage reruns.

## Stages (each saves `runs/<dataset>/checkpoints/NN_<stage>.h5ad` and writes `reports/<stage>.md` + `.json`)
0 fetch, 1 assemble and integrity checks, 2 QC and filtering, 3 normalisation and variable genes, 4 PCA / Harmony / neighbours / UMAP,
5 clustering with a stability-chosen resolution, 6 cell type annotation (marker panel + optional reference, confidence tiers),
7 condition comparison (design gate, abundance, pseudobulk differential expression with null and consistency checks, within-type state change, pre-ranked pathways; set `condition`/`sample` in the stage parameters, `None` skips it),
8 biological interpretation and validation (pathways, gene programs, expectations, held-out-sample signature test, robustness and sensitivity; uses stage 7's differential expression when present).

Every stage ends with PASS / WARN / FAIL, the reason, a table of the checks that ran (with what each one tests) and a "Not run" list. A FAIL stops the run
until it is fixed or accepted with a written reason (the `ACCEPT` dict in the notebook). After each stage `runs/<dataset>/RUN_SUMMARY.md` (and `run_summary.json`)
is rewritten: a quality check table at the top (Good / Usable with caveats / Accepted despite a failed check / Not reliable / Not applicable / Not done, per stage and overall, with the command to go back), verdict per stage, all check tables, parameters, versions, timing and memory, input SHA-256, figures with one-line captions, folder map, how to resume.

## Inputs you edit (copy `inputs/<dataset>/*` into `runs/<dataset>/`)
- `marker_panel.csv`  cell_type, positive, negative, source (genes separated by `;`); editing it reruns annotation
- `expectations.csv`  known biology the data must reproduce, including negative controls (kinds: gene_up, set_up, fraction_range, set_enriched)
- `gene_sets/*.gmt`   optional; with `gene_sets_kegg: "hsa"` the KEGG human pathways are fetched once from rest.kegg.jp (network needed)
Panel, gene-set and expectation file contents are hashed into the stage parameters.

## Layout
- `scrna_checkpoints.ipynb`  one section per stage
- `scrna_workflow/ckpt.py`   checkpoint engine (save, skip current stages, stale detection, FAIL gate, accept, invalidate, run summary)
- `scrna_workflow/captions.py`  one-line figure captions written beside every figure set (`figure_captions.json`)
- `scrna_workflow/loaders.py`, `adapters.py`  chunked sparse readers and one adapter per dataset
- `scrna_workflow/stages.py`  stage functions; `scrna_workflow/scrna_*_kernel.py` are copies of the published skills' kernel.py
- `tests/`  pytest suites; `run_nb.py`, `run_stress.py` run the notebook cells without Jupyter and a memory-guarded large-data run

Environment: python 3.11, scanpy 1.10, anndata 0.10, numpy<2, pandas<3, scipy, scikit-learn, statsmodels, matplotlib, python-igraph, harmonypy, psutil, pytest.

## Known limits
- Stage state tracks parameters and input checkpoint hashes, not code. After editing a stage function, call `wf.invalidate("<stage>")` (a FAILED stage reruns by itself).
- Recorded peak memory is the process-lifetime peak, so it is stage-specific only in a fresh process.
- Run the tests inside one process (`pytest.main([...])`) on Windows: separate Python subprocesses crashed in scikit-learn's thread pool on the development machine.
  `scanpy.tl.leiden` was about 1000 times slower on that machine, so clustering calls igraph directly with a seeded random generator.
- On the development machine importing torch (needed by harmonypy) failed in some kernels (DLL load error); the Harmony test is the only one affected.
- The notebook was executed cell by cell in plain Python, not opened in Jupyter.
- Stage 8 pathway analysis is over-representation of marker genes (and pre-ranked enrichment once stage 7 exists), not a measure of pathway activity. KEGG is fetched; Reactome could not be downloaded automatically.
- Pseudobulk differential expression uses pydeseq2 single-process (`n_cpus=1`); its default worker pool crashed on Windows. Pathways: pre-ranked enrichment uses 1,500 permutations per set-size group.
- Stage 7 cannot separate a treatment effect from a library effect when each condition is its own library (Kang PBMC); it warns instead of failing.
- Large dataset (GSE120575): stages 3 and later have not been run on it (compute was limited).
- Tests need testdata/GSM2230757_human1_umifm_counts.csv.gz for the full-data integrity test; the pipeline and workflow suites do not.

## Repository layout
- `scrna_checkpoints.ipynb`, `run_nb.py`, `run_stress.py`  the notebook and two runners that execute its cells without Jupyter
- `scrna_workflow/`  the engine (`ckpt.py`), adapters, loaders, captions, stage functions and the eight kernel files
- `skills/`  a copy of the nine published Claude skills (`SKILL.md` and `kernel.py` each): the master rule `scrna-seq-workflow` and one per checkpoint
- `tests/`, `conftest.py`, `pytest.ini`  run `pytest` from this folder (about 5 minutes; 149 passed and 4 skipped on the development machine). A few tests use the public file `testdata/GSM2230757_human1_umifm_counts.csv.gz`; `python fetch_testdata.py` downloads it, and those tests are skipped while it is absent
- `inputs/<dataset>/`  editable marker panels and expectation files; copy them into `runs/<dataset>/` before the first run
- `runs/`  created by the notebook (checkpoints, reports, figures, `RUN_SUMMARY.md`); ignored by git

## Setup
`conda env create -f environment.yml`, then `conda activate scrna-repro`, then `set NB_DATASET=kang` (Windows) or `export NB_DATASET=kang` and run `python run_nb.py` or open the notebook.
The environment file records the versions the workflow was developed and tested with.

## Keeping the skills and this folder in step
The published skills are the source new Claude conversations use. `skills/*/kernel.py` and `scrna_workflow/scrna_*_kernel.py` are copies of the same code; if you change one, change the other (or republish the skill).
