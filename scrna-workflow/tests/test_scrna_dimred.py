"""Tests for scrna-dimensionality-reduction (checkpoint 4). Synthetic data with known structure for the logic; real counts for the pipeline."""
import os
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, pytest
import scrna_features_kernel as F
import scrna_dimred_kernel as K

PATH = "testdata/GSM2230757_human1_umifm_counts.csv.gz"

def status(recs, name):
    return [r["status"] for r in recs if r["check"] == name][0]

def prepared(n_cells=600, n_genes=1500, groups=3, batch_shift=0.0, seed=0, depth_effect=False, n_top=400):
    """Counts with `groups` real cell groups and 3 batches (optional batch shift on 30% of genes), normalized with HVGs selected."""
    rng = np.random.default_rng(seed)
    g = np.repeat(np.arange(groups), n_cells // groups)
    batch = np.array([f"B{i % 3}" for i in range(len(g))])
    base = rng.gamma(0.4, 1.0, n_genes) + 0.05
    mu = np.tile(base, (len(g), 1))
    for k in range(groups):
        mu[g == k, k * 120:(k + 1) * 120] *= 8
    if batch_shift:
        shifted = np.arange(n_genes - 450, n_genes)
        for b, f in (("B1", batch_shift), ("B2", 1 / batch_shift)):
            mu[np.ix_(batch == b, shifted)] *= f
    depth = rng.gamma(6, 1 / 6, (len(g), 1))
    if depth_effect:
        depth = depth * np.where(g[:, None] == 0, 3.0, 1.0)
    X = rng.poisson(mu * depth * 2).astype(np.float32)
    a = ad.AnnData(sp.csr_matrix(X))
    a.var_names = [f"G{i}" for i in range(n_genes)]
    a.obs_names = [f"c{i}" for i in range(len(g))]
    a.obs["batch"], a.obs["group"] = batch, g.astype(str)
    a.obs["n_genes"] = np.diff(a.X.indptr)
    a.obs["total"] = np.asarray(a.X.sum(1)).ravel()
    F.feat_run(a, "counts", n_top_genes=n_top)
    return a

def test_pcs_chosen_from_data_and_outputs_exist(tmp_path):
    a = prepared()
    a, recs = K.dr_run(a, fig_dir=str(tmp_path))
    d = a.uns["dimred"]
    assert 5 <= d["n_pcs"] <= d["nmax"] and d["n_null"] >= 2          # 3 real groups: at least 2 informative PCs
    assert a.obsm["X_pca"].shape[1] == d["nmax"] and a.obsm["X_umap"].shape == (a.n_obs, 2) and np.isfinite(a.obsm["X_umap"]).all()
    assert a.varm["PCs"].shape == (a.n_vars, d["nmax"]) and "connectivities" in a.obsp
    assert status(recs, "neighbor graph fragments") == "PASS" and status(recs, "UMAP finite") == "PASS"
    assert status(recs, "variance threshold reached") == "PASS" and status(recs, "chosen PCs beat shuffled data") == "PASS"
    assert all(r["status"] != "FAIL" for r in recs)

def test_threshold_helper_is_monotone_and_exact():
    eig = np.array([10, 5, 2, 1, 0.5]); null = np.full(5, 0.6); ratio = np.array([.4, .3, .1, .05, .02])
    n = lambda T, b: K.dr_pcs_for_threshold(eig, null, ratio, T, b)[0]
    assert [n(.5, "signal"), n(.9, "signal"), n(.99, "signal")] == [1, 3, 4]       # signal = eig - null = 9.4, 4.4, 1.4, 0.4, 0
    assert [n(.5, "total"), n(.75, "total"), n(.85, "total")] == [2, 3, 4]      # cumulative .4 .7 .8 .85 .87
    assert K.dr_pcs_for_threshold(eig, null, ratio, .99, "total")[1:2] == (False,)  # 0.97 total never reaches 0.99
    assert K.dr_pcs_for_threshold(eig, eig + 1, ratio, .9, "signal")[0] == 0         # nothing beats the null
    late = np.array([10, 5, 0.5, 0.5, 3.0])                                           # a later PC above the null must not count after a failure
    assert K.dr_pcs_for_threshold(late, np.full(5, 0.6), ratio, .99, "signal")[0] == 2

def test_higher_threshold_never_fewer_pcs():
    ns = []
    for T in (0.5, 0.8, 0.95):
        a = prepared(groups=6, n_cells=900, seed=2)
        K.dr_run(a, variance_threshold=T, min_pcs=2); ns.append(a.uns["dimred"]["n_pcs"])
    assert ns == sorted(ns) and ns[0] < ns[-1]

def test_more_structure_gives_more_pcs():
    few = prepared(groups=2, seed=1); K.dr_run(few)
    many = prepared(groups=6, seed=1, n_cells=900); K.dr_run(many)
    assert many.uns["dimred"]["n_null"] > few.uns["dimred"]["n_null"]

def test_pure_noise_falls_back_to_minimum_with_warning():
    rng = np.random.default_rng(5)
    a = ad.AnnData(sp.csr_matrix(rng.poisson(1.0, (500, 1200)).astype(np.float32)))
    a.var_names = [f"G{i}" for i in range(1200)]
    F.feat_run(a, "counts", n_top_genes=300)
    a, recs = K.dr_run(a, min_pcs=5)
    assert a.uns["dimred"]["n_pcs"] == 5 and status(recs, "number of PCs chosen from the data") == "WARN"

def test_threshold_not_reached_warns():
    a = prepared(groups=6, n_cells=900, seed=2)
    _, recs = K.dr_run(a, max_pcs=4, min_pcs=2, variance_basis="total", variance_threshold=0.99)
    assert status(recs, "variance threshold reached") == "WARN"

def test_total_basis_can_keep_pcs_that_do_not_beat_noise():
    a = prepared(groups=2, seed=1)
    _, recs = K.dr_run(a, variance_basis="total", variance_threshold=0.3)
    assert status(recs, "chosen PCs beat shuffled data") in ("PASS", "WARN")

def test_same_seed_gives_identical_result():
    a1, a2 = prepared(seed=3), prepared(seed=3)
    K.dr_run(a1, seed=7); K.dr_run(a2, seed=7)
    assert a1.uns["dimred"]["n_pcs"] == a2.uns["dimred"]["n_pcs"]
    assert np.allclose(a1.obsm["X_pca"], a2.obsm["X_pca"]) and np.allclose(a1.obsm["X_umap"], a2.obsm["X_umap"])

def test_depth_driven_pc_is_flagged():
    a = prepared(depth_effect=True, seed=4)
    _, recs = K.dr_run(a, tech_var_limit=0.1)
    rec = [r for r in recs if r["check"] == "PCs driven by technical covariates"][0]
    assert rec["status"] == "WARN" and "PCs" in rec["detail"]

def test_batch_effect_detected_and_corrected():
    pytest.importorskip("harmonypy")
    a = prepared(batch_shift=3.0, seed=6, n_cells=900)
    a, recs = K.dr_run(a, batch_key="batch", correct_batch="auto")
    d = a.uns["dimred"]
    assert d["corrected"] and "X_pca_harmony" in a.obsm and a.uns["dimred"]["rep"] == "X_pca_harmony"
    assert d["share_after"] < d["share_before"] and d["mix_after"] >= d["mix_before"]
    assert status(recs, "batch correction effect") in ("PASS", "WARN")

def test_no_correction_when_not_needed_or_declined():
    clean = prepared(seed=8)
    a, recs = K.dr_run(clean, batch_key="batch", correct_batch="auto")
    assert not a.uns["dimred"]["corrected"] and status(recs, "batch structure in PCs") == "PASS"
    shifted = prepared(batch_shift=3.0, seed=6, n_cells=900)
    b, recs = K.dr_run(shifted, batch_key="batch", correct_batch=False)
    assert not b.uns["dimred"]["corrected"] and status(recs, "batch structure in PCs") == "WARN"

def test_forced_correction_without_batch_column_fails():
    _, recs = K.dr_run(prepared(seed=9), correct_batch=True)
    assert status(recs, "batch correction effect") == "FAIL"

def test_missing_features_is_refused():
    a = ad.AnnData(sp.csr_matrix(np.ones((10, 10), dtype=np.float32)))
    _, recs = K.dr_run(a)
    assert recs[0]["status"] == "FAIL"

def test_figures_and_legends_with_computed_numbers(tmp_path):
    a = prepared(seed=10)
    a, _ = K.dr_run(a, batch_key="batch", correct_batch=False, fig_dir=str(tmp_path))
    for f in ("pca_variance.png", "pca_covariates.png", "pca_scatter.png", "umap.png", "figure_legends.md"):
        assert (tmp_path / f).exists() and (tmp_path / f).stat().st_size > 1000
    text = (tmp_path / "figure_legends.md").read_text()
    n = a.uns["dimred"]["n_pcs"]
    assert f"{n} PCs are used downstream" in text and "not quantitative" in text and "batch (batch)" in text
    assert "UMAP" in text and "dashed line" in text

def test_report_has_reason_and_not_run(tmp_path):
    a = prepared(seed=11)
    _, recs = K.dr_run(a)
    overall = K.write_dr_report(recs, str(tmp_path / "r"))
    md = (tmp_path / "r.md").read_text()
    assert overall in ("PASS", "WARN") and "**Why:**" in md and "- batch structure in PCs:" in md and "- batch correction effect:" in md

def test_real_counts_pipeline(tmp_path):
    cols = pd.read_csv(PATH, nrows=1, index_col=0).columns
    df = pd.read_csv(PATH, index_col=0, dtype={c: np.int32 for c in cols[2:]})
    genes = list(cols[2:])
    a = ad.AnnData(sp.csr_matrix(df[genes].to_numpy(dtype=np.float32)))
    a.var_names = genes
    a = a[:, np.asarray(a.X.getnnz(axis=0) >= 3)].copy()
    a.obs["n_genes"] = np.diff(a.X.indptr); a.obs["total"] = np.asarray(a.X.sum(1)).ravel()
    F.feat_run(a, "counts", n_top_genes=2000)
    a, recs = K.dr_run(a, fig_dir=str(tmp_path))
    d = a.uns["dimred"]
    assert 5 <= d["n_pcs"] <= 50 and all(r["status"] != "FAIL" for r in recs)
    assert (tmp_path / "umap.png").exists()
