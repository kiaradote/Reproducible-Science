"""Tests for scrna-normalization-features (checkpoint 3). Synthetic data for the logic, a real count matrix (GSE84133 donor 1, 400-cell subsample) for the counts path."""
import numpy as np, pandas as pd, scipy.sparse as sp, pytest, anndata as ad
import scrna_features_kernel as K

PATH = "testdata/GSM2230757_human1_umifm_counts.csv.gz"

def status(recs, name):
    return [r["status"] for r in recs if r["check"] == name][0]

def synth(n_cells=300, n_genes=1500, groups=3, seed=0, depth_by_sample=None, names=None):
    """Counts with group-specific genes (so there is real variable-gene structure) and per-cell depth differences."""
    rng = np.random.default_rng(seed)
    base = rng.gamma(0.4, 1.0, n_genes) + 0.02
    g = np.repeat(np.arange(groups), n_cells // groups)
    mu = np.tile(base, (len(g), 1))
    for k in range(groups):
        mu[g == k, k * 100:(k + 1) * 100] *= 10
    depth = rng.gamma(6, 1 / 6, (len(g), 1))
    sample = np.array([f"S{i % 3}" for i in range(len(g))])
    if depth_by_sample:
        for s, f in depth_by_sample.items():
            depth[sample == s] *= f
    X = rng.poisson(mu * depth * 2).astype(np.float32)
    a = ad.AnnData(sp.csr_matrix(X))
    a.var_names = names or [f"G{i}" for i in range(n_genes)]
    a.obs_names = [f"c{i}" for i in range(len(g))]
    a.obs["sample"] = sample
    return a

def to_log2tpm(a):
    X = a.X.toarray()
    tpm = X / X.sum(1, keepdims=True) * 1e6
    b = a.copy(); b.X = sp.csr_matrix(np.log2(tpm + 1).astype(np.float32)); return b

# ---------- counts path ----------
def test_counts_path_normalizes_and_keeps_counts():
    a = synth(); raw = a.X.copy()
    a, recs = K.feat_run(a, "counts", n_top_genes=300, sample_col="sample")
    assert "counts" in a.layers and (a.layers["counts"] != raw).nnz == 0
    assert a.uns["log1p"]["base"] is None and a.uns["feature_selection"]["flavor"] == "seurat_v3"
    assert status(recs, "per-cell totals comparable") == "PASS" and status(recs, "sparsity preserved") == "PASS"
    assert status(recs, "HVG count") == "PASS" and int(a.var["highly_variable"].sum()) == 300
    assert all(r["status"] != "FAIL" for r in recs)

def test_real_counts_subsample():
    cols = pd.read_csv(PATH, nrows=1, index_col=0).columns
    df = pd.read_csv(PATH, index_col=0, dtype={c: np.int32 for c in cols[2:]}, nrows=1937)
    genes = list(cols[2:])
    idx = np.arange(0, df.shape[0], max(1, df.shape[0] // 400))
    a = ad.AnnData(sp.csr_matrix(df[genes].to_numpy(dtype=np.float32)[idx]))
    a.var_names = genes; a.obs_names = list(df.index[idx])
    a = a[:, np.asarray(a.X.getnnz(axis=0) >= 3)].copy()
    a, recs = K.feat_run(a, "counts", n_top_genes=2000)
    assert a.uns["feature_selection"]["flavor"] == "seurat_v3" and int(a.var["highly_variable"].sum()) == 2000
    assert all(r["status"] != "FAIL" for r in recs), recs

# ---------- log2(TPM+1): the log base must reach the variable-gene step ----------
def test_log2tpm_leaves_x_unchanged_and_records_base():
    a = to_log2tpm(synth()); X0 = a.X.copy()
    a, recs = K.feat_run(a, "log2tpm", n_top_genes=300, sample_col="sample")
    assert (a.X != X0).nnz == 0 and a.uns["log1p"]["base"] == 2.0 and a.uns["feature_selection"]["flavor"] == "seurat"
    assert status(recs, "per-cell totals comparable") == "PASS"

def test_log2_base_gives_same_genes_as_explicit_natural_log():
    a = to_log2tpm(synth(seed=1))
    ours, _ = K.feat_run(a.copy(), "log2tpm", n_top_genes=300)
    nat = a.copy(); nat.X = nat.X * np.log(2); K.feat_normalize(nat, "log1p_cp10k")
    K.feat_select(nat, 300)
    assert (ours.var["highly_variable"].to_numpy() == nat.var["highly_variable"].to_numpy()).all()

def test_linear_tpm_is_logged():
    a = synth(); X = a.X.toarray(); tpm = X / X.sum(1, keepdims=True) * 1e6
    a.X = sp.csr_matrix(tpm.astype(np.float32))
    a, recs = K.feat_run(a, "linear_tpm", n_top_genes=300)
    assert "input" in a.layers and status(recs, "input values kept") == "PASS"

def test_scaled_data_refused_and_untouched():
    a = synth(); X = a.X.toarray(); z = (X - X.mean(0)) / (X.std(0) + 1e-9); a.X = sp.csr_matrix(z.astype(np.float32))
    before = a.X.copy()
    a, recs = K.feat_run(a, "scaled_or_unknown", n_top_genes=300)
    assert len(recs) == 1 and recs[0]["status"] == "FAIL" and (a.X != before).nnz == 0 and "highly_variable" not in a.var

# ---------- checks that should warn ----------
def test_depth_spread_warns():
    a = synth(depth_by_sample={"S0": 12.0})
    _, recs = K.feat_run(a, "counts", n_top_genes=300, sample_col="sample")
    assert status(recs, "depth spread across samples") == "WARN"

def test_dominant_gene_cells_warn():
    a = synth(); X = a.X.toarray(); X[:20, 0] = X[:20].sum(1) * 6; a.X = sp.csr_matrix(X)
    _, recs = K.feat_run(a, "counts", n_top_genes=300)
    assert status(recs, "dominant-gene cells") == "WARN"

def test_technical_and_lowdetection_genes_warn():
    names = [f"RPL{i}" if i < 100 else f"G{i}" for i in range(1500)]
    a = synth(names=names)
    K.feat_normalize(a, "counts")
    a.var["highly_variable"] = False; a.var.iloc[:100, a.var.columns.get_loc("highly_variable")] = True
    recs = K.feat_check_selection(a, "seurat_v3", "counts", n_top_genes=100)
    assert status(recs, "technical genes among HVGs") == "WARN"
    # genes seen in under 1% of cells
    rare = synth(); rare.X = sp.csr_matrix(rare.X.toarray()); X = rare.X.toarray(); X[:, 1000:1100] = 0; X[0, 1000:1100] = 1; rare.X = sp.csr_matrix(X)
    K.feat_normalize(rare, "counts"); rare.var["highly_variable"] = False
    rare.var.iloc[1000:1100, rare.var.columns.get_loc("highly_variable")] = True
    assert status(K.feat_check_selection(rare, "seurat_v3", "counts", n_top_genes=100), "low-detection HVGs") == "WARN"

def test_unstable_selection_warns_on_pure_noise():
    rng = np.random.default_rng(3)
    a = ad.AnnData(sp.csr_matrix(rng.poisson(1.0, (400, 1500)).astype(np.float32)))
    a.var_names = [f"G{i}" for i in range(1500)]
    a, recs = K.feat_run(a, "counts", n_top_genes=300)
    assert status(recs, "HVG stability (two halves)") == "WARN"

def test_structured_data_is_stable():
    a = synth(n_cells=600)
    _, recs = K.feat_run(a, "counts", n_top_genes=300)
    assert status(recs, "HVG stability (two halves)") == "PASS"

def test_batch_key_adds_robustness_check():
    a = synth(n_cells=600)
    a, recs = K.feat_run(a, "counts", n_top_genes=300, batch_key="sample")
    assert "highly_variable_nbatches" in a.var and status(recs, "batch robustness") in ("PASS", "WARN")
    _, recs2 = K.feat_run(synth(n_cells=600), "counts", n_top_genes=300)
    assert not [r for r in recs2 if r["check"] == "batch robustness"]

# ---------- report ----------
def test_report_reason_and_not_run(tmp_path):
    a = synth(); a, recs = K.feat_run(a, "counts", n_top_genes=300)          # no sample_col, no batch_key
    overall = K.write_feat_report(recs, str(tmp_path / "r"))
    md = (tmp_path / "r.md").read_text()
    assert overall in ("PASS", "WARN") and "**Why:**" in md and "**Not run:**" in md
    assert "- depth spread across samples:" in md and "- batch robustness:" in md
    bad = K.write_feat_report(K.feat_run(synth(), "scaled_or_unknown")[1], str(tmp_path / "s"))
    assert bad == "FAIL" and "method matches scale" in (tmp_path / "s.md").read_text()
