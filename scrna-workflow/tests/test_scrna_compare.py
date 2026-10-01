"""Tests for scrna-condition-comparison (checkpoint 7). Synthetic paired data with planted changes in share, in gene expression and in state."""
import os, json
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, pytest
import scrna_compare_kernel as C

GENES = 300

def make_cmp(n_donors=6, per=40, seed=0, shares_stim=(0.4, 0.25, 0.35), effect=3.0, paired=True, batch="mixed", donor_sd=0.15):
    """Types A, B, C. Stim: genes 0-14 are multiplied by `effect` in A only; the share of B falls and C rises. Donors differ in baseline."""
    rng = np.random.default_rng(seed)
    rows, obs = [], []
    shares = {"ctrl": (0.4, 0.4, 0.2), "stim": shares_stim}
    for d in range(n_donors):
        base = np.exp(rng.normal(0, donor_sd, GENES))
        for cond in ("ctrl", "stim"):
            sh = np.array(shares[cond]) * np.exp(rng.normal(0, 0.1, 3))
            sh = sh / sh.sum()
            n_tot = per * 3
            ty = rng.choice(3, size=n_tot, p=sh)
            for t in ty:
                mu = base * 1.0
                mu[20 + 10 * t:30 + 10 * t] *= 6.0
                if cond == "stim" and t == 0:
                    mu[:15] *= effect
                rows.append(rng.poisson(mu))
                sample = f"d{d}" if paired else f"d{d}_{cond}"
                lib = {"mixed": f"lib{d % 2}", "split": "libA" if cond == "ctrl" else "libB"}[batch]
                obs.append((sample, cond, "ABC"[t], lib, f"{'ABC'[t]}_{cond}"))
    X = np.array(rows, dtype=np.float32)
    a = ad.AnnData(sp.csr_matrix(X))
    a.var_names = [f"G{i}" for i in range(GENES)]
    a.obs_names = [f"c{i}" for i in range(len(obs))]
    a.obs["sample"], a.obs["condition"], a.obs["cell_type"], a.obs["library"], a.obs["cluster"] = zip(*obs)
    for c in ("sample", "condition", "cell_type", "library", "cluster"):
        a.obs[c] = a.obs[c].astype(str)
    a.layers["counts"] = a.X.copy()
    lib_size = np.asarray(a.X.sum(axis=1)).ravel()
    a.X = sp.csr_matrix(np.log1p(X / lib_size[:, None] * 1e4))
    Z = a.X.toarray()
    Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-6)
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    a.obsm["X_pca"] = (U[:, :20] * S[:20]).astype(np.float32)
    return a

@pytest.fixture(scope="module")
def data():
    return make_cmp()

def status(recs, name):
    return [r["status"] for r in recs if r["check"] == name][0]

def test_bh_and_ols_match_reference():
    assert np.allclose(C.cmp_bh([0.01, 0.04, 0.03, 0.2]), [0.04, 0.0533333, 0.0533333, 0.2], atol=1e-4)
    assert np.isnan(C.cmp_bh([0.01, np.nan])[1])
    rng = np.random.default_rng(1)
    units = pd.DataFrame({"sample": list("aabbccdd"), "treated": [0, 1] * 4})
    X = C.cmp_matrix(units, True)
    Y = rng.normal(size=(8, 5))
    Y[:, 0] += 8 * units["treated"].to_numpy()
    b, t, p = C.cmp_ols(Y, X)
    ref = np.linalg.lstsq(X, Y, rcond=None)[0][1]
    assert np.allclose(b, ref) and p[0] < 0.01 and (p[1:] > 0.02).all()

def test_design_passes_for_a_paired_replicated_experiment(data):
    recs, des = C.cmp_design(data, "condition", "sample", batch_col="library")
    assert status(recs, "replicates") == "PASS" and status(recs, "pairing") == "PASS" and status(recs, "condition vs batch") == "PASS"
    assert des["paired"] and des["ref"] == "ctrl" and des["test"] == "stim" and set(des["testable"]) == {"A", "B", "C"}

def test_design_fails_without_replicates_and_warns_with_two(data):
    one = data[data.obs["sample"] == "d0"].copy()
    recs, _ = C.cmp_design(one, "condition", "sample")
    assert status(recs, "replicates") == "FAIL"
    two = data[data.obs["sample"].isin(["d0", "d1"])].copy()
    recs, _ = C.cmp_design(two, "condition", "sample")
    assert status(recs, "replicates") == "WARN"

def test_condition_confounded_with_library_is_flagged():
    a = make_cmp(batch="split")
    recs, des = C.cmp_design(a, "condition", "sample", batch_col="library")
    assert status(recs, "condition vs batch") == "WARN" and des["batch_confounded"]

def test_unpaired_design_is_detected():
    a = make_cmp(paired=False)
    recs, des = C.cmp_design(a, "condition", "sample")
    assert not des["paired"] and status(recs, "pairing") == "INFO"

def test_cell_types_and_clusters_that_follow_the_condition_are_flagged(data):
    a = data.copy()
    a.obs["cluster"] = a.obs["condition"] + "_" + a.obs["cell_type"]
    a.obs.loc[a.obs["cell_type"] == "C", "cell_type"] = np.where(a.obs.loc[a.obs["cell_type"] == "C", "condition"] == "stim", "C", "Cx")
    recs, des = C.cmp_design(a, "condition", "sample")
    assert status(recs, "clusters follow condition") == "WARN" and status(recs, "cell types follow condition") == "WARN"
    assert set(des["followers"]) == {"C", "Cx"} and {"C", "Cx"} <= set(des["excluded"]) | set(des["testable"])

def test_wrong_inputs_fail_cleanly(data):
    assert C.cmp_design(data, "nope", "sample")[0][0]["status"] == "FAIL"
    a = data.copy()
    a.obs["condition"] = np.where(np.arange(a.n_obs) % 3 == 0, "x", a.obs["condition"])
    assert C.cmp_design(a, "condition", "sample")[0][0]["status"] == "FAIL"
    assert C.cmp_design(a, "condition", "sample", groups=("ctrl", "stim"))[0][0]["check"] == "groups"

def test_abundance_finds_planted_shares_and_ignores_unchanged_types(data):
    _, des = C.cmp_design(data, "condition", "sample")
    ab = C.cmp_abundance(data, des).set_index("type")
    assert ab.loc["C", "log2_ratio"] > 0.5 and ab.loc["B", "log2_ratio"] < -0.3
    assert ab.loc["C", "fdr"] < 0.05 and ab.loc["B", "fdr"] < 0.05

def test_abundance_does_not_invent_changes_when_shares_are_equal():
    a = make_cmp(shares_stim=(0.4, 0.4, 0.2), seed=3)
    _, des = C.cmp_design(a, "condition", "sample")
    ab = C.cmp_abundance(a, des)
    assert (ab["p_quasibinomial"] > 1e-4).all() and not ab["agree"].any()          # regression: an unscaled dispersion once gave p of 1e-248

def test_pseudobulk_sums_counts_and_drops_small_units(data):
    _, des = C.cmp_design(data, "condition", "sample")
    m = (data.obs["cell_type"] == "A").to_numpy()
    mat, units, kind = C.cmp_pseudobulk(data, m, des)
    assert kind == "counts" and len(units) == 12
    u = units.iloc[0]
    sel = m & (data.obs["sample"] == u["sample"]).to_numpy() & (data.obs["condition"] == u["condition"]).to_numpy()
    assert np.allclose(mat.iloc[0].to_numpy(), np.asarray(data.layers["counts"][sel].sum(axis=0)).ravel()) and u["n_cells"] == sel.sum()
    d2 = data.copy()
    del d2.layers["counts"]
    mat2, units2, kind2 = C.cmp_pseudobulk(d2, m, des)
    assert kind2 == "mean log expression" and mat2.to_numpy().max() < 20

@pytest.mark.parametrize("method", ["deseq2", "ols"])
def test_de_finds_planted_genes_with_the_right_sign_and_only_in_the_right_type(data, method):
    _, des = C.cmp_design(data, "condition", "sample")
    res = {}
    for t in ("A", "B"):
        mat, units, kind = C.cmp_pseudobulk(data, (data.obs["cell_type"] == t).to_numpy(), des)
        res[t] = C.cmp_de_type(mat, units, des, kind, method=method)
    sa = res["A"][res["A"]["significant"]]
    planted = {f"G{i}" for i in range(15)}
    assert len(planted & set(sa["gene"])) >= 12 and (sa[sa["gene"].isin(planted)]["log2fc"] > 0.8).all()
    assert res["B"]["significant"].sum() <= 3
    assert res["A"]["method"].iloc[0].startswith("pydeseq2" if method == "deseq2" else "linear model")

def test_label_swap_null_and_sample_consistency(data):
    _, des = C.cmp_design(data, "condition", "sample")
    mat, units, kind = C.cmp_pseudobulk(data, (data.obs["cell_type"] == "A").to_numpy(), des)
    real, nulls = C.cmp_null_swaps(mat, units, des, kind, n_swaps=10)
    assert real > 10 and len(nulls) >= 5 and np.median(nulls) <= 3
    d = C.cmp_de_type(mat, units, des, kind, method="ols")
    sg = d[d["significant"]]
    assert C.cmp_consistency(mat, units, des, kind, list(sg["gene"]), sg["log2fc"].to_numpy()) >= 0.8

def test_state_change_found_in_the_changed_type_and_not_in_the_unchanged_one(data):
    _, des = C.cmp_design(data, "condition", "sample")
    a = C.cmp_state(data, des, "A", n_pcs=20)
    b = C.cmp_state(data, des, "B", n_pcs=20)
    assert a["median_auc"] > 0.9 and a["label"] == "clear shift" and a["n_folds"] == 6
    assert b["median_auc"] < 0.7 and b["label"] != "clear shift"

def test_enrichment_recovers_the_planted_set(data):
    _, des = C.cmp_design(data, "condition", "sample")
    mat, units, kind = C.cmp_pseudobulk(data, (data.obs["cell_type"] == "A").to_numpy(), des)
    d = C.cmp_de_type(mat, units, des, kind, method="ols")
    sets = {"planted": {f"G{i}" for i in range(15)}, "other": {f"G{i}" for i in range(100, 130)}, "other2": {f"G{i}" for i in range(150, 185)}, "o3": {f"G{i}" for i in range(200, 240)},
            "o4": {f"G{i}" for i in range(240, 280)}, "o5": {f"G{i}" for i in range(50, 90)}}
    e = C.cmp_enrich(d, sets, {f"G{i}" for i in range(GENES)}, n_perm=600)
    assert "planted" in set(e["term"]) and e.set_index("term").loc["planted", "nes"] > 0

def test_full_run_writes_tables_figures_and_a_question_first_report(data, tmp_path):
    sets = {"planted": {f"G{i}" for i in range(15)}, "other": {f"G{i}" for i in range(100, 130)}, "other2": {f"G{i}" for i in range(150, 185)}, "o3": {f"G{i}" for i in range(200, 240)},
            "o4": {f"G{i}" for i in range(240, 280)}, "o5": {f"G{i}" for i in range(50, 90)}}
    a, recs, det = C.cmp_run(data.copy(), "condition", "sample", batch_col="library", gene_sets=sets, method="ols", n_swaps=8, fig_dir=str(tmp_path))
    for c in ("replicates", "abundance", "differential expression", "null check", "sample consistency", "state change", "pathway enrichment"):
        assert any(r["check"] == c for r in recs), c
    assert not any(r["status"] == "FAIL" for r in recs)
    for f in ("cmp_abundance.png", "cmp_de.png", "cmp_state.png", "figure_captions.json", "figure_legends.md", "abundance.csv", "differential_expression.csv", "state_change.csv"):
        assert os.path.getsize(tmp_path / f) > 0
    assert "comparison_de" in a.uns and {"type", "gene", "stat", "padj"} <= set(a.uns["comparison_de"].columns) and a.uns["comparison"]["treatment"] == "stim"
    caps = json.load(open(tmp_path / "figure_captions.json"))
    assert "cmp_abundance.png" in caps and "cmp_de.png" in caps
    overall = C.write_cmp_report(recs, det, [{"stage": "qc", "kind": "WARN", "text": "x"}], str(tmp_path / "rep"), question="What does stim change?")
    text = (tmp_path / "rep.md").read_text(encoding="utf8")
    assert overall in ("PASS", "WARN") and "**Question:** What does stim change?" in text and "## 1. Is the design able to answer the question?" in text
    assert "Which genes change within each cell type?" in text and "**Not run:**" in text and "transcription factor activity" in text and "Caveats carried" in text

def test_run_stops_at_the_design_gate_when_there_are_no_replicates(data, tmp_path):
    one = data[data.obs["sample"] == "d0"].copy()
    a, recs, det = C.cmp_run(one, "condition", "sample", fig_dir=str(tmp_path))
    assert any(r["status"] == "FAIL" for r in recs) and "differential expression" not in [r["check"] for r in recs] and "comparison_de" not in a.uns
    assert C.write_cmp_report(recs, det, None, str(tmp_path / "r")) == "FAIL"
