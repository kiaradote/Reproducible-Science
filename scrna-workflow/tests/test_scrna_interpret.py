"""Tests for scrna-interpretation-validation (checkpoint 8). Synthetic data with planted structure for the logic; real runs for the pipeline."""
import os, json
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, pytest
import scrna_interpret_kernel as I
from test_scrna_annotate import make, PANEL, status

BG = {f"G{i}" for i in range(200)}
SETS = {"setA": {f"G{i}" for i in range(0, 20)}, "setB": {f"G{i}" for i in range(20, 40)}, "setA_copy": {f"G{i}" for i in range(0, 19)}, "setC": {f"G{i}" for i in range(100, 120)}}

def test_bh_matches_known_values():
    q = I.int_bh([0.01, 0.04, 0.03, 0.2])
    assert np.allclose(q, [0.04, 0.0533333, 0.0533333, 0.2], atol=1e-4)

def test_ora_finds_planted_set_and_uses_the_given_background():
    genes = [f"G{i}" for i in range(0, 15)] + ["G150", "G151"]
    o = I.int_ora(genes, SETS, BG)
    assert o.iloc[0]["term"] in ("setA", "setA_copy") and o.iloc[0]["fdr"] < 1e-6 and "setB" not in set(o["term"])
    narrow = I.int_ora(genes, SETS, set(BG) - {f"G{i}" for i in range(60, 200)})        # a smaller background changes the odds
    assert narrow.iloc[0]["p"] != o.iloc[0]["p"]

def test_null_genes_are_not_enriched():
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in rng.choice(np.arange(40, 200), 30, replace=False)]
    o = I.int_ora(genes, SETS, BG)
    assert (o["fdr"] >= 0.05).all() if len(o) else True

def test_trim_removes_redundant_terms():
    genes = [f"G{i}" for i in range(0, 15)]
    o = I.int_ora(genes, SETS, BG)
    sig = o[o["fdr"] < 0.05]
    assert len(sig) >= 2 and len(I.int_trim(sig)) == 1

def test_gsea_positive_when_set_is_at_the_top_and_null_otherwise():
    rng = np.random.default_rng(1)
    s = pd.Series(rng.normal(0, 1, 200), index=[f"G{i}" for i in range(200)])
    s[[f"G{i}" for i in range(0, 20)]] += 3.0
    g = I.int_gsea(s, {k: v for k, v in SETS.items() if k in ("setA", "setC")}, n_perm=200, min_overlap=10)
    a, c = g.set_index("term").loc["setA"], g.set_index("term").loc["setC"]
    assert a["nes"] > 1.5 and a["p"] < 0.05 and c["p"] > 0.05

def test_group_genes_find_planted_markers():
    a = make()
    gl, df = I.int_group_genes(a, a.obs["cluster"].astype(str).to_numpy(), groups=["0", "1"], top_n=10)
    assert set(gl["0"][:5]) == {"G0", "G1", "G2", "G3", "G4"} and set(gl["1"][:5]) == {"G10", "G11", "G12", "G13", "G14"}

def test_set_scores_are_higher_where_the_set_is_expressed():
    a = make()
    sc = I.int_score_sets(a, {"s": {"G0", "G1", "G2", "G3", "G4"}})
    m = (a.obs["cluster"].astype(str) == "0").to_numpy()
    assert sc["s"][m].mean() - sc["s"][~m].mean() > 0.5

def test_nmf_recovers_stable_programs_and_their_cell_groups():
    a = make()
    nm = I.int_nmf(a, k_values=(3, 4), n_restarts=2, max_cells=600)
    assert nm["k"] in (3, 4) and nm["table"]["stability"].min() > 0.5
    top = set(g for p in nm["loadings"].index for g in nm["loadings"].loc[p].sort_values(ascending=False).index[:5])
    assert {"G0", "G10", "G20"} & top

def test_expectations_met_not_met_and_not_testable():
    a = make()
    a.obs["cell_type"] = a.obs["cluster"].astype(str).map({"0": "T1", "1": "T2", "2": "T3"}).fillna("other")
    t = pd.DataFrame([
        dict(expectation="T1 markers up in T1", kind="gene_up", column="cell_type", group="T1", genes="G0;G1;G2"),
        dict(expectation="T2 markers up in T1 (wrong)", kind="gene_up", column="cell_type", group="T1", genes="G10;G11;G12"),
        dict(expectation="negative control holds", kind="gene_up", column="cell_type", group="T1", genes="G10;G11", direction="down"),
        dict(expectation="share is plausible", kind="fraction_range", column="cell_type", group="T1", low=0.1, high=0.3),
        dict(expectation="missing genes", kind="gene_up", column="cell_type", group="T1", genes="NOPE1;NOPE2"),
        dict(expectation="missing column", kind="gene_up", column="nothing", group="x", genes="G0"),
    ])
    r = {x["expectation"]: x["status"] for x in I.int_expectations(a, t)}
    assert r["T1 markers up in T1"] == "met" and r["T2 markers up in T1 (wrong)"] == "not met" and r["negative control holds"] == "met"
    assert r["share is plausible"] == "met" and r["missing genes"] == "not testable" and r["missing column"] == "not testable"

def _loso_data(effect=2.0, n_units=4, cells=40, seed=0):
    rng = np.random.default_rng(seed)
    rows, lab, unit = [], [], []
    for cond in ("ctrl", "stim"):
        for u in range(n_units):
            mu = np.full((cells, 30), 1.0)
            if cond == "stim":
                mu[:, :10] += effect
            rows.append(rng.normal(mu, 1.0))
            lab += [cond] * cells
            unit += [f"{cond}_{u}"] * cells
    X = np.vstack(rows).astype(np.float32)
    a = ad.AnnData(sp.csr_matrix(X))
    a.var_names = [f"G{i}" for i in range(30)]
    a.obs["condition"], a.obs["unit"] = lab, unit
    return a

def test_signature_generalises_to_unseen_samples_only_when_real():
    r = I.int_signature_loso(_loso_data(2.0), [f"G{i}" for i in range(10)], "condition", "unit", "stim", n_perm=100)
    assert r["auc_units"] == 1.0 and r["p"] < 0.05 and r["n_units"] == 8
    n = I.int_signature_loso(_loso_data(0.0, seed=3), [f"G{i}" for i in range(10)], "condition", "unit", "stim", n_perm=100)
    assert n["auc_units"] < 0.95 and n["p"] > 0.05

def test_signature_refuses_mixed_or_tiny_samples():
    a = _loso_data()
    a.obs["unit"] = "one"
    assert "error" in I.int_signature_loso(a, [f"G{i}" for i in range(10)], "condition", "unit", "stim")
    b = _loso_data(n_units=2)
    r = I.int_signature_loso(b, [f"G{i}" for i in range(10)], "condition", "unit", "stim", min_cells=100)
    assert "error" in r

def annotated(seed=0, clean=False):
    from test_scrna_annotate import A as ANN
    a = make(seed=seed)
    if clean:                                   # three clearly separated types, no mixed or noise clusters
        a = a[a.obs["cluster"].isin(["0", "1", "2"])].copy()
        a.obs["cluster"] = a.obs["cluster"].astype(str).astype("category")
    a.layers["counts"] = sp.csr_matrix(np.expm1(a.X.toarray()).round().astype(np.float32))
    a.var["highly_variable"] = True
    a, _ = ANN.ann_run(a, PANEL)
    a.uns["annotation_panel"] = PANEL
    a.uns["dimred"] = {"n_pcs": 5}
    a.uns["clustering"] = {"resolution": 0.5}
    return a

def test_perturbation_calls_clear_structure_robust():
    a = annotated(clean=True)
    specs = [{"name": "seed 1", "seed": 1}, {"name": "80% of cells", "frac": 0.8, "seed": 11}, {"name": "3000 variable genes", "n_hvg": 150}]
    p = I.int_perturb(a, PANEL, specs=specs, max_cells=600)
    assert p["cluster_label"] == "robust" and all(p["summary"][t]["label"] == "robust" for t in ("T1", "T2", "T3"))
    assert (p["runs"]["error"] == "").all()

def test_perturbation_flags_labels_that_were_not_reproducible():
    a = annotated()
    rng = np.random.default_rng(0)
    a.obs["cell_type"] = pd.Categorical(rng.choice(["T1", "T2", "T3"], a.n_obs))      # names unrelated to the data
    a.obs["cluster"] = pd.Categorical(rng.choice(["0", "1", "2"], a.n_obs))
    p = I.int_perturb(a, PANEL, specs=[{"name": "seed 1", "seed": 1}, {"name": "80%", "frac": 0.8, "seed": 2}], max_cells=600)
    assert all(v["label"] in ("unstable", "sensitive") for v in p["summary"].values()) and p["cluster_label"] != "robust"

def test_collect_caveats_reads_reports(tmp_path):
    (tmp_path / "qc.json").write_text(json.dumps({"checks": [{"check": "x", "status": "WARN", "detail": "y"}, {"check": "z", "status": "PASS", "detail": ""}], "not_run": ["batch mixing"]}))
    c = I.int_collect_caveats(str(tmp_path))
    assert {x["kind"] for x in c} == {"WARN", "not run"} and any("batch mixing" in x["text"] for x in c)

def test_collect_caveats_skips_own_stage(tmp_path):
    (tmp_path / "interpret.json").write_text(json.dumps({"checks": [{"check": "old", "status": "WARN", "detail": "from a previous run"}], "not_run": ["x"]}))
    assert I.int_collect_caveats(str(tmp_path)) == []

def test_read_gmt_ignores_header(tmp_path):
    f = tmp_path / "s.gmt"
    f.write_text("# source\nsetA\tdesc\tg1\tG2\tg3\n")
    assert I.int_read_gmt(str(f)) == {"setA": {"G1", "G2", "G3"}}

def test_full_run_writes_outputs_report_and_findings(tmp_path):
    a = annotated()
    exp = pd.DataFrame([dict(expectation="T1 genes up in T1", kind="gene_up", column="cell_type", group="T1", genes="G0;G1;G2")])
    a, recs, det = I.int_run(a, gene_sets=SETS, expectations=exp, max_cells=600, nmf_k=(3, 4), perturb_specs=[{"name": "seed 1", "seed": 1}, {"name": "80%", "frac": 0.8, "seed": 2}],
                             caveats=[{"stage": "qc", "kind": "WARN", "text": "mito: none"}], fig_dir=str(tmp_path))
    for c in ("gene sets", "pathway analysis", "gene programs", "expectations", "robustness: clusters", "robustness: cell-type calls", "caveats carried"):
        assert any(r["check"] == c for r in recs), c
    assert det["findings"] and "X_programs" in a.obsm and "interpretation" in a.uns
    for f in ("int_robustness.png", "int_programs.png", "figure_legends.md", "figure_captions.json", "findings.csv"):
        assert os.path.getsize(tmp_path / f) > 0
    caps = json.load(open(tmp_path / "figure_captions.json"))
    assert set(caps) <= {"int_pathways.png", "int_programs.png", "int_robustness.png", "int_validation.png"} and "int_robustness.png" in caps
    overall = I.write_int_report(recs, det["findings"], det["caveats"], prefix=str(tmp_path / "rep"))
    text = (tmp_path / "rep.md").read_text(encoding="utf-8")
    assert overall in ("PASS", "WARN") and "Findings, their evidence" in text and "Caveats carried" in text and "Not run" in text and "pre-ranked enrichment" in text.split("Not run")[1]

def test_missing_cell_type_column_fails():
    a = make()
    _, recs, _ = I.int_run(a, perturb=False)
    assert recs[0]["status"] == "FAIL"
