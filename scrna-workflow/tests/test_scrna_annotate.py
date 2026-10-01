"""Tests for scrna-cell-type-annotation (checkpoint 6). Synthetic clusters with known marker programs for the logic; real pancreas data for the pipeline."""
import os, json
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, pytest
import scrna_annotate_kernel as A

RAW = "handoff/cellguide_pancreas_raw.json"
CK = "runs/baron/checkpoints/05_cluster.h5ad"

def status(recs, name):
    return [r["status"] for r in recs if r["check"] == name][0]

def make(n_per=60, n_genes=200, seed=0, with_umap=True):
    """Clusters: 0 type T1 markers G0-4; 1 T2 markers G10-14; 2 T3 markers G20-24; 3 half T1 / half T2 (mixed);
    4 no marker program (noise); 5 T1 markers plus G30-31, which are negative markers of T1."""
    rng = np.random.default_rng(seed)
    prog = {0: [range(0, 5)], 1: [range(10, 15)], 2: [range(20, 25)], 3: [range(0, 3), range(10, 13)], 4: [], 5: [range(0, 5), range(30, 32)]}
    labs, rows = [], []
    for c, groups in prog.items():
        mu = np.full((n_per, n_genes), 0.05)
        for g in groups:
            mu[:, list(g)] = 3.0
        rows.append(rng.poisson(mu) * rng.uniform(0.8, 1.2, (n_per, n_genes)))
        labs += [str(c)] * n_per
    X = np.log1p(np.vstack(rows)).astype(np.float32)
    a = ad.AnnData(sp.csr_matrix(X))
    a.var_names = [f"G{i}" for i in range(n_genes)]
    a.obs_names = [f"c{i}" for i in range(len(labs))]
    a.obs["cluster"] = pd.Categorical(labs)
    if with_umap:
        cen = rng.normal(0, 6, (6, 2))
        a.obsm["X_umap"] = np.vstack([cen[int(c)] + rng.normal(0, 0.6, 2) for c in labs]).astype(np.float32)
    return a

PANEL = pd.DataFrame({"cell_type": ["T1", "T2", "T3"],
                      "positive": ["G0;G1;G2;G3;G4", "G10;G11;G12;G13;G14", "G20;G21;G22;G23;G24"],
                      "negative": ["G30;G31", "", ""], "source": ["test", "test", "test"]})

def row(a, c):
    t = a.uns["annotation"]
    return t[t["cluster"] == str(c)].iloc[0]

def test_clear_clusters_named_with_right_tier():
    a, recs = A.ann_run(make(), PANEL)
    assert [row(a, c)["label"] for c in (0, 1, 2)] == ["T1", "T2", "T3"]
    assert all(row(a, c)["tier"] == "confident" for c in (0, 1, 2))
    assert {"cell_type", "annotation_tier", "cell_type_candidates"} <= set(a.obs) and "annotation" in a.uns
    assert all(r["status"] != "FAIL" for r in recs)

def test_mixed_cluster_is_ambiguous_and_lists_both_candidates():
    a, _ = A.ann_run(make(), PANEL)
    r = row(a, 3)
    assert r["tier"] == "ambiguous" and r["label"] == "ambiguous" and set(r["candidates"].split(" / ")) == {"T1", "T2"}

def test_cluster_without_markers_is_unassigned():
    a, _ = A.ann_run(make(), PANEL)
    assert row(a, 4)["tier"] == "unassigned" and row(a, 4)["label"] == "unassigned"

def test_negative_markers_on_make_a_cluster_ambiguous():
    a, recs = A.ann_run(make(), PANEL)
    assert row(a, 5)["tier"] == "ambiguous" and row(a, 5)["negatives_on"] == 2
    assert status(recs, "negative marker conflicts") == "WARN"

def test_obs_labels_follow_clusters():
    a, _ = A.ann_run(make(), PANEL)
    assert set(a.obs.loc[a.obs["cluster"] == "1", "cell_type"]) == {"T2"}
    assert set(a.obs.loc[a.obs["cluster"] == "4", "annotation_tier"]) == {"unassigned"}

def test_gene_names_match_case_insensitively():
    p = PANEL.copy(); p["positive"] = p["positive"].str.lower(); p["negative"] = p["negative"].str.lower()
    a, _ = A.ann_run(make(), p)
    assert row(a, 0)["label"] == "T1"

def test_panel_genes_missing_leaves_type_out_or_fails():
    p = pd.concat([PANEL, pd.DataFrame({"cell_type": ["T9"], "positive": ["NOPE1;NOPE2"], "negative": [""], "source": ["test"]})], ignore_index=True)
    a, recs = A.ann_run(make(), p)
    assert status(recs, "marker panel") == "WARN" and "T9" in [r["detail"] for r in recs if r["check"] == "marker panel"][0]
    bad = PANEL.assign(positive=["X1;X2", "X3;X4", "X5;X6"])
    _, recs = A.ann_run(make(), bad)
    assert status(recs, "marker panel") == "FAIL"

def test_single_gene_panel_is_flagged_weak():
    p = PANEL.copy(); p.loc[2, "positive"] = "G20"
    _, recs = A.ann_run(make(), p)
    assert status(recs, "panel strength") == "WARN"

def test_missing_cluster_column_fails():
    a = make(); del a.obs["cluster"]
    _, recs = A.ann_run(a, PANEL)
    assert recs[0]["status"] == "FAIL"

def test_call_rules_are_pure_and_ordered():
    idx = ["0", "1", "2", "3", "4"]
    score = pd.DataFrame({"A": [0.9, 0.6, 0.5, 0.2, 0.9], "B": [0.3, 0.5, 0.48, 0.1, 0.3]}, index=idx)
    pdet = pd.DataFrame({"A": [1.0] * 5, "B": [1.0] * 5}, index=idx)
    nneg = pd.DataFrame({"A": [0, 0, 0, 0, 1], "B": [0] * 5}, index=idx)
    sizes = pd.Series([100] * 5, index=idx)
    t = A.ann_call(score, pdet, nneg, None, sizes).set_index("cluster")
    assert list(t["tier"]) == ["confident", "probable", "ambiguous", "unassigned", "ambiguous"]
    t2 = A.ann_call(score.iloc[:1], pdet.iloc[:1], nneg.iloc[:1], {"0": ("B", 0.9)}, sizes).iloc[0]       # reference disagrees: confident -> probable
    assert t2["tier"] == "probable" and "reference says B" in t2["reason"]
    t3 = A.ann_call(score.iloc[:1], pdet.iloc[:1], nneg.iloc[:1], None, pd.Series([5], index=["0"])).iloc[0]   # tiny cluster capped
    assert t3["tier"] == "probable" and "5 cells" in t3["reason"]

def test_suggest_negatives_only_fills_empty_ones():
    out, auto = A.ann_suggest_negatives(PANEL, per_type=2)
    assert auto == ["T2", "T3"] and out.loc[0, "negative"] == "G30;G31"
    assert out.loc[1, "negative"].split(";")[:2] == ["G0", "G1"]

def test_match_label_and_map():
    types = ["pancreatic beta cell", "T cell"]
    assert A.ann_match_label("beta", types) == "pancreatic beta cell" and A.ann_match_label("t_cell", types) is None
    assert A.ann_match_label("t_cell", types, {"t_cell": "T cell"}) == "T cell"

def _ref(shuffle_names=False, counts=False):
    a = make(seed=1, with_umap=False)
    names = {"0": "T1", "1": "T2", "2": "T3", "3": "mix", "4": "noise", "5": "T1b"}
    if shuffle_names:
        names = {"0": "T2", "1": "T1", "2": "T3", "3": "mix", "4": "noise", "5": "T1b"}
    a.obs["cell_type"] = a.obs["cluster"].astype(str).map(names).values
    if counts:
        a.X = sp.csr_matrix(np.expm1(a.X.toarray()).round())
    return a

def test_reference_agrees_even_on_a_different_scale():
    a, recs = A.ann_run(make(), PANEL, reference=_ref(counts=True))
    r = [x for x in recs if x["check"] == "reference agreement"][0]
    assert r["status"] == "PASS" and row(a, 0)["reference_label"] == "T1"

def test_reference_disagreement_is_reported_and_downgrades():
    a, recs = A.ann_run(make(), PANEL, reference=_ref(shuffle_names=True))
    assert status(recs, "reference agreement") == "WARN"
    assert row(a, 1)["label"] == "T2" and row(a, 1)["tier"] == "probable" and "reference says T1" in row(a, 1)["reason"]   # cluster 1 has a unique reference match

def test_reference_too_few_shared_genes_warns():
    r = _ref(); r.var_names = [f"X{i}" for i in range(r.n_vars)]
    _, recs = A.ann_run(make(), PANEL, reference=r)
    assert status(recs, "reference agreement") == "WARN"

def test_cellguide_cleaning_filters_species_duplicates_and_shared_genes():
    raw = {"x": {"id": "CL:1", "canonical": [{"symbol": "Aaa"}, {"symbol": "Aaa"}],
                 "computational": [{"symbol": "Mm1", "organism": "Mus musculus", "marker_score": 9}, {"symbol": "H1", "organism": "Homo sapiens", "marker_score": 2},
                                   {"symbol": "H2", "organism": "Homo sapiens", "marker_score": 5}, {"symbol": "MT-RNR1", "organism": "Homo sapiens", "marker_score": 99},
                                   {"symbol": "SHARED", "organism": "Homo sapiens", "marker_score": 50}]},
           "y": {"id": "CL:2", "canonical": [], "computational": [{"symbol": "SHARED", "organism": "Homo sapiens", "marker_score": 3}]},
           "z": {"id": "CL:3", "canonical": [], "computational": [{"symbol": "SHARED", "organism": "Homo sapiens", "marker_score": 3}]},
           "bad": {"error": "no match"}}
    p = A.ann_clean_cellguide(raw, max_types_sharing=2).set_index("cell_type")
    assert p.loc["x", "positive"] == "AAA;H2;H1" and p.loc["y", "positive"] == "" and "bad" not in p.index

def test_outputs_legends_and_report(tmp_path):
    a, recs = A.ann_run(make(), PANEL, fig_dir=str(tmp_path))
    for f in ("annotation_scores.png", "annotation_umap.png", "annotation_markers.png", "annotation_table.csv", "marker_panel_used.csv", "figure_legends.md"):
        assert os.path.getsize(tmp_path / f) > 0
    assert "3 of 6 clusters are named" in (tmp_path / "figure_legends.md").read_text() or "clusters are named" in (tmp_path / "figure_legends.md").read_text()
    overall = A.write_ann_report(recs, prefix=str(tmp_path / "ann"))
    text = (tmp_path / "ann.md").read_text()
    assert overall in ("PASS", "WARN") and "Not run" in text and "pretrained annotation tools" in text and "reference agreement" in text.split("Not run")[1]

@pytest.mark.skipif(not (os.path.exists(RAW) and os.path.exists(CK)), reason="pancreas checkpoint or CellGuide fixture not present")
def test_real_pancreas_calls_match_authors_for_every_named_cluster():
    raw = json.load(open(RAW))
    a = ad.read_h5ad(CK)
    a, recs = A.ann_run(a, A.ann_clean_cellguide(raw))
    truth = pd.crosstab(a.obs["cluster"], a.obs["assigned_cluster"]).idxmax(axis=1)
    key = {"pancreatic beta cell": "beta", "pancreatic alpha cell": "alpha", "pancreatic ductal cell": "ductal", "endothelial cell": "endothelial",
           "mast cell": "mast", "pancreatic stellate cell": "quiescent_stellate"}
    t = a.uns["annotation"]
    named = t[t["tier"].isin(["confident", "probable"])]
    assert len(named) >= 6 and all(key[r.label] == truth[r.cluster] for r in named.itertuples())
    assert all(r["status"] != "FAIL" for r in recs)

def test_sparse_droplet_data_still_names_clusters():
    """Thin the counts so even the best cluster detects its markers in under 25% of cells (10x-like depth): the relative presence rule must still work."""
    a = make()
    rng = np.random.default_rng(5)
    X = a.X.toarray()
    counts = np.expm1(X).round()
    thin = rng.binomial(counts.astype(int), 0.08)
    a.X = sp.csr_matrix(np.log1p(thin).astype(np.float32))
    det = (a.X.toarray()[a.obs["cluster"].astype(str).to_numpy() == "0"][:, :5] > 0).mean()
    assert det < 0.25                                           # sparse_droplet: the strict 25% rule would find nothing
    a, recs = A.ann_run(a, PANEL)
    assert row(a, 0)["label"] == "T1" and row(a, 1)["label"] == "T2" and row(a, 2)["label"] == "T3"
    assert row(a, 4)["tier"] in ("unassigned", "ambiguous")      # noise cluster is still not named

def test_undetectable_panel_genes_are_set_aside_not_counted_against_the_type():
    p = PANEL.copy()
    p.loc[0, "positive"] = "G0;G1;G2;G3;G4;G40;G41;G42;G43;G44;G45"           # six extra genes that are never expressed (dropout)
    a0 = make()
    m_ = a0.X.tolil(); m_[:, 40:46] = 0; a0.X = m_.tocsr()
    a, recs = A.ann_run(a0, p)
    assert row(a, 0)["label"] == "T1" and row(a, 0)["tier"] == "confident"
    r = [x for x in recs if x["check"] == "sparse panel genes"][0]
    assert r["status"] == "INFO" and "G40" in r["detail"] and "T1" in r["detail"]

def test_one_usable_gene_is_not_enough_to_name_a_cluster():
    p = PANEL.copy()
    p.loc[2, "positive"] = "G20;G40"                                  # G40 is never expressed, leaving a single usable gene for T3
    a0 = make()
    m_ = a0.X.tolil(); m_[:, 40:41] = 0; a0.X = m_.tocsr()
    a, _ = A.ann_run(a0, p)
    assert row(a, 2)["label"] != "T3"

def test_type_with_too_few_usable_markers_is_reported_and_never_named():
    p = PANEL.copy()
    p.loc[2, "positive"] = "G20;G21;G40;G41;G42;G43"                  # only two of six genes are ever detected: fewer than min(3, 6)
    a0 = make()
    m_ = a0.X.tolil(); m_[:, 40:44] = 0; a0.X = m_.tocsr()
    a, recs = A.ann_run(a0, p)
    assert row(a, 2)["label"] != "T3" and status(recs, "enough usable markers") == "WARN"
