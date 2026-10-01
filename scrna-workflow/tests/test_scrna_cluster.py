"""Tests for scrna-clustering (checkpoint 5). Synthetic data with known groups for the logic; real counts for the pipeline."""
import os
import numpy as np, pytest
import scrna_dimred_kernel as K
import scrna_cluster_kernel as C
from test_scrna_dimred import prepared, status, PATH

_cache = {}

def embedded(groups=4, n_cells=800, batch_shift=0.0, seed=0, **kw):
    key = (groups, n_cells, batch_shift, seed, tuple(sorted(kw.items())))
    if key not in _cache:
        a = prepared(n_cells=n_cells, groups=groups, batch_shift=batch_shift, seed=seed)
        a, _ = K.dr_run(a, **kw)
        _cache[key] = a
    return _cache[key].copy()

def ari(a, b):
    from sklearn.metrics import adjusted_rand_score
    return adjusted_rand_score(a, b)

def test_clear_groups_recovered_and_stable(tmp_path):
    a = embedded(groups=4)
    a, recs = C.cl_run(a, fig_dir=str(tmp_path))
    assert a.obs["cluster"].nunique() == 4 and ari(a.obs["group"], a.obs["cluster"]) > 0.95
    assert status(recs, "chosen resolution is stable") == "PASS" and status(recs, "no tiny clusters") == "PASS"
    assert status(recs, "no dominant cluster") == "PASS" and status(recs, "cluster markers") == "PASS"
    assert status(recs, "neighboring clusters distinguishable") == "PASS"
    assert all(r["status"] != "FAIL" for r in recs)
    for k in ("cluster", "leiden_r%s" % a.uns["clustering"]["resolution"]):
        assert k in a.obs
    assert {"clustering", "cluster_resolution_scan", "cluster_markers"} <= set(a.uns)

def test_outputs_and_legends_exist(tmp_path):
    a = embedded(groups=3)
    a, recs = C.cl_run(a, fig_dir=str(tmp_path))
    for f in ("cluster_resolution.png", "cluster_umap.png", "cluster_sizes.png", "cluster_markers.png", "cluster_markers.csv", "resolution_scan.csv", "figure_legends.md"):
        assert os.path.getsize(tmp_path / f) > 0
    leg = (tmp_path / "figure_legends.md").read_text()
    assert str(a.uns["clustering"]["n_clusters"]) in leg and str(a.uns["clustering"]["resolution"]) in leg

def test_same_seed_same_clusters(tmp_path):
    a1, _ = C.cl_run(embedded(groups=3), seed=3)
    a2, _ = C.cl_run(embedded(groups=3), seed=3)
    assert (a1.obs["cluster"].values == a2.obs["cluster"].values).all()

def test_tiny_cluster_is_flagged():
    a = embedded(groups=4)
    keep = np.ones(a.n_obs, bool)
    idx = np.where(a.obs["group"].values == "3")[0]
    keep[idx[6:]] = False                      # leave only 6 cells of one real group
    b = a[keep].copy()
    b, recs = C.cl_run(b)
    assert status(recs, "no tiny clusters") == "WARN"

def test_dominant_cluster_is_flagged():
    a = embedded(groups=3)
    keep = np.ones(a.n_obs, bool)
    for g in ("1", "2"):
        idx = np.where(a.obs["group"].values == g)[0]
        keep[idx[40:]] = False                 # group 0 now holds ~75% of cells
    b = a[keep].copy()
    b, recs = C.cl_run(b)
    assert status(recs, "no dominant cluster") == "WARN"

def test_no_structure_is_not_called_stable():
    """Pure noise: forcing a fine resolution must not pass the stability or distinguishability checks."""
    a = embedded(groups=1)
    a, recs = C.cl_run(a, resolution=2.0)
    bad = [status(recs, "chosen resolution is stable"), status(recs, "neighboring clusters distinguishable"), status(recs, "per-cluster stability")]
    assert "WARN" in bad or "FAIL" in bad

def test_fixed_resolution_is_used():
    a, recs = C.cl_run(embedded(groups=3), resolution=0.5)
    assert a.uns["clustering"]["resolution"] == 0.5

def test_pick_resolution_rules():
    t = [dict(resolution=0.1, n_clusters=3, stability=0.95, min_size=50),
         dict(resolution=0.5, n_clusters=6, stability=0.85, min_size=40),
         dict(resolution=1.0, n_clusters=12, stability=0.60, min_size=5)]
    assert C.cl_pick_resolution(t) == (0.5, True)                                   # finest that qualifies
    t2 = [dict(r, min_size=3) for r in t]
    assert C.cl_pick_resolution(t2) == (0.1, False)                                 # none qualifies: most stable
    t3 = [dict(resolution=0.1, n_clusters=1, stability=1.0, min_size=800), dict(resolution=0.5, n_clusters=4, stability=0.7, min_size=30)]
    assert C.cl_pick_resolution(t3) == (0.5, False)                                 # one cluster is never a result

def test_best_jaccard_and_rank_labels():
    ref = np.array(list("aaabbb")); same = np.array(list("xxxyyy")); half = np.array(list("xxyyyy"))
    assert np.allclose(C.cl_best_jaccard(ref, same, ["a", "b"]), [1.0, 1.0])
    j = C.cl_best_jaccard(ref, half, ["a", "b"])
    assert j[0] == pytest.approx(2 / 3) and j[1] == pytest.approx(3 / 4)
    lab = C.cl_rank_labels(np.array(["7", "7", "7", "2", "2", "9"]))
    assert list(lab) == ["0", "0", "0", "1", "1", "2"]                              # largest cluster is 0

def test_batch_mixing_flags_single_batch_cluster():
    labels = np.array(["0"] * 40 + ["1"] * 40 + ["2"] * 40)
    batch = np.array(["A", "B"] * 20 + ["A"] * 40 + ["A", "B"] * 20)
    out = C.cl_batch_mixing(labels, batch)
    assert [o[0] for o in out] == ["1"]

def test_batch_driven_clusters_are_flagged():
    a = embedded(groups=2, batch_shift=6.0, seed=1)
    a, recs = C.cl_run(a, batch_key="batch", resolution=1.0)
    assert any(r["check"].startswith("batch") for r in recs)

def test_harmony_representation_is_used_when_present():
    a = embedded(groups=3)
    a.obsm["X_pca_harmony"] = a.obsm["X_pca"].copy()
    a.uns["dimred"] = dict(a.uns["dimred"], rep="X_pca_harmony")
    a, recs = C.cl_run(a)
    d = [r["detail"] for r in recs if r["check"] == "representation used"][0]
    assert "X_pca_harmony" in d and "batch-corrected" in d

def test_report_lists_checks_and_not_run(tmp_path):
    a, recs = C.cl_run(embedded(groups=3))
    overall = C.write_cl_report(recs, prefix=str(tmp_path / "cl"))
    text = open(str(tmp_path / "cl") + ".md", encoding="utf8").read()
    assert overall in ("PASS", "WARN") and "Not run" in text and "what it tests" in text.lower() and "Why" in text

@pytest.mark.skipif(not os.path.exists(PATH), reason="real data not present")
def test_real_pancreas_clusters_match_known_cell_types(tmp_path):
    import scanpy as sc
    a = sc.read_h5ad("runs/baron/checkpoints/04_dimred.h5ad") if os.path.exists("runs/baron/checkpoints/04_dimred.h5ad") else None
    if a is None:
        pytest.skip("checkpoint 4 output not present")
    a, recs = C.cl_run(a, fig_dir=str(tmp_path))
    assert ari(a.obs["assigned_cluster"], a.obs["cluster"]) > 0.75
    assert all(r["status"] != "FAIL" for r in recs)
