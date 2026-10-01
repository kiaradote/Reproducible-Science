"""Unit tests for the checkpoint engine (ckpt.py) and the chunked loader (loaders.py)."""
import gzip, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "scrna_workflow"))
import numpy as np, pytest
import ckpt, loaders

# ---------- loader ----------
def write_gz(path, text):
    with gzip.open(path, "wt", newline="") as f:
        f.write(text)

SF = "\tc1\tc2\tc3\n\tS1\tS1\tS2\ng1\t1\t0\t2\t\ng2\t0\t0\t3\t\ng3\t4\t5\t0\t\n"

def test_loader_sade_feldman_layout(tmp_path):
    p = tmp_path / "m.txt.gz"; write_gz(p, SF)
    X, ids, genes, labels, lay = loaders.read_genes_by_cells_text(str(p), has_label_row=True, chunk_rows=1)
    assert ids == ["c1", "c2", "c3"] and genes == ["g1", "g2", "g3"] and labels == ["S1", "S1", "S2"]
    assert X.shape == (3, 3) and X.format == "csr"
    assert X.toarray().tolist() == [[1, 0, 4], [0, 0, 5], [2, 3, 0]]       # cells x genes
    assert lay["trailing_fields_after_cells"] == 1 and lay["header_has_blank_gene_column"]

def test_loader_chunk_size_does_not_change_result(tmp_path):
    p = tmp_path / "m.txt.gz"; write_gz(p, SF)
    a = loaders.read_genes_by_cells_text(str(p), has_label_row=True, chunk_rows=1)[0]
    b = loaders.read_genes_by_cells_text(str(p), has_label_row=True, chunk_rows=1000)[0]
    assert (a != b).nnz == 0

def test_loader_rejects_extra_values(tmp_path):
    p = tmp_path / "m.txt.gz"; write_gz(p, "\tc1\tc2\ng1\t1\t2\t9\n")
    with pytest.raises(ValueError):
        loaders.read_genes_by_cells_text(str(p))

def test_loader_rejects_short_rows(tmp_path):
    p = tmp_path / "m.txt.gz"; write_gz(p, "\tc1\tc2\tc3\ng1\t1\t2\n")
    with pytest.raises(ValueError):
        loaders.read_genes_by_cells_text(str(p))

# ---------- engine ----------
def make(tmp_path, params, calls, fail_stage=None):
    wf = ckpt.Workflow(tmp_path / "run", params, ["a", "b", "c"], kinds={"a": "json", "b": "json", "c": "json"})
    def fn_for(name):
        def fn(prev, p, g):
            calls.append(name)
            return {"name": name, "x": p.get("x"), "prev": prev}, [{"check": name, "status": "PASS", "detail": ""}]
        return fn
    def writer_for(name):
        def w(records, prefix):
            open(prefix + ".md", "w").write("report")
            return "FAIL" if name == fail_stage else "PASS"
        return w
    return wf, fn_for, writer_for

def run_all(wf, fn_for, writer_for):
    for s in ("a", "b", "c"):
        wf.run(s, fn_for(s), writer_for(s))

def test_rerun_skips_current_stages(tmp_path):
    calls = []; wf, f, w = make(tmp_path, {"a": {"x": 1}, "b": {"x": 1}, "c": {"x": 1}}, calls)
    run_all(wf, f, w); assert calls == ["a", "b", "c"]
    run_all(wf, f, w); assert calls == ["a", "b", "c"]                  # nothing recomputed
    assert set(wf.status()["state"]) == {"current"}

def test_changed_parameter_marks_only_downstream_stale(tmp_path):
    calls = []; params = {"a": {"x": 1}, "b": {"x": 1}, "c": {"x": 1}}
    wf, f, w = make(tmp_path, params, calls); run_all(wf, f, w)
    params["b"]["x"] = 2
    st = wf.status()["state"].to_dict()
    assert st == {"a": "current", "b": "stale", "c": "stale"}
    calls.clear(); run_all(wf, f, w)
    assert calls == ["b", "c"]                                          # restart from b; a is loaded, not recomputed

def test_restart_state_survives_a_new_session(tmp_path):
    calls = []; params = {"a": {"x": 1}, "b": {"x": 1}, "c": {"x": 1}}
    wf, f, w = make(tmp_path, params, calls); run_all(wf, f, w)
    calls.clear()
    wf2, f2, w2 = make(tmp_path, params, calls)                         # fresh object, same folder
    run_all(wf2, f2, w2); assert calls == []

def test_fail_gate_blocks_downstream_until_accepted(tmp_path):
    calls = []; params = {"a": {"x": 1}, "b": {"x": 1}, "c": {"x": 1}}
    wf, f, w = make(tmp_path, params, calls, fail_stage="b")
    wf.run("a", f("a"), w("a"))
    with pytest.raises(ckpt.StageFailed):
        wf.run("b", f("b"), w("b"))
    assert wf.ckpt_path("b").exists()                                   # checkpoint kept for inspection
    with pytest.raises(ckpt.StageFailed, match="upstream stage 'b' is failed"):
        wf.run("c", f("c"), w("c"))
    wf.accept("b", "no mitochondrial genes in this dataset; mito filter dropped")
    wf.run("c", f("c"), w("c"))
    assert wf.status().loc["b", "accepted"]

def test_invalidate_removes_stage_and_downstream(tmp_path):
    calls = []; params = {"a": {"x": 1}, "b": {"x": 1}, "c": {"x": 1}}
    wf, f, w = make(tmp_path, params, calls); run_all(wf, f, w)
    wf.invalidate("b")
    assert wf.status()["state"].to_dict() == {"a": "current", "b": "not_run", "c": "not_run"}
    assert wf.ckpt_path("a").exists() and not wf.ckpt_path("b").exists() and not wf.ckpt_path("c").exists()

def test_modified_checkpoint_file_makes_downstream_stale(tmp_path):
    calls = []; params = {"a": {"x": 1}, "b": {"x": 1}, "c": {"x": 1}}
    wf, f, w = make(tmp_path, params, calls); run_all(wf, f, w)
    p = wf.ckpt_path("a"); wf.m["a"]["output_hash"] = "0" * 64; wf._save_manifest()   # simulate upstream content change
    assert wf.status()["state"].to_dict()["b"] == "stale"

def test_cannot_run_stage_before_upstream(tmp_path):
    calls = []; wf, f, w = make(tmp_path, {"a": {}, "b": {}, "c": {}}, calls)
    with pytest.raises(ckpt.StageFailed, match="Run it first"):
        wf.run("b", f("b"), w("b"))
