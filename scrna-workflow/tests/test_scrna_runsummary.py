"""Tests for the run summary written by the checkpoint engine."""
import json
import pytest
from ckpt import Workflow, StageFailed
import scrna_cluster_kernel as K

def mk(tmp_path, fail=False):
    def s1(prev, p, g):
        return {"files": [{"name": "m", "path": "x", "bytes": 5, "sha256": "ab" * 32}]}, [K.cl_record("cluster count", "PASS", "9 clusters")]
    def s2(prev, p, g):
        return {"x": 1}, [K.cl_record("no tiny clusters", "FAIL" if fail else "WARN", "below 10 cells: 8 (8 cells)")]
    wf = Workflow(tmp_path, {"global": {"species": "human"}, "s2": {"a": 1}}, ["fetch", "s2", "s3"], kinds={"fetch": "json", "s2": "json", "s3": "json"})
    return wf, s1, s2

def test_summary_is_written_after_each_stage_and_explains_checks(tmp_path):
    wf, s1, s2 = mk(tmp_path)
    wf.run("fetch", s1, K.write_cl_report)
    first = (tmp_path / "RUN_SUMMARY.md").read_text(encoding="utf-8")
    assert "1 of 3" in first and "not run" in first            # rebuilt after the first stage, later stages shown as not run
    wf.run("s2", s2, K.write_cl_report)
    md = (tmp_path / "RUN_SUMMARY.md").read_text(encoding="utf-8")
    assert "no cluster is too small to trust" in md             # the 'what it tests' text from the stage report is included
    assert "below 10 cells: 8" in md and "Checks that did not run" in md and "Going back and resuming" in md
    assert "`" + "ab" * 32 + "`" in md and '"a": 1' in md       # input hash and stage parameters recorded
    js = json.load(open(tmp_path / "run_summary.json"))
    assert js["worst"] == "WARN" and js["totals"]["WARN"] == 1 and [r["stage"] for r in js["stages"]] == ["fetch", "s2"]

def test_failed_and_accepted_stage_are_recorded(tmp_path):
    wf, s1, s2 = mk(tmp_path, fail=True)
    wf.run("fetch", s1, K.write_cl_report)
    with pytest.raises(StageFailed):
        wf.run("s2", s2, K.write_cl_report)
    md = (tmp_path / "RUN_SUMMARY.md").read_text(encoding="utf-8")        # written even though the stage failed
    assert "FAIL" in md and "failed" in md
    wf.accept("s2", "not applicable to this dataset")
    assert "accepted: not applicable to this dataset" in wf.summary()

def test_summary_error_does_not_break_a_stage(tmp_path, monkeypatch):
    wf, s1, s2 = mk(tmp_path)
    monkeypatch.setattr(Workflow, "summary", lambda self, write=True: 1 / 0)
    wf.run("fetch", s1, K.write_cl_report)
    assert wf.state("fetch")[0] == "current"

def test_figures_and_legends_are_listed(tmp_path):
    wf, s1, s2 = mk(tmp_path)
    (tmp_path / "figures" / "cluster").mkdir(parents=True)
    (tmp_path / "figures" / "cluster" / "umap.png").write_bytes(b"x")
    (tmp_path / "figures" / "cluster" / "figure_legends.md").write_text("legend")
    wf.run("fetch", s1, K.write_cl_report)
    md = (tmp_path / "RUN_SUMMARY.md").read_text(encoding="utf-8")
    assert "![umap](figures/cluster/umap.png)" in md and "figures/cluster/figure_legends.md" in md

def test_quality_check_rates_each_stage_in_plain_words(tmp_path):
    wf, s1, s2 = mk(tmp_path)
    wf.run("fetch", s1, K.write_cl_report)
    md = (tmp_path / "RUN_SUMMARY.md").read_text(encoding="utf-8")
    assert "## Quality check" in md and "Incomplete: some stages have not been run" in md and "| Good |" in md and "| Not done |" in md
    wf.run("s2", s2, K.write_cl_report)
    md = (tmp_path / "RUN_SUMMARY.md").read_text(encoding="utf-8")
    assert "Usable, with 1 caveat(s)" in md and "wf.invalidate(\"s2\")" in md and md.index("## Quality check") < md.index("## Stages at a glance")

def test_quality_check_marks_failed_and_accepted_stages(tmp_path):
    wf, s1, s2 = mk(tmp_path, fail=True)
    wf.run("fetch", s1, K.write_cl_report)
    with pytest.raises(StageFailed):
        wf.run("s2", s2, K.write_cl_report)
    md = (tmp_path / "RUN_SUMMARY.md").read_text(encoding="utf-8")
    assert "Not reliable yet" in md
    wf.accept("s2", "not applicable to this dataset")
    md = wf.summary()
    assert "Accepted despite a failed check" in md and "not applicable to this dataset" in md and "Not reliable yet" not in md

def test_stage_without_scored_checks_is_not_applicable(tmp_path):
    wf = Workflow(tmp_path, {"global": {}}, ["fetch"], kinds={"fetch": "json"})
    wf.run("fetch", lambda prev, p, g: ({"x": 1}, [K.cl_record("nothing to compare", "INFO", "not requested")]), K.write_cl_report)
    md = (tmp_path / "RUN_SUMMARY.md").read_text(encoding="utf-8")
    assert "| Not applicable |" in md and "1 not applicable" in md and "Overall: Good" in md
