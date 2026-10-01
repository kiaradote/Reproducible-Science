"""Tests for scrna-qc-filtering (checkpoint 2: mitochondrial gene set, per-filter and per-sample cell loss).
Real data: GSE84133 / GSM2230757 (human islets donor 1, UMI counts). Fault injection: contaminated or missing mito gene sets, a sample losing most of its cells."""
import numpy as np, pandas as pd, pytest
import scrna_qc_kernel as K

PATH = "testdata/GSM2230757_human1_umifm_counts.csv.gz"

def status(recs, name):
    return [r["status"] for r in recs if r["check"] == name][0]

@pytest.fixture(scope="module")
def data():
    cols = pd.read_csv(PATH, nrows=1, index_col=0).columns
    return None, None, [c for c in cols if c not in ("barcode", "assigned_cluster")]

def test_real_data_has_no_mito_genes(data):
    # Real finding: this matrix contains no MT- genes (and metallothioneins MT1A etc. must not be mistaken for them)
    recs, mito = K.check_mito(data[2], "human")
    assert mito == [] and recs[0]["status"] == "FAIL"

def test_mito_prefix_rules():
    names = ["MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "MT-CYB", "MT-CO1", "MT-CO2",
             "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-RNR1", "MTRNR2L1", "MTRNR2L8", "MT1A", "MTHFD1", "ACTB"]
    recs, mito = K.check_mito(names, "human")
    assert "MTRNR2L1" not in mito and "MT1A" not in mito and len(mito) == 14
    assert status(recs, "nuclear MTRNR2L pseudogenes") == "WARN"
    recs, mito = K.check_mito(["ENSG00000198888", "ENSG00000198763"], "human")
    assert recs[0]["status"] == "FAIL"
    recs, mito = K.check_mito(["mt-Nd1", "Actb"], "mouse")
    assert mito == ["mt-Nd1"] and recs[0]["status"] == "WARN"      # missing 12 protein-coding genes

def test_filter_report_flags_sample_loss():
    obs = pd.DataFrame({"sample": ["A"] * 350 + ["B"] * 50}, index=[f"c{i}" for i in range(400)])
    mask = pd.Series([False] * 350 + [True] * 30 + [False] * 20, index=obs.index)   # 7.5% overall, 60% of sample B
    r = K.qc_filter_report(obs, {"pct_mt>20": mask}, sample_col="sample")
    assert status(r, "per-sample loss") == "WARN"
    assert status(r, "total cells removed") == "PASS"
    mask2 = pd.Series([True] * 100 + [False] * 300, index=obs.index)
    assert status(K.qc_filter_report(obs, {"x": mask2}), "total cells removed") == "WARN"

def test_report_written(tmp_path, data):
    recs, _ = K.check_mito(data[2], "human")
    overall = K.write_qc_report(recs, str(tmp_path / "rep"))
    assert overall == "FAIL" and (tmp_path / "rep.md").exists() and (tmp_path / "rep.json").exists()

def test_report_explains_verdict_and_lists_checks(tmp_path, data):
    recs, _ = K.check_mito(data[2], "human")
    obs = pd.DataFrame({"s": ["A"] * 10}, index=[f"c{i}" for i in range(10)])
    recs += K.qc_filter_report(obs, {"x": pd.Series([False] * 10, index=obs.index)})
    overall = K.write_qc_report(recs, str(tmp_path / "r"))
    md = (tmp_path / "r.md").read_text()
    assert overall == "FAIL" and "**Why:** 1 of 3 checks failed: mito genes found" in md
    assert "**Not run:**" in md and "- per-sample loss:" in md

def test_mito_genes_listed_but_without_counts_is_flagged():
    zero = K.check_mito_signal(np.zeros(500), 20)
    assert zero["status"] == "WARN" and "removed upstream" in zero["detail"]
    ok = K.check_mito_signal(np.random.default_rng(0).uniform(0.5, 10, 500), 20)
    assert ok["status"] == "PASS"
    assert K.check_mito_signal([np.nan] * 5)["status"] == "WARN"
    assert "mito signal present" in K.qc_catalog()
