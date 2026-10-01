"""Tests for scrna-data-integrity (checkpoint 1: files, alignment, values, expected numbers).
Real data: GSE84133 / GSM2230757 (human islets donor 1, UMI counts, 1937 cells).
Fault injection: each defect below mirrors a failure seen in the Sade-Feldman example or a common scRNA-seq pitfall."""
import gzip, json, shutil
import numpy as np, pandas as pd, pytest
import scrna_integrity_kernel as K

PATH = "testdata/GSM2230757_human1_umifm_counts.csv.gz"

def status(recs, name):
    return [r["status"] for r in recs if r["check"] == name][0]

@pytest.fixture(scope="module")
def data():
    """Full table read as int32 (low memory); X is a 400-cell subsample so the injected-fault tests stay light.
    Labels, IDs and shape expectations come from the full 1937-cell table."""
    df = pd.read_csv(PATH, index_col=0, dtype={c: np.int32 for c in pd.read_csv(PATH, nrows=1, index_col=0).columns[2:]})
    genes = [c for c in df.columns if c not in ("barcode", "assigned_cluster")]
    full_shape = (df.shape[0], len(genes))
    lab, ids = df["assigned_cluster"], list(df.index)
    X = df[genes].to_numpy()[:: max(1, df.shape[0] // 400)].astype(np.float32)
    del df
    return X, ids, genes, lab, full_shape

# ---------- real data: expected behaviour ----------
def test_real_file_ok():
    r = K.check_file(PATH)
    assert status(r, "gzip stream") == "PASS" and status(r, "file sha256") == "INFO"

def test_real_counts_detected(data):
    X, ids, genes, lab = data[:4]
    recs, det = K.check_values(X, declared="counts")
    assert det == "counts" and status(recs, "scale matches declaration") == "PASS"

def test_real_alignment_ok(data):
    X, ids, genes, lab = data[:4]
    r = K.check_alignment(ids, ids, matrix_labels=lab.to_numpy(), meta_labels=lab)
    assert all(x["status"] == "PASS" for x in r)

def test_real_expected_counts(data):
    X, ids, genes, lab = data[:4]
    assert K.check_expected("n beta cells", int((lab == "beta").sum()), 872)["status"] == "PASS"
    assert K.check_expected("shape", data[4], (1937, 20125))["status"] == "PASS"


# ---------- injected faults ----------
def test_shifted_header_fails_order_but_passes_set_check(data):
    X, ids, genes, lab = data[:4]
    shifted = ids[1:] + ids[:1]                      # same set, off by one position
    r = K.check_alignment(shifted, ids, matrix_labels=lab.to_numpy(), meta_labels=lab)
    assert status(r, "cell ID sets match metadata") == "PASS"
    assert status(r, "column order vs metadata") == "FAIL"

def test_missing_label_row_warns(data):
    r = K.check_alignment(data[1], data[1])
    assert status(r, "column order vs metadata") == "WARN"

def test_duplicate_ids_fail(data):
    ids = data[1][:]; ids[5] = ids[4]
    assert status(K.check_alignment(ids, data[1]), "unique cell IDs") == "FAIL"

def test_nan_fails(data):
    X = data[0].astype(float).copy(); X[3, 7] = np.nan
    recs, _ = K.check_values(X)
    assert status(recs, "NaN/inf values") == "FAIL"

def test_normalized_declared_as_counts_fails(data):
    X = data[0]
    norm = np.log1p(X / X.sum(1, keepdims=True) * 1e4)
    recs, det = K.check_values(norm, declared="counts")
    assert det == "log1p_cp10k" and status(recs, "scale matches declaration") == "FAIL"

def test_log2tpm_detected_like_sade_feldman(data):
    X = data[0]
    tpm = X / X.sum(1, keepdims=True) * 1e6
    _, det = K.check_values(np.log2(tpm + 1))
    assert det == "log2tpm"
    _, det = K.check_values(tpm)
    assert det == "linear_tpm"

def test_scaled_data_not_misread(data):
    X = data[0].astype(float)
    z = (X - X.mean(0)) / (X.std(0) + 1e-9)
    _, det = K.check_values(z)
    assert det == "scaled_or_unknown"

def test_truncated_gzip_fails(tmp_path):
    bad = tmp_path / "trunc.csv.gz"
    raw = open(PATH, "rb").read()
    bad.write_bytes(raw[: len(raw) // 2])
    r = K.check_file(str(bad))
    assert status(r, "gzip stream") == "FAIL"

def test_wrong_checksum_fails():
    r = K.check_file(PATH, expected_sha256="0" * 64)
    assert status(r, "file sha256") == "FAIL"


def test_report_written(tmp_path, data):
    recs = K.check_file(PATH) + [K.check_expected("shape", (1, 2), (1, 3))]
    overall = K.write_integrity_report(recs, str(tmp_path / "rep"))
    assert overall == "FAIL" and (tmp_path / "rep.md").exists() and (tmp_path / "rep.json").exists()

def test_report_explains_verdict_and_lists_checks(tmp_path):
    recs = K.check_file(PATH, expected_bytes=1) + [K.check_expected("shape", (1, 2), (1, 2))]
    overall = K.write_integrity_report(recs, str(tmp_path / "r"))
    md = (tmp_path / "r.md").read_text()
    why = md.split("**Why:**")[1].split("\n")[0]
    assert overall == "FAIL" and why.startswith(" 1 of") and "file size" in why
    assert "| what it tests |" in md
    assert "**Not run:**" in md and "- unique cell IDs:" in md and "- gzip stream:" not in md
    js = json.load(open(tmp_path / "r.json"))
    assert "unique cell IDs" in js["not_run"] and js["counts"]["FAIL"] == 1

def test_report_pass_lists_everything_and_no_not_run(tmp_path, data):
    X, ids, genes, lab = data[:4]
    recs = (K.check_file(PATH) + K.check_alignment(ids, ids, lab.to_numpy(), lab)
            + K.check_values(X, "counts")[0] + [K.check_expected("shape", data[4], data[4])])
    overall = K.write_integrity_report(recs, str(tmp_path / "r"))
    md = (tmp_path / "r.md").read_text()
    assert overall == "PASS" and "all " in md.split("**Why:**")[1].split("\n")[0] and "**Not run:**" not in md

def test_sparse_input_same_result(data):
    import scipy.sparse as sp
    X = data[0]
    for name, M in {"counts": X, "log2tpm": np.log2(X / X.sum(1, keepdims=True) * 1e6 + 1)}.items():
        _, d_dense = K.check_values(M)
        _, d_sparse = K.check_values(sp.csr_matrix(M))
        assert d_dense == d_sparse
    assert d_sparse == "log2tpm"
