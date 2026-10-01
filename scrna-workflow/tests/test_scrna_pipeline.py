"""Pipeline tests on a synthetic dataset in the GSE120575 file format (genes x cells TPM with a label row,
trailing tabs, and a GEO sample sheet with comment lines). Good files pass; a shifted header must stop stage 1."""
import gzip, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "scrna_workflow"))
import numpy as np, pytest
import ckpt, stages

MITO = ["MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "MT-CYB", "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8"]

def make_files(raw, shift_header=False, seed=0, suffix=False):
    rng = np.random.default_rng(seed)
    genes = MITO + [f"G{i}" for i in range(87)]
    samples = {"Pre_P1": 5, "Post_P1": 4, "Pre_P2": 6, "Post_P2": 5}
    cells, labels = [], []
    for s, n in samples.items():
        for j in range(n):
            cells.append(f"{s}_c{j}"); labels.append(s + "_T_enriched" if (suffix and s == "Post_P2" and j < 3) else s)
    base = rng.gamma(0.7, 1.0, len(genes))
    cols = []
    for _ in cells:
        p = base * rng.gamma(8, 0.125, len(genes)); p /= p.sum()
        cols.append(np.log2(p * 1e6 + 1))
    M = np.array(cols).T                                              # genes x cells
    header_names = cells[1:] + cells[:1] if shift_header else cells
    os.makedirs(raw, exist_ok=True)
    with gzip.open(os.path.join(raw, "m.txt.gz"), "wt", newline="") as f:
        f.write("\t" + "\t".join(header_names) + "\n")
        f.write("\t" + "\t".join(labels) + "\n")
        for g, row in zip(genes, M):
            f.write(g + "\t" + "\t".join(f"{v:.4f}" for v in row) + "\t\n")
    labels_meta = [l.replace("_T_enriched", "") for l in labels]
    resp = {"Pre_P1": "Responder", "Post_P1": "Responder", "Pre_P2": "Non-responder", "Post_P2": "Non-responder"}
    with gzip.open(os.path.join(raw, "meta.txt.gz"), "wt", newline="") as f:
        for i in range(19):
            f.write(f"# comment line {i}\t\t\t\t\t\t\t\n")
        f.write("Sample name\ttitle\tsource name\torganism\tcharacteristics: patinet ID (Pre=baseline; Post= on treatment)\tcharacteristics: response\tcharacteristics: therapy\tmolecule\n")
        for k, (c, l) in enumerate(zip(cells, labels_meta)):
            f.write(f"Sample {k+1}\t{c}\tmelanoma\tHomo sapiens\t{l}\t{resp[l]}\tanti-PD1\t\n")
    return len(cells), len(genes)

def build(tmp_path, shift_header=False, suffix=False):
    raw = str(tmp_path / "raw")
    n_cells, n_genes = make_files(raw, shift_header, suffix=suffix)
    params = {
        "global": {"raw_dir": raw, "run_dir": str(tmp_path / "run"), "species": "human"},
        "fetch": {"files": [{"name": "matrix", "filename": "m.txt.gz", "url": "file:///unused"},
                            {"name": "meta", "filename": "meta.txt.gz", "url": "file:///unused"}]},
        "assemble": {"adapter": "sade_feldman", "declared_scale": "log2tpm", "chunk_rows": 30,
                     "expected": {"shape": [n_cells, n_genes], "n_unique": {"patient_timepoint": 4}}},
        "qc": {"min_genes": 20, "max_genes": None, "max_pct_mt": 40, "min_cells": 3, "sample_col": "patient_timepoint"},
    }
    wf = ckpt.Workflow(tmp_path / "run", params, ["fetch", "assemble", "qc"], kinds={"fetch": "json"})
    return wf

def run3(wf):
    wf.run("fetch", stages.fetch, stages.write_fetch)
    wf.run("assemble", stages.assemble, stages.write_assemble)
    return wf.run("qc", stages.qc, stages.write_qc)

def test_good_files_pass_all_three_stages(tmp_path):
    wf = build(tmp_path)
    ad = run3(wf)
    st = wf.status()
    assert st["state"].tolist() == ["current"] * 3 and st["report"].tolist() == ["PASS", "PASS", "PASS"]
    assert ad.n_obs == 20 and (tmp_path / "run" / "figures" / "qc_before_filtering.png").exists()
    md = (tmp_path / "run" / "reports" / "assemble.md").read_text()
    assert "column order vs metadata" in md and "1 columns" not in md and "**Not run:**" not in md

def test_shifted_header_stops_at_stage_1(tmp_path):
    wf = build(tmp_path, shift_header=True)
    wf.run("fetch", stages.fetch, stages.write_fetch)
    with pytest.raises(ckpt.StageFailed):
        wf.run("assemble", stages.assemble, stages.write_assemble)
    md = (tmp_path / "run" / "reports" / "assemble.md").read_text()
    row = [l for l in md.splitlines() if "column order vs metadata" in l and l.startswith("| FAIL")]
    assert row, md
    assert "cell ID sets match metadata" in md and "| PASS | cell ID sets match metadata" in md   # set check alone would have passed
    with pytest.raises(ckpt.StageFailed, match="upstream stage 'assemble' is failed"):
        wf.run("qc", stages.qc, stages.write_qc)

def test_changing_qc_parameter_recomputes_only_qc(tmp_path):
    wf = build(tmp_path); run3(wf)
    wf.params["qc"]["min_genes"] = 25
    assert wf.status()["state"].tolist() == ["current", "current", "stale"]
    m_before = wf.m["assemble"]["time_utc"], wf.m["assemble"]["output_hash"]
    wf.run("assemble", stages.assemble, stages.write_assemble)         # loads checkpoint
    wf.run("qc", stages.qc, stages.write_qc)
    assert (wf.m["assemble"]["time_utc"], wf.m["assemble"]["output_hash"]) == m_before
    assert wf.status()["state"].tolist() == ["current"] * 3

def test_wrong_declared_scale_fails_stage_1(tmp_path):
    wf = build(tmp_path)
    wf.params["assemble"]["declared_scale"] = "counts"
    wf.run("fetch", stages.fetch, stages.write_fetch)
    with pytest.raises(ckpt.StageFailed):
        wf.run("assemble", stages.assemble, stages.write_assemble)
    assert "declared counts, detected log2tpm" in (tmp_path / "run" / "reports" / "assemble.md").read_text()

def test_replacing_raw_file_makes_everything_stale(tmp_path):
    wf = build(tmp_path); run3(wf)
    make_files(wf.params["global"]["raw_dir"], seed=1)                  # different content, same names
    wf.run("fetch", stages.fetch, stages.write_fetch, force=True)       # fetch re-checks the files
    assert wf.status()["state"].tolist() == ["current", "stale", "stale"]

def test_sort_fraction_suffix_is_separated_not_failed(tmp_path):
    wf = build(tmp_path, suffix=True)
    ad = run3(wf)
    md = (tmp_path / "run" / "reports" / "assemble.md").read_text()
    assert "| PASS | column order vs metadata" in md and "3 cells carry a sort-fraction suffix" in md
    assert (ad.obs["sort_fraction"] == "T_enriched").sum() == 3 and (ad.obs["sort_fraction"] == "unsorted").sum() == 17
