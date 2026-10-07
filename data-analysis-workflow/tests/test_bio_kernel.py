"""Tests for the generic biological-interpretation helpers."""
import numpy as np, pandas as pd, pytest
import bio_kernel as B

BG = {f"G{i}" for i in range(300)}
SETS = {"A": {f"G{i}" for i in range(0, 25)}, "A_copy": {f"G{i}" for i in range(0, 24)}, "B": {f"G{i}" for i in range(50, 80)}, "C": {f"G{i}" for i in range(120, 150)},
        "D": {f"G{i}" for i in range(180, 215)}, "E": {f"G{i}" for i in range(220, 260)}}

def test_bh_known_values():
    assert np.allclose(B.bio_bh([0.01, 0.04, 0.03, 0.2]), [0.04, 0.0533333, 0.0533333, 0.2], atol=1e-4)

def test_enrichment_finds_the_planted_set_and_trims_its_near_copy():
    genes = [f"G{i}" for i in range(0, 18)] + ["G200", "G201"]
    o = B.bio_ora(genes, SETS, BG)
    assert o.iloc[0]["term"] in ("A", "A_copy") and o.iloc[0]["fdr"] < 1e-6 and "B" not in set(o["term"])
    assert len(B.bio_trim(o[o["fdr"] < 0.05])) == 1

def test_enrichment_uses_the_background_given():
    genes = [f"G{i}" for i in range(0, 12)]
    wide = B.bio_ora(genes, SETS, BG).iloc[0]["p"]
    narrow = B.bio_ora(genes, SETS, {f"G{i}" for i in range(0, 60)}).iloc[0]["p"]
    assert narrow > wide                                                  # a smaller background makes the same overlap less surprising

def test_random_genes_are_not_enriched():
    rng = np.random.default_rng(0)
    o = B.bio_ora([f"G{i}" for i in rng.choice(np.arange(260, 300), 25, replace=False)], SETS, BG)
    assert len(o) == 0 or (o["fdr"] >= 0.05).all()

def test_preranked_enrichment_has_power_across_many_sets_and_no_false_alarm():
    rng = np.random.default_rng(1)
    s = pd.Series(rng.normal(0, 1, 300), index=[f"G{i}" for i in range(300)])
    s[[f"G{i}" for i in range(0, 25)]] += 3.0
    many = {**SETS, **{f"X{k}": {f"G{i}" for i in rng.choice(300, 30, replace=False)} for k in range(120)}}
    g = B.bio_gsea(s, many, n_perm=800, min_overlap=10).set_index("term")
    assert g.loc["A", "fdr"] < 0.05 and g.loc["A", "nes"] > 1.5 and (g.drop(index=["A", "A_copy"])["fdr"] < 0.05).sum() <= 3

def test_gmt_files_are_read_and_headers_ignored(tmp_path):
    f = tmp_path / "s.gmt"
    f.write_text("# KEGG downloaded today\nsetA\tdesc\tg1\tG2\tg3\n")
    assert B.bio_read_gmt(str(f)) == {"setA": {"G1", "G2", "G3"}}

def test_generic_terms_are_flagged():
    df = pd.DataFrame({"term": ["Ribosome [hsa03010]", "Influenza A [hsa05164]", "Insulin secretion [hsa04911]"]})
    fl = B.bio_flag_terms(df)["flag"].tolist()
    assert fl[0] == "housekeeping machinery" and fl[1].startswith("disease or infection") and fl[2] == ""

def result_table():
    rng = np.random.default_rng(2)
    r = pd.DataFrame({"feature": [f"G{i}" for i in range(300)], "effect": rng.normal(0, 0.2, 300), "padj": rng.uniform(0.2, 1, 300)})
    r.loc[:19, "effect"] = 1.5; r.loc[:19, "padj"] = 0.001
    r.loc[20:39, "effect"] = -1.2; r.loc[20:39, "padj"] = 0.002
    return r

def test_expectations_met_not_met_and_not_testable():
    r = result_table()
    ex = [{"name": "marker set rises", "features": [f"G{i}" for i in range(10)], "direction": "up"},
          {"name": "repressed set falls", "features": [f"G{i}" for i in range(20, 30)], "direction": "down"},
          {"name": "wrong way", "features": [f"G{i}" for i in range(10)], "direction": "down"},
          {"name": "negative control is quiet", "features": [f"G{i}" for i in range(100, 120)], "direction": "none"},
          {"name": "negative control that is not quiet", "features": [f"G{i}" for i in range(0, 20)], "direction": "none"},
          {"name": "features absent", "features": ["NOPE1", "NOPE2", "G1"], "direction": "up"}]
    s = {e["expectation"]: e["status"] for e in B.bio_expectations(r, ex)}
    assert s["marker set rises"] == "met" and s["repressed set falls"] == "met" and s["wrong way"] == "not met"
    assert s["negative control is quiet"] == "met" and s["negative control that is not quiet"] == "not met" and s["features absent"] == "not testable"

def test_agreement_with_another_study():
    a = result_table()
    same = B.bio_concordance(a, a.copy())
    assert same["same_direction"] == 1.0 and same["overlap"] == 40 and same["overlap_p"] < 1e-10 and same["effect_correlation"] > 0.99
    flipped = a.copy(); flipped["effect"] = -flipped["effect"]
    assert B.bio_concordance(a, flipped)["same_direction"] == 0.0
    rng = np.random.default_rng(5)
    other = a.copy(); other["effect"] = rng.normal(0, 1, len(a)); other["padj"] = rng.uniform(0, 1, len(a))
    r = B.bio_concordance(a, other)
    assert r["overlap_p"] > 0.01 and abs(r["effect_correlation"]) < 0.3
    assert B.bio_concordance(a.iloc[:5], a.iloc[:5])["status"] == "not testable"

def test_ranked_series_is_upper_case_and_deduplicated():
    r = pd.DataFrame({"feature": ["g1", "G1", "g2", "g3"], "stat": [2.0, 1.0, np.nan, -1.0]})
    s = B.bio_ranked(r)
    assert list(s.index) == ["G1", "G3"] and s["G1"] == 2.0
