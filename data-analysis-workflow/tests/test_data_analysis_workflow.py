"""Tests for the shared helpers of the general data-analysis workflow."""
import json
import pytest
import da_kernel as D

def rec(status, name="t", why="because the rule says so"):
    return D.da_record(name, status, f"what {name} checks", f"found for {name}", why)

def test_a_record_must_explain_itself_and_use_a_known_status():
    with pytest.raises(ValueError):
        D.da_record("x", "PASS", "checks something", "3 found", "   ")
    with pytest.raises(ValueError):
        D.da_record("x", "MAYBE", "checks something", "3 found", "reason")
    r = D.da_record("x", "WARN", "checks something", "3 found", "above the limit of 2", group="g")
    assert r["group"] == "g" and r["why"] == "above the limit of 2"

def test_overall_is_the_worst_status_and_info_never_counts():
    assert D.da_overall([rec("PASS"), rec("INFO")]) == "PASS"
    assert D.da_overall([rec("PASS"), rec("WARN"), rec("INFO")]) == "WARN"
    assert D.da_overall([rec("WARN"), rec("FAIL")]) == "FAIL"
    assert D.da_overall([]) == "PASS"

def test_ratings_cover_every_case_in_plain_words():
    assert D.da_rating([rec("PASS"), rec("INFO")])[0] == "Good"
    assert D.da_rating([rec("PASS"), rec("WARN"), rec("WARN", "u")])[0] == "Usable, with 2 caveat(s)"
    assert D.da_rating([rec("FAIL")])[0] == "Not reliable until fixed"
    r, note = D.da_rating([rec("FAIL")], accepted_reason="no mitochondrial genes in this assay")
    assert r == "Accepted despite a failed check" and "no mitochondrial genes" in note
    assert D.da_rating([rec("INFO")])[0] == "Not applicable"
    assert D.da_rating(None)[0] == "Not done" and D.da_rating([rec("PASS")], ran=False)[0] == "Not done"

def test_card_shows_every_problem_a_limited_number_of_passes_and_the_count_of_the_rest():
    recs = [rec("FAIL", "f1"), rec("WARN", "w1")] + [rec("PASS", f"p{i}") for i in range(9)] + [rec("INFO", "i1"), rec("INFO", "i2"), rec("INFO", "i3")]
    card = D.da_card("Checkpoint 2: QC", "which measurements are damaged?", recs, "usable with care", not_run=["replicate concordance (one replicate)"],
                     caveats=["batch not balanced"], go_back="outputs/02_qc.csv", max_pass=5)
    assert card.startswith("**Checkpoint 2: QC: FAIL, not reliable until fixed**")
    assert "found for f1" in card and "found for w1" in card and "found for p4" in card and "found for p5" not in card
    assert "| What it checked | What we found | Result and why |" in card
    assert "FAIL: because the rule says so" in card
    assert "5 more test(s) not shown (4 passed, 1 informational)" in card
    assert "*Conclusion:* usable with care" in card and "*Did not run:* replicate concordance" in card and "*Carried caveats:* batch not balanced" in card and "*Go back:* outputs/02_qc.csv" in card
    assert card.index("f1") < card.index("w1") < card.index("p0")

def test_report_files_list_what_did_not_run_and_the_reasons(tmp_path):
    catalog = {"duplicate samples": "no two samples are near copies", "batch balance": "each batch holds every group", "outliers": "no sample is far from the rest"}
    recs = [D.da_record("duplicate samples", "PASS", "no two samples are near copies", "none above 0.999", "limit is 0.999"),
            D.da_record("outliers", "WARN", "no sample is far from the rest", "2 of 24 far away", "we expect none")]
    overall = D.da_write_report(recs, catalog, str(tmp_path / "qc"), "QC report", asks="which samples are damaged?", conclusion="keep but watch", caveats=["two weak samples kept"])
    assert overall == "WARN"
    js = json.load(open(tmp_path / "qc.json"))
    assert js["not_run"] == ["batch balance"] and js["rating"] == "Usable, with 1 caveat(s)" and js["checks"][1]["why"] == "we expect none"
    md = (tmp_path / "qc.md").read_text(encoding="utf8")
    assert "**Quality rating:** Usable, with 1 caveat(s)" in md and "- batch balance: each batch holds every group" in md and "| WARN | outliers |" in md and "two weak samples kept" in md

def test_run_summary_rates_the_run_by_its_worst_checkpoint(tmp_path):
    base = [{"stage": s, "rating": "Good", "note": "ok", "go_back": f"outputs/{s}.csv"} for s in D.CHECKPOINTS]
    assert D.da_run_summary(base, str(tmp_path / "a")).startswith("Good")
    base[1]["rating"] = "Usable, with 2 caveat(s)"
    assert D.da_run_summary(base, str(tmp_path / "b")).startswith("Usable")
    base[4]["rating"] = "Not applicable"
    assert D.da_run_summary(base, str(tmp_path / "c")).startswith("Usable")
    base[5]["rating"] = "Not done"
    assert D.da_run_summary(base, str(tmp_path / "d")).startswith("Incomplete")
    base[2]["rating"] = "Not reliable until fixed"
    assert D.da_run_summary(base, str(tmp_path / "e")).startswith("Not reliable")
    md = (tmp_path / "e.md").read_text(encoding="utf8")
    assert md.index("## Quality check") < md.index("| checkpoint |") and "3. Normalization and features | Not reliable until fixed" in md

def test_restart_logic_notices_every_kind_of_change(tmp_path):
    out = tmp_path / "02_qc.csv"
    out.write_text("a,b\n1,2\n")
    params = {"min_features": 200, "max_pct": 20}
    man = {"qc": {"params_hash": D.da_hash_params(params), "input_hash": "abc", "output": str(out), "output_hash": D.da_hash_file(str(out))}}
    assert D.da_is_current(man, "qc", params, "abc") == (True, "current")
    assert D.da_is_current(man, "normalize", params, "abc") == (False, "never run")
    assert D.da_is_current(man, "qc", {**params, "max_pct": 25}, "abc") == (False, "parameters changed")
    assert D.da_is_current(man, "qc", params, "different") == (False, "input changed")
    out.write_text("a,b\n9,9\n")
    assert D.da_is_current(man, "qc", params, "abc") == (False, "output file changed since it was written")
    out.unlink()
    assert D.da_is_current(man, "qc", params, "abc") == (False, "output file missing")
    assert D.da_hash_params({"a": 1, "b": 2}) == D.da_hash_params({"b": 2, "a": 1})
