"""The worked example doubles as a test that the general rules behave sensibly on data with known planted effects."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "examples"))
import proteomics_example as PE
import da_kernel as D

def test_example_finds_the_planted_effects_and_rates_each_checkpoint():
    cards, entries, ratings = PE.run()
    assert sorted(cards) == [1, 2, 3, 4, 5, 6]
    assert all(c.startswith("**Checkpoint") and "| What it checked | What we found | Result and why |" in c for c in cards.values())
    rate = {e["stage"]: e["rating"] for e in entries}
    assert rate["compare"] == "Good" and rate["interpret"] == "Good" and rate["intake"].startswith("Usable") and rate["qc"].startswith("Usable")
    assert "0 of 20 null proteins called" in cards[5] and "of 40 planted found" in cards[5]
    assert "S07" in cards[2] and "S10, S11" in cards[1]                         # the weak sample and the duplicate pair are caught
    assert "batch explains" in cards[4] and "WARN" in cards[4]
    assert D.da_run_summary(entries, str(Path(__file__).parent / "_tmp_summary")).startswith("Usable")
    (Path(__file__).parent / "_tmp_summary.md").unlink()
