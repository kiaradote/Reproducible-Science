
def claim_record(claim_id, claim, location, required_evidence, reproduction="not attempted",
                 method_assessment="not assessed", conclusion="Not assessable", remaining_uncertainty="", quote=""):
    """One row of the claim ledger (the six-part structure). location = page/section/figure/table.
    required_evidence = dict with keys data, labels, formula, threshold (use 'n/a' where truly not needed).
    conclusion in Verified / Contradicted / Supported concern / Not assessable.
    Verified and Contradicted require a performed reproduction; the function refuses them otherwise."""
    verdicts = ("Verified", "Contradicted", "Supported concern", "Not assessable")
    if conclusion not in verdicts:
        raise ValueError("conclusion must be one of " + str(verdicts))
    need = ("data", "labels", "formula", "threshold")
    if not isinstance(required_evidence, dict) or any(k not in required_evidence for k in need):
        raise ValueError("required_evidence must be a dict with keys " + str(need))
    if conclusion in ("Verified", "Contradicted") and str(reproduction).lower().startswith("not attempted"):
        raise ValueError("cannot assign " + conclusion + " without a performed reproduction")
    if conclusion == "Not assessable" and not remaining_uncertainty:
        raise ValueError("Not assessable needs remaining_uncertainty: say what evidence is missing")
    return {"claim_id": claim_id, "claim": claim, "quote": quote, "location": location,
            "required_evidence": required_evidence, "reproduction": reproduction,
            "method_assessment": method_assessment, "conclusion": conclusion,
            "remaining_uncertainty": remaining_uncertainty}


def ledger_to_csv(records, path):
    """Flatten the claim ledger to CSV."""
    import csv
    cols = ["claim_id", "claim", "quote", "location", "data_needed", "labels_needed", "formula_needed",
            "threshold_needed", "reproduction", "method_assessment", "conclusion", "remaining_uncertainty"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in records:
            e = r["required_evidence"]
            w.writerow({"claim_id": r["claim_id"], "claim": r["claim"], "quote": r["quote"],
                        "location": r["location"], "data_needed": e["data"], "labels_needed": e["labels"],
                        "formula_needed": e["formula"], "threshold_needed": e["threshold"],
                        "reproduction": r["reproduction"], "method_assessment": r["method_assessment"],
                        "conclusion": r["conclusion"], "remaining_uncertainty": r["remaining_uncertainty"]})
    return path
