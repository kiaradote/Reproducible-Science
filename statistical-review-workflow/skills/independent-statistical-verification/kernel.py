
def repro_compare(reported, recomputed, decimals=None, rel_tol=0.02):
    """Classify a reproduced number: 'match' (equal at reported precision), 'minor' (within rel_tol),
    or 'major'. Pass decimals = number of decimals the paper reported. Record the verdict with both values."""
    if decimals is not None and round(recomputed, decimals) == round(reported, decimals):
        return {"reported": reported, "recomputed": recomputed, "class": "match"}
    denom = abs(reported) if reported != 0 else 1.0
    rel = abs(recomputed - reported) / denom
    return {"reported": reported, "recomputed": recomputed, "rel_diff": rel,
            "class": "minor" if rel <= rel_tol else "major"}


def sim_record(null_hypothesis, assumptions, randomization, n_iter, seed, statistic, result,
               interval, limitations):
    """Document a simulation or sensitivity analysis. Every field is required and non-empty.
    A null-model result describes what chance alone produces under the stated null; it is NOT an estimate of
    the share of truly affected units unless that inference is separately justified."""
    fields = {"null_hypothesis": null_hypothesis, "assumptions": assumptions, "randomization": randomization,
              "n_iter": n_iter, "seed": seed, "statistic": statistic, "result": result,
              "interval": interval, "limitations": limitations}
    missing = [k for k, v in fields.items() if v is None or v == ""]
    if missing:
        raise ValueError("sim_record missing: " + ", ".join(missing))
    fields["interpretation_guard"] = ("Null-model output is a chance expectation, not a direct estimate "
                                      "of how many units are genuinely affected.")
    return fields
