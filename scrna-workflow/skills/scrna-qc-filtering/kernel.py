import json

def qc_record(check, status, detail, group=None):
    """One result row. status is PASS, WARN, FAIL or INFO (same format as the integrity checkpoint, so reports can be merged)."""
    r = {"check": check, "status": status, "detail": detail}
    if group is not None:
        r["group"] = group
    return r

def check_mito(var_names, species="human"):
    """Validate the mitochondrial gene set. Returns (records, clean_mito_gene_list).
    Uses the MT-/mt- prefix only: MTRNR2L* are nuclear pseudogenes and are excluded."""
    pre = "MT-" if species == "human" else "mt-"
    names = list(var_names)
    mito = [g for g in names if g.startswith(pre)]
    pseudo = [g for g in names if g.upper().startswith("MTRNR2L")]
    out = []
    if not mito:
        out.append(qc_record("mito genes found", "FAIL", f"no '{pre}' genes: Ensembl IDs, wrong species, or genes filtered?"))
        return out, mito
    core = {"ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYB", "CO1", "CO2", "CO3", "ATP6", "ATP8"}
    missing = sorted(core - {g[len(pre):].upper() for g in mito})
    out.append(qc_record("mito genes found", "PASS" if not missing else "WARN",
                                f"{len(mito)} genes; missing protein-coding: {missing}" if missing else f"{len(mito)} genes incl. all 13 protein-coding"))
    if pseudo:
        out.append(qc_record("nuclear MTRNR2L pseudogenes", "WARN",
                                    f"{len(pseudo)} present and EXCLUDED; a 'MT-'+'MTRNR' prefix rule would count them as mitochondrial"))
    return out, mito

def check_mito_signal(pct_mt, max_pct_mt=20.0, min_share_nonzero=0.01):
    """The mitochondrial gene set exists, but does it carry signal? If almost no cell has any mitochondrial counts the matrix was filtered
    upstream, and the mitochondrial filter cannot find stressed or dying cells. pct_mt: per-cell percent (array-like, may contain NaN)."""
    import numpy as np
    v = np.asarray(pct_mt, dtype=float)
    v = v[~np.isnan(v)]
    if v.size == 0:
        return qc_record("mito signal present", "WARN", "no mitochondrial percentage could be computed")
    share = float((v > 0).mean())
    if share < min_share_nonzero:
        return qc_record("mito signal present", "WARN",
                         f"only {100 * share:.2f}% of cells have any mitochondrial counts although the genes are listed: mitochondrial reads were probably removed upstream, "
                         f"so the {max_pct_mt:g}% filter removes nothing and cannot find stressed or dying cells")
    return qc_record("mito signal present", "PASS", f"{100 * share:.0f}% of cells have mitochondrial counts (median {np.median(v):.1f}%, 95th percentile {np.percentile(v, 95):.1f}%)")

def qc_filter_report(obs, fail_masks, sample_col=None, max_total_frac=0.2, max_sample_frac=0.5):
    """Per-filter removal counts, combined loss, and per-sample loss. fail_masks: {name: boolean Series, True = cell removed}."""
    import pandas as pd
    out = []
    n = len(obs)
    any_fail = pd.Series(False, index=obs.index)
    for name, mask in fail_masks.items():
        mask = pd.Series(mask, index=obs.index).astype(bool)
        any_fail |= mask
        out.append(qc_record(f"filter: {name}", "INFO", f"removes {int(mask.sum())} of {n} cells ({100 * mask.mean():.1f}%)", "filters"))
    tot = float(any_fail.mean())
    out.append(qc_record("total cells removed", "PASS" if tot <= max_total_frac else "WARN",
                                f"{int(any_fail.sum())} of {n} ({100 * tot:.1f}%); limit {100 * max_total_frac:.0f}%"))
    if sample_col is not None:
        frac = any_fail.groupby(obs[sample_col], observed=True).mean()
        worst = frac.sort_values(ascending=False).head(3)
        flag = bool((frac > max_sample_frac).any())
        out.append(qc_record("per-sample loss", "WARN" if flag else "PASS",
                                    "worst: " + ", ".join(f"{k} {100 * v:.0f}%" for k, v in worst.items()) + f"; limit {100 * max_sample_frac:.0f}%"))
    return out

def qc_catalog():
    """Checks this checkpoint can run: name (or group) -> what it tests."""
    return {
        "mito genes found": "a valid mitochondrial gene set exists in the matrix",
        "mito signal present": "the mitochondrial genes carry counts, so the mitochondrial filter can actually detect stressed cells",
        "filters": "cells removed by each QC filter",
        "total cells removed": "combined cell loss is within the limit",
        "per-sample loss": "no sample loses more than the limit",
    }

def write_qc_report(records, prefix="qc_report"):
    """Write <prefix>.json and <prefix>.md; return the overall status (worst of PASS/WARN/FAIL).
    The report states WHY the verdict was reached, lists every check that ran with what it tests,
    and lists catalogued checks that did not run."""
    cat = qc_catalog()
    rank = {"PASS": 0, "INFO": 0, "WARN": 1, "FAIL": 2}
    overall = ["PASS", "WARN", "FAIL"][max((rank[r["status"]] for r in records), default=0)]
    count = {s: sum(r["status"] == s for r in records) for s in ("PASS", "WARN", "FAIL", "INFO")}
    seen = [r.get("group") or r["check"] for r in records]
    not_run = {k: v for k, v in cat.items() if k not in seen}
    def short(r):
        return r["check"] + " (" + r["detail"][:110].replace("|", "/") + ")"
    fails = [short(r) for r in records if r["status"] == "FAIL"]
    warns = [short(r) for r in records if r["status"] == "WARN"]
    if fails:
        why = f"{len(fails)} of {len(records)} checks failed: " + "; ".join(fails)
        if warns:
            why += f". Also {len(warns)} warning(s): " + "; ".join(warns)
    elif warns:
        why = f"{len(warns)} warning(s), no failures: " + "; ".join(warns)
    else:
        why = f"all {count['PASS']} scored checks passed, none failed or warned ({count['INFO']} informational)"
    if not_run:
        why += f". {len(not_run)} catalogued check(s) did not run."
    with open(prefix + ".json", "w") as f:
        json.dump({"overall": overall, "why": why, "counts": count, "not_run": list(not_run), "checks": records}, f, indent=1)
    lines = [f"# QC report: {overall}", "", f"**Why:** {why}", "",
             f"**Checks run ({len(records)}):** {count['PASS']} pass, {count['WARN']} warn, {count['FAIL']} fail, {count['INFO']} info", "",
             "| status | check | what it tests | result |", "|---|---|---|---|"]
    for r in records:
        what = cat.get(r.get("group") or r["check"], "additional check")
        lines.append(f"| {r['status']} | {r['check']} | {what} | {r['detail'].replace('|', '/')} |")
    if not_run:
        lines += ["", "**Not run:**", ""] + [f"- {k}: {v}" for k, v in not_run.items()]
    with open(prefix + ".md", "w") as f:
        f.write("\n".join(lines) + "\n")
    return overall
