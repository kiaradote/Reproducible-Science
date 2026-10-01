import json

def integrity_record(check, status, detail, group=None):
    """One result row. status is PASS, WARN, FAIL or INFO. group (optional) maps the row to a catalogued check."""
    r = {"check": check, "status": status, "detail": detail}
    if group is not None:
        r["group"] = group
    return r

def check_file(path, expected_sha256=None, expected_bytes=None):
    """Size, SHA-256 and (for .gz) a full decompression pass to catch truncated downloads."""
    import os, gzip, hashlib
    out = []
    size = os.path.getsize(path)
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    digest = h.hexdigest()
    if expected_bytes is None:
        out.append(integrity_record("file size", "INFO", f"{size} bytes; no reference given"))
    else:
        out.append(integrity_record("file size", "PASS" if size == expected_bytes else "FAIL",
                                    f"{size} bytes" if size == expected_bytes else f"{size} bytes, expected {expected_bytes}"))
    if expected_sha256 is not None:
        ok = digest == expected_sha256
        out.append(integrity_record("file sha256", "PASS" if ok else "FAIL", digest if ok else f"{digest} != expected {expected_sha256}"))
    else:
        out.append(integrity_record("file sha256", "INFO", f"{digest} ({size} bytes); no reference given, record this in DATA_SOURCES.md"))
    if path.endswith(".gz"):
        try:
            with gzip.open(path, "rb") as g:
                n = sum(len(c) for c in iter(lambda: g.read(1 << 22), b""))
            out.append(integrity_record("gzip stream", "PASS", f"decompresses fully ({n} bytes)"))
        except Exception as e:
            out.append(integrity_record("gzip stream", "FAIL", f"{type(e).__name__}: {e}"))
    return out

def check_expected(name, observed, expected):
    """Compare any observed value (shape, counts dict, n samples) with a published/expected one."""
    ok = observed == expected
    return integrity_record(name, "PASS" if ok else "FAIL", f"{observed}" if ok else f"observed {observed}, expected {expected}", "published numbers")

def check_alignment(matrix_ids, meta_ids, matrix_labels=None, meta_labels=None):
    """Are the IDs attached to matrix columns the right ones, in the right order?
    matrix_ids: IDs as attached to matrix columns. meta_ids: metadata index.
    matrix_labels: labels read from the matrix file's own label row (same column order).
    meta_labels: pandas Series (index = meta_ids) holding the same label. Any disagreement FAILS:
    a one-position shift only disagrees at sample boundaries, so a % agreement is not acceptable."""
    import pandas as pd
    out = []
    m = list(matrix_ids)
    dup = len(m) - len(set(m))
    out.append(integrity_record("unique cell IDs", "PASS" if dup == 0 else "FAIL", f"{dup} duplicated IDs"))
    only_m, only_meta = set(m) - set(meta_ids), set(meta_ids) - set(m)
    out.append(integrity_record("cell ID sets match metadata", "PASS" if not only_m and not only_meta else "FAIL",
                                f"{len(only_m)} only in matrix, {len(only_meta)} only in metadata"))
    if matrix_labels is None or meta_labels is None:
        out.append(integrity_record("column order vs metadata", "WARN",
                                    "no label row supplied: set equality cannot detect a shifted header; supply matrix_labels and meta_labels"))
        return out
    exp = pd.Series(meta_labels).reindex(m).astype(str).to_numpy()
    got = pd.Series(list(matrix_labels)).astype(str).to_numpy()
    if len(exp) != len(got):
        out.append(integrity_record("column order vs metadata", "FAIL", f"{len(got)} labels for {len(exp)} columns"))
        return out
    bad = int((exp != got).sum())
    detail = f"{len(got)} columns agree in order" if bad == 0 else f"{bad} of {len(got)} columns disagree (shifted or misordered header?)"
    out.append(integrity_record("column order vs metadata", "PASS" if bad == 0 else "FAIL", detail))
    return out

def check_values(X, declared=None, n_cells=2000, seed=0):
    """Detect what the matrix holds. Needs ALL genes (sum tests are invalid after HVG subsetting or scaling).
    Works on a random sample of cells and only on stored (non-zero) values, so it is safe for large sparse matrices.
    declared: one of 'counts', 'log2tpm', 'log1p_cp10k', 'linear_tpm'. Returns (records, detected)."""
    import numpy as np
    import scipy.sparse as sp
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(X.shape[0], size=min(n_cells, X.shape[0]), replace=False))
    A = sp.csr_matrix(X[idx])
    d = np.asarray(A.data, dtype=np.float64)
    out = []
    bad = int(np.isnan(d).sum() + np.isinf(d).sum())
    out.append(integrity_record("NaN/inf values", "PASS" if bad == 0 else "FAIL", f"{bad} in {len(idx)} sampled cells"))
    neg = int((d < 0).sum())
    out.append(integrity_record("negative values", "PASS" if neg == 0 else "WARN", f"{neg} (scaled data legitimately has negatives)"))
    d = np.nan_to_num(d)
    is_int = bool(np.all(d == np.round(d)))
    maxv = float(d.max()) if d.size else 0.0
    islog = maxv < 50          # log-scale values never reach 50; avoids overflow in the 2**x tests

    def med(f):
        B = A.copy()
        B.data = f(d)
        return float(np.median(np.asarray(B.sum(axis=1)).ravel()))

    lo, hi = 0.8e6, 1.2e6
    if neg > 0:
        detected = "scaled_or_unknown"
    elif is_int:
        detected = "counts"
    elif islog and lo < med(lambda v: 2.0 ** v - 1) < hi:
        detected = "log2tpm"
    elif islog and 0.8e4 < med(np.expm1) < 1.2e4:
        detected = "log1p_cp10k"
    elif lo < med(lambda v: v) < hi:
        detected = "linear_tpm"
    else:
        detected = "scaled_or_unknown"
    out.append(integrity_record("detected scale", "INFO", f"{detected}; max {maxv:.2f}, integer-valued {is_int}"))
    if declared is not None:
        out.append(integrity_record("scale matches declaration", "PASS" if declared == detected else "FAIL",
                                    f"declared {declared}, detected {detected}"))
    return out, detected

def integrity_catalog():
    """Checks this checkpoint can run: name (or group) -> what it tests."""
    return {
        "file size": "file size equals the expected size",
        "file sha256": "file checksum equals the reference, or is recorded",
        "gzip stream": "gzip file decompresses fully (not truncated)",
        "unique cell IDs": "no duplicated cell IDs",
        "cell ID sets match metadata": "matrix and metadata describe the same cells",
        "column order vs metadata": "matrix columns are in the metadata's order (label-row test)",
        "NaN/inf values": "no missing or infinite values",
        "negative values": "no negative values (scaled data legitimately has some)",
        "detected scale": "what the matrix holds: counts, TPM, log-normalized or scaled",
        "scale matches declaration": "detected scale equals the scale the source states",
        "published numbers": "shape and label counts equal the published numbers",
    }

def write_integrity_report(records, prefix="integrity_report", scope=None):
    """Write <prefix>.json and <prefix>.md; return the overall status (worst of PASS/WARN/FAIL).
    The report states WHY the verdict was reached, lists every check that ran with what it tests,
    and lists catalogued checks that did not run. scope (optional list of catalogue names) limits the Not run list
    to the checks this stage is expected to run."""
    cat = integrity_catalog()
    if scope is not None:
        cat = {a: b for a, b in cat.items() if a in scope}
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
    lines = [f"# Data integrity report: {overall}", "", f"**Why:** {why}", "",
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
