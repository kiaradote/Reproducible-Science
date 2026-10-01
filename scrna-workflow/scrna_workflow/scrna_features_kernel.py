import json

def feat_record(check, status, detail, group=None):
    """One result row. status is PASS, WARN, FAIL or INFO."""
    r = {"check": check, "status": status, "detail": detail}
    if group is not None:
        r["group"] = group
    return r

def feat_catalog():
    """Checks this checkpoint can run: name -> what it tests."""
    return {
        "method matches scale": "the normalization fits what the matrix holds (scaled or unknown data is refused)",
        "per-cell totals comparable": "after normalization every cell has about the same linear-scale total",
        "finite and non-negative values": "no NaN, infinite or negative values after the transform",
        "sparsity preserved": "the matrix was not densified or padded",
        "input values kept": "the original values survive (counts layer, or X unchanged)",
        "depth spread across samples": "median depth does not differ wildly between samples",
        "dominant-gene cells": "few cells are dominated by a single gene",
        "HVG count": "the requested number of unique variable genes was selected",
        "technical genes among HVGs": "mitochondrial, ribosomal and haemoglobin genes do not dominate the selection",
        "low-detection HVGs": "selected genes are detected in enough cells",
        "HVG stability (two halves)": "two random halves of the cells select similar genes",
        "batch robustness": "selected genes are variable in many batches, not one",
        "mean expression coverage": "selected genes are not concentrated among the lowest-expressed genes",
    }

def linear_stats(X, base=None):
    """Per-cell linear-scale total and largest single-gene share, from log-scale values (base None = natural log1p)."""
    import numpy as np
    import scipy.sparse as sp
    tot = np.zeros(X.shape[0])
    top = np.zeros(X.shape[0])
    for i in range(0, X.shape[0], 2000):
        B = sp.csr_matrix(X[i:i + 2000])
        d = B.data.astype(np.float64)
        d = np.expm1(d) if base is None else base ** d - 1
        B = sp.csr_matrix((d, B.indices, B.indptr), shape=B.shape)
        tot[i:i + 2000] = np.asarray(B.sum(axis=1)).ravel()
        top[i:i + 2000] = B.max(axis=1).toarray().ravel() / np.maximum(tot[i:i + 2000], 1e-12)
    return tot, top

def feat_normalize(adata, scale, target_sum=10000.0):
    """Normalize in place according to the detected scale. Returns (records, ok, nnz_before).
    counts: CP(target_sum) + log1p, raw counts kept in layers['counts']. linear_tpm: log1p, input kept in layers['input'].
    log2tpm and log1p_cp10k are already normalized: X is left unchanged and the log base is recorded in uns['log1p']."""
    import scipy.sparse as sp
    import scanpy as sc
    X = sp.csr_matrix(adata.X)
    adata.X = X
    nnz0 = int(X.nnz)
    if scale == "counts":
        adata.layers["counts"] = X.copy()
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        method = f"normalize to {int(target_sum)} per cell, then natural log1p; raw counts kept in layers['counts']"
    elif scale == "linear_tpm":
        adata.layers["input"] = X.copy()
        sc.pp.log1p(adata)
        method = "natural log1p of TPM; input kept in layers['input']"
    elif scale == "log2tpm":
        adata.uns["log1p"] = {"base": 2.0}
        method = "none needed (already log2(TPM+1)); log base 2 recorded so variable-gene selection uses the right scale"
    elif scale == "log1p_cp10k":
        adata.uns["log1p"] = {"base": None}
        method = "none needed (already natural log1p of per-cell-normalized values)"
    else:
        return [feat_record("method matches scale", "FAIL", f"scale '{scale}' cannot be normalized here: scaled or unrecognized values")], False, nnz0
    adata.uns["normalization"] = {"scale_in": scale, "method": method}
    return [feat_record("method matches scale", "PASS", f"{scale}: {method}")], True, nnz0

def feat_check_normalization(adata, nnz0, sample_col=None, max_cv=0.1, max_depth_ratio=5.0, max_dominant_frac=0.01):
    """Did normalization make cells comparable? Thresholds are defaults, not literature values."""
    import numpy as np
    import pandas as pd
    out = []
    base = adata.uns.get("log1p", {}).get("base")
    tot, top = linear_stats(adata.X, base)
    cv = float(tot.std() / max(tot.mean(), 1e-12))
    out.append(feat_record("per-cell totals comparable", "PASS" if cv < max_cv else "WARN",
                           f"linear-scale total median {np.median(tot):.0f}, CV {cv:.3f} (limit {max_cv})"))
    d = np.asarray(adata.X.data, dtype=np.float64)
    nbad, nneg = int((~np.isfinite(d)).sum()), int((d < 0).sum())
    out.append(feat_record("finite and non-negative values", "PASS" if nbad + nneg == 0 else "FAIL", f"{nbad} NaN/inf, {nneg} negative"))
    out.append(feat_record("sparsity preserved", "PASS" if adata.X.nnz == nnz0 else "WARN", f"{adata.X.nnz} stored values, {nnz0} before"))
    kept = [k for k in ("counts", "input") if k in adata.layers]
    if kept:
        L = adata.layers[kept[0]]
        ok = L.nnz == nnz0 and (kept[0] != "counts" or bool(np.all(L.data == np.round(L.data))))
        out.append(feat_record("input values kept", "PASS" if ok else "FAIL", f"layers['{kept[0]}'] holds {L.nnz} stored values" + (", integer-valued" if kept[0] == "counts" else "")))
    else:
        out.append(feat_record("input values kept", "PASS", "no transformation applied: X still holds the input values"))
    out.append(feat_record("dominant-gene cells", "PASS" if float((top > 0.5).mean()) <= max_dominant_frac else "WARN",
                           f"{int((top > 0.5).sum())} of {len(top)} cells have one gene above 50% of their total (limit {100 * max_dominant_frac:.0f}% of cells)"))
    if sample_col is not None:
        if "counts" in adata.layers:
            depth, what = np.asarray(adata.layers["counts"].sum(axis=1)).ravel(), "total counts"
        else:
            depth, what = np.diff(adata.X.indptr), "genes detected (depth is not recoverable from normalized values)"
        med = pd.Series(depth).groupby(adata.obs[sample_col].to_numpy()).median()
        ratio = float(med.max() / max(med.min(), 1e-12))
        out.append(feat_record("depth spread across samples", "PASS" if ratio <= max_depth_ratio else "WARN",
                               f"median {what} differs {ratio:.1f}-fold between samples (limit {max_depth_ratio}-fold)"))
    return out

def feat_select(adata, n_top_genes=2000, batch_key=None):
    """Select highly variable genes in place. Counts: seurat_v3 on layers['counts']. Otherwise: seurat on log-scale X
    (scanpy reads the log base from uns['log1p']). Returns (flavor, layer)."""
    import scanpy as sc
    flavor, layer = ("seurat_v3", "counts") if "counts" in adata.layers else ("seurat", None)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor, layer=layer, batch_key=batch_key, subset=False, inplace=True)
    adata.uns["feature_selection"] = {"flavor": flavor, "n_top_genes": n_top_genes, "batch_key": str(batch_key)}
    return flavor, layer

def feat_check_selection(adata, flavor, layer, n_top_genes=2000, species="human", batch_key=None, seed=0, max_cells=5000,
                         max_tech_frac=0.1, max_lowdet_frac=0.1, min_jaccard=0.5):
    """Are the selected genes informative? Thresholds are defaults, not literature values."""
    import re
    import numpy as np
    import scanpy as sc
    out = []
    hv = adata.var["highly_variable"].to_numpy()
    n_sel = int(hv.sum())
    uniq = adata.var_names.is_unique
    st = "FAIL" if (n_sel == 0 or not uniq) else ("PASS" if n_sel == n_top_genes else "WARN")
    out.append(feat_record("HVG count", st, f"{n_sel} selected, {n_top_genes} requested, gene names unique: {uniq}"))
    names = adata.var_names[hv]
    pats = {"mitochondrial": r"^mt-", "ribosomal": r"^rp[sl]\d", "haemoglobin": r"^hb[ab]"}
    counts = {k: int(sum(bool(re.match(p, g, re.I)) for g in names)) for k, p in pats.items()}
    frac = sum(counts.values()) / max(n_sel, 1)
    out.append(feat_record("technical genes among HVGs", "PASS" if frac <= max_tech_frac else "WARN",
                           f"{sum(counts.values())} of {n_sel} ({100 * frac:.1f}%): " + ", ".join(f"{k} {v}" for k, v in counts.items()) + f"; limit {100 * max_tech_frac:.0f}%"))
    det = adata.X.getnnz(axis=0) / adata.n_obs
    low = float((det[hv] < 0.01).mean()) if n_sel else 0.0
    out.append(feat_record("low-detection HVGs", "PASS" if low <= max_lowdet_frac else "WARN",
                           f"{100 * low:.1f}% of selected genes are detected in under 1% of cells; limit {100 * max_lowdet_frac:.0f}%"))
    means = np.asarray(adata.X.mean(axis=0)).ravel()
    cut = np.percentile(means[means > 0], 10) if (means > 0).any() else 0
    lowmean = float((means[hv] <= cut).mean()) if n_sel else 0.0
    out.append(feat_record("mean expression coverage", "PASS" if lowmean <= 0.3 else "WARN",
                           f"{100 * lowmean:.1f}% of selected genes are in the lowest-expressed 10% of genes; limit 30%"))
    try:
        rng = np.random.default_rng(seed)
        perm = rng.permutation(adata.n_obs)
        half = min(adata.n_obs // 2, max_cells)
        sets = []
        for idx in (perm[:half], perm[half:2 * half]):
            df = sc.pp.highly_variable_genes(adata[np.sort(idx)], n_top_genes=n_top_genes, flavor=flavor, layer=layer, batch_key=batch_key, inplace=False)
            sets.append(set(df.index[df["highly_variable"]]))
        jac = len(sets[0] & sets[1]) / max(len(sets[0] | sets[1]), 1)
        out.append(feat_record("HVG stability (two halves)", "PASS" if jac >= min_jaccard else "WARN",
                               f"overlap (Jaccard) {jac:.2f} between selections on two random halves of {half} cells; limit {min_jaccard}"))
    except Exception as e:
        out.append(feat_record("HVG stability (two halves)", "WARN", f"could not run: {type(e).__name__}: {str(e)[:120]}"))
    if batch_key is not None and "highly_variable_nbatches" in adata.var:
        nb = int(adata.obs[batch_key].nunique())
        need = max(1, int(np.ceil(nb / 2)))
        robust = float((adata.var["highly_variable_nbatches"].to_numpy()[hv] >= need).mean()) if n_sel else 0.0
        out.append(feat_record("batch robustness", "PASS" if robust >= 0.5 else "WARN",
                               f"{100 * robust:.1f}% of selected genes are variable in at least {need} of {nb} batches; limit 50%"))
    return out

def feat_plot(adata, path):
    """Mean versus normalized variability, selected genes highlighted."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ycol = next((c for c in ("dispersions_norm", "variances_norm") if c in adata.var), None)
    fig, ax = plt.subplots(figsize=(5, 4))
    if ycol is not None:
        hv = adata.var["highly_variable"].to_numpy()
        ax.scatter(adata.var["means"][~hv], adata.var[ycol][~hv], s=3, c="0.7", label="other genes")
        ax.scatter(adata.var["means"][hv], adata.var[ycol][hv], s=3, c="tab:red", label="selected")
        ax.set_xscale("symlog", linthresh=0.01)
        ax.set_xlabel("mean expression")
        ax.set_ylabel(ycol.replace("_", " "))
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)

def feat_run(adata, scale, species="human", n_top_genes=2000, target_sum=10000.0, batch_key=None, sample_col=None, seed=0, plot_path=None):
    """Normalize, check, select variable genes, check. Returns (adata, records). Stops after the first step if scale is refused."""
    recs, ok, nnz0 = feat_normalize(adata, scale, target_sum)
    if not ok:
        return adata, recs
    recs += feat_check_normalization(adata, nnz0, sample_col)
    flavor, layer = feat_select(adata, n_top_genes, batch_key)
    recs += feat_check_selection(adata, flavor, layer, n_top_genes, species, batch_key, seed)
    if plot_path:
        feat_plot(adata, plot_path)
    return adata, recs

def write_feat_report(records, prefix="features_report", scope=None):
    """Write <prefix>.json and <prefix>.md; return the overall status (worst of PASS/WARN/FAIL).
    The report states WHY the verdict was reached, lists every check that ran with what it tests,
    and lists catalogued checks that did not run."""
    cat = feat_catalog()
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
    lines = [f"# Normalization and feature selection report: {overall}", "", f"**Why:** {why}", "",
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
