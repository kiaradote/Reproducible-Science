"""Checkpoint 7 helpers: control vs treatment comparison. Prefix cmp_. Needs numpy, scipy, pandas, statsmodels, scikit-learn, matplotlib; pydeseq2 is used when installed."""
import json
import os

def cmp_record(check, status, detail, group=None):
    r = {"check": check, "status": status, "detail": detail}
    if group:
        r["group"] = group
    return r

def cmp_catalog():
    """Every check this checkpoint can run, with what it tests. Checks that did not run are listed as 'Not run'."""
    return {
        "groups": "exactly two conditions are compared and both contain cells",
        "replicates": "each condition has enough independent samples (donors, patients, animals) for a test; cells of one sample are not replicates",
        "pairing": "samples measured in both conditions are compared with themselves (paired design)",
        "condition vs batch": "the condition is not perfectly confounded with a batch or library, which would make a treatment effect indistinguishable from a technical one",
        "clusters follow condition": "clusters are not almost entirely made of one condition, which would mean the treatment itself drives the clustering",
        "cell types follow condition": "a cell-type name is not almost entirely used in one condition only",
        "testable cell types": "each compared cell type has enough cells from enough samples in both conditions",
        "unlabelled cells": "the share of cells without a cell-type name is small",
        "abundance": "the share of each cell type per sample is compared between conditions by two different models that agree",
        "differential expression": "genes are compared between conditions within each cell type on sample-level (pseudobulk) totals, not on single cells",
        "null check": "swapping the condition labels within samples gives (almost) no discoveries",
        "sample consistency": "differentially expressed genes move in the same direction in most individual samples",
        "state change": "cells of one type can be told apart by condition in samples the classifier never saw",
        "pathway enrichment": "ranked differential-expression results are tested against gene sets with a pre-ranked method",
        "transcription factor activity": "regulator activity inferred from a transcription-factor network (not implemented in this version: needs the OmniPath network and the decoupler package)",
        "cell-cell communication": "ligand-receptor signalling between cell types (not implemented in this version: needs a ligand-receptor resource)",
        "spatial association": "cell types or genes that co-vary in tissue position (needs coordinates in obsm['spatial'])",
        "gene modules": "co-expressed gene modules that change with the condition (use the gene programs of checkpoint 8 for the descriptive version)",
    }

def cmp_bh(p):
    import numpy as np
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    q = p[ok]
    n = len(q)
    if n == 0:
        return out
    o = np.argsort(q)
    v = q[o] * n / (np.arange(n) + 1)
    v = np.minimum.accumulate(v[::-1])[::-1]
    r = np.empty(n)
    r[o] = np.minimum(v, 1.0)
    out[ok] = r
    return out

def cmp_named_types(adata, type_col="cell_type"):
    import pandas as pd
    lab = adata.obs[type_col].astype(str)
    return [t for t in lab.value_counts().index if t not in ("ambiguous", "unassigned")]

def cmp_design(adata, condition_col, sample_col, type_col="cell_type", batch_col=None, cluster_col="cluster", groups=None, min_units=3, min_cells=10, min_pure=0.9):
    """The first gate. Returns (records, design). design holds ref, test, paired, units table, testable types, excluded types with reasons, and the types that follow the condition.
    FAIL means the comparison cannot be done (no two groups, fewer than two samples in a group, no testable cell type)."""
    import numpy as np
    import pandas as pd
    recs, des = [], {}
    obs = adata.obs
    for c in (condition_col, sample_col, type_col):
        if c not in obs:
            return [cmp_record("groups", "FAIL", f"obs['{c}'] not found: name the condition, sample and cell-type columns")], des
    cond = obs[condition_col].astype(str)
    levels = sorted(cond.unique()) if groups is None else [str(g) for g in groups]
    if len(levels) != 2 or any(l not in set(cond) for l in levels):
        return [cmp_record("groups", "FAIL", f"need exactly two conditions, found {sorted(cond.unique())}: pass groups=(control, treatment)")], des
    ref, test = levels
    keep = cond.isin(levels).to_numpy()
    o = obs.loc[keep, [sample_col, condition_col, type_col] + ([cluster_col] if cluster_col in obs else []) + ([batch_col] if batch_col and batch_col in obs else [])].copy()
    o[condition_col] = o[condition_col].astype(str)
    o[sample_col] = o[sample_col].astype(str)
    o[type_col] = o[type_col].astype(str)
    n_ref, n_test = int((o[condition_col] == ref).sum()), int((o[condition_col] == test).sum())
    recs.append(cmp_record("groups", "INFO", f"control = {ref} ({n_ref:,} cells), treatment = {test} ({n_test:,} cells); {int((~keep).sum()):,} cells of other conditions are left out"))
    u = o.groupby([sample_col, condition_col]).size().unstack(fill_value=0).reindex(columns=[ref, test], fill_value=0)
    ok = u >= min_cells
    n_u_ref, n_u_test = int(ok[ref].sum()), int(ok[test].sum())
    both = [s for s in u.index if ok.loc[s, ref] and ok.loc[s, test]]
    lowest = min(n_u_ref, n_u_test)
    if lowest < 2:
        recs.append(cmp_record("replicates", "FAIL", f"{n_u_ref} sample(s) in {ref} and {n_u_test} in {test} with at least {min_cells} cells: a test needs at least 2 per condition (cells of one sample are not replicates)"))
    elif lowest < min_units:
        recs.append(cmp_record("replicates", "WARN", f"{n_u_ref} sample(s) in {ref} and {n_u_test} in {test}: fewer than {min_units} per condition, so p-values are fragile and only large effects can be found"))
    else:
        recs.append(cmp_record("replicates", "PASS", f"{n_u_ref} sample(s) in {ref} and {n_u_test} in {test} with at least {min_cells} cells each"))
    paired = len(both) >= min_units and len(both) == len(u.index)
    if paired:
        recs.append(cmp_record("pairing", "PASS", f"all {len(both)} samples were measured in both conditions: each sample is compared with itself (sample is a blocking factor)"))
    elif both:
        recs.append(cmp_record("pairing", "WARN", f"{len(both)} of {len(u)} samples appear in both conditions: treated as unpaired, which loses power"))
    else:
        recs.append(cmp_record("pairing", "INFO", "no sample appears in both conditions: unpaired comparison between different samples"))
    batch_confounded = False
    if batch_col and batch_col in o:
        bt = o.groupby([batch_col, condition_col]).size().unstack(fill_value=0)
        single = int(((bt > 0).sum(axis=1) == 1).sum())
        batch_confounded = single == len(bt) and len(bt) >= 1
        if batch_confounded:
            recs.append(cmp_record("condition vs batch", "WARN", f"every {batch_col} holds one condition only ({len(bt)} {batch_col} values): a treatment effect cannot be separated from a {batch_col} effect; "
                                   "pairing by sample removes donor differences but not this. Treat results as treatment-or-library effects"))
        elif single:
            recs.append(cmp_record("condition vs batch", "WARN", f"{single} of {len(bt)} {batch_col} values hold one condition only"))
        else:
            recs.append(cmp_record("condition vs batch", "PASS", f"every {batch_col} contains both conditions"))
    if cluster_col in o:
        ct = o.groupby([cluster_col, condition_col]).size().unstack(fill_value=0)
        share = ct.max(axis=1) / ct.sum(axis=1)
        pure = share[share >= min_pure]
        frac = len(pure) / max(len(share), 1)
        recs.append(cmp_record("clusters follow condition", "WARN" if frac > 0.5 else "PASS",
                               f"{len(pure)} of {len(share)} clusters have {min_pure:.0%} or more of their cells from one condition" + ("; the treatment itself drives the clustering, so cell types and the condition are partly the same thing" if frac > 0.5 else "")))
    tt = o.groupby([type_col, condition_col]).size().unstack(fill_value=0).reindex(columns=[ref, test], fill_value=0)
    tt = tt[~tt.index.isin(["ambiguous", "unassigned"])]
    tshare = tt.max(axis=1) / tt.sum(axis=1)
    followers = [t for t in tt.index if tshare[t] >= min_pure and tt.loc[t].sum() >= 3 * min_cells]
    if followers:
        recs.append(cmp_record("cell types follow condition", "WARN", "; ".join(f"{t} ({tt.loc[t, ref]} {ref} / {tt.loc[t, test]} {test})" for t in followers) +
                               ": the name may have moved with the treatment instead of the cells, so these types are not comparable between conditions"))
    else:
        recs.append(cmp_record("cell types follow condition", "PASS", f"no cell type has {min_pure:.0%} or more of its cells in one condition"))
    cs = o.groupby([sample_col, condition_col, type_col]).size()
    testable, excluded = [], {}
    for t in tt.index:
        c = cs.xs(t, level=type_col) if t in cs.index.get_level_values(type_col) else pd.Series(dtype=int)
        piv = c.unstack(condition_col).reindex(index=u.index, columns=[ref, test]).fillna(0) if len(c) else pd.DataFrame(0, index=u.index, columns=[ref, test])
        okc = piv >= min_cells
        nb = int((okc[ref] & okc[test]).sum()) if paired else min(int(okc[ref].sum()), int(okc[test].sum()))
        if nb >= min_units:
            testable.append(t)
        else:
            excluded[t] = f"only {int(okc[ref].sum())} {ref} and {int(okc[test].sum())} {test} samples have {min_cells}+ cells" + (f" ({nb} paired)" if paired else "")
    if not testable:
        recs.append(cmp_record("testable cell types", "FAIL", "no cell type has enough cells in enough samples in both conditions: " + "; ".join(f"{k}: {v}" for k, v in list(excluded.items())[:4])))
    else:
        recs.append(cmp_record("testable cell types", "PASS" if not excluded else "WARN", f"{len(testable)} of {len(tt)} cell types can be compared (at least {min_units} samples with {min_cells}+ cells" + (" in both conditions)" if paired else " per condition)") +
                               (". Not compared: " + "; ".join(f"{k} ({v})" for k, v in excluded.items()) if excluded else "")))
    unl = float(o[type_col].isin(["ambiguous", "unassigned"]).mean())
    recs.append(cmp_record("unlabelled cells", "INFO" if unl < 0.2 else "WARN", f"{unl:.1%} of the compared cells are ambiguous or unassigned and are left out of per-type results" + ("" if unl < 0.2 else "; above 20%")))
    des = {"ref": ref, "test": test, "paired": bool(paired), "units": u, "testable": testable, "excluded": excluded, "followers": followers, "batch_confounded": bool(batch_confounded),
           "n_cells": int(keep.sum()), "condition_col": condition_col, "sample_col": sample_col, "type_col": type_col, "min_units": min_units, "min_cells": min_cells, "type_table": tt}
    return recs, des

def cmp_matrix(units, paired):
    """Design matrix for sample-level models: intercept, treatment (1 = treated), and sample dummies when paired."""
    import numpy as np
    cols = [np.ones(len(units)), units["treated"].to_numpy(float)]
    if paired:
        for s in sorted(units["sample"].unique())[1:]:
            cols.append((units["sample"] == s).to_numpy(float))
    return np.column_stack(cols)

def cmp_ols(Y, X, col=1):
    """Ordinary least squares for every column of Y (units x genes) with design X. Returns effect, t, p of column col."""
    import numpy as np
    from scipy import stats
    n, p = X.shape
    df = n - p
    if df < 1:
        z = np.full(Y.shape[1], np.nan)
        return z, z, z
    XtXi = np.linalg.pinv(X.T @ X)
    B = XtXi @ X.T @ Y
    R = Y - X @ B
    s2 = (R ** 2).sum(axis=0) / df
    se = np.sqrt(np.maximum(s2 * XtXi[col, col], 1e-300))
    t = B[col] / se
    return B[col], t, 2 * stats.t.sf(np.abs(t), df)

def cmp_abundance(adata, des, min_type_cells=30):
    """Share of each named cell type per sample, compared between conditions. Model 1: quasi-binomial regression on the cell counts (treatment + sample when paired).
    Model 2: ordinary regression on centred log-ratio shares. Returns a DataFrame: type, share_ref, share_test, log2_ratio, p_quasibinomial, fdr, p_clr, fdr_clr, agree, caution."""
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    from scipy import stats
    obs = adata.obs
    sc_, cc_, tc_ = des["sample_col"], des["condition_col"], des["type_col"]
    o = obs[[sc_, cc_, tc_]].astype(str)
    o = o[o[cc_].isin([des["ref"], des["test"]])]
    tot = o.groupby([sc_, cc_]).size().rename("n").reset_index()
    tot = tot[tot["n"] >= des["min_cells"]].rename(columns={sc_: "sample", cc_: "condition"})
    tot["treated"] = (tot["condition"] == des["test"]).astype(int)
    paired = des["paired"]
    if paired:
        good = set(tot.groupby("sample").size()[lambda x: x == 2].index)
        tot = tot[tot["sample"].isin(good)]
    tot = tot.sort_values(["sample", "condition"]).reset_index(drop=True)
    types = [t for t in cmp_named_types(adata, tc_) if (o[tc_] == t).sum() >= min_type_cells]
    cnt = o.groupby([sc_, cc_, tc_]).size().unstack(fill_value=0).reindex(columns=types, fill_value=0)
    K = np.array([[cnt.loc[(r.sample, r.condition), t] if (r.sample, r.condition) in cnt.index else 0 for t in types] for r in tot.itertuples()], dtype=float)
    n = tot["n"].to_numpy(float)
    X = cmp_matrix(tot, paired)
    lk = np.log(K + 0.5)
    clr = lk - lk.mean(axis=1, keepdims=True)
    rows = []
    for j, t in enumerate(types):
        k = K[:, j]
        pr = k / n
        r = {"type": t, "share_ref": float(pr[tot["treated"] == 0].mean()), "share_test": float(pr[tot["treated"] == 1].mean()), "n_samples": int(tot["sample"].nunique())}
        r["log2_ratio"] = float(np.log2((r["share_test"] + 1e-4) / (r["share_ref"] + 1e-4)))
        try:
            fit = sm.GLM(np.column_stack([k, n - k]), X, family=sm.families.Binomial()).fit()
            dfr = max(int(fit.df_resid), 1)
            phi = max(float(fit.pearson_chi2) / dfr, 1.0)      # quasi-binomial dispersion, never below the binomial value
            tval = float(fit.params[1] / (fit.bse[1] * np.sqrt(phi)))
            r["p_quasibinomial"] = float(2 * stats.t.sf(abs(tval), dfr))
            r["dispersion"] = phi
        except Exception:
            r["p_quasibinomial"] = np.nan
        _, _, pc = cmp_ols(clr[:, [j]], X)
        r["p_clr"] = float(pc[0])
        rows.append(r)
    df = pd.DataFrame(rows)
    if len(df):
        df["fdr"] = cmp_bh(df["p_quasibinomial"].to_numpy())
        df["fdr_clr"] = cmp_bh(df["p_clr"].to_numpy())
        df["agree"] = (df["fdr"] < 0.05) & (df["fdr_clr"] < 0.05)
        df["caution"] = df["type"].isin(des["followers"])
        df = df.sort_values("p_quasibinomial").reset_index(drop=True)
    return df

def cmp_pseudobulk(adata, type_mask, des, layer="counts"):
    """Sum counts (or average log expression when no counts exist) per sample and condition for the cells of one type. Returns (matrix DataFrame units x genes, units DataFrame, kind)."""
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    use_counts = layer in adata.layers if layer else False
    X = adata.layers[layer] if use_counts else adata.X
    obs = adata.obs
    cond = obs[des["condition_col"]].astype(str).to_numpy()
    keep = type_mask & np.isin(cond, [des["ref"], des["test"]])
    idx = np.where(keep)[0]
    key = np.array([f"{s}||{c}" for s, c in zip(obs[des["sample_col"]].astype(str).to_numpy()[idx], cond[idx])])
    uniq, inv = np.unique(key, return_inverse=True)
    M = sp.csr_matrix((np.ones(len(idx)), (inv, np.arange(len(idx)))), shape=(len(uniq), len(idx)))
    sub = X[idx]
    S = M @ (sub if sp.issparse(sub) else sp.csr_matrix(sub))
    S = S.toarray()
    ncell = np.bincount(inv, minlength=len(uniq))
    if not use_counts:
        S = S / np.maximum(ncell[:, None], 1)
    units = pd.DataFrame({"unit": uniq, "sample": [u.split("||")[0] for u in uniq], "condition": [u.split("||")[1] for u in uniq], "n_cells": ncell})
    units["treated"] = (units["condition"] == des["test"]).astype(int)
    ok = (units["n_cells"] >= des["min_cells"]).to_numpy()
    units = units[ok].reset_index(drop=True)
    mat = pd.DataFrame(S[ok], index=units["unit"], columns=[str(g).upper() for g in adata.var_names])
    if des["paired"]:
        have = units.groupby("sample")["treated"].nunique()
        good = set(have[have == 2].index)
        sel = units["sample"].isin(good).to_numpy()
        units, mat = units[sel].reset_index(drop=True), mat[sel]
    return mat, units, ("counts" if use_counts else "mean log expression")

def cmp_logcpm(mat):
    import numpy as np
    m = mat.to_numpy(float)
    return np.log2(m / np.maximum(m.sum(axis=1, keepdims=True), 1) * 1e6 + 1)

def cmp_de_type(mat, units, des, kind, method="auto", min_lfc=0.5, fdr=0.05):
    """Differential expression of one cell type on sample-level totals. method 'deseq2' (pydeseq2, negative binomial, counts only), 'ols' (log-CPM regression), or 'auto'
    (deseq2 for counts when importable, otherwise ols). Returns (DataFrame gene, base_mean, log2fc, stat, pvalue, padj, method)."""
    import numpy as np
    import pandas as pd
    paired = des["paired"]
    if kind == "counts":
        tot = mat.sum(axis=0).to_numpy()
        nz = (mat.to_numpy() > 0).sum(axis=0)
        keep = (tot >= 10) & (nz >= des["min_units"])
    else:
        keep = (mat.to_numpy().max(axis=0) > 0)
    m = mat.loc[:, keep]
    used = None
    if kind == "counts" and method in ("auto", "deseq2"):
        try:
            import warnings
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
            from pydeseq2.default_inference import DefaultInference
            md = pd.DataFrame({"sample": ["s" + str(s) for s in units["sample"]], "condition": units["condition"].astype(str).to_numpy()}, index=units["unit"])
            inf = DefaultInference(n_cpus=1)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                dds = DeseqDataSet(counts=m.round().astype(int), metadata=md, design="~sample + condition" if paired else "~condition", inference=inf, quiet=True)
                dds.deseq2()
                st = DeseqStats(dds, contrast=["condition", des["test"], des["ref"]], inference=inf, quiet=True)
                st.summary()
            r = st.results_df
            out = pd.DataFrame({"gene": r.index, "base_mean": r["baseMean"].to_numpy(), "log2fc": r["log2FoldChange"].to_numpy(), "stat": r["stat"].to_numpy(),
                                "pvalue": r["pvalue"].to_numpy(), "padj": r["padj"].to_numpy()})
            used = "pydeseq2 (negative binomial, sample-level counts)"
        except Exception as e:
            if method == "deseq2":
                raise
            used = None
    if used is None:
        Y = cmp_logcpm(m) if kind == "counts" else m.to_numpy(float) / np.log(2)
        X = cmp_matrix(units, paired)
        b, t, p = cmp_ols(Y, X)
        out = pd.DataFrame({"gene": m.columns, "base_mean": m.mean(axis=0).to_numpy(), "log2fc": b, "stat": t, "pvalue": p})
        out["padj"] = cmp_bh(out["pvalue"].to_numpy())
        used = "linear model on log-CPM per sample" if kind == "counts" else "linear model on mean log expression per sample (no counts available)"
    out["padj"] = out["padj"].astype(float)
    out["method"] = used
    out["significant"] = (out["padj"] < fdr) & (out["log2fc"].abs() >= min_lfc)
    return out.sort_values("pvalue").reset_index(drop=True)

def cmp_consistency(mat, units, des, kind, genes, lfc):
    """For significant genes: share of individual samples (pairs) whose treated-minus-control difference has the same sign as the overall effect; median over genes."""
    import numpy as np
    if not des["paired"] or len(genes) == 0:
        return float("nan")
    Y = cmp_logcpm(mat) if kind == "counts" else mat.to_numpy(float)
    cols = [mat.columns.get_loc(g) for g in genes]
    diffs = []
    for s, grp in units.groupby("sample"):
        a = grp.index[grp["treated"] == 1]
        b = grp.index[grp["treated"] == 0]
        if len(a) == 1 and len(b) == 1:
            diffs.append(Y[a[0], cols] - Y[b[0], cols])
    if len(diffs) < 2:
        return float("nan")
    D = np.array(diffs)
    same = (np.sign(D) == np.sign(np.asarray(lfc))[None, :]).mean(axis=0)
    return float(np.median(same))

def cmp_null_swaps(mat, units, des, kind, genes_mask=None, n_swaps=20, fdr=0.05, seed=0):
    """Swap the condition labels (within samples when paired, among samples otherwise) and count discoveries of the log-CPM regression; returns (real count, list of null counts)."""
    import numpy as np
    rng = np.random.default_rng(seed)
    Y = cmp_logcpm(mat) if kind == "counts" else mat.to_numpy(float)
    keep = Y.std(axis=0) > 0
    Y = Y[:, keep]
    paired = des["paired"]
    u = units.copy()

    def count(tr):
        uu = u.copy()
        uu["treated"] = tr
        _, _, p = cmp_ols(Y, cmp_matrix(uu, paired))
        return int((cmp_bh(p) < fdr).sum())

    real = count(u["treated"].to_numpy())
    nulls = []
    samples = list(u["sample"].unique())
    for _ in range(n_swaps):
        if paired:
            flip = {s: rng.random() < 0.5 for s in samples}
            tr = np.array([(1 - t) if flip[s] else t for s, t in zip(u["sample"], u["treated"])])
            if (tr == u["treated"].to_numpy()).all() or (tr != u["treated"].to_numpy()).all():
                continue
        else:
            if len(u) < 6:
                return real, []
            tr = rng.permutation(u["treated"].to_numpy())
            if (tr == u["treated"].to_numpy()).all() or (tr == 1 - u["treated"].to_numpy()).all():
                continue
        nulls.append(count(tr))
    return real, nulls

def cmp_state(adata, des, ctype, emb_key="X_pca", max_cells=2000, n_pcs=30, seed=0, n_null=3):
    """Can cells of one type be told apart by condition in held-out samples? Logistic regression on principal components; the held-out unit is one sample (paired) or one control and one
    treated sample (unpaired). Returns dict: auc per held-out unit, median, minimum, and the same with condition labels shuffled within samples."""
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    obs = adata.obs
    cond = obs[des["condition_col"]].astype(str).to_numpy()
    samp = obs[des["sample_col"]].astype(str).to_numpy()
    m = (obs[des["type_col"]].astype(str).to_numpy() == ctype) & np.isin(cond, [des["ref"], des["test"]])
    idx = np.where(m)[0]
    if len(idx) > max_cells:
        idx = np.sort(rng.choice(idx, max_cells, replace=False))
    E = np.asarray(adata.obsm[emb_key])[idx][:, :n_pcs]
    y = (cond[idx] == des["test"]).astype(int)
    s = samp[idx]
    ok_s = [x for x in np.unique(s) if (s == x).sum() >= des["min_cells"]]
    folds = []
    if des["paired"]:
        for x in ok_s:
            if len(set(y[s == x])) == 2:
                folds.append(np.isin(s, [x]))
    else:
        r_s = [x for x in ok_s if set(y[s == x]) == {0}]
        t_s = [x for x in ok_s if set(y[s == x]) == {1}]
        for _ in range(min(10, len(r_s) * len(t_s))):
            folds.append(np.isin(s, [r_s[rng.integers(len(r_s))], t_s[rng.integers(len(t_s))]]))
    def run(yy):
        out = []
        for te in folds:
            tr = ~te
            if len(set(yy[tr])) < 2 or len(set(yy[te])) < 2:
                continue
            mu, sd = E[tr].mean(0), E[tr].std(0) + 1e-9
            clf = LogisticRegression(max_iter=300, class_weight="balanced", C=0.5).fit((E[tr] - mu) / sd, yy[tr])
            out.append(float(roc_auc_score(yy[te], clf.decision_function((E[te] - mu) / sd))))
        return out
    auc = run(y)
    nulls = []
    for _ in range(n_null):
        yy = y.copy()
        for x in np.unique(s):
            w = np.where(s == x)[0]
            yy[w] = rng.permutation(y[w]) if des["paired"] else yy[w]
        if not des["paired"]:
            yy = rng.permutation(y)
        nulls += run(yy)
    res = {"type": ctype, "n_cells": int(len(idx)), "auc": auc, "median_auc": float(np.median(auc)) if auc else float("nan"), "min_auc": float(np.min(auc)) if auc else float("nan"),
           "null_median_auc": float(np.median(nulls)) if nulls else float("nan"), "n_folds": len(auc)}
    med, mn = res["median_auc"], res["min_auc"]
    res["label"] = "not testable" if not auc else ("clear shift" if (med >= 0.8 and mn >= 0.65) else ("weak shift" if med >= 0.65 else "no clear shift"))
    return res

def cmp_props(adata, des):
    """Long table of the share of each named cell type per sample and condition (all cells of the sample as denominator)."""
    import pandas as pd
    o = adata.obs[[des["sample_col"], des["condition_col"], des["type_col"]]].astype(str)
    o = o[o[des["condition_col"]].isin([des["ref"], des["test"]])]
    n = o.groupby([des["sample_col"], des["condition_col"]]).size().rename("n")
    k = o.groupby([des["sample_col"], des["condition_col"], des["type_col"]]).size().rename("k").reset_index()
    k = k.merge(n.reset_index(), on=[des["sample_col"], des["condition_col"]])
    k = k[k["n"] >= des["min_cells"]]
    k["share"] = k["k"] / k["n"]
    return k.rename(columns={des["sample_col"]: "sample", des["condition_col"]: "condition", des["type_col"]: "type"})

def cmp_enrich(de, sets, background, n_perm=1500, fdr=0.05, seed=0):
    """Pre-ranked enrichment of one cell type's differential expression (statistic = Wald or t statistic, positive = higher in treated)."""
    import pandas as pd
    try:
        import scrna_interpret_kernel as K7
        gsea, filt = K7.int_gsea, K7.int_filter_sets
    except ImportError:
        g = globals()
        gsea, filt = g["int_gsea"], g["int_filter_sets"]
    r = de.dropna(subset=["stat"]).drop_duplicates("gene").set_index("gene")["stat"]
    fs = filt(sets, background & set(r.index))
    if len(fs) < 5 or len(r) < 200:
        return pd.DataFrame(columns=["term", "size", "es", "nes", "p", "fdr"])
    out = gsea(r, fs, n_perm=n_perm, seed=seed)
    if not len(out):
        return out
    out = out[out["fdr"] < fdr]
    return out.assign(_a=out["nes"].abs()).sort_values(["fdr", "_a"], ascending=[True, False]).drop(columns="_a").reset_index(drop=True)      # ties at the smallest p are ordered by effect size

def cmp_summaries(des, ab, de, state, enr, fdr=0.05):
    """Plain sentences answering each question from the computed tables."""
    import numpy as np
    A = {}
    if ab is not None and len(ab):
        both = ab[ab["agree"]]
        one = ab[((ab["fdr"] < fdr) | (ab["fdr_clr"] < fdr)) & ~ab["agree"]]
        def fmt(r):
            return f"{r.type} ({'more' if r.log2_ratio > 0 else 'fewer'}: {r.share_ref:.1%} to {r.share_test:.1%})" + (" [name follows the condition]" if r.caution else "")
        A["abundance"] = (f"{len(both)} of {len(ab)} cell types change in share between {des['ref']} and {des['test']} in both models" + (": " + "; ".join(fmt(r) for r in both.itertuples()) if len(both) else "") +
                          (f". {len(one)} more change in one model only: " + "; ".join(one["type"]) if len(one) else "") + ".")
    if de:
        tested = {t: d for t, d in de.items()}
        sig = {t: d[d["significant"]] for t, d in tested.items()}
        with_any = [t for t, d in sig.items() if len(d)]
        parts = []
        for t in with_any:
            d = sig[t]
            up, dn = d[d["log2fc"] > 0], d[d["log2fc"] < 0]
            parts.append(f"{t}: {len(up)} up, {len(dn)} down (top up {', '.join(up.sort_values('padj')['gene'].head(3)) or 'none'})")
        A["de"] = (f"{len(with_any)} of {len(tested)} compared cell types have genes that differ between {des['ref']} and {des['test']} (adjusted p below {fdr} and at least a 1.4-fold change)" + (": " + "; ".join(parts) if parts else "") + ".")
    if state:
        lab = {}
        for t, r in state.items():
            lab.setdefault(r["label"], []).append(t)
        A["state"] = "; ".join(f"{k}: {', '.join(v)}" for k, v in lab.items()) + "."
    if enr:
        parts = [f"{t}: {', '.join(e['term'].head(2))}" for t, e in enr.items() if len(e)]
        A["pathways"] = ("Pathways moving with the treatment (pre-ranked, adjusted p below 0.05): " + "; ".join(parts) + ".") if parts else "No pathway passes the threshold in any compared cell type."
    return A

def cmp_run(adata, condition_col, sample_col, type_col="cell_type", batch_col=None, groups=None, gene_sets=None, counts_layer="counts", min_units=3, min_cells=10,
            fdr=0.05, min_lfc=0.5, method="auto", state=True, n_swaps=20, fig_dir=None, seed=0):
    """Control vs treatment. Order: design gate; abundance; per-type pseudobulk differential expression with a label-swap null and a per-sample consistency check;
    within-type state change on held-out samples; pre-ranked pathway enrichment. Returns (adata, records, details); writes uns['comparison'] and uns['comparison_de']."""
    import numpy as np
    import pandas as pd
    recs, des = cmp_design(adata, condition_col, sample_col, type_col, batch_col, groups=groups, min_units=min_units, min_cells=min_cells)
    det = {"design": des}
    if not des or any(r["status"] == "FAIL" for r in recs):
        return adata, recs, det
    ab = cmp_abundance(adata, des)
    det["abundance"] = ab
    if len(ab):
        both = ab[ab["agree"]]
        flag = [t for t in both["type"] if t in des["followers"]]
        recs.append(cmp_record("abundance", "WARN" if flag else "PASS",
                               f"{len(both)} of {len(ab)} cell types change in share in both models (quasi-binomial and log-ratio, FDR below {fdr}); {int(((ab['fdr'] < fdr) ^ (ab['fdr_clr'] < fdr)).sum())} in one model only" +
                               (f"; {', '.join(flag)} follow the condition, so a change in share may be a change of name" if flag else "")))
    de, mats = {}, {}
    ct = adata.obs[type_col].astype(str).to_numpy()
    kinds = set()
    for t in des["testable"]:
        mat, units, kind = cmp_pseudobulk(adata, ct == t, des, counts_layer)
        if units["treated"].nunique() < 2 or units.groupby("treated").size().min() < 2:
            continue
        det.setdefault("units", {})[t] = (int((units["treated"] == 0).sum()), int((units["treated"] == 1).sum()))
        d = cmp_de_type(mat, units, des, kind, method=method, min_lfc=min_lfc, fdr=fdr)
        d["type"] = t
        de[t], mats[t] = d, (mat, units, kind)
        kinds.add(kind)
    det["de"] = de
    if de:
        n_sig = {t: int(d["significant"].sum()) for t, d in de.items()}
        meth = sorted(set(m for d in de.values() for m in d["method"].unique()))
        fallback = ("counts" in kinds) and any(not m.startswith("pydeseq2") for m in meth) and method == "auto"
        recs.append(cmp_record("differential expression", "WARN" if fallback else "PASS",
                               f"{sum(v > 0 for v in n_sig.values())} of {len(de)} compared cell types have significant genes ({', '.join(f'{t} {n}' for t, n in n_sig.items())}); "
                               f"method: {'; '.join(meth)}" + ("; the negative-binomial fit failed for at least one type and the simpler log-CPM model was used" if fallback else "")))
        real, nulls = 0, []
        for t, (mat, units, kind) in mats.items():
            r_, n_ = cmp_null_swaps(mat, units, des, kind, n_swaps=n_swaps, fdr=fdr, seed=seed)
            real += r_
            nulls.append(n_)
        flat = [np.sum(x) for x in zip(*[n for n in nulls if n])] if any(nulls) else []
        if flat:
            med = float(np.median(flat))
            recs.append(cmp_record("null check", "PASS" if med <= max(5, 0.05 * real) else "WARN",
                                   f"with condition labels swapped {'within samples' if des['paired'] else 'among samples'} the log-CPM model finds a median of {med:.0f} genes across all compared types ({len(flat)} swaps) against {real} with the real labels"))
        cons = {}
        for t, (mat, units, kind) in mats.items():
            sg = de[t][de[t]["significant"]]
            sg = sg[sg["gene"].isin(mat.columns)]
            c = cmp_consistency(mat, units, des, kind, list(sg["gene"]), sg["log2fc"].to_numpy())
            if c == c:
                cons[t] = (c, len(sg))
        if cons:
            tot = sum(n for _, n in cons.values())
            wmed = float(np.median([c for c, _ in cons.values()]))
            recs.append(cmp_record("sample consistency", "PASS" if wmed >= 0.75 else "WARN",
                                   f"significant genes keep their direction in a median of {wmed:.0%} of individual samples ({tot} genes in {len(cons)} cell types; 75% or more expected)"))
        det["consistency"] = cons
        allde = pd.concat([d.assign(significant=d["significant"].astype(int)) for d in de.values()], ignore_index=True)
        adata.uns["comparison_de"] = allde[["type", "gene", "base_mean", "log2fc", "stat", "pvalue", "padj", "significant", "method"]].astype({"type": str, "gene": str, "method": str})
    st = {}
    if state and "X_pca" in adata.obsm:
        for t in des["testable"]:
            r = cmp_state(adata, des, t, seed=seed)
            st[t] = r
        det["state"] = st
        usable = {t: r for t, r in st.items() if r["n_folds"]}
        if usable:
            bad = [t for t, r in usable.items() if r["null_median_auc"] > 0.65]
            lab = {}
            for t, r in usable.items():
                lab.setdefault(r["label"], []).append(t)
            recs.append(cmp_record("state change", "WARN" if bad else "PASS",
                                   "; ".join(f"{k}: {', '.join(v)}" for k, v in lab.items()) + f" (held-out AUC: clear shift median 0.8+, weak 0.65+)" +
                                   (f"; shuffled-label AUC above 0.65 in {', '.join(bad)}, so the classifier finds structure without the true labels" if bad else "")))
    enr = {}
    if gene_sets and de:
        bg = set(str(g).upper() for g in adata.var_names)
        for t, d in de.items():
            e = cmp_enrich(d, gene_sets, bg, fdr=fdr, seed=seed)
            enr[t] = e
        det["enrichment"] = enr
        recs.append(cmp_record("pathway enrichment", "PASS", f"{sum(len(e) > 0 for e in enr.values())} of {len(enr)} cell types have pathways below FDR {fdr}; " +
                               "; ".join(f"{t} {len(e)}" for t, e in enr.items())))
    det["answers"] = cmp_summaries(des, ab, de, st, enr, fdr)
    adata.uns["comparison"] = {"condition_col": condition_col, "sample_col": sample_col, "reference": des["ref"], "treatment": des["test"], "paired": bool(des["paired"]),
                               "types_compared": list(de), "min_lfc": float(min_lfc), "fdr": float(fdr)}
    det["props"] = cmp_props(adata, des)
    if fig_dir:
        cmp_figures(fig_dir, det)
    return adata, recs, det

def cmp_captions(det):
    des = det["design"]
    return {
        "cmp_abundance.png": f"Share of each cell type in each sample, {des['ref']} to {des['test']}, one line per sample (point = one sample); panel titles give the adjusted p-value of the count model.",
        "cmp_de.png": f"Left: number of genes higher (right) and lower (left) in {des['test']} than {des['ref']} per cell type at the stated thresholds. Right: volcano plot of the cell type with the most such genes (x = log2 fold change, y = -log10 adjusted p).",
        "cmp_state.png": f"Area under the ROC curve for telling {des['ref']} from {des['test']} cells of the same type in held-out samples (one dot per held-out sample); crosses are the same test with shuffled labels.",
        "cmp_pathways.png": "Normalised enrichment score of the strongest pathways (columns) in each cell type's ranked differential expression (rows); red = higher with treatment, blue = lower.",
    }

def cmp_legends(det, fdr=0.05, min_lfc=0.5):
    des = det["design"]
    L = {
        "cmp_abundance": (f"Each panel is one cell type. Each line joins the share of all cells that belong to the type in one sample under {des['ref']} and under {des['test']}; lines going the same way in most samples mean a consistent change. "
                          "Shares add up to one within a sample, so when one type grows the others shrink; the counts are tested with a quasi-binomial model and the log-ratio of shares with a second model, and a type is called changed only when both agree. "
                          "Types whose name follows the condition are marked: their share may change because the label changed, not the cells."),
        "cmp_de": (f"Genes are compared within each cell type after adding up the counts of its cells in each sample, so that the number of independent observations is the number of samples and not the number of cells (single-cell tests treat thousands of cells as replicates and find far too many genes). "
                   f"A gene counts when its adjusted p-value is below {fdr} and its change is at least {2 ** min_lfc:.1f}-fold."),
        "cmp_state": ("A classifier is trained on cells from all samples but one and asked to tell the conditions apart in the held-out sample, so a high value means the cells of that type look different under treatment in a sample the model never saw; 0.5 is chance. "
                      + ("Control and treated cells come from different libraries here, so a shift can also be a library effect." if des.get("batch_confounded") else "")),
        "cmp_pathways": ("Genes of each cell type are ranked by their test statistic and the ranking is compared with each gene set; the normalised enrichment score is positive when the set's genes sit near the top (higher with treatment). "
                         "Only pathways below the stated false-discovery rate are shown."),
    }
    return L

def cmp_figures(fig_dir, det, fdr=0.05, top_label=8):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(fig_dir, exist_ok=True)
    des, ab, de, st = det["design"], det.get("abundance"), det.get("de") or {}, det.get("state") or {}
    c_ref, c_test = "#4C72B0", "#C44E52"
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    if ab is not None and len(ab):
        props = det["props"]
        types = list(ab["type"])
        ncol = min(4, len(types))
        nrow = int(np.ceil(len(types) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol, 2.4 * nrow), squeeze=False)
        for ax, r in zip(axes.ravel(), ab.itertuples()):
            sub = props[props["type"] == r.type]
            piv = sub.pivot(index="sample", columns="condition", values="share").reindex(columns=[des["ref"], des["test"]]).fillna(0) * 100
            for s_, row in piv.iterrows():
                ax.plot([0, 1], row.to_numpy(), color="#888888", lw=0.8, marker="o", ms=3)
            ax.set_xticks([0, 1])
            ax.set_xticklabels([des["ref"], des["test"]])
            ax.set_title(f"{r.type}\nFDR {r.fdr:.2g}" + (" (name follows condition)" if r.caution else ""), fontsize=8)
            ax.set_ylabel("% of cells")
        for ax in axes.ravel()[len(types):]:
            ax.axis("off")
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "cmp_abundance.png"), dpi=170)
        plt.close(fig)
    if de:
        types = list(de)
        up = [int(((d["log2fc"] > 0) & d["significant"]).sum()) for d in de.values()]
        dn = [int(((d["log2fc"] < 0) & d["significant"]).sum()) for d in de.values()]
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 0.45 * len(types) + 2.5), gridspec_kw={"width_ratios": [1, 1.1]})
        y = np.arange(len(types))
        a1.barh(y, up, color=c_test)
        a1.barh(y, [-v for v in dn], color=c_ref)
        a1.set_yticks(y)
        a1.set_yticklabels(types)
        a1.axvline(0, color="k", lw=0.6)
        a1.set_xlabel(f"genes lower (left) / higher (right) in {des['test']}")
        best = types[int(np.argmax(np.array(up) + np.array(dn)))]
        d = de[best].dropna(subset=["pvalue"])
        x, yv = d["log2fc"].to_numpy(), -np.log10(np.maximum(d["pvalue"].to_numpy(), 1e-300))
        a2.scatter(x, yv, s=4, c=np.where(d["significant"], np.where(d["log2fc"] > 0, c_test, c_ref), "#BBBBBB"), linewidths=0)
        for r in d[d["significant"]].sort_values("padj").head(top_label).itertuples():
            a2.annotate(r.gene, (r.log2fc, -np.log10(max(r.pvalue, 1e-300))), fontsize=7)
        a2.set_xlabel("log2 fold change")
        a2.set_ylabel("-log10 p")
        a2.set_title(best, fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "cmp_de.png"), dpi=170)
        plt.close(fig)
    use = {t: r for t, r in st.items() if r["n_folds"]}
    if use:
        fig, ax = plt.subplots(figsize=(7, 0.45 * len(use) + 1.8))
        for i, (t, r) in enumerate(use.items()):
            ax.scatter(r["auc"], np.full(len(r["auc"]), i), color="#555555", s=18, zorder=3)
            ax.scatter([r["null_median_auc"]], [i], marker="x", color="#C44E52", s=30, zorder=3)
        ax.axvline(0.5, color="#999999", lw=0.8, ls=":")
        ax.axvline(0.8, color="#999999", lw=0.8, ls="--")
        ax.set_yticks(range(len(use)))
        ax.set_yticklabels([f"{t} ({r['label']})" for t, r in use.items()])
        ax.set_xlim(0.3, 1.02)
        ax.set_xlabel("held-out sample AUC (cross = shuffled labels)")
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "cmp_state.png"), dpi=170)
        plt.close(fig)
    enr = {t: e for t, e in (det.get("enrichment") or {}).items() if len(e)}
    if enr:
        terms = []
        for t, e in enr.items():
            terms += [x for x in e["term"].head(3) if x not in terms]
        M = np.full((len(enr), len(terms)), np.nan)
        for i, (t, e) in enumerate(enr.items()):
            for r in e.itertuples():
                if r.term in terms:
                    M[i, terms.index(r.term)] = r.nes
        fig, ax = plt.subplots(figsize=(max(5, 0.45 * len(terms) + 3), 0.45 * len(enr) + 3.5))
        lim = np.nanmax(np.abs(M)) if np.isfinite(M).any() else 1
        im = ax.imshow(M, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
        ax.set_yticks(range(len(enr)))
        ax.set_yticklabels(list(enr))
        ax.set_xticks(range(len(terms)))
        ax.set_xticklabels([x[:48] for x in terms], rotation=60, ha="right", fontsize=7)
        fig.colorbar(im, ax=ax, shrink=0.6, label="normalised enrichment score")
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "cmp_pathways.png"), dpi=170)
        plt.close(fig)
    caps = cmp_captions(det)
    with open(os.path.join(fig_dir, "figure_captions.json"), "w", encoding="utf8") as f:
        json.dump({k: v for k, v in caps.items() if os.path.exists(os.path.join(fig_dir, k))}, f, indent=1)
    leg = cmp_legends(det)
    with open(os.path.join(fig_dir, "figure_legends.md"), "w", encoding="utf8") as f:
        for k, v in leg.items():
            if os.path.exists(os.path.join(fig_dir, k + ".png")):
                f.write(f"**{k}** ({k}.png): {v}\n\n")
    if ab is not None:
        ab.to_csv(os.path.join(fig_dir, "abundance.csv"), index=False)
    if de:
        import pandas as pd
        pd.concat(de.values(), ignore_index=True).to_csv(os.path.join(fig_dir, "differential_expression.csv"), index=False)
    if st:
        import pandas as pd
        pd.DataFrame([{k: v for k, v in r.items() if k != "auc"} | {"auc": ";".join(f"{a:.3f}" for a in r["auc"])} for r in st.values()]).to_csv(os.path.join(fig_dir, "state_change.csv"), index=False)
    if enr:
        import pandas as pd
        pd.concat([e.assign(type=t) for t, e in enr.items()], ignore_index=True).to_csv(os.path.join(fig_dir, "pathway_enrichment.csv"), index=False)

def write_cmp_report(records, det, caveats=None, prefix="comparison_report", question=None, scope=None):
    """Write <prefix>.json and <prefix>.md and return the overall status. Question first, then the answer to each sub-question with its evidence, the design gate, the checks with what
    they test, the checks that did not run, and the caveats carried from earlier checkpoints."""
    cat = cmp_catalog()
    rank = {"PASS": 0, "INFO": 0, "WARN": 1, "FAIL": 2}
    overall = ["PASS", "WARN", "FAIL"][max((rank[r["status"]] for r in records), default=0)]
    count = {s: sum(r["status"] == s for r in records) for s in ("PASS", "WARN", "FAIL", "INFO")}
    seen = [r.get("group") or r["check"] for r in records]
    not_run = {k: v for k, v in cat.items() if k not in seen}
    short = lambda r: r["check"] + " (" + r["detail"][:110].replace("|", "/") + ")"
    fails, warns = [short(r) for r in records if r["status"] == "FAIL"], [short(r) for r in records if r["status"] == "WARN"]
    if fails:
        why = f"{len(fails)} of {len(records)} checks failed: " + "; ".join(fails)
    elif warns:
        why = f"{len(warns)} warning(s), no failures: " + "; ".join(warns)
    else:
        why = f"all {count['PASS']} scored checks passed ({count['INFO']} informational)"
    if not_run:
        why += f". {len(not_run)} catalogued check(s) did not run."
    des = det.get("design") or {}
    A = det.get("answers", {})
    with open(prefix + ".json", "w") as f:
        json.dump({"overall": overall, "why": why, "counts": count, "not_run": list(not_run), "checks": records, "answers": A}, f, indent=1, default=str)
    L = [f"# Condition comparison report: {overall}", "", f"**Why:** {why}", ""]
    if des:
        q = question or f"What changes in these cells when {des['test']} is compared with {des['ref']}?"
        L += [f"**Question:** {q}", ""]
        L += ["**Answer in short** (read the design checks first; they decide how far the answer can be trusted):", ""]
        for k, h in (("abundance", "Which cell types change in abundance?"), ("de", "Which genes change within each cell type?"), ("state", "Do cells shift state within a type?"), ("pathways", "Which pathways move?")):
            if k in A:
                L.append(f"- *{h}* {A[k]}")
        L.append("")
        L += [f"**Checks run ({len(records)}):** {count['PASS']} pass, {count['WARN']} warn, {count['FAIL']} fail, {count['INFO']} info", ""]
        L += ["## 1. Is the design able to answer the question?", "", "| status | check | what it tests | result |", "|---|---|---|---|"]
        for r in records:
            if r["check"] in ("groups", "replicates", "pairing", "condition vs batch", "clusters follow condition", "cell types follow condition", "testable cell types", "unlabelled cells"):
                L.append(f"| {r['status']} | {r['check']} | {cat.get(r['check'], '')} | {r['detail'].replace('|', '/')} |")
        L.append("")
        ab = det.get("abundance")
        if ab is not None and len(ab):
            L += ["## 2. Which cell types change in abundance?", "", f"**Result:** {A.get('abundance', '')}", "",
                  "| cell type | share " + des["ref"] + " | share " + des["test"] + " | log2 ratio | FDR (counts) | FDR (log-ratio) | both agree | caution |", "|---|---|---|---|---|---|---|---|"]
            for r in ab.itertuples():
                L.append(f"| {r.type} | {r.share_ref:.1%} | {r.share_test:.1%} | {r.log2_ratio:+.2f} | {r.fdr:.2g} | {r.fdr_clr:.2g} | {'yes' if r.agree else 'no'} | {'name follows condition' if r.caution else ''} |")
            L.append("")
        de = det.get("de")
        if de:
            L += ["## 3. Which genes change within each cell type?", "", f"**Result:** {A.get('de', '')}", "",
                  "| cell type | samples (control/treated) | genes tested | higher in treated | lower in treated | strongest higher | strongest lower | method |", "|---|---|---|---|---|---|---|---|"]
            for t, d in de.items():
                sg = d[d["significant"]]
                up, dn = sg[sg["log2fc"] > 0].sort_values("padj"), sg[sg["log2fc"] < 0].sort_values("padj")
                L.append(f"| {t} | {'/'.join(str(x) for x in det.get('units', {}).get(t, ('', '')))} | {len(d):,} | {len(up)} | {len(dn)} | {', '.join(up['gene'].head(4))} | {', '.join(dn['gene'].head(4))} | {d['method'].iloc[0]} |")
            L.append("")
        st = {t: r for t, r in (det.get("state") or {}).items() if r["n_folds"]}
        if st:
            L += ["## 4. Do cells shift state within a type?", "", f"**Result:** {A.get('state', '')}", "", "| cell type | cells used | held-out samples | median AUC | lowest AUC | shuffled-label AUC | label |", "|---|---|---|---|---|---|---|"]
            for t, r in st.items():
                L.append(f"| {t} | {r['n_cells']:,} | {r['n_folds']} | {r['median_auc']:.2f} | {r['min_auc']:.2f} | {r['null_median_auc']:.2f} | {r['label']} |")
            L.append("")
        enr = det.get("enrichment")
        if enr:
            L += ["## 5. Which pathways move with the treatment?", "", f"**Result:** {A.get('pathways', '')}", "", "| cell type | pathway | normalised enrichment | FDR |", "|---|---|---|---|"]
            for t, e in enr.items():
                for r in e.head(3).itertuples():
                    L.append(f"| {t} | {r.term} | {r.nes:+.2f} | {r.fdr:.2g} |")
            L.append("")
    L += ["## Checks", "", "| status | check | what it tests | result |", "|---|---|---|---|"]
    for r in records:
        L.append(f"| {r['status']} | {r['check']} | {cat.get(r.get('group') or r['check'], 'additional check')} | {r['detail'].replace('|', '/')} |")
    if not_run:
        L += ["", "**Not run:**", ""] + [f"- {k}: {v}" for k, v in not_run.items()]
    if caveats:
        L += ["", "## Caveats carried from earlier checkpoints", "", "These limit every conclusion above.", ""]
        for kind in ("FAIL", "WARN", "not run"):
            for c in [c for c in caveats if c["kind"] == kind]:
                L.append(f"- [{kind}] {c['stage']}: {c['text']}")
    with open(prefix + ".md", "w", encoding="utf8") as f:
        f.write("\n".join(L) + "\n")
    return overall
