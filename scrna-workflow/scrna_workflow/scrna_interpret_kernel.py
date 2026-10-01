"""Checkpoint 8 helpers: biological interpretation, validation, robustness and sensitivity (functions and imports only)."""
import json

LABELS = ("robust", "sensitive", "unstable")

def int_record(check, status, detail, group=None):
    r = {"check": check, "status": status, "detail": detail}
    if group:
        r["group"] = group
    return r

def int_catalog():
    """Every check this checkpoint can run, with what it tests. Checks that did not run are listed as 'Not run'."""
    return {
        "gene sets": "at least one gene-set collection is loaded and usable against the genes detected in this dataset",
        "background": "the genes tested for enrichment are compared with the genes actually detected, not with the whole genome",
        "pathway analysis": "cell-type marker genes are tested for over-representation of gene sets with multiple-testing correction",
        "pre-ranked enrichment": "a ranked differential-expression result (checkpoint 7) is tested against the gene sets",
        "gene programs": "data-driven expression programs are reproducible across random restarts",
        "expectations": "known biology listed as positive and negative controls is reproduced by the data",
        "signature validation": "a gene signature separates two groups in samples it was not trained on",
        "robustness: clusters": "cluster membership is stable when seeds, cells, variable genes, PCs and resolution change",
        "robustness: cell-type calls": "each cell-type name survives the same changes",
        "robustness: pathways": "the enriched pathways are found again after the same changes",
        "caveats carried": "warnings and unrun checks from earlier checkpoints are carried into the conclusions",
    }

def int_split(x):
    if x is None or (isinstance(x, float) and x != x):
        return []
    if isinstance(x, (list, tuple, set)):
        return [str(g).strip() for g in x if str(g).strip()]
    return [g.strip() for g in str(x).replace(",", ";").split(";") if g.strip()]

def int_fetch_kegg_gmt(path, species="hsa"):
    """Download KEGG pathway membership (gene symbols) through the KEGG REST service and save it as a GMT file.
    Records the source and date in a header line starting with '#'. KEGG is free for academic use; check its terms before other uses."""
    import time
    import urllib.request

    def get(p):
        with urllib.request.urlopen("https://rest.kegg.jp/" + p, timeout=90) as r:
            return r.read().decode("utf8")

    names = dict(l.split("\t") for l in get(f"list/pathway/{species}").strip().splitlines())
    sym = {}
    for l in get(f"list/{species}").strip().splitlines():
        p = l.split("\t")
        if len(p) >= 4:
            sym[p[0]] = p[3].split(";")[0].split(",")[0].strip()
    sets = {}
    for l in get(f"link/{species}/pathway").strip().splitlines():
        pw, g = l.split("\t")
        s = sym.get(g)
        if s and " " not in s:
            sets.setdefault(pw.replace("path:", ""), set()).add(s.upper())
    with open(path, "w", encoding="utf8") as f:
        f.write(f"# KEGG {species} pathways via rest.kegg.jp, downloaded {time.strftime('%Y-%m-%d', time.gmtime())}\n")
        for k, v in sets.items():
            f.write(names.get(k, k).split(" - ")[0].replace("\t", " ") + " [" + k + "]\t" + k + "\t" + "\t".join(sorted(v)) + "\n")
    return path

def int_read_gmt(path):
    """GMT file -> {set name: set of upper-case gene symbols}. Lines starting with '#' are ignored."""
    out = {}
    for line in open(path, encoding="utf8"):
        if not line.strip() or line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            out[p[0]] = set(g.upper() for g in p[2:] if g)
    return out

def int_background(adata, min_cells=3):
    """Genes detected (above zero) in at least min_cells cells, upper-case."""
    import numpy as np
    import scipy.sparse as sp
    X = adata.X
    n = np.asarray((X > 0).sum(axis=0)).ravel() if sp.issparse(X) else (np.asarray(X) > 0).sum(axis=0)
    return set(str(g).upper() for g in adata.var_names[np.asarray(n) >= min_cells])

def int_filter_sets(sets, background, min_size=10, max_size=500):
    """Restrict every set to the background and keep sets with min_size to max_size genes left."""
    out = {}
    for k, v in sets.items():
        w = v & background
        if min_size <= len(w) <= max_size:
            out[k] = w
    return out

def int_bh(p):
    import numpy as np
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out

def int_ora(genes, sets, background, min_overlap=3):
    """Over-representation: hypergeometric test of `genes` against each set, with the detected genes as background (sets already restricted to it).
    Returns a DataFrame term, size, overlap, fold_enrichment, p, fdr, genes sorted by p."""
    import numpy as np
    import pandas as pd
    from scipy.stats import hypergeom
    g = set(x.upper() for x in genes) & background
    N, n = len(background), len(g)
    rows = []
    for t, s in sets.items():
        ov = g & s
        if len(ov) < min_overlap:
            continue
        p = float(hypergeom.sf(len(ov) - 1, N, len(s), n))
        rows.append({"term": t, "size": len(s), "overlap": len(ov), "n_input": n, "fold_enrichment": (len(ov) / n) / (len(s) / N), "p": p, "genes": ";".join(sorted(ov))})
    df = pd.DataFrame(rows, columns=["term", "size", "overlap", "n_input", "fold_enrichment", "p", "genes"])
    if len(df):
        df["fdr"] = int_bh(df["p"].to_numpy())
        df = df.sort_values("p").reset_index(drop=True)
    else:
        df["fdr"] = []
    return df

def int_trim(df, max_jaccard=0.5):
    """Drop a significant term whose overlapping genes largely repeat those of a better term already kept."""
    keep, seen = [], []
    for r in df.itertuples():
        s = set(r.genes.split(";"))
        if all(len(s & t) / max(len(s | t), 1) < max_jaccard for t in seen):
            keep.append(r.Index)
            seen.append(s)
    return df.loc[keep].reset_index(drop=True)

def int_group_genes(adata, labels, groups=None, top_n=100, min_pct=0.1, min_diff=0.25, min_cells=10):
    """Genes that distinguish each group from the rest, ranked by a t-like score on the normalised expression
    (mean difference over pooled standard error). A gene is kept if its mean is higher by min_diff and it is detected in min_pct of the group's cells.
    Returns {group: [genes in rank order]} and a DataFrame of scores. Fast and deterministic."""
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    X = sp.csr_matrix(adata.X)
    labels = np.asarray(labels).astype(str)
    groups = list(groups) if groups is not None else sorted(np.unique(labels))
    n = X.shape[0]
    mean_all = np.asarray(X.mean(axis=0)).ravel()
    sq = X.multiply(X)
    msq_all = np.asarray(sq.mean(axis=0)).ravel()
    out, rows = {}, []
    for g in groups:
        m = labels == g
        k = int(m.sum())
        if k < min_cells or k == n:
            continue
        Xi = X[m]
        mi = np.asarray(Xi.mean(axis=0)).ravel()
        det = np.asarray((Xi > 0).mean(axis=0)).ravel()
        mo = (mean_all * n - mi * k) / (n - k)
        vi = np.maximum(np.asarray(Xi.multiply(Xi).mean(axis=0)).ravel() - mi ** 2, 1e-9)
        vo = np.maximum((msq_all * n - np.asarray(Xi.multiply(Xi).sum(axis=0)).ravel()) / (n - k) - mo ** 2, 1e-9)
        t = (mi - mo) / np.sqrt(vi / k + vo / (n - k))
        ok = (mi - mo >= min_diff) & (det >= min_pct)
        idx = np.where(ok)[0]
        idx = idx[np.argsort(-t[idx])][:top_n]
        out[g] = [str(adata.var_names[i]).upper() for i in idx]
        rows += [{"group": g, "gene": str(adata.var_names[i]).upper(), "score": float(t[i]), "mean_diff": float(mi[i] - mo[i]), "detected": float(det[i])} for i in idx]
    return out, pd.DataFrame(rows, columns=["group", "gene", "score", "mean_diff", "detected"])

def int_gsea(ranked, sets, n_perm=1500, min_overlap=10, seed=0, chunk=250):
    """Pre-ranked gene-set enrichment. ranked: Series of scores indexed by upper-case gene (higher = more up). Enrichment score is the running-sum
    statistic weighted by |score|; the null comes from random gene sets of the same size (sizes grouped on a logarithmic grid), drawn n_perm times per size group.
    The p-value is the share of null scores at least as extreme in absolute value, so its smallest value is 1/(n_perm+1): with hundreds of sets, n_perm of 1000 or more
    is needed before any set can reach FDR 0.05. Returns term, size, es, nes, p, fdr."""
    import numpy as np
    import pandas as pd
    rng = np.random.default_rng(seed)
    r = ranked.sort_values(ascending=False)
    genes = np.array(r.index)
    w = np.abs(r.to_numpy()).astype(float)
    N = len(genes)
    pos = {g: i for i, g in enumerate(genes)}

    def es_of(idx):
        hit = np.zeros(N)
        hit[idx] = w[idx]
        run = np.cumsum(hit / max(hit.sum(), 1e-12) - (~np.isin(np.arange(N), idx)) / max(N - len(idx), 1))
        return run[np.argmax(np.abs(run))]

    def null_of(k):
        out = []
        done = 0
        while done < n_perm:
            c = min(chunk, n_perm - done)
            sel = np.argpartition(rng.random((c, N)), k - 1, axis=1)[:, :k]
            mask = np.zeros((c, N))
            np.put_along_axis(mask, sel, 1.0, axis=1)
            hit = mask * w[None, :]
            run = np.cumsum(hit / np.maximum(hit.sum(axis=1, keepdims=True), 1e-12) - (1.0 - mask) / max(N - k, 1), axis=1)
            out.append(run[np.arange(c), np.argmax(np.abs(run), axis=1)])
            done += c
        return np.concatenate(out)

    null_cache = {}
    rows = []
    for t, s in sets.items():
        idx = np.array([pos[g] for g in s if g in pos])
        if len(idx) < min_overlap:
            continue
        key = int(round(10 ** (round(np.log10(len(idx)) * 8) / 8)))
        key = min(max(key, 2), N - 1)
        if key not in null_cache:
            null_cache[key] = np.abs(null_of(key))
        nul = null_cache[key]
        es = es_of(idx)
        p = (np.sum(nul >= abs(es)) + 1) / (len(nul) + 1)
        rows.append({"term": t, "size": len(idx), "es": float(es), "nes": float(es / max(np.mean(nul), 1e-12)), "p": float(p)})
    df = pd.DataFrame(rows, columns=["term", "size", "es", "nes", "p"])
    if len(df):
        df["fdr"] = int_bh(df["p"].to_numpy())
        df = df.sort_values("p").reset_index(drop=True)
    else:
        df["fdr"] = []
    return df

def int_score_sets(adata, sets, n_control=50, n_bins=25, seed=0):
    """Per-cell score for each gene set: mean expression of the set minus the mean of expression-matched control genes
    (controls drawn from the same expression bin as each set gene). Returns a DataFrame cells x sets."""
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    rng = np.random.default_rng(seed)
    X = adata.X.tocsc() if sp.issparse(adata.X) else sp.csc_matrix(adata.X)
    up = pd.Series(range(adata.n_vars), index=[str(v).upper() for v in adata.var_names])
    up = up[~up.index.duplicated()]
    mean = np.asarray(X.mean(axis=0)).ravel()
    bins = pd.qcut(pd.Series(mean).rank(method="first"), n_bins, labels=False).to_numpy()
    out = {}
    for name, genes in sets.items():
        gi = np.array([up[g] for g in genes if g in up.index])
        if len(gi) < 2:
            continue
        ctrl = []
        for b in np.unique(bins[gi]):
            pool = np.where(bins == b)[0]
            pool = pool[~np.isin(pool, gi)]
            ctrl.append(rng.choice(pool, size=min(n_control * int((bins[gi] == b).sum()), len(pool)), replace=False))
        ctrl = np.concatenate(ctrl)
        s = np.asarray(X[:, gi].mean(axis=1)).ravel()
        c = np.asarray(X[:, ctrl].mean(axis=1)).ravel()
        out[name] = s - c
    return pd.DataFrame(out, index=adata.obs_names)

def int_nmf(adata, k_values=(5, 8, 12), n_restarts=3, max_cells=3000, min_stability=0.8, seed=0):
    """Data-driven gene programs by non-negative matrix factorisation of the variable genes. For every k, several restarts are compared
    (cosine similarity of best-matching programs); the largest k whose restarts agree (mean similarity at least min_stability) is used.
    Returns dict: k, table (k, stability, error), loadings (programs x genes DataFrame), usage (cells x programs DataFrame), qualified."""
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import warnings
    from sklearn.decomposition import NMF
    from sklearn.exceptions import ConvergenceWarning
    warnings.simplefilter('ignore', ConvergenceWarning)
    rng = np.random.default_rng(seed)
    hv = adata.var["highly_variable"].to_numpy() if "highly_variable" in adata.var else np.ones(adata.n_vars, bool)
    genes = np.array([str(g) for g in adata.var_names])[hv]
    X = adata.X[:, np.where(hv)[0]]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    X = np.clip(X, 0, None)
    idx = rng.choice(X.shape[0], size=min(max_cells, X.shape[0]), replace=False)
    Xs = X[idx]
    table, fits = [], {}
    for k in k_values:
        hs, err = [], []
        for r in range(n_restarts):
            m = NMF(n_components=k, init="random", random_state=seed + r, max_iter=250, tol=1e-3)
            m.fit(Xs)
            hs.append(m.components_ / np.maximum(np.linalg.norm(m.components_, axis=1, keepdims=True), 1e-12))
            err.append(m.reconstruction_err_)
            fits[(k, r)] = m
        sims = []
        for a in range(1, n_restarts):
            S = hs[0] @ hs[a].T
            sims.append(S.max(axis=1).mean())
        table.append({"k": k, "stability": float(np.mean(sims)) if sims else np.nan, "error": float(np.mean(err))})
    tab = pd.DataFrame(table)
    ok = tab[tab["stability"] >= min_stability]
    kbest = int(ok["k"].max()) if len(ok) else int(tab.sort_values("stability", ascending=False)["k"].iloc[0])
    m = fits[(kbest, 0)]
    U = m.transform(np.clip(X, 0, None))
    names = [f"program_{i + 1}" for i in range(kbest)]
    return {"k": kbest, "table": tab, "qualified": bool(len(ok)),
            "loadings": pd.DataFrame(m.components_, index=names, columns=genes),
            "usage": pd.DataFrame(U, index=adata.obs_names, columns=names)}

def int_expectations(adata, table, sets=None, ora_hits=None, min_diff=0.25, alpha=0.01):
    """Check known biology listed in a table. Columns: expectation, kind, column, group, genes or set, optional direction (up/down), against, min_fraction, low, high.
    kinds: gene_up (listed genes higher in group than in `against`, default the rest), set_up (set score higher), fraction_range (share of cells in group between low and high),
    set_enriched (the set is significant among the markers of group; needs ora_hits {group: set of terms}).
    Returns a list of dicts: expectation, status (met, not met, not testable), detail."""
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from scipy.stats import mannwhitneyu
    up = pd.Series(range(adata.n_vars), index=[str(v).upper() for v in adata.var_names])
    up = up[~up.index.duplicated()]
    out = []
    for r in table.to_dict("records"):
        name, kind = str(r.get("expectation", "")), str(r.get("kind", ""))
        col, grp = str(r.get("column", "")), str(r.get("group", ""))
        res = {"expectation": name, "kind": kind, "column": col, "group": grp}
        try:
            if col not in adata.obs:
                raise KeyError(f"column {col} not found")
            lab = adata.obs[col].astype(str).to_numpy()
            inn = lab == grp
            if not inn.any():
                raise KeyError(f"no cells with {col} = {grp}")
            against = str(r.get("against", "") or "")
            outm = (lab == against) if against and against != "nan" else ~inn
            direction = str(r.get("direction", "up") or "up")
            if kind == "fraction_range":
                f = float(inn.mean())
                lo, hi = float(r.get("low", 0)), float(r.get("high", 1))
                res.update(status="met" if lo <= f <= hi else "not met", detail=f"{100 * f:.1f}% of cells (expected {100 * lo:.0f} to {100 * hi:.0f}%)")
            elif kind in ("gene_up", "set_up"):
                if kind == "gene_up":
                    genes = [g.upper() for g in int_split(r.get("genes"))]
                    gi = [up[g] for g in genes if g in up.index]
                    missing = [g for g in genes if g not in up.index]
                    if len(gi) < 1:
                        raise KeyError("none of the listed genes are in the data")
                    M = adata.X[:, gi]
                    M = M.toarray() if sp.issparse(M) else np.asarray(M)
                    cols = [genes[i] for i in range(len(genes)) if genes[i] in up.index]
                else:
                    sname = str(r.get("set"))
                    if sets is None or sname not in sets:
                        raise KeyError(f"gene set {sname} not available")
                    M = int_score_sets(adata, {sname: sets[sname]}).to_numpy()
                    cols, missing = [sname], []
                sign = 1.0 if direction != "down" else -1.0
                ok, parts = 0, []
                for j, c in enumerate(cols):
                    a, b = M[inn, j], M[outm, j]
                    d = float(a.mean() - b.mean())
                    p = float(mannwhitneyu(a, b, alternative="greater" if sign > 0 else "less").pvalue) if len(a) > 2 and len(b) > 2 else 1.0
                    good = sign * d >= min_diff * (1 if kind == "gene_up" else 0.4) and p < alpha
                    ok += good
                    parts.append(f"{c} {'+' if d >= 0 else ''}{d:.2f}{'*' if good else ''}")
                need = float(r.get("min_fraction", 0.6) or 0.6)
                frac = ok / len(cols)
                res.update(status="met" if frac >= need else "not met",
                           detail=f"{ok} of {len(cols)} {'genes' if kind == 'gene_up' else 'set'} {'higher' if sign > 0 else 'lower'} in {grp} (mean difference; * = passes): " + ", ".join(parts[:8])
                           + (f"; not in data: {', '.join(missing)}" if missing else ""))
            elif kind == "set_enriched":
                sname = str(r.get("set"))
                if ora_hits is None or grp not in ora_hits:
                    raise KeyError("no pathway results for this group")
                hit = [t for t in ora_hits[grp] if sname.lower() in t.lower()]
                res.update(status="met" if hit else "not met", detail=("significant: " + "; ".join(hit[:3])) if hit else f"no significant term matching '{sname}'")
            else:
                raise KeyError(f"unknown kind '{kind}'")
        except Exception as e:
            res.update(status="not testable", detail=str(e)[:200])
        out.append(res)
    return out

def int_signature_loso(adata, genes, label_col, unit_col, positive, cells=None, n_perm=200, seed=0, min_cells=10):
    """Does a gene signature separate two groups in samples it never saw? Leave-one-sample-out logistic regression on the signature genes.
    Every unit (sample) must belong to one group. Returns per-unit mean predicted probability of `positive`, the unit-level AUC, a permutation p
    and the pooled cell-level AUC. cells: optional boolean mask (for example one cell type)."""
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    up = pd.Series(range(adata.n_vars), index=[str(v).upper() for v in adata.var_names])
    up = up[~up.index.duplicated()]
    gi = [up[g.upper()] for g in genes if g.upper() in up.index]
    if len(gi) < 3:
        return {"error": f"only {len(gi)} signature genes in the data"}
    m = np.ones(adata.n_obs, bool) if cells is None else np.asarray(cells)
    X = adata.X[m][:, gi]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    X = (X - X.mean(axis=0)) / np.maximum(X.std(axis=0), 1e-6)
    y = (adata.obs[label_col].astype(str).to_numpy()[m] == str(positive)).astype(int)
    u = adata.obs[unit_col].astype(str).to_numpy()[m]
    counts = pd.Series(u).value_counts()
    small = set(counts[counts < min_cells].index)
    if small:
        keep = ~np.isin(u, list(small))
        X, y, u = X[keep], y[keep], u[keep]
    units = sorted(np.unique(u))
    if not units:
        return {"error": f"no sample has at least {min_cells} cells"}
    ulab = {k: int(round(y[u == k].mean())) for k in units}
    mixed = [k for k in units if 0.0 < y[u == k].mean() < 1.0]
    if mixed:
        return {"error": "samples with both groups: " + ", ".join(mixed[:5])}
    if len(set(ulab.values())) < 2 or min(sum(v == 0 for v in ulab.values()), sum(v == 1 for v in ulab.values())) < 2:
        return {"error": f"need at least 2 samples with {min_cells} or more cells in each group" + (f" ({len(small)} small samples dropped)" if small else "")}
    prob = np.zeros(len(y))
    for k in units:
        te = u == k
        clf = LogisticRegression(C=0.5, max_iter=300)
        clf.fit(X[~te], y[~te])
        prob[te] = clf.predict_proba(X[te])[:, 1]
    unit_p = pd.Series({k: float(prob[u == k].mean()) for k in units})
    ul = pd.Series(ulab)
    auc = float(roc_auc_score(ul.loc[unit_p.index], unit_p))
    null = []
    vals = ul.to_numpy()
    for _ in range(n_perm):
        sh = rng.permutation(vals)
        null.append(roc_auc_score(sh, unit_p.loc[ul.index].to_numpy()) if len(set(sh)) > 1 else 0.5)
    p = (np.sum(np.array(null) >= auc) + 1) / (n_perm + 1)
    return {"unit_prob": unit_p, "unit_label": ul, "auc_units": auc, "p": float(p), "auc_cells": float(roc_auc_score(y, prob)), "n_genes": len(gi), "n_units": len(units), "n_dropped": len(small)}

def int_collect_caveats(report_dir, exclude=("interpret",)):
    """Warnings, failures and unrun checks recorded by earlier checkpoints (reports/*.json; stages named in exclude are skipped). Returns a list of dicts: stage, kind, text."""
    import os
    out = []
    if not os.path.isdir(report_dir):
        return out
    for fn in sorted(os.listdir(report_dir)):
        if not fn.endswith(".json"):
            continue
        stage = fn[:-5]
        if stage in exclude:                 # this checkpoint's own earlier report is not an earlier checkpoint
            continue
        try:
            d = json.load(open(os.path.join(report_dir, fn), encoding="utf8"))
        except Exception:
            continue
        for c in d.get("checks", []):
            if c["status"] in ("WARN", "FAIL"):
                out.append({"stage": stage, "kind": c["status"], "text": c["check"] + ": " + c["detail"][:200]})
        for nr in d.get("not_run", []):
            out.append({"stage": stage, "kind": "not run", "text": nr})
    return out

def int_perturb(adata, panel, sets=None, specs=None, max_cells=5000, seed=0, n_neighbors=15, fdr=0.05):
    """Sensitivity analysis. Re-run the variable-gene, PCA, clustering and annotation steps on a capped subsample of cells under changed settings and compare with the
    baseline labels in obs['cluster'] and obs['cell_type']. Returns dict: runs (DataFrame), agreement (types x runs), pathway (types x runs, Jaccard of
    significant terms), summary (type -> label robust/sensitive/unstable with median agreement), cluster_label."""
    import types
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import anndata as ad
    import scanpy as sc
    from sklearn.metrics import adjusted_rand_score
    from sklearn.utils.extmath import randomized_svd
    try:
        import scrna_cluster_kernel as K5
        import scrna_annotate_kernel as K6
    except ImportError:
        g = globals()
        K5 = types.SimpleNamespace(cl_graph=g["cl_graph"], cl_leiden=g["cl_leiden"])
        K6 = types.SimpleNamespace(ann_run=g["ann_run"])
    rng = np.random.default_rng(seed)
    base_idx = np.arange(adata.n_obs)
    if adata.n_obs > max_cells:
        base_idx = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
    A = adata[base_idx]
    d = adata.uns.get("dimred", {})
    n_hvg0 = int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var else 2000
    n_pcs0 = int(d.get("n_pcs", 20))
    res0 = float(adata.uns.get("clustering", {}).get("resolution", 0.8))
    base_cl = A.obs["cluster"].astype(str).to_numpy()
    base_ct = A.obs["cell_type"].astype(str).to_numpy()
    named = [t for t in sorted(set(base_ct)) if t not in ("ambiguous", "unassigned")]
    if specs is None:
        specs = [{"name": "seed 1", "seed": 1}, {"name": "seed 2", "seed": 2}, {"name": "80% of cells (a)", "frac": 0.8, "seed": 11},
                 {"name": "80% of cells (b)", "frac": 0.8, "seed": 12}, {"name": "1000 variable genes", "n_hvg": 1000},
                 {"name": "3000 variable genes", "n_hvg": 3000}, {"name": "3 fewer PCs", "n_pcs": max(n_pcs0 - 3, 3)},
                 {"name": "5 more PCs", "n_pcs": n_pcs0 + 5}, {"name": "resolution x0.6", "resolution": round(res0 * 0.6, 3)},
                 {"name": "resolution x1.5", "resolution": round(res0 * 1.5, 3)}]
    bg = int_background(A) if sets is not None else None
    fsets = int_filter_sets(sets, bg) if sets is not None else None

    def term_sets(labels, adata_sub):
        gl, _ = int_group_genes(adata_sub, labels, groups=[t for t in named if (labels == t).sum() >= 10], top_n=100)
        res = {}
        for t, gs in gl.items():
            o = int_ora(gs, fsets, bg)
            o = o[o["fdr"] < fdr]
            res[t] = set(int_trim(o)["term"]) if len(o) else set()
        return res

    base_terms = term_sets(base_ct, A) if fsets else {}
    counts = adata.layers["counts"] if "counts" in adata.layers else None
    rows, agree, pw = [], {}, {}
    for sp_ in specs:
        r = np.random.default_rng(seed * 1000 + int(sp_.get("seed", 0)))
        keep = np.arange(A.n_obs)
        if sp_.get("frac", 1.0) < 1.0:
            keep = np.sort(r.choice(A.n_obs, size=int(sp_["frac"] * A.n_obs), replace=False))
        sub = A[keep]
        n_hvg, n_pcs, res = int(sp_.get("n_hvg", n_hvg0)), int(sp_.get("n_pcs", n_pcs0)), float(sp_.get("resolution", res0))
        try:
            n_hvg = min(n_hvg, adata.n_vars - 1)
            try:
                if counts is not None:
                    tmp = ad.AnnData(sp.csr_matrix(counts[base_idx][keep]), var=pd.DataFrame(index=adata.var_names))
                    sc.pp.highly_variable_genes(tmp, n_top_genes=n_hvg, flavor="seurat_v3")
                else:
                    tmp = ad.AnnData(sub.X.copy(), var=pd.DataFrame(index=adata.var_names))
                    sc.pp.highly_variable_genes(tmp, n_top_genes=n_hvg, flavor="seurat")
                hv = tmp.var["highly_variable"].to_numpy()
            except Exception:                    # the loess fit can fail on very small or degenerate data: fall back to ranking genes by variance
                Xs = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
                hv = np.zeros(adata.n_vars, dtype=bool)
                hv[np.argsort(-Xs.var(0))[:n_hvg]] = True
            Z = sub.X[:, np.where(hv)[0]]
            Z = Z.toarray() if sp.issparse(Z) else np.asarray(Z)
            Z = np.clip((Z - Z.mean(0)) / np.maximum(Z.std(0), 1e-6), -10, 10)
            U, S, Vt = randomized_svd(Z, n_components=min(n_pcs, Z.shape[1] - 1), n_iter=5, random_state=int(sp_.get("seed", 0)))
            emb = U * S
            gr = K5.cl_graph(emb, n_neighbors, int(sp_.get("seed", 0)))
            lab = K5.cl_leiden(gr, res, int(sp_.get("seed", 0)))
            s2 = ad.AnnData(sub.X, var=pd.DataFrame(index=adata.var_names))
            s2.obs["cluster"] = pd.Categorical(lab)
            s2, _ = K6.ann_run(s2, panel, suggest_negatives=False)
            ct = s2.obs["cell_type"].astype(str).to_numpy()
            ari = float(adjusted_rand_score(base_cl[keep], lab))
            bct = base_ct[keep]
            col = {}
            for t in named:
                m = bct == t
                col[t] = float((ct[m] == t).mean()) if m.sum() >= 5 else np.nan
            nm = np.isin(bct, named)
            rows.append({"run": sp_["name"], "ari_clusters": ari, "n_clusters": int(len(set(lab))), "named_agreement": float((ct[nm] == bct[nm]).mean()) if nm.any() else np.nan, "error": ""})
            agree[sp_["name"]] = col
            if fsets:
                nt = term_sets(ct, sub)
                pw[sp_["name"]] = {t: (len(base_terms.get(t, set()) & nt.get(t, set())) / max(len(base_terms.get(t, set()) | nt.get(t, set())), 1)) if base_terms.get(t) else np.nan for t in named}
        except Exception as e:
            rows.append({"run": sp_["name"], "ari_clusters": np.nan, "n_clusters": np.nan, "named_agreement": np.nan, "error": str(e)[:150]})
    runs = pd.DataFrame(rows)
    ag = pd.DataFrame(agree).reindex(named)
    pwd = pd.DataFrame(pw).reindex(named) if pw else pd.DataFrame(index=named)
    summ = {}
    for t in named:
        v = ag.loc[t].dropna().to_numpy()
        if len(v) == 0:
            summ[t] = {"label": "not tested", "median": np.nan, "min": np.nan}
            continue
        med, mn = float(np.median(v)), float(v.min())
        summ[t] = {"label": "robust" if (med >= 0.9 and mn >= 0.75) else ("sensitive" if med >= 0.7 else "unstable"), "median": med, "min": mn}
    ari_med = float(np.nanmedian(runs["ari_clusters"])) if runs["ari_clusters"].notna().any() else np.nan
    cl_label = "robust" if ari_med >= 0.8 else ("sensitive" if ari_med >= 0.6 else "unstable")
    return {"runs": runs, "agreement": ag, "pathway": pwd, "summary": summ, "cluster_label": cl_label, "ari_median": ari_med, "n_cells": int(A.n_obs),
            "baseline_terms": {t: sorted(v) for t, v in base_terms.items()}}

def int_legends(info):
    """Plain-language legends built from the computed numbers; one-line captions are produced separately."""
    return {
        "int_pathways": (f"Each row is a cell type and each column a pathway that its marker genes share more than chance allows ({info['n_terms']} shown). Dot size is how many marker genes fall in the "
                         f"pathway; colour is the strength of evidence after correcting for the number of pathways tested (brighter = stronger). The comparison is against the {info['n_bg']:,} genes "
                         f"actually detected here, not the whole genome. A pathway name is a hint about what the genes have in common, not proof that the pathway is active."),
        "int_programs": (f"Each row is a group of genes found by the data itself that rise and fall together ({info['n_programs']} programs); each column is a cell type. Colour shows how much of the "
                         f"cell type's activity is spent on that program. The gene names beside each row are its strongest genes. Programs were checked for reproducibility across random restarts "
                         f"(similarity {info['prog_stability']:.2f})."),
        "int_robustness": ("Each setting that could reasonably have been chosen differently was changed one at a time and the analysis re-run on the same cells. Left: how much the cluster assignment "
                           "agrees with the original (1 = identical). Right: for each cell type, the share of its cells that still get the same name. Rows that stay bright across all settings are robust; "
                           "rows that fade depend on the choice."),
        "int_validation": ("Left: each known fact listed in advance as a check, shown as met, not met or not testable. Right: how well a gene signature separates the two groups in samples it was not "
                           "trained on (area under the curve; 0.5 = no better than chance, 1 = perfect)."),
    }

def int_captions(info):
    return {
        "int_pathways.png": f"Dot plot of pathways (columns) enriched among the marker genes of each cell type (rows); dot size = marker genes in the pathway, colour = -log10 adjusted p-value.",
        "int_programs.png": f"Heat map of the usage of {info['n_programs']} data-driven gene programs (rows) in each cell type (columns), with each program's top genes.",
        "int_robustness.png": f"Left: adjusted Rand index of the clusters against the baseline for each changed setting. Right: share of each cell type's cells keeping the same name under each setting.",
        "int_validation.png": "Left: outcome of each listed expectation. Right: leave-one-sample-out area under the curve of the gene signature.",
    }

def int_findings(adata, ctype, ora, prog, perturb, expect, loso, ranked_res=None):
    """List of dicts: finding, evidence, robustness, caveats. Built only from computed results."""
    import numpy as np
    out = []
    tab = adata.uns.get("annotation")
    summ = perturb["summary"] if perturb else {}
    for t in sorted(set(adata.obs[ctype].astype(str)) - {"ambiguous", "unassigned"}):
        n = int((adata.obs[ctype].astype(str) == t).sum())
        ev, cav = [f"{n:,} cells"], []
        if tab is not None:
            sub = tab[tab["label"] == t]
            ev.append("clusters " + ", ".join(str(c) for c in sub["cluster"]) + " (" + ", ".join(sorted(set(sub["tier"]))) + ")")
            for r in sub.itertuples():
                if r.tier != "confident":
                    cav.append(f"cluster {r.cluster}: {r.reason}")
        s = summ.get(t, {})
        rob = s.get("label", "not tested")
        if s and s.get("median") == s.get("median"):
            ev.append(f"same name in a median {100 * s['median']:.0f}% of cells across settings (worst {100 * s['min']:.0f}%)")
        out.append({"finding": f"{t} cells are present", "evidence": "; ".join(ev), "robustness": rob, "caveats": "; ".join(cav) or "none recorded"})
    for t, df in (ora or {}).items():
        for r in df.head(2).itertuples():
            pr = None
            if perturb is not None and len(perturb.get("pathway", [])) and t in perturb["pathway"].index:
                v = perturb["pathway"].loc[t].dropna()
                pr = float(np.median(v)) if len(v) else None
            rob = "not tested" if pr is None else ("robust" if pr >= 0.6 else ("sensitive" if pr >= 0.3 else "unstable"))
            cav = ["enrichment describes shared genes, not pathway activity"]
            low = r.term.lower()
            if any(w in low for w in ("ribosom", "oxidative phosph", "spliceosome", "proteasome", "mitochond")):
                cav.append("housekeeping-type term: differences often reflect cell size, translation or RNA quality rather than identity")
            if "[hsa05" in low or "[hsa049" in low:
                cav.append("KEGG disease or infection pathway: it contains generic immune genes, so the name does not mean the disease is present")
            out.append({"finding": f"{t}: marker genes are enriched for {r.term}", "evidence": f"{r.overlap} of {r.size} pathway genes among the markers, fold enrichment {r.fold_enrichment:.1f}, adjusted p {r.fdr:.1e}",
                        "robustness": rob + ("" if pr is None else f" (term set overlap {pr:.2f} across settings)"), "caveats": "; ".join(cav)})
    if prog is not None:
        for p in prog.get("dominant", []):
            out.append({"finding": f"{p['program']} is used mostly by {p['type']}", "evidence": f"{100 * p['share']:.0f}% of its usage; top genes {', '.join(p['genes'][:5])}",
                        "robustness": "reproducible across restarts" if prog.get("qualified") else "restarts disagree", "caveats": "programs are descriptive and not tied to a pathway unless tested"})
    for e in expect or []:
        out.append({"finding": f"Expectation: {e['expectation']}", "evidence": e["detail"], "robustness": "n/a", "caveats": "none" if e["status"] == "met" else e["status"]})
    for name, r in (loso or {}).items():
        if "error" in r:
            continue
        out.append({"finding": f"{name} separates the groups in unseen samples", "evidence": f"unit-level AUC {r['auc_units']:.2f} over {r['n_units']} samples (permutation p {r['p']:.3f}); cell-level AUC {r['auc_cells']:.2f}",
                    "robustness": "held-out samples", "caveats": "samples are few; AUC is limited by the number of samples, not cells" + (f"; {r['n_dropped']} sample(s) with too few cells dropped" if r.get("n_dropped") else "")})
    return out

def int_figures(fig_dir, info, ora, prog, perturb, expect, loso, ctype_order):
    """Four presentation figures. Returns the list of file names written."""
    import os
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rc = {"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
          "axes.spines.top": False, "axes.spines.right": False, "axes.titleweight": "normal", "axes.titlelocation": "left", "legend.frameon": False}
    os.makedirs(fig_dir, exist_ok=True)
    files = []
    with matplotlib.rc_context(rc):
        rows = []
        for t in ctype_order:
            df = (ora or {}).get(t)
            if df is not None:
                for r in df.head(3).itertuples():
                    rows.append((t, r.term, r.overlap, -np.log10(max(r.fdr, 1e-300))))
        if rows:
            terms = list(dict.fromkeys(r[1] for r in rows))
            types = [t for t in ctype_order if any(r[0] == t for r in rows)]
            fig, ax = plt.subplots(figsize=(0.32 * len(terms) + 2.2, 0.32 * len(types) + 3.4))
            xs = [terms.index(r[1]) for r in rows]
            ys = [types.index(r[0]) for r in rows]
            sc_ = ax.scatter(xs, ys, s=[20 + 12 * r[2] for r in rows], c=[r[3] for r in rows], cmap="viridis", linewidths=0)
            ax.set_xticks(range(len(terms)))
            ax.set_xticklabels([t[:48] for t in terms], rotation=70, ha="right")
            ax.set_yticks(range(len(types)))
            ax.set_yticklabels(types)
            ax.set_ylim(len(types) - 0.5, -0.5)
            ax.margins(0.05)
            ax.set_title("pathways enriched among each cell type's markers")
            fig.colorbar(sc_, ax=ax, label="-log10 adjusted p", shrink=0.6)
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, "int_pathways.png"), dpi=300)
            plt.close(fig)
            files.append("int_pathways.png")
        if prog is not None:
            U = prog["by_type"]
            fig, ax = plt.subplots(figsize=(0.55 * U.shape[1] + 3.6, 0.32 * U.shape[0] + 1.6))
            im = ax.imshow(U.to_numpy(), cmap="viridis", aspect="auto", vmin=0)
            ax.set_xticks(range(U.shape[1]))
            ax.set_xticklabels(U.columns, rotation=60, ha="right")
            ax.set_yticks(range(U.shape[0]))
            ax.set_yticklabels([f"{p}: {', '.join(prog['top_genes'][p][:4])}" for p in U.index], fontsize=6.5)
            ax.set_title("gene program usage by cell type")
            fig.colorbar(im, ax=ax, label="share of the cell type's program usage", shrink=0.6)
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, "int_programs.png"), dpi=300)
            plt.close(fig)
            files.append("int_programs.png")
        if perturb is not None and len(perturb["runs"]):
            runs, ag = perturb["runs"], perturb["agreement"]
            fig, axes = plt.subplots(1, 2, figsize=(8.2, 0.3 * max(len(runs), len(ag)) + 2.2), gridspec_kw={"width_ratios": [1, 1.6]})
            axes[0].barh(range(len(runs)), runs["ari_clusters"].fillna(0), color="0.55")
            axes[0].axvline(0.8, color="0.6", ls=":", lw=1)
            axes[0].set_yticks(range(len(runs)))
            axes[0].set_yticklabels(runs["run"])
            axes[0].set_xlim(0, 1.02)
            axes[0].set_xlabel("cluster agreement (ARI)")
            axes[0].invert_yaxis()
            if ag.shape[0]:
                im = axes[1].imshow(ag.to_numpy(dtype=float), cmap="viridis", vmin=0, vmax=1, aspect="auto")
                axes[1].set_xticks(range(ag.shape[1]))
                axes[1].set_xticklabels(ag.columns, rotation=60, ha="right")
                axes[1].set_yticks(range(ag.shape[0]))
                axes[1].set_yticklabels([f"{t} ({perturb['summary'][t]['label']})" for t in ag.index])
                fig.colorbar(im, ax=axes[1], label="cells keeping their name", shrink=0.7)
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, "int_robustness.png"), dpi=300)
            plt.close(fig)
            files.append("int_robustness.png")
        aucs = {k: v for k, v in (loso or {}).items() if "error" not in v}
        if expect or aucs:
            fig, axes = plt.subplots(1, 2, figsize=(8.4, 0.3 * max(len(expect or []), len(aucs), 3) + 1.8), gridspec_kw={"width_ratios": [1.6, 1]})
            col = {"met": "tab:green", "not met": "tab:red", "not testable": "0.6"}
            for i, e in enumerate(expect or []):
                axes[0].scatter([0], [i], color=col[e["status"]], s=40)
                axes[0].text(0.08, i, f"{e['expectation'][:60]} ({e['status']})", va="center", fontsize=7)
            axes[0].set_xlim(-0.1, 2.2)
            axes[0].set_ylim(len(expect or []) - 0.5, -0.5)
            axes[0].axis("off")
            axes[0].set_title("expectations")
            if aucs:
                names = list(aucs)
                axes[1].barh(range(len(names)), [aucs[n]["auc_units"] for n in names], color="0.55")
                axes[1].axvline(0.5, color="0.6", ls=":", lw=1)
                axes[1].set_yticks(range(len(names)))
                axes[1].set_yticklabels(names)
                axes[1].set_xlim(0, 1.02)
                axes[1].set_xlabel("held-out sample AUC")
                axes[1].invert_yaxis()
            else:
                axes[1].axis("off")
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, "int_validation.png"), dpi=300)
            plt.close(fig)
            files.append("int_validation.png")
    return files

def int_run(adata, gene_sets=None, panel=None, expectations=None, ranked=None, signature=None, cell_type_key="cell_type", perturb=True, perturb_specs=None,
            max_cells=5000, caveats=None, nmf_k=(5, 8, 12), top_n=100, fdr=0.05, min_group_cells=20, seed=0, fig_dir=None):
    """Interpret the annotated object: pathway over-representation of cell-type markers, data-driven gene programs, expectations, signature validation on held-out samples,
    sensitivity analysis and a findings table. gene_sets: {name: set of symbols}; expectations: DataFrame (see int_expectations); ranked: {group: Series} from a
    differential-expression result; signature: {name, genes, label_col, positive, unit_col or sample_col, by_type}. Returns (adata, records, details)."""
    import os
    import numpy as np
    import pandas as pd
    recs, det = [], {}
    if cell_type_key not in adata.obs:
        return adata, [int_record("gene sets", "FAIL", f"obs['{cell_type_key}'] not found: run checkpoint 6 first")], det
    lab = adata.obs[cell_type_key].astype(str).to_numpy()
    sizes = pd.Series(lab).value_counts()
    ctypes = [t for t in sizes.index if t not in ("ambiguous", "unassigned") and sizes[t] >= min_group_cells]
    bg = int_background(adata)
    fsets = None
    if gene_sets:
        fsets = int_filter_sets(gene_sets, bg)
        recs.append(int_record("gene sets", "PASS" if len(fsets) >= 10 else "WARN", f"{len(fsets)} of {len(gene_sets)} sets keep 10 to 500 genes among the {len(bg):,} detected genes"))
    recs.append(int_record("background", "INFO", f"{len(bg):,} genes detected in at least 3 cells are the background for every enrichment test"))
    glists, gdf = int_group_genes(adata, lab, groups=ctypes, top_n=top_n)
    det["marker_genes"] = gdf
    ora, ora_all = {}, []
    if fsets:
        for t, gs in glists.items():
            o = int_ora(gs, fsets, bg)
            o["cell_type"] = t
            ora_all.append(o)
            sig = o[o["fdr"] < fdr] if len(o) else o
            ora[t] = int_trim(sig) if len(sig) else sig
        det["ora"] = pd.concat(ora_all) if ora_all else pd.DataFrame()
        n_sig = sum(len(v) for v in ora.values())
        with_hits = [t for t, v in ora.items() if len(v)]
        recs.append(int_record("pathway analysis", "PASS" if with_hits else "WARN",
                               (f"{n_sig} non-redundant pathways below FDR {fdr} in {len(with_hits)} of {len(ora)} cell types" if with_hits else
                                f"no pathway below FDR {fdr} for any of {len(ora)} cell types; marker lists may be short or the collection may not cover this biology") +
                               f"; markers per type: median {int(np.median([len(v) for v in glists.values()])) if glists else 0}"))
    if ranked:
        rows = []
        for g, s in ranked.items():
            if fsets:
                x = int_gsea(s, fsets, seed=seed)
                x["group"] = g
                rows.append(x)
        if rows:
            det["gsea"] = pd.concat(rows)
            recs.append(int_record("pre-ranked enrichment", "PASS", f"{int((det['gsea']['fdr'] < fdr).sum())} set(s) below FDR {fdr} across {len(rows)} ranked lists"))
    prog = None
    try:
        nm = int_nmf(adata, k_values=nmf_k, seed=seed)
        U = nm["usage"]
        by = pd.DataFrame({t: U[lab == t].mean() for t in ctypes}) if ctypes else pd.DataFrame(index=U.columns)
        shares = by.div(by.sum(axis=1).replace(0, np.nan), axis=0)
        dom = []
        for p in by.index:
            if len(by.columns) and shares.loc[p].max() >= 0.5:
                dom.append({"program": p, "type": shares.loc[p].idxmax(), "share": float(shares.loc[p].max()),
                            "genes": list(nm["loadings"].loc[p].sort_values(ascending=False).index[:8])})
        prog = {"k": nm["k"], "table": nm["table"], "qualified": nm["qualified"], "by_type": shares.fillna(0), "dominant": dom, "usage": U, "loadings": nm["loadings"],
                "top_genes": {p: list(nm["loadings"].loc[p].sort_values(ascending=False).index[:6]) for p in nm["loadings"].index}}
        stab = float(nm["table"].loc[nm["table"]["k"] == nm["k"], "stability"].iloc[0])
        recs.append(int_record("gene programs", "PASS" if nm["qualified"] else "WARN",
                               f"{nm['k']} programs; restart similarity {stab:.2f} (target 0.8); {len(dom)} used mostly (50% or more) by one cell type"))
        adata.obsm["X_programs"] = U.to_numpy().astype(np.float32)
        adata.uns["program_names"] = list(U.columns)
    except Exception as e:
        recs.append(int_record("gene programs", "WARN", "programs could not be computed: " + str(e)[:150]))
    expect_res = None
    if expectations is not None and len(expectations):
        expect_res = int_expectations(adata, expectations, sets=fsets, ora_hits={t: list(v["term"]) for t, v in ora.items()} if ora else None)
        n = pd.Series([e["status"] for e in expect_res]).value_counts()
        bad = [e["expectation"] for e in expect_res if e["status"] == "not met"]
        recs.append(int_record("expectations", "WARN" if bad else "PASS", f"{n.get('met', 0)} met, {n.get('not met', 0)} not met, {n.get('not testable', 0)} not testable" +
                               (": not met: " + "; ".join(bad[:5]) if bad else "")))
        det["expectations"] = pd.DataFrame(expect_res)
    loso = None
    if signature:
        loso = {}
        sg = signature
        unit = sg.get("unit_col")
        if not unit:
            adata.obs["_unit"] = adata.obs[sg["sample_col"]].astype(str) + "_" + adata.obs[sg["label_col"]].astype(str)
            unit = "_unit"
        for nm_ in (ctypes if sg.get("by_type") else [None]):
            cells = None if nm_ is None else (lab == nm_)
            r = int_signature_loso(adata, sg["genes"], sg["label_col"], unit, sg["positive"], cells=cells, seed=seed)
            loso[f"{sg.get('name', 'signature')}" + ("" if nm_ is None else f" in {nm_}")] = r
        if "_unit" in adata.obs:
            del adata.obs["_unit"]
        ok = {k: v for k, v in loso.items() if "error" not in v}
        if ok:
            worst = min(v["auc_units"] for v in ok.values())
            skipped = [k for k, v in loso.items() if "error" in v]
            recs.append(int_record("signature validation", "PASS" if worst >= 0.8 else "WARN",
                                   f"unit-level AUC {worst:.2f} to {max(v['auc_units'] for v in ok.values()):.2f} across {len(ok)} held-out analyses (lowest must reach 0.8); permutation p up to {max(v['p'] for v in ok.values()):.3f}" +
                                   (f"; not testable: {', '.join(skipped)} (too few samples with enough cells in each group)" if skipped else "")))
        else:
            recs.append(int_record("signature validation", "WARN", "not testable: " + "; ".join(sorted(set(v.get("error", "") for v in loso.values())))[:200]))
    pert = None
    if perturb:
        pnl = panel if panel is not None else adata.uns.get("annotation_panel")
        if pnl is None:
            recs.append(int_record("robustness: clusters", "WARN", "no marker panel available for the sensitivity analysis"))
        else:
            pert = int_perturb(adata, pnl, sets=gene_sets, specs=perturb_specs, max_cells=max_cells, seed=seed)
            fails = pert["runs"][pert["runs"]["error"] != ""]
            recs.append(int_record("robustness: clusters", "PASS" if pert["cluster_label"] == "robust" else "WARN",
                                   f"median adjusted Rand index {pert['ari_median']:.2f} across {len(pert['runs'])} changed settings on {pert['n_cells']:,} cells ({pert['cluster_label']})" +
                                   (f"; {len(fails)} setting(s) failed" if len(fails) else "")))
            lab_counts = pd.Series([v["label"] for v in pert["summary"].values()]).value_counts()
            weak = [t for t, v in pert["summary"].items() if v["label"] != "robust"]
            recs.append(int_record("robustness: cell-type calls", "PASS" if not weak else "WARN",
                                   f"{lab_counts.get('robust', 0)} robust, {lab_counts.get('sensitive', 0)} sensitive, {lab_counts.get('unstable', 0)} unstable of {len(pert['summary'])} named cell types" +
                                   (": " + "; ".join(f"{t} ({pert['summary'][t]['label']}, median {100 * pert['summary'][t]['median']:.0f}%)" for t in weak[:6]) if weak else "")))
            if len(pert["pathway"]) and pert["pathway"].notna().any().any():
                med = float(np.nanmedian(pert["pathway"].to_numpy(dtype=float)))
                recs.append(int_record("robustness: pathways", "PASS" if med >= 0.6 else "WARN",
                                       f"median overlap of significant pathway sets with the baseline {med:.2f} across settings (0.6 or more = robust)"))
            det["perturb"] = pert
    cav = caveats or []
    recs.append(int_record("caveats carried", "INFO", f"{sum(c['kind'] == 'WARN' for c in cav)} warning(s), {sum(c['kind'] == 'FAIL' for c in cav)} failure(s) and {sum(c['kind'] == 'not run' for c in cav)} unrun check(s) from earlier checkpoints"))
    findings = int_findings(adata, cell_type_key, ora, prog, pert, expect_res, loso)
    det.update(findings=findings, ora_trimmed=ora, programs=prog, loso=loso, caveats=cav)
    adata.uns["interpretation"] = {"n_types": len(ctypes), "n_findings": len(findings), "n_background": len(bg), "cluster_robustness": (pert or {}).get("cluster_label", "not tested")}
    if fig_dir is not None:
        os.makedirs(fig_dir, exist_ok=True)
        pd.DataFrame(findings).to_csv(os.path.join(fig_dir, "findings.csv"), index=False)
        if "ora" in det and len(det["ora"]):
            det["ora"].to_csv(os.path.join(fig_dir, "pathway_enrichment.csv"), index=False)
        if prog is not None:
            prog["loadings"].T.to_csv(os.path.join(fig_dir, "program_loadings.csv"))
            prog["by_type"].to_csv(os.path.join(fig_dir, "program_usage_by_type.csv"))
        if pert is not None:
            pert["runs"].to_csv(os.path.join(fig_dir, "robustness_runs.csv"), index=False)
            pert["agreement"].to_csv(os.path.join(fig_dir, "robustness_cell_type_agreement.csv"))
        if expect_res:
            pd.DataFrame(expect_res).to_csv(os.path.join(fig_dir, "expectations.csv"), index=False)
        info = {"n_terms": sum(min(len(v), 3) for v in ora.values()), "n_bg": len(bg), "n_programs": prog["k"] if prog else 0, "prog_stability": float(prog["table"]["stability"].max()) if prog else 0.0}
        files = int_figures(fig_dir, info, ora, prog, pert, expect_res, loso, ctypes)
        leg, cap = int_legends(info), int_captions(info)
        from_name = {"int_pathways": "int_pathways.png", "int_programs": "int_programs.png", "int_robustness": "int_robustness.png", "int_validation": "int_validation.png"}
        with open(os.path.join(fig_dir, "figure_legends.md"), "w", encoding="utf8") as f:
            f.write("# Figure legends\n\n" + "\n\n".join(f"**{k.replace('_', ' ')}** ({from_name[k]}): {v}" for k, v in leg.items() if from_name[k] in files) + "\n")
        with open(os.path.join(fig_dir, "figure_captions.json"), "w", encoding="utf8") as f:
            json.dump({k: v for k, v in cap.items() if k in files}, f, indent=1)
    return adata, recs, det

def write_int_report(records, findings=None, caveats=None, prefix="interpretation_report", scope=None):
    """Write <prefix>.json and <prefix>.md and return the overall status. Opens with the verdict and why, then the findings with evidence and robustness,
    the checks with what they test, the checks that did not run, and the caveats carried from earlier checkpoints."""
    cat = int_catalog()
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
        json.dump({"overall": overall, "why": why, "counts": count, "not_run": list(not_run), "checks": records, "findings": findings or []}, f, indent=1, default=str)
    L = [f"# Biological interpretation and validation report: {overall}", "", f"**Why:** {why}", "",
         f"**Checks run ({len(records)}):** {count['PASS']} pass, {count['WARN']} warn, {count['FAIL']} fail, {count['INFO']} info", ""]
    if findings:
        good = [x for x in findings if x["robustness"].startswith(("robust", "reproducible", "held-out"))]
        L += ["## Findings, their evidence and how far to trust them", "",
              f"{len(findings)} findings; {len(good)} carry a robust or validated label. A finding labelled sensitive or unstable depends on a choice that could reasonably have been made differently.", "",
              "| finding | evidence | robustness | caveats |", "|---|---|---|---|"]
        for x in findings:
            L.append("| " + " | ".join(str(x[k]).replace("|", "/") for k in ("finding", "evidence", "robustness", "caveats")) + " |")
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
    with open(prefix + ".md", "w") as f:
        f.write("\n".join(L) + "\n")
    return overall
