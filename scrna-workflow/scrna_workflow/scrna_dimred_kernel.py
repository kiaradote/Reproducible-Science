import json

def dr_record(check, status, detail, group=None):
    """One result row. status is PASS, WARN, FAIL or INFO."""
    r = {"check": check, "status": status, "detail": detail}
    if group is not None:
        r["group"] = group
    return r

def dr_catalog():
    """Checks this checkpoint can run: name -> what it tests."""
    return {
        "PCA on selected genes": "PCA ran on the scaled selected genes and gave finite results",
        "number of PCs chosen from the data": "the PC count is the smallest number capturing a set share of the real structure in this dataset",
        "variance threshold reached": "the chosen share of variation is reached within the PCs computed",
        "chosen PCs beat shuffled data": "every PC kept rises above what shuffled gene values would give",
        "PCs driven by technical covariates": "retained PCs are not mainly depth, genes-detected or mitochondrial effects",
        "batch structure in PCs": "how much of the retained PC variation is explained by the batch column",
        "batch correction effect": "correction lowered batch structure and improved mixing",
        "neighbor graph fragments": "no cells are stranded in tiny disconnected pieces of the neighbor graph",
        "UMAP finite": "UMAP coordinates are finite for every cell",
    }

def dr_scale(adata, max_value=10.0):
    """Dense, scaled (zero mean, unit variance, clipped) matrix of the selected genes. Returns (matrix, gene mask)."""
    import numpy as np
    import scipy.sparse as sp
    hv = adata.var["highly_variable"].to_numpy()
    X = adata.X[:, np.where(hv)[0]]
    X = (X.toarray() if sp.issparse(X) else np.asarray(X)).astype(np.float32)
    sd = X.std(axis=0, ddof=1)
    sd[sd == 0] = 1
    X = np.clip((X - X.mean(axis=0)) / sd, -max_value, max_value)
    X -= X.mean(axis=0)
    return X, hv

def dr_svd(X, n, seed):
    from sklearn.utils.extmath import randomized_svd
    return randomized_svd(X, n_components=n, n_iter=5, random_state=seed)

def dr_choose_pcs(X, eig, n_perm=10, seed=0):
    """Data-driven PC count. Parallel analysis: each gene is shuffled across cells (destroying structure), PCA is repeated n_perm
    times, and a PC is kept while its variance exceeds the largest shuffled variance at the same rank. The elbow (point farthest
    from the line joining the first and last variance) is returned as a cross-check. Returns (n_null, n_elbow, null_max)."""
    import numpy as np
    rng = np.random.default_rng(seed)
    nmax = len(eig)
    null = np.zeros((n_perm, nmax))
    for p in range(n_perm):
        S = dr_svd(rng.permuted(X, axis=0), nmax, seed + 1 + p)[1]
        null[p] = S ** 2 / (X.shape[0] - 1)
    null_max = null.max(axis=0)
    above = eig > null_max
    n_null = nmax if above.all() else int(np.argmin(above))
    y = (eig - eig[-1]) / max(eig[0] - eig[-1], 1e-12)
    x = np.arange(nmax) / max(nmax - 1, 1)
    n_elbow = int(np.argmax((1 - x) - y)) + 1
    return n_null, n_elbow, null_max

def dr_pcs_for_threshold(eig, null_max, ratio, threshold=0.9, basis="signal"):
    """Smallest number of PCs whose cumulative share reaches the threshold. basis 'signal': share of the variation that rises above
    the shuffled-data null (noise floor removed); 'total': share of all variation. Returns (n, reached, cumulative curve)."""
    import numpy as np
    if basis == "signal":
        above = eig > null_max
        lead = len(eig) if above.all() else int(np.argmin(above))      # only the leading run of PCs that beat the null counts
        w = np.clip(eig - null_max, 0, None)
        w[lead:] = 0
        if w.sum() <= 0:
            return 0, True, np.zeros(len(eig))
        cum = np.cumsum(w) / w.sum()
    else:
        cum = np.cumsum(ratio)
    hit = np.where(cum >= threshold)[0]
    return (int(hit[0]) + 1 if len(hit) else len(eig)), bool(len(hit)), cum

def dr_pc_corr(pcs, obs, cols):
    """Correlation of each PC with technical covariates present in obs (depth and genes-detected on a log scale)."""
    import numpy as np
    names, rows = [], []
    Pz = (pcs - pcs.mean(0)) / np.maximum(pcs.std(0), 1e-12)
    for c in cols:
        if c not in obs:
            continue
        v = obs[c].to_numpy(dtype=float)
        if not np.isfinite(v).all() or v.std() == 0:
            continue
        v = np.log1p(v) if c in ("total", "n_genes") else v
        z = (v - v.mean()) / v.std()
        names.append(c)
        rows.append(Pz.T @ z / len(z))
    return names, (np.array(rows) if rows else np.zeros((0, pcs.shape[1])))

def dr_batch_r2(emb, labels):
    """Share of each component's variance explained by the batch label (one-way ANOVA R squared)."""
    import numpy as np
    codes = np.unique(np.asarray(labels).astype(str), return_inverse=True)[1]
    counts = np.bincount(codes).astype(float)
    r2 = np.zeros(emb.shape[1])
    for j in range(emb.shape[1]):
        col = emb[:, j].astype(float)
        means = np.bincount(codes, weights=col) / counts
        ssb = (counts * (means - col.mean()) ** 2).sum()
        sst = ((col - col.mean()) ** 2).sum()
        r2[j] = ssb / max(sst, 1e-12)
    return r2

def dr_batch_share(emb, labels):
    """Variance-weighted mean batch R squared over the components of an embedding."""
    import numpy as np
    w = emb.var(axis=0)
    return float((dr_batch_r2(emb, labels) * w).sum() / max(w.sum(), 1e-12))

def dr_mixing(emb, labels, k=30, max_cells=3000, seed=0):
    """Batch mixing in the neighborhood of each cell: mean entropy of neighbor batch labels divided by the entropy of the global
    batch proportions. 1 means neighbors look like the whole dataset; near 0 means neighbors share the cell's batch."""
    import numpy as np
    from sklearn.neighbors import NearestNeighbors
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(labels), size=min(max_cells, len(labels)), replace=False)
    lab = np.unique(np.asarray(labels).astype(str)[idx], return_inverse=True)[1]
    E = np.asarray(emb)[idx]
    nb = NearestNeighbors(n_neighbors=min(k + 1, len(idx))).fit(E).kneighbors(E, return_distance=False)[:, 1:]
    nl = int(lab.max()) + 1
    p = np.stack([np.bincount(lab[row], minlength=nl) for row in nb]).astype(float)
    p /= p.sum(axis=1, keepdims=True)
    h = -(np.where(p > 0, p * np.log(np.where(p > 0, p, 1)), 0)).sum(axis=1)
    g = np.bincount(lab, minlength=nl) / len(lab)
    hg = -(g[g > 0] * np.log(g[g > 0])).sum()
    return float(h.mean() / hg) if hg > 0 else float("nan")

def dr_harmony(P, labels, batch_key, seed=0):
    """Harmony batch correction of a PC matrix (cells x PCs). Returns a corrected matrix of the same shape."""
    import numpy as np
    import pandas as pd
    import harmonypy
    meta = pd.DataFrame({batch_key: np.asarray(labels).astype(str)})
    try:
        ho = harmonypy.run_harmony(P, meta, batch_key, random_state=seed)
    except TypeError:
        ho = harmonypy.run_harmony(P, meta, batch_key)
    Z = np.asarray(ho.Z_corr)
    if Z.shape == (P.shape[1], P.shape[0]):
        Z = Z.T
    return Z.astype(np.float32)

def dr_legends(info):
    """Plain-language figure legends built from the computed numbers. No biological interpretation is made."""
    n, g = info["n_pcs"], info["n_genes_selected"]
    thr = 100 * info["threshold"]
    sig = info["basis"] == "signal"
    rule = (f"the smallest number that together capture {thr:.0f}% of the real structure (the variation that rises above the dashed line)" if sig
            else f"the smallest number that together capture {thr:.0f}% of all variation")
    why = ""
    if info["n_null"] < info["min_pcs"]:
        why = f" Fewer than {info['min_pcs']} PCs rose above that line, so the minimum of {info['min_pcs']} is used."
    elif not info["reached"]:
        why = f" The {thr:.0f}% target was not reached within the {info['nmax']} PCs computed, so all {info['nmax']} are used."
    bk = info["batch_key"]
    out = {
        "pca_variance": (f"Each bar is one principal component (PC): a direction along which cells differ from one another, ordered from the "
                         f"largest source of variation (PC1, {100 * info['var_ratio'][0]:.1f}% of the total) downward. The dashed line shows how much "
                         f"variation a PC would capture if gene values were shuffled between cells, which destroys any real structure; "
                         f"{info['n_null']} PCs rise above it.{why} {n} PCs are used downstream: {rule}. Together they hold "
                         f"{100 * info['cum_at_n']:.0f}% of all variation in the {g} selected genes. The right panel shows the running total of "
                         f"{'the real structure' if sig else 'all variation'}; the red line marks the PCs used."),
        "pca_covariates": (f"Each column is a PC and each row is a measurement of how a cell was captured ({', '.join(info['cov_names']) or 'none available'}"
                           f"{'; batch' if info['has_batch'] else ''}). Brighter squares mean the PC rises and falls with that measurement, so it may reflect "
                           f"measurement rather than cell biology. {len(info['tech_pcs'])} of the {n} PCs, holding {100 * info['tech_var_share']:.0f}% of the retained "
                           f"variation, have a link stronger than {info['tech_r_limit']}."),
        "pca_scatter": (f"Each dot is a cell placed by its first two PCs; cells near each other have similar expression of the selected genes. Colour shows "
                        f"{info['color_what']}. " + (f"Batch explained {100 * info['share_before']:.0f}% of the variation in the retained PCs before correction"
                        + (f" and {100 * info['share_after']:.0f}% after batch correction (right panel)." if info["corrected"] else " (no correction applied).")
                           if info["has_batch"] else "No batch column was given.") + " The axes have no units."),
        "umap": (f"UMAP squeezes the similarity between cells, measured on the first {n} PCs"
                 f"{' after batch correction' if info['corrected'] else ''}, into two dimensions for viewing. Each dot is a cell and nearby cells have similar "
                 f"expression. Distances between separated groups and the sizes of groups are not quantitative, so UMAP is for viewing and not for measuring; "
                 f"all analysis uses the PCs. Colour shows {info['color_what']}" + (" (left) and total signal per cell on a log scale (right)." if info["color_kind"] == "batch" else ".")),
    }
    return out

def dr_axes_arrows(ax, xlabel, ylabel):
    """Embedding axes: no ticks, a small corner arrow pair names the axes."""
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    arrow = dict(arrowstyle="->", lw=0.8)
    ax.annotate("", xy=(0.14, 0.0), xytext=(0.0, 0.0), xycoords="axes fraction", arrowprops=arrow)
    ax.annotate("", xy=(0.0, 0.14), xytext=(0.0, 0.0), xycoords="axes fraction", arrowprops=arrow)
    ax.text(0.07, -0.03, xlabel, transform=ax.transAxes, ha="center", va="top", fontsize=7)
    ax.text(-0.03, 0.07, ylabel, transform=ax.transAxes, ha="right", va="center", rotation=90, fontsize=7)

def dr_figures(adata, info, fig_dir, batch_key=None):
    """Four presentation figures: variance by PC, PC links to technical measurements, PCA scatter, UMAP. Returns file names."""
    import os
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    rc = {"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 8,
          "axes.spines.top": False, "axes.spines.right": False, "xtick.direction": "out", "ytick.direction": "out",
          "axes.titleweight": "normal", "axes.titlelocation": "left", "legend.frameon": False}
    os.makedirs(fig_dir, exist_ok=True)
    files = []
    ratio, null = np.asarray(info["var_ratio"]), np.asarray(info["null_ratio"])
    n = info["n_pcs"]
    depth = None
    for c in ("total", "n_genes"):
        if c in adata.obs:
            depth = np.log1p(adata.obs[c].to_numpy(dtype=float))
            break
    lab = None
    if batch_key is not None and batch_key in adata.obs and 2 <= adata.obs[batch_key].nunique() <= 12:
        lab = adata.obs[batch_key].astype(str).to_numpy()
    cats = sorted(set(lab)) if lab is not None else []
    cmap = plt.get_cmap("tab10" if len(cats) <= 10 else "tab20")
    pt = 3 if adata.n_obs < 20000 else 1.5
    with matplotlib.rc_context(rc):
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
        x = np.arange(1, len(ratio) + 1)
        axes[0].bar(x, 100 * ratio, color="0.55", width=0.8)
        axes[0].plot(x, 100 * null, color="k", ls="--", lw=1, label="shuffled genes")
        axes[0].axvline(n + 0.5, color="tab:red", lw=1)
        axes[0].text(n + 1, 100 * ratio[0] * 0.9, f"{n} PCs used", color="tab:red", fontsize=8, va="top")
        axes[0].set_xlabel("principal component")
        axes[0].set_ylabel("variation captured (%)")
        axes[0].set_title(f"{info['n_null']} of {len(ratio)} components stand out from shuffled data")
        axes[0].legend(loc="center right")
        curve = np.asarray(info["cum_curve"])
        axes[1].plot(x, 100 * curve, color="0.3", marker="o", ms=2.5, lw=1)
        axes[1].axhline(100 * info["threshold"], color="0.6", ls=":", lw=1)
        axes[1].axvline(n + 0.5, color="tab:red", lw=1)
        axes[1].set_xlabel("principal component")
        if info["basis"] == "signal":
            axes[1].set_ylabel("cumulative real structure (%)")
            axes[1].set_title(f"{n} PCs capture {100 * curve[n - 1]:.0f}% of real structure")
        else:
            axes[1].set_ylabel("cumulative variation (%)")
            axes[1].set_title(f"{100 * info['cum_at_n']:.0f}% of variation in the first {n}")
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "pca_variance.png"), dpi=300)
        plt.close(fig)
        files.append("pca_variance.png")
        names, r = info["cov_names"], np.asarray(info["cov_matrix"])
        rows, vals = list(names), np.abs(r) if r.size else np.zeros((0, n))
        if info["has_batch"]:
            rows.append("batch")
            vals = np.vstack([vals, np.asarray(info["batch_r2"])[None, :]]) if vals.size else np.asarray(info["batch_r2"])[None, :]
        k = min(n, 20)
        if rows:
            fig, ax = plt.subplots(figsize=(0.32 * k + 2.2, 0.5 * len(rows) + 1.3))
            im = ax.imshow(vals[:, :k], vmin=0, vmax=1, cmap="viridis", aspect="auto")
            for i in range(vals.shape[0]):
                for j in range(k):
                    ax.text(j, i, f"{vals[i, j]:.1f}".lstrip("0") if vals[i, j] < 1 else "1", ha="center", va="center", fontsize=6.5,
                            color="white" if vals[i, j] < 0.6 else "black")
            ax.set_xticks(range(k))
            ax.set_xticklabels(range(1, k + 1))
            ax.set_yticks(range(len(rows)))
            ax.set_yticklabels(rows)
            ax.set_xlabel("principal component")
            ax.set_title(f"{len(info['tech_pcs'])} of {n} components track a technical measurement")
            fig.colorbar(im, ax=ax, label="link strength", shrink=0.8)
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, "pca_covariates.png"), dpi=300)
            plt.close(fig)
            files.append("pca_covariates.png")

        def paint(ax, xy, title, xl, yl):
            if lab is not None:
                for i, c in enumerate(cats):
                    m = lab == c
                    ax.scatter(xy[m, 0], xy[m, 1], s=pt, color=cmap(i), rasterized=True, linewidths=0)
                ax.legend(handles=[Line2D([], [], marker="o", ls="", color=cmap(i), label=c) for i, c in enumerate(cats)], loc="center left",
                          bbox_to_anchor=(1.0, 0.5), markerscale=1.2)
            elif depth is not None:
                s_ = ax.scatter(xy[:, 0], xy[:, 1], s=pt, c=depth, cmap="viridis", rasterized=True, linewidths=0)
                ax.figure.colorbar(s_, ax=ax, label="log total signal per cell", shrink=0.7)
            else:
                ax.scatter(xy[:, 0], xy[:, 1], s=pt, color="0.4", rasterized=True, linewidths=0)
            ax.set_title(title)
            dr_axes_arrows(ax, xl, yl)

        two = info["corrected"]
        fig, axes = plt.subplots(1, 2 if two else 1, figsize=(7.6 if two else 4.2, 3.4))
        axes = np.atleast_1d(axes)
        paint(axes[0], adata.obsm["X_pca"][:, :2], "before batch correction" if two else "first two components", "PC1", "PC2")
        if two:
            paint(axes[1], adata.obsm["X_pca_harmony"][:, :2], "after batch correction", "PC1", "PC2")
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "pca_scatter.png"), dpi=300)
        plt.close(fig)
        files.append("pca_scatter.png")
        have_both = lab is not None and depth is not None
        fig, axes = plt.subplots(1, 2 if have_both else 1, figsize=(7.6 if have_both else 4.2, 3.4))
        axes = np.atleast_1d(axes)
        paint(axes[0], adata.obsm["X_umap"], "cells coloured by batch" if lab is not None else "cells coloured by total signal", "UMAP1", "UMAP2")
        if have_both:
            s_ = axes[1].scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1], s=pt, c=depth, cmap="viridis", rasterized=True, linewidths=0)
            fig.colorbar(s_, ax=axes[1], label="log total signal per cell", shrink=0.7)
            axes[1].set_title("cells coloured by total signal")
            dr_axes_arrows(axes[1], "UMAP1", "UMAP2")
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "umap.png"), dpi=300)
        plt.close(fig)
        files.append("umap.png")
    return files

def dr_run(adata, batch_key=None, correct_batch="auto", n_neighbors=15, max_pcs=50, min_pcs=5, n_perm=10, seed=0, technical_cols=None,
           tech_r_limit=0.5, tech_var_limit=0.3, batch_share_limit=0.2, fig_dir=None, variance_threshold=0.9, variance_basis="signal"):
    """Scale the selected genes, run PCA, choose the number of PCs from the data, examine the PCs, optionally correct batch (Harmony),
    build the neighbor graph and UMAP, draw figures with plain-language legends. Returns (adata, records). Thresholds are defaults.
    PC count: the smallest number of PCs capturing variance_threshold of the variation that rises above shuffled data (basis 'signal'),
    or of all variation (basis 'total'). The shuffled-null count and the elbow are reported as cross-checks."""
    import numpy as np
    import scanpy as sc
    from scipy.sparse.csgraph import connected_components
    if "highly_variable" not in adata.var:
        return adata, [dr_record("PCA on selected genes", "FAIL", "no highly_variable column in var: run checkpoint 3 first")]
    cols = list(technical_cols) if technical_cols is not None else ["n_genes", "total", "pct_mt"]
    X, hv = dr_scale(adata)
    n_cells, n_sel = X.shape
    nmax = int(min(max_pcs, n_cells - 1, n_sel - 1))
    U, S, Vt = dr_svd(X, nmax, seed)
    eig = S ** 2 / (n_cells - 1)
    ratio = eig / float(X.var(axis=0, ddof=1).sum())
    pcs = (U * S).astype(np.float32)
    recs = []
    if not np.isfinite(pcs).all():
        return adata, [dr_record("PCA on selected genes", "FAIL", "non-finite values in the PCA result")]
    recs.append(dr_record("PCA on selected genes", "PASS", f"{n_sel} selected genes, {n_cells} cells, {nmax} components computed"))
    n_null, n_elbow, null_max = dr_choose_pcs(X, eig, n_perm, seed)
    del X
    n_thr, reached, cum = dr_pcs_for_threshold(eig, null_max, ratio, variance_threshold, variance_basis)
    n_pcs = int(min(max(min_pcs, n_thr), nmax))
    low = n_null < min_pcs
    what = "real structure (variation above shuffled data)" if variance_basis == "signal" else "total variation"
    recs.append(dr_record("number of PCs chosen from the data", "WARN" if low else "INFO",
                          f"{n_pcs} PCs used: the smallest number capturing {100 * variance_threshold:.0f}% of the {what}; {n_null} PCs beat shuffled data ({n_perm} shuffles), "
                          f"elbow at {n_elbow}" + (f"; fewer than {min_pcs} stand out, minimum used" if low else "")))
    recs.append(dr_record("variance threshold reached", "PASS" if reached else "WARN",
                          f"{100 * float(cum[n_pcs - 1]):.0f}% of the {what} captured by {n_pcs} PCs (target {100 * variance_threshold:.0f}%)"
                          + ("" if reached else f"; not reached within {nmax} PCs, raise max_pcs")))
    beyond = max(0, n_pcs - max(n_null, min_pcs))
    recs.append(dr_record("chosen PCs beat shuffled data", "PASS" if beyond == 0 else "WARN",
                          f"{n_pcs - beyond} of {n_pcs} PCs rise above shuffled data" + ("" if beyond == 0 else f"; {beyond} do not, consider the 'signal' basis")))
    load = np.zeros((adata.n_vars, nmax), dtype=np.float32)
    load[hv] = Vt.T
    adata.obsm["X_pca"], adata.varm["PCs"] = pcs, load
    adata.uns["pca"] = {"variance": eig.astype(np.float32), "variance_ratio": ratio.astype(np.float32)}
    names, r = dr_pc_corr(pcs[:, :n_pcs], adata.obs, cols)
    strong = (np.abs(r) > tech_r_limit).any(axis=0) if r.size else np.zeros(n_pcs, dtype=bool)
    share = float(eig[:n_pcs][strong].sum() / eig[:n_pcs].sum())
    tech_pcs = [int(i) + 1 for i in np.where(strong)[0]]
    if names:
        which = {nm: [int(i) + 1 for i in np.where(np.abs(r[k]) > tech_r_limit)[0]] for k, nm in enumerate(names)}
        recs.append(dr_record("PCs driven by technical covariates", "PASS" if share <= tech_var_limit else "WARN",
                              f"{len(tech_pcs)} of {n_pcs} PCs ({100 * share:.0f}% of retained variation) correlate above {tech_r_limit} with a measurement"
                              + ("; " + "; ".join(f"{k}: PCs {v}" for k, v in which.items() if v) if any(which.values()) else "") + f"; limit {100 * tech_var_limit:.0f}%"))
    else:
        recs.append(dr_record("PCs driven by technical covariates", "INFO", "no usable technical covariates in obs (n_genes, total, pct_mt)"))
    labels = adata.obs[batch_key].astype(str).to_numpy() if batch_key is not None and batch_key in adata.obs else None
    has_batch = labels is not None and len(set(labels)) >= 2
    share_b = share_a = mix_b = mix_a = float("nan")
    corrected, rep, batch_r2 = False, "X_pca", np.zeros(n_pcs)
    if has_batch:
        batch_r2 = dr_batch_r2(pcs[:, :n_pcs], labels)
        share_b = dr_batch_share(pcs[:, :n_pcs], labels)
        mix_b = dr_mixing(pcs[:, :n_pcs], labels, seed=seed)
        want = correct_batch is True or (correct_batch == "auto" and share_b > batch_share_limit)
        if want:
            try:
                Z = dr_harmony(pcs[:, :n_pcs], labels, batch_key, seed)
                adata.obsm["X_pca_harmony"] = Z
                share_a, mix_a = dr_batch_share(Z, labels), dr_mixing(Z, labels, seed=seed)
                corrected, rep = True, "X_pca_harmony"
            except ImportError:
                recs.append(dr_record("batch correction effect", "FAIL" if correct_batch is True else "WARN", "harmonypy is not installed; correction skipped"))
            except Exception as e:
                recs.append(dr_record("batch correction effect", "WARN", f"correction failed: {type(e).__name__}: {str(e)[:100]}"))
        stt = "PASS" if share_b <= batch_share_limit else ("INFO" if corrected else "WARN")
        recs.append(dr_record("batch structure in PCs", stt, f"batch '{batch_key}' ({len(set(labels))} groups) explains {100 * share_b:.0f}% of retained PC variation; "
                              f"neighbor mixing {mix_b:.2f}; limit {100 * batch_share_limit:.0f}%" + ("; correction applied" if corrected else "")))
        if corrected:
            good = share_a < share_b and mix_a >= mix_b and share_a <= batch_share_limit
            recs.append(dr_record("batch correction effect", "PASS" if good else "WARN",
                                  f"batch variation {100 * share_b:.0f}% to {100 * share_a:.0f}%, neighbor mixing {mix_b:.2f} to {mix_a:.2f}. Correction cannot tell batch "
                                  "from real biology: check that known groups stay separate"))
    elif correct_batch is True:
        recs.append(dr_record("batch correction effect", "FAIL", "correct_batch=True but no usable batch column with at least 2 groups"))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep, random_state=seed)
    ncomp, comp = connected_components(adata.obsp["connectivities"], directed=False)
    sizes = np.bincount(comp)
    stranded = int(sizes[sizes < 10].sum())
    recs.append(dr_record("neighbor graph fragments", "PASS" if stranded <= 0.01 * adata.n_obs else "WARN",
                          f"{ncomp} connected piece(s), largest holds {100 * sizes.max() / adata.n_obs:.1f}% of cells; {stranded} cells ({100 * stranded / adata.n_obs:.1f}%) "
                          "sit in pieces under 10 cells. Well separated cell groups can form their own large pieces, which is not a fault"))
    sc.tl.umap(adata, random_state=seed)
    ok = bool(np.isfinite(adata.obsm["X_umap"]).all())
    recs.append(dr_record("UMAP finite", "PASS" if ok else "FAIL", "all coordinates finite" if ok else "non-finite UMAP coordinates"))
    depth_cols = [c for c in ("total", "n_genes") if c in adata.obs]
    lab_ok = has_batch and 2 <= len(set(labels)) <= 12
    info = {"n_cells": int(n_cells), "n_genes_selected": int(n_sel), "n_pcs": n_pcs, "n_null": int(n_null), "n_elbow": int(n_elbow), "nmax": nmax,
            "min_pcs": int(min_pcs), "var_ratio": [float(v) for v in ratio], "null_ratio": [float(v) for v in null_max / float(eig.sum() / ratio.sum())],
            "cum_at_n": float(ratio[:n_pcs].sum()), "threshold": float(variance_threshold), "basis": variance_basis, "reached": bool(reached), "cum_curve": [float(v) for v in cum], "cov_names": names, "cov_matrix": [[float(v) for v in row[:n_pcs]] for row in r], "tech_pcs": tech_pcs,
            "tech_var_share": share, "tech_r_limit": tech_r_limit, "has_batch": bool(has_batch), "batch_key": batch_key, "batch_r2": [float(v) for v in batch_r2],
            "share_before": share_b, "share_after": share_a, "mix_before": mix_b, "mix_after": mix_a, "corrected": corrected, "rep": rep,
            "color_kind": "batch" if lab_ok else ("depth" if depth_cols else "none"),
            "color_what": (f"batch ({batch_key})" if lab_ok else ("total signal per cell on a log scale" if depth_cols else "nothing (single colour)"))}
    adata.uns["dimred"] = {k: v for k, v in info.items() if not isinstance(v, (list, dict))}
    if fig_dir is not None:
        import os
        files = dr_figures(adata, info, fig_dir, batch_key)
        legends = dr_legends(info)
        os.makedirs(fig_dir, exist_ok=True)
        with open(os.path.join(fig_dir, "figure_legends.md"), "w", encoding="utf-8") as f:
            f.write("# Figure legends\n\n" + "\n\n".join(f"**{k.replace('_', ' ')}** ({k}.png): {v}" for k, v in legends.items() if k + ".png" in files) + "\n")
    return adata, recs

def write_dr_report(records, prefix="dimred_report", scope=None):
    """Write <prefix>.json and <prefix>.md; return the overall status (worst of PASS/WARN/FAIL).
    The report states WHY the verdict was reached, lists every check that ran with what it tests,
    and lists catalogued checks that did not run."""
    cat = dr_catalog()
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
    lines = [f"# Dimensionality reduction report: {overall}", "", f"**Why:** {why}", "",
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
