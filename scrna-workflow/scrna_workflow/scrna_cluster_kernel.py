import json

def cl_record(check, status, detail, group=None):
    """One result row. status is PASS, WARN, FAIL or INFO."""
    r = {"check": check, "status": status, "detail": detail}
    if group is not None:
        r["group"] = group
    return r

def cl_catalog():
    """Checks this checkpoint can run: name -> what it tests."""
    return {
        "representation used": "clusters are found on the same (batch-corrected if needed) PCs chosen at the previous checkpoint",
        "resolution chosen from the data": "the resolution is the finest one whose clusters stay stable when cells are resampled and whose neighboring clusters still differ by enough genes",
        "chosen resolution is stable": "clusters at the chosen resolution are reproduced when 20% of cells are left out",
        "cluster count": "the number of clusters is plausible for the number of cells",
        "no tiny clusters": "no cluster is too small to trust",
        "no dominant cluster": "no single cluster holds most of the cells",
        "per-cluster stability": "every cluster is recovered when cells are resampled",
        "neighboring clusters distinguishable": "clusters next to each other differ by enough genes to be kept apart",
        "cluster markers": "every cluster has marker genes that define it",
        "batch mixing within clusters": "no cluster is made almost entirely of one batch",
    }

def cl_rank_labels(lab):
    """Relabel clusters 0, 1, 2... from the largest to the smallest."""
    import numpy as np
    import pandas as pd
    order = pd.Series(lab).value_counts().index
    m = {old: str(i) for i, old in enumerate(order)}
    return np.array([m[x] for x in lab])

def cl_graph(emb, n_neighbors=15, seed=0):
    """Neighbor graph on an embedding, held in a light AnnData (no expression matrix)."""
    import numpy as np
    import scipy.sparse as sp
    import anndata as ad
    import scanpy as sc
    a = ad.AnnData(sp.csr_matrix((emb.shape[0], 1), dtype=np.float32))
    a.obsm["X_rep"] = np.asarray(emb, dtype=np.float32)
    sc.pp.neighbors(a, n_neighbors=n_neighbors, use_rep="X_rep", random_state=seed)
    return a

def cl_leiden(a, resolution, seed=0):
    """Leiden (modularity) on the graph held in `a`, run directly in igraph with Python's seeded random generator.
    (scanpy's wrapper routes igraph through numpy's generator, which raises thousands of silent errors on Windows and is ~1000x slower.)"""
    import random
    import numpy as np
    import igraph
    from scanpy._utils import get_igraph_from_adjacency
    G = get_igraph_from_adjacency(a.obsp["connectivities"], directed=False)
    igraph.set_random_number_generator(random)
    random.seed(int(seed))
    part = G.community_leiden(objective_function="modularity", weights=G.es["weight"], resolution=float(resolution), n_iterations=2)
    return cl_rank_labels(np.array([str(x) for x in part.membership]))

def cl_best_jaccard(ref, sub, names):
    """For each reference cluster, the Jaccard overlap with its best-matching cluster in a resampled clustering (NaN if absent)."""
    import numpy as np
    import pandas as pd
    ri = pd.Categorical(ref, categories=names).codes
    sn = np.unique(sub)
    si = pd.Categorical(sub, categories=sn).codes
    M = np.zeros((len(names), len(sn)))
    np.add.at(M, (ri, si), 1)
    union = M.sum(1, keepdims=True) + M.sum(0, keepdims=True) - M
    J = np.where(union > 0, M / np.maximum(union, 1), 0.0)
    out = J.max(axis=1)
    out[M.sum(1) == 0] = np.nan
    return out

def cl_stability(emb, full, n_neighbors=15, n_boot=5, frac=0.8, seed=0):
    """Resample cells n_boot times (frac of cells), rebuild the graph, recluster at every resolution and compare with the full-data
    clusters on the shared cells. Returns {resolution: {'ari': [...], 'jac': array(n_boot x n_clusters)}}."""
    import numpy as np
    from sklearn.metrics import adjusted_rand_score
    rng = np.random.default_rng(seed)
    n = emb.shape[0]
    out = {r: {"ari": [], "jac": []} for r in full}
    for b in range(n_boot):
        idx = np.sort(rng.choice(n, size=int(frac * n), replace=False))
        g = cl_graph(emb[idx], n_neighbors, seed + 1 + b)
        for r, lab in full.items():
            sub = cl_leiden(g, r, seed + 1 + b)
            ref = lab[idx]
            names = [str(i) for i in range(len(np.unique(lab)))]
            out[r]["ari"].append(float(adjusted_rand_score(ref, sub)))
            out[r]["jac"].append(cl_best_jaccard(ref, sub, names))
    for r in out:
        out[r]["jac"] = np.array(out[r]["jac"])
    return out

def cl_silhouette(emb, labels, seed=0, max_cells=5000):
    import numpy as np
    from sklearn.metrics import silhouette_score
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(labels), size=min(max_cells, len(labels)), replace=False)
    if len(np.unique(labels[idx])) < 2:
        return float("nan")
    return float(silhouette_score(emb[idx], labels[idx]))

def cl_pick_resolution(table, min_stability=0.8, min_size=10):
    """Finest resolution that is stable (median ARI at least min_stability), has at least 2 clusters, no cluster below min_size and,
    where it was tested (`distinct` key), neighboring clusters that differ from one another. If none qualifies, the most stable
    resolution with at least 2 clusters. Returns (resolution, qualified)."""
    ok = [t for t in table if t["n_clusters"] >= 2 and t["stability"] >= min_stability and t["min_size"] >= min_size and t.get("distinct", True) is not False]
    if ok:
        return max(ok, key=lambda t: t["resolution"])["resolution"], True
    cand = [t for t in table if t["n_clusters"] >= 2] or table
    return max(cand, key=lambda t: t["stability"])["resolution"], False

def cl_nearest_pairs(emb, labels):
    """Each cluster paired with its nearest other cluster (centroid distance). Returns sorted unique pairs."""
    import numpy as np
    names = sorted(np.unique(labels), key=int)
    cent = np.stack([emb[labels == c].mean(axis=0) for c in names])
    d = np.linalg.norm(cent[:, None] - cent[None], axis=2)
    np.fill_diagonal(d, np.inf)
    return sorted({tuple(sorted((names[i], names[int(d[i].argmin())]), key=int)) for i in range(len(names))}, key=lambda p: (int(p[0]), int(p[1])))

def cl_pair_de(adata, labels, pairs, gene_mask=None, max_cells=500, log2fc=0.5, padj=0.01, seed=0):
    """Number of genes separating each cluster pair (Wilcoxon, |log2 fold change| and adjusted p thresholds). None if a cluster has under 3 cells."""
    import numpy as np
    import pandas as pd
    import anndata as ad
    import scanpy as sc
    rng = np.random.default_rng(seed)
    X = adata.X if gene_mask is None else adata.X[:, np.where(gene_mask)[0]]
    out = {}
    for a, b in pairs:
        parts = []
        for c in (a, b):
            ii = np.where(labels == c)[0]
            parts.append(rng.choice(ii, size=min(max_cells, len(ii)), replace=False))
        if min(len(parts[0]), len(parts[1])) < 3:
            out[(a, b)] = None
            continue
        idx = np.concatenate(parts)
        sub = ad.AnnData(X[idx])
        sub.obs["g"] = pd.Categorical(labels[idx])
        sc.tl.rank_genes_groups(sub, "g", groups=[a], reference=b, method="wilcoxon", n_genes=sub.n_vars, use_raw=False)
        df = sc.get.rank_genes_groups_df(sub, group=a)
        out[(a, b)] = int(((df["logfoldchanges"].abs() >= log2fc) & (df["pvals_adj"] < padj)).sum())
    return out

def cl_markers(adata, labels, min_pct=0.25, log2fc=0.5, padj=0.01, max_cells=500, seed=0):
    """Marker genes of each cluster against all other cells (Wilcoxon). Only genes detected in at least min_pct of one cluster's cells are tested,
    on up to max_cells cells per cluster. A marker has log2 fold change and adjusted p beyond the thresholds and detection in min_pct of its cluster.
    Returns (table, detection fraction matrix clusters x genes, cluster names)."""
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import anndata as ad
    import scanpy as sc
    X = sp.csr_matrix(adata.X)
    names = sorted(np.unique(labels), key=int)
    codes = pd.Categorical(labels, categories=names).codes
    sizes = np.bincount(codes).astype(float)
    M = sp.csr_matrix((np.ones(len(labels)), (codes, np.arange(len(labels)))), shape=(len(names), len(labels)))
    Xb = X.copy()
    Xb.data = np.ones_like(Xb.data)
    det = np.asarray((M @ Xb).todense()) / sizes[:, None]
    valid = [i for i in range(len(names)) if sizes[i] >= 3]
    keep = np.where(det[valid].max(axis=0) >= min_pct)[0]
    rng = np.random.default_rng(seed)
    idx = np.concatenate([rng.choice(np.where(codes == i)[0], size=int(min(max_cells, sizes[i])), replace=False) for i in valid])
    sub = ad.AnnData(X[idx][:, keep])
    sub.var_names = adata.var_names[keep]
    sub.obs["g"] = pd.Categorical([names[c] for c in codes[idx]], categories=[names[i] for i in valid])
    sc.tl.rank_genes_groups(sub, "g", method="wilcoxon", n_genes=sub.n_vars, use_raw=False)
    pos = pd.Index(sub.var_names)
    rows = []
    for i in valid:
        df = sc.get.rank_genes_groups_df(sub, group=names[i])
        df["cluster"] = names[i]
        df["pct_in_cluster"] = det[i, keep][pos.get_indexer(df["names"])]
        rows.append(df)
    tab = pd.concat(rows, ignore_index=True).rename(columns={"names": "gene", "logfoldchanges": "log2fc"})
    tab["is_marker"] = (tab["log2fc"] >= log2fc) & (tab["pvals_adj"] < padj) & (tab["pct_in_cluster"] >= min_pct)
    return tab, det, names

def cl_batch_mixing(labels, batch, dominant=0.9, min_cells=10):
    """Clusters (of at least min_cells cells) in which one batch supplies at least `dominant` of the cells although that batch is not
    the bulk of the dataset. Returns a list of (cluster, batch, share)."""
    import numpy as np
    import pandas as pd
    t = pd.crosstab(pd.Series(labels, name="c"), pd.Series(np.asarray(batch).astype(str), name="b"))
    overall = t.sum(0) / t.values.sum()
    out = []
    for c, row in t.iterrows():
        n = row.sum()
        if n < min_cells:
            continue
        top = row.idxmax()
        if row[top] / n >= dominant and overall[top] < dominant:
            out.append((str(c), str(top), float(row[top] / n)))
    return sorted(out, key=lambda x: int(x[0]))

def cl_legends(info):
    """Plain-language figure legends built from the computed numbers. No cell-type names are given."""
    k, r = info["n_clusters"], info["resolution"]
    return {
        "cluster_resolution": (f"Each point is one clustering resolution, the setting that controls how finely cells are split (higher gives more, smaller clusters). "
                               f"Top: number of clusters. Middle: stability, meaning how similar the clusters stay when 20% of cells are left out and the analysis is repeated "
                               f"({info['n_boot']} repeats; 1 means identical; line = median); the dotted line is the {info['min_stability']:.1f} target. Bottom: silhouette, how well a "
                               f"cell fits its own cluster compared with the nearest other one (higher is better separated). The red line marks the resolution used ({r}), the finest "
                               f"one that stays stable with neighboring clusters that differ" + ("." if info["qualified"] else ", but none reached the target, so the most stable one is used.")),
        "cluster_umap": (f"Each dot is a cell and each colour a cluster ({k} clusters, resolution {r}); numbers mark cluster centres (very small clusters are drawn with larger outlined dots). Clusters were found on "
                         f"{info['rep_text']}, not on this picture. UMAP is for viewing: distances between separated groups and group sizes are not quantitative."),
        "cluster_sizes": (f"Left: number of cells in each cluster (red = fewer than {info['min_cells']} cells, too few to trust). Right: how consistently each cluster is recovered "
                          f"when cells are resampled (1 = the same cells always end up together; below {info['min_jaccard']:.1f} the cluster is unreliable)."),
        "cluster_markers": (f"Each row is a cluster and each column a gene that is higher in that cluster than in the other cells (up to {info['top_m']} per cluster). Dot size is the share "
                            f"of the cluster's cells in which the gene is detected; colour is average expression scaled for each gene (dark = lowest, bright = highest across clusters). "
                            f"Markers show how clusters differ; they are not cell-type names."),
    }

def cl_figures(adata, info, fig_dir, tab, det, names, resolutions, table):
    """Four presentation figures: resolution scan, UMAP by cluster, cluster sizes and stability, marker dot plot. Returns file names."""
    import os
    import numpy as np
    import scipy.sparse as sp
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rc = {"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 8,
          "axes.spines.top": False, "axes.spines.right": False, "xtick.direction": "out", "ytick.direction": "out",
          "axes.titleweight": "normal", "axes.titlelocation": "left", "legend.frameon": False}
    os.makedirs(fig_dir, exist_ok=True)
    files = []
    lab = adata.obs["cluster"].astype(str).to_numpy()
    k = len(names)
    sizes = np.array([(lab == c).sum() for c in names])
    with matplotlib.rc_context(rc):
        fig, axes = plt.subplots(3, 1, figsize=(5.2, 5.0), sharex=True)
        x = [t["resolution"] for t in table]
        axes[0].plot(x, [t["n_clusters"] for t in table], color="0.3", marker="o", ms=3)
        axes[0].set_ylabel("clusters")
        axes[1].plot(x, [t["stability"] for t in table], color="0.3", marker="o", ms=3)
        axes[1].fill_between(x, [t["stability_min"] for t in table], [t["stability"] for t in table], color="0.8", lw=0)
        axes[1].axhline(info["min_stability"], color="0.6", ls=":", lw=1)
        axes[1].set_ylabel("stability (ARI)")
        axes[1].set_ylim(0, 1.05)
        axes[2].plot(x, [t["silhouette"] for t in table], color="0.3", marker="o", ms=3)
        axes[2].set_ylabel("silhouette")
        axes[2].set_xlabel("clustering resolution")
        for ax in axes:
            ax.axvline(info["resolution"], color="tab:red", lw=1)
        axes[0].set_title(f"{k} clusters at resolution {info['resolution']}, the finest stable setting" if info["qualified"] else f"{k} clusters at resolution {info['resolution']}, the most stable setting")
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "cluster_resolution.png"), dpi=300)
        plt.close(fig)
        files.append("cluster_resolution.png")
        cmap = plt.get_cmap("tab20")
        um = adata.obsm["X_umap"]
        codes = np.array([names.index(c) for c in lab])
        fig, ax = plt.subplots(figsize=(4.6, 4.0))
        ax.scatter(um[:, 0], um[:, 1], s=3 if len(lab) < 20000 else 1.5, c=[cmap(i % 20) for i in codes], rasterized=True, linewidths=0)
        span = np.ptp(um, axis=0)
        for i, c in enumerate(names):
            sel = lab == c
            m = np.median(um[sel], axis=0)
            if sel.sum() < max(30, 0.01 * len(lab)):
                # a very small cluster is hidden by its own label: enlarge its dots and put the number beside them
                ax.scatter(um[sel, 0], um[sel, 1], s=14, c=[cmap(i % 20)], edgecolors="black", linewidths=0.4, zorder=3)
                ax.text(m[0], m[1] + 0.05 * span[1], c, fontsize=8, ha="center", va="bottom", weight="bold")
            else:
                ax.text(m[0], m[1], c, fontsize=8, ha="center", va="center", weight="bold", bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.7))
        ax.set_title(f"{k} clusters, labelled by number")
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        arrow = dict(arrowstyle="->", lw=0.8)
        ax.annotate("", xy=(0.14, 0.0), xytext=(0.0, 0.0), xycoords="axes fraction", arrowprops=arrow)
        ax.annotate("", xy=(0.0, 0.14), xytext=(0.0, 0.0), xycoords="axes fraction", arrowprops=arrow)
        ax.text(0.07, -0.03, "UMAP1", transform=ax.transAxes, ha="center", va="top", fontsize=7)
        ax.text(-0.03, 0.07, "UMAP2", transform=ax.transAxes, ha="right", va="center", rotation=90, fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "cluster_umap.png"), dpi=300)
        plt.close(fig)
        files.append("cluster_umap.png")
        jac = np.asarray(info["cluster_jaccard"])
        fig, axes = plt.subplots(1, 2, figsize=(7.0, 0.22 * k + 1.6), sharey=True)
        y = np.arange(k)
        axes[0].barh(y, sizes, color=["tab:red" if s < info["min_cells"] else "0.55" for s in sizes])
        axes[0].set_xlabel("cells")
        axes[0].set_yticks(y)
        axes[0].set_yticklabels(names)
        axes[0].set_ylabel("cluster")
        axes[0].invert_yaxis()
        axes[0].set_title("cluster sizes")
        axes[1].barh(y, np.nan_to_num(jac), color="0.55")
        axes[1].axvline(info["min_jaccard"], color="0.6", ls=":", lw=1)
        axes[1].set_xlim(0, 1.02)
        axes[1].set_xlabel("recovery under resampling")
        axes[1].set_title("how reliably each cluster is found")
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "cluster_sizes.png"), dpi=300)
        plt.close(fig)
        files.append("cluster_sizes.png")
        top = tab[tab["is_marker"]].sort_values(["cluster", "scores"], ascending=[True, False]).groupby("cluster", sort=False).head(info["top_m"])
        genes = list(dict.fromkeys(top.sort_values("cluster", key=lambda s: s.astype(int))["gene"]))[:60]
        if genes:
            X = sp.csr_matrix(adata.X)
            gi = adata.var_names.get_indexer(genes)
            M = sp.csr_matrix((np.ones(len(lab)), (codes, np.arange(len(lab)))), shape=(k, len(lab)))
            mean = np.asarray((M @ X[:, gi]).todense()) / sizes[:, None]
            frac = det[:, gi]
            lo, hi = mean.min(axis=0), mean.max(axis=0)
            col = (mean - lo) / np.maximum(hi - lo, 1e-9)
            fig, ax = plt.subplots(figsize=(0.24 * len(genes) + 1.8, 0.26 * k + 1.5))
            gx, cy = np.meshgrid(np.arange(len(genes)), np.arange(k))
            sc_ = ax.scatter(gx.ravel(), cy.ravel(), s=4 + 70 * frac.ravel(), c=col.ravel(), cmap="viridis", vmin=0, vmax=1, linewidths=0)
            ax.set_xticks(range(len(genes)))
            ax.set_xticklabels(genes, rotation=90, fontstyle="italic", fontsize=6.5)
            ax.set_yticks(range(k))
            ax.set_yticklabels(names)
            ax.set_ylabel("cluster")
            ax.set_ylim(k - 0.5, -0.5)
            ax.margins(0.03)
            ax.set_title("genes that distinguish each cluster")
            fig.colorbar(sc_, ax=ax, label="scaled mean expression", shrink=0.6)
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, "cluster_markers.png"), dpi=300)
            plt.close(fig)
            files.append("cluster_markers.png")
    return files

def cl_run(adata, resolution=None, resolutions=None, n_neighbors=15, n_boot=5, frac=0.8, seed=0, batch_key=None, min_stability=0.8,
           min_cluster_cells=10, min_cluster_frac=0.005, max_cluster_frac=0.5, min_cluster_jaccard=0.6, min_markers=5, marker_log2fc=0.5,
           marker_padj=0.01, marker_min_pct=0.25, min_pair_de=10, max_cells=500, top_m=3, fig_dir=None):
    """Cluster on the representation chosen at checkpoint 4, scan resolutions, choose the finest stable one (or use `resolution`),
    check sizes, per-cluster stability, distinguishability of neighboring clusters, markers and batch mixing, and draw figures with legends.
    Returns (adata, records). Thresholds are defaults, not literature values."""
    import os
    import numpy as np
    import pandas as pd
    res = sorted(set(list(resolutions) if resolutions is not None else [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]) | ({float(resolution)} if resolution is not None else set()))
    d = adata.uns.get("dimred", {})
    rep, npc = d.get("rep", "X_pca"), int(d.get("n_pcs", 30))
    if rep not in adata.obsm:
        return adata, [cl_record("representation used", "FAIL", f"{rep} not found in obsm: run checkpoint 4 first")]
    npc = min(npc, adata.obsm[rep].shape[1])
    emb = np.asarray(adata.obsm[rep][:, :npc])
    rep_text = f"the first {npc} {'batch-corrected ' if rep == 'X_pca_harmony' else ''}principal components"
    recs = [cl_record("representation used", "INFO", f"{rep}, first {npc} components" + (" (batch-corrected)" if rep == "X_pca_harmony" else "") + f"; neighbor graph rebuilt with {n_neighbors} neighbors")]
    n = adata.n_obs
    g = cl_graph(emb, n_neighbors, seed)
    full = {r: cl_leiden(g, r, seed) for r in res}
    stab = cl_stability(emb, full, n_neighbors, n_boot, frac, seed)
    min_size_ok = max(min_cluster_cells, int(np.ceil(min_cluster_frac * n)))
    table = []
    for r in res:
        sz = pd.Series(full[r]).value_counts()
        table.append({"resolution": r, "n_clusters": int(len(sz)), "min_size": int(sz.min()), "max_frac": float(sz.max() / n), "silhouette": cl_silhouette(emb, full[r], seed),
                      "stability": float(np.median(stab[r]["ari"])), "stability_min": float(np.min(stab[r]["ari"]))})
    hv = adata.var["highly_variable"].to_numpy() if "highly_variable" in adata.var else None
    if resolution is None:
        # from the finest stable setting downwards, keep the first one whose neighboring clusters differ by enough genes
        for t in sorted(table, key=lambda t: -t["resolution"]):
            if t["n_clusters"] >= 2 and t["stability"] >= min_stability and t["min_size"] >= min_size_ok:
                lab_t = full[t["resolution"]]
                de_t = [v for v in cl_pair_de(adata, lab_t, cl_nearest_pairs(emb, lab_t), hv, max_cells, marker_log2fc, marker_padj, seed).values() if v is not None]
                t["min_pair_genes"] = int(min(de_t)) if de_t else None
                t["distinct"] = (not de_t) or min(de_t) >= min_pair_de
                if t["distinct"]:
                    break
        chosen, qualified = cl_pick_resolution(table, min_stability, min_size_ok)
    else:
        chosen = float(resolution)
        qualified = next(t for t in table if t["resolution"] == chosen)["stability"] >= min_stability
    row = next(t for t in table if t["resolution"] == chosen)
    labels = full[chosen]
    names = [str(i) for i in range(row["n_clusters"])]
    k = len(names)
    how = "set by the user" if resolution is not None else ("the finest resolution that is stable and whose neighboring clusters differ" if qualified else "no resolution reached the stability target; the most stable one")
    recs.append(cl_record("resolution chosen from the data", "INFO" if (qualified or resolution is not None) else "WARN",
                          f"resolution {chosen} ({how}): {k} clusters, stability {row['stability']:.2f}, silhouette {row['silhouette']:.2f}; scanned {res[0]} to {res[-1]}"))
    recs.append(cl_record("chosen resolution is stable", "PASS" if row["stability"] >= min_stability else "WARN",
                          f"median adjusted Rand index {row['stability']:.2f} (lowest {row['stability_min']:.2f}) over {n_boot} resamples keeping {int(100 * frac)}% of cells; target {min_stability}"))
    max_k = max(2, min(60, n // 20))
    recs.append(cl_record("cluster count", "PASS" if 2 <= k <= max_k else "WARN", f"{k} clusters for {n} cells (plausible range 2 to {max_k})"))
    sizes = pd.Series(labels).value_counts()
    tiny = [f"{c} ({sizes[c]} cells)" for c in sorted(sizes.index, key=int) if sizes[c] < min_size_ok]
    recs.append(cl_record("no tiny clusters", "PASS" if not tiny else "WARN", "smallest cluster has " + f"{int(sizes.min())} cells; limit {min_size_ok}" if not tiny else f"below {min_size_ok} cells: " + ", ".join(tiny)))
    big = float(sizes.max() / n)
    recs.append(cl_record("no dominant cluster", "PASS" if big <= max_cluster_frac else "WARN", f"largest cluster {sizes.idxmax()} holds {100 * big:.0f}% of cells; limit {100 * max_cluster_frac:.0f}%"))
    jac = np.nanmean(stab[chosen]["jac"], axis=0)
    weak = [f"{names[i]} ({jac[i]:.2f})" for i in range(k) if not jac[i] >= min_cluster_jaccard]
    recs.append(cl_record("per-cluster stability", "PASS" if not weak else "WARN", f"lowest recovery {np.nanmin(jac):.2f}; limit {min_cluster_jaccard}" if not weak else f"recovered under {min_cluster_jaccard}: " + ", ".join(weak)))
    if k >= 2:
        pairs = cl_nearest_pairs(emb, labels)
        de = cl_pair_de(adata, labels, pairs, hv, max_cells, marker_log2fc, marker_padj, seed)
        close = [f"{a}-{b} ({v} genes)" for (a, b), v in de.items() if v is not None and v < min_pair_de]
        recs.append(cl_record("neighboring clusters distinguishable", "PASS" if not close else "WARN",
                              f"each cluster's nearest neighbor differs by at least {min_pair_de} genes (|log2FC| at least {marker_log2fc}, adjusted p under {marker_padj}; selected genes)" if not close
                              else "few separating genes between: " + ", ".join(close) + "; consider merging"))
        tab, det, mnames = cl_markers(adata, labels, marker_min_pct, marker_log2fc, marker_padj, max_cells, seed)
        cnt = tab[tab["is_marker"]].groupby("cluster").size().reindex(names, fill_value=0)
        few = [f"{c} ({cnt[c]})" for c in names if cnt[c] < min_markers and c in set(tab["cluster"])]
        recs.append(cl_record("cluster markers", "PASS" if not few else "WARN",
                              (f"every cluster has at least {min_markers} markers; median {int(cnt.median())} per cluster" if not few else f"fewer than {min_markers} markers: " + ", ".join(few))
                              + f". Tested on up to {max_cells} cells per cluster and genes detected in at least {int(100 * marker_min_pct)}% of a cluster's cells"))
    else:
        tab, det = pd.DataFrame(columns=["cluster", "gene", "scores", "log2fc", "pvals_adj", "pct_in_cluster", "is_marker"]), np.zeros((1, adata.n_vars))
        recs.append(cl_record("neighboring clusters distinguishable", "INFO", "only one cluster: nothing to compare"))
        recs.append(cl_record("cluster markers", "INFO", "only one cluster: markers are not defined"))
    if batch_key is not None and batch_key in adata.obs:
        flagged = cl_batch_mixing(labels, adata.obs[batch_key].to_numpy(), min_cells=min_size_ok)
        recs.append(cl_record("batch mixing within clusters", "PASS" if not flagged else "WARN",
                              f"no cluster is 90% or more one batch ('{batch_key}')" if not flagged else "single-batch clusters: " + ", ".join(f"{c} ({int(100 * s)}% {b})" for c, b, s in flagged)))
    adata.obs["cluster"] = pd.Categorical(labels, categories=names)
    for r in res:
        adata.obs[f"leiden_r{r}"] = pd.Categorical(full[r])
    adata.uns["clustering"] = {"resolution": float(chosen), "n_clusters": int(k), "rep": rep, "n_pcs": int(npc), "stability": float(row["stability"]), "silhouette": float(row["silhouette"]),
                               "n_neighbors": int(n_neighbors), "n_boot": int(n_boot), "qualified": bool(qualified)}
    scan = pd.DataFrame(table)
    scan["min_pair_genes"] = scan["min_pair_genes"].astype(float) if "min_pair_genes" in scan else np.nan
    scan["distinct"] = scan["distinct"].map({True: "yes", False: "no"}).fillna("not tested") if "distinct" in scan else "not tested"
    adata.uns["cluster_resolution_scan"] = scan
    top50 = tab[tab["is_marker"]].sort_values(["cluster", "scores"], ascending=[True, False]).groupby("cluster", sort=False).head(50)
    adata.uns["cluster_markers"] = top50.reset_index(drop=True)
    if fig_dir is not None:
        os.makedirs(fig_dir, exist_ok=True)
        tab[tab["is_marker"]].sort_values(["cluster", "scores"], ascending=[True, False]).to_csv(os.path.join(fig_dir, "cluster_markers.csv"), index=False)
        scan.to_csv(os.path.join(fig_dir, "resolution_scan.csv"), index=False)
        info = {"n_clusters": k, "resolution": chosen, "n_boot": n_boot, "min_stability": min_stability, "qualified": bool(qualified), "rep_text": rep_text, "min_cells": min_size_ok,
                "min_jaccard": min_cluster_jaccard, "top_m": top_m, "cluster_jaccard": [float(v) for v in jac]}
        files = cl_figures(adata, info, fig_dir, tab, det, names, res, table) if "X_umap" in adata.obsm else []
        leg = cl_legends(info)
        with open(os.path.join(fig_dir, "figure_legends.md"), "w", encoding="utf-8") as f:
            f.write("# Figure legends\n\n" + "\n\n".join(f"**{a.replace('_', ' ')}** ({a}.png): {b}" for a, b in leg.items() if a + ".png" in files) + "\n")
    return adata, recs

def write_cl_report(records, prefix="clustering_report", scope=None):
    """Write <prefix>.json and <prefix>.md; return the overall status (worst of PASS/WARN/FAIL).
    The report states WHY the verdict was reached, lists every check that ran with what it tests,
    and lists catalogued checks that did not run."""
    cat = cl_catalog()
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
    lines = [f"# Clustering report: {overall}", "", f"**Why:** {why}", "",
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
