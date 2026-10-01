"""Checkpoint 6 helpers: marker-based cell type annotation with confidence tiers (functions and imports only)."""
import json

TIERS = ("confident", "probable", "ambiguous", "unassigned")

def ann_record(check, status, detail, group=None):
    r = {"check": check, "status": status, "detail": detail}
    if group:
        r["group"] = group
    return r

def ann_catalog():
    """Every check this checkpoint can run, with what it tests. Checks that did not run are listed as 'Not run'."""
    return {
        "marker panel": "each cell type in the panel has enough of its positive marker genes present in the data to be scored",
        "panel strength": "every cell type has at least 2 positive markers, so no call rests on a single gene",
        "enough usable markers": "every panel type keeps at least 3 detectable positive markers (or all it was given, if fewer), so a call never rests on one or two noisy genes",
        "sparse panel genes": "panel genes too rarely detected in this dataset (dropout) to count as evidence are set aside and listed",
        "negative markers": "genes that should be off in a cell type are defined (given, or suggested from the other types)",
        "annotation tiers": "how many clusters could be named with confidence, and how many are ambiguous or unassigned",
        "cells named": "most cells sit in clusters that could be named (confident or probable)",
        "negative marker conflicts": "no named cluster switches on genes that should be off in its assigned type",
        "data-driven marker support": "the cluster's own top marker genes (found without the panel) include genes from its assigned type's panel",
        "reference agreement": "an independent labelled reference dataset gives the same name as the marker panel",
        "same type in several clusters": "clusters sharing a name are reported, since they may be subtypes, states or over-splitting",
        "panel types not found": "cell types in the panel that no cluster was assigned to",
        "pretrained annotation tools": "automated annotation with pretrained models (for example CellTypist) was compared",
    }

def ann_split(x):
    if x is None:
        return []
    if isinstance(x, float) and x != x:
        return []
    if isinstance(x, (list, tuple, set)):
        return [str(g).strip() for g in x if str(g).strip()]
    return [g.strip() for g in str(x).replace(",", ";").split(";") if g.strip()]

def ann_read_panel(panel):
    """Panel as a DataFrame with columns cell_type, positive, negative (genes separated by ';'), optional source. Accepts a path or a DataFrame."""
    import pandas as pd
    df = pd.read_csv(panel) if isinstance(panel, str) else panel.copy()
    miss = {"cell_type", "positive"} - set(df.columns)
    if miss:
        raise ValueError(f"panel needs columns cell_type and positive; missing {sorted(miss)}")
    if "negative" not in df:
        df["negative"] = ""
    if "source" not in df:
        df["source"] = ""
    df["cell_type"] = df["cell_type"].astype(str)
    df["positive"] = [";".join(ann_split(v)) for v in df["positive"]]
    df["negative"] = [";".join(ann_split(v)) for v in df["negative"]]
    return df[["cell_type", "positive", "negative", "source"]].reset_index(drop=True)

def ann_clean_cellguide(raw, organism="Homo sapiens", max_per_type=8, max_types_sharing=2):
    """Turn CellGuide marker results into a panel table.
    raw: {requested name: {"id", "name", "canonical": [gene dicts], "computational": [gene dicts]}} as returned by the connector.
    Canonical (literature) markers come first, then data-derived ones for `organism` by marker score. Genes listed for more than
    `max_types_sharing` of the requested types are dropped as non-specific; mitochondrial and ribosomal genes are dropped."""
    import pandas as pd
    cand = {}
    for name, v in raw.items():
        if not isinstance(v, dict) or "error" in v:
            continue
        can = [g["symbol"].upper() for g in v.get("canonical", []) if g.get("symbol")]
        com = sorted([g for g in v.get("computational", []) if g.get("symbol") and g.get("organism") == organism],
                     key=lambda g: -float(g.get("marker_score") or 0))
        genes = list(dict.fromkeys(can + [g["symbol"].upper() for g in com]))
        genes = [g for g in genes if not g.startswith(("MT-", "RPL", "RPS", "MRPL", "MRPS"))]
        cand[name] = genes
    share = {}
    for genes in cand.values():
        for g in genes:
            share[g] = share.get(g, 0) + 1
    rows = []
    for name, genes in cand.items():
        keep = [g for g in genes if share[g] <= max_types_sharing][:max_per_type]
        rows.append({"cell_type": name, "positive": ";".join(keep), "negative": "",
                     "source": f"CellGuide {raw[name].get('id', '')} canonical + {organism} computational markers; genes shared by more than {max_types_sharing} types dropped"})
    return pd.DataFrame(rows, columns=["cell_type", "positive", "negative", "source"])

def ann_suggest_negatives(panel, per_type=3):
    """Negative markers for types that have none: the first `per_type` positive genes of the other types (not positive for this type)."""
    import pandas as pd
    out = panel.copy()
    pos = {r.cell_type: ann_split(r.positive) for r in panel.itertuples()}
    auto = []
    for i, r in enumerate(panel.itertuples()):
        if ann_split(r.negative):
            continue
        mine = set(pos[r.cell_type])
        neg = []
        for t, genes in pos.items():
            if t == r.cell_type:
                continue
            neg += [g for g in genes if g not in mine][:per_type]
        out.at[i, "negative"] = ";".join(dict.fromkeys(neg))
        auto.append(r.cell_type)
    return out, auto

def ann_gene_stats(adata, genes, labels, names):
    """Mean expression and detected fraction of `genes` in each cluster. Gene names match case-insensitively. Returns (means, detect) DataFrames."""
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    up = pd.Series(range(adata.n_vars), index=[str(v).upper() for v in adata.var_names])
    up = up[~up.index.duplicated()]
    genes = [g for g in dict.fromkeys(g.upper() for g in genes) if g in up.index]
    gi = up[genes].to_numpy()
    X = adata.X[:, gi]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    codes = pd.Categorical(labels, categories=names).codes
    means = np.zeros((len(names), len(genes)))
    det = np.zeros_like(means)
    for k in range(len(names)):
        sub = X[codes == k]
        if len(sub):
            means[k] = sub.mean(axis=0)
            det[k] = (sub > 0).mean(axis=0)
    return pd.DataFrame(means, index=names, columns=genes), pd.DataFrame(det, index=names, columns=genes)

def ann_score_table(panel, means, detect, min_detect=0.25, neg_on=0.5, min_gene_detect=0.05):
    """Score every cluster against every panel type.
    score = mean over the type's positive genes of the gene's cluster mean rescaled to 0 (lowest cluster) .. 1 (highest cluster);
    Positive genes detected in under min_gene_detect of the cells of every cluster are set aside (dropout makes them uninformative); a type left with fewer than min(3, genes given) usable genes (and never fewer than 2) scores 0.
    pos_detect = share of the remaining positive genes detected in the cluster (at least min_detect of its cells, or half the rate of the gene's best cluster when it is never that common, never under 10%);
    neg_on = negative genes that are on in the cluster (rescaled mean at least neg_on and detected in at least min_detect)."""
    import numpy as np
    import pandas as pd
    rng = (means.max() - means.min()).replace(0, np.nan)
    resc = ((means - means.min()) / rng).fillna(0.0)
    # a gene counts as detected in a cluster if it reaches min_detect of the cells, or, when it is never that common (sparse droplet data),
    # at least half the rate of its best cluster; never below 10%, since a marker found in fewer cells cannot define a cluster
    cut = np.maximum(np.minimum(min_detect, 0.5 * detect.max()), 0.10)
    on_det = detect.ge(cut, axis=1)
    score, pdet, nneg, cover, neg_genes = {}, {}, {}, {}, {}
    for r in panel.itertuples():
        pos_all = [g.upper() for g in ann_split(r.positive)]
        neg_all = [g.upper() for g in ann_split(r.negative)]
        pos = [g for g in pos_all if g in resc.columns and detect[g].max() >= min_gene_detect]    # genes seen in under 5% of cells everywhere carry no evidence
        neg = [g for g in neg_all if g in resc.columns and g not in pos_all]
        cover[r.cell_type] = len(pos) / len(pos_all) if pos_all else 0.0
        if len(pos) >= min(3, len(pos_all)) and len(pos) >= 2:   # too few usable genes (after dropout) cannot name a cluster
            score[r.cell_type] = resc[pos].mean(axis=1)
            pdet[r.cell_type] = on_det[pos].mean(axis=1)
        else:
            score[r.cell_type] = pd.Series(0.0, index=means.index)
            pdet[r.cell_type] = pd.Series(0.0, index=means.index)
        on = ((resc[neg] >= neg_on) & on_det[neg]) if neg else pd.DataFrame(index=means.index)
        nneg[r.cell_type] = on.sum(axis=1) if neg else pd.Series(0, index=means.index)
        neg_genes[r.cell_type] = {c: [g for g in on.columns if on.loc[c, g]] for c in means.index} if neg else {c: [] for c in means.index}
    return pd.DataFrame(score), pd.DataFrame(pdet), pd.DataFrame(nneg), pd.Series(cover), neg_genes

def ann_reference_labels(adata, ref, label_key, n_informative=500, max_per_label=500, min_label_cells=10, min_shared=200, seed=0):
    """Name every cell by its rank correlation with each reference cell type's average profile (a simplified SingleR).
    Uses genes present in both datasets (query selected genes when marked), then the n_informative genes whose average differs most between
    reference types; within each cell, genes are ranked, so the comparison does not depend on the two datasets' normalisation.
    Returns (labels array, best correlation array, info dict) or (None, None, info with 'error')."""
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from scipy.stats import rankdata
    rng = np.random.default_rng(seed)
    qn = pd.Series(range(adata.n_vars), index=[str(v).upper() for v in adata.var_names])
    qn = qn[~qn.index.duplicated()]
    rn = pd.Series(range(ref.n_vars), index=[str(v).upper() for v in ref.var_names])
    rn = rn[~rn.index.duplicated()]
    shared = [g for g in qn.index if g in rn.index]
    if "highly_variable" in adata.var:
        hv = set(str(v).upper() for v in adata.var_names[adata.var["highly_variable"].to_numpy()])
        shared_hv = [g for g in shared if g in hv]
        if len(shared_hv) >= min_shared:
            shared = shared_hv
    if len(shared) < min_shared:
        return None, None, {"error": f"only {len(shared)} genes shared with the reference (need {min_shared})"}
    lab = ref.obs[label_key].astype(str)
    keep = lab.value_counts()
    keep = keep[keep >= min_label_cells].index
    idx = np.concatenate([rng.choice(np.where(lab.to_numpy() == k)[0], size=min(max_per_label, int((lab == k).sum())), replace=False) for k in keep])
    R = ref.X[idx][:, rn[shared].to_numpy()]
    R = R.toarray() if sp.issparse(R) else np.asarray(R)
    if R.max() > 30:   # raw counts: log-normalise
        full = ref.X[idx]
        tot = np.asarray(full.sum(axis=1)).ravel()
        R = np.log1p(R / np.maximum(tot, 1)[:, None] * 1e4)
    rl = lab.to_numpy()[idx]
    cent = np.stack([R[rl == k].mean(axis=0) for k in keep])
    top = np.argsort(-cent.var(axis=0))[:min(n_informative, cent.shape[1])]
    cr = np.stack([rankdata(c[top]) for c in cent])
    cr = (cr - cr.mean(1, keepdims=True)) / np.maximum(cr.std(1, keepdims=True), 1e-9)
    cols = qn[[shared[i] for i in top]].to_numpy()
    out_lab = np.empty(adata.n_obs, dtype=object)
    out_cor = np.zeros(adata.n_obs)
    for s in range(0, adata.n_obs, 2000):
        Q = adata.X[s:s + 2000][:, cols]
        Q = Q.toarray() if sp.issparse(Q) else np.asarray(Q)
        Qr = np.apply_along_axis(rankdata, 1, Q)
        Qr = (Qr - Qr.mean(1, keepdims=True)) / np.maximum(Qr.std(1, keepdims=True), 1e-9)
        C = Qr @ cr.T / Qr.shape[1]
        out_lab[s:s + 2000] = np.array(list(keep))[C.argmax(1)]
        out_cor[s:s + 2000] = C.max(1)
    return out_lab, out_cor, {"n_shared": len(shared), "n_informative": len(top), "labels": list(keep)}

def ann_match_label(ref_label, panel_types, label_map=None):
    """Panel type that a reference label corresponds to: explicit map, then exact (case-insensitive) match, then substring either way."""
    if label_map and ref_label in label_map:
        return label_map[ref_label]
    r = ref_label.lower()
    for t in panel_types:
        if t.lower() == r:
            return t
    for t in panel_types:
        if r in t.lower() or t.lower() in r:
            return t
    return None

def ann_call(score, pdet, nneg, ref=None, sizes=None, min_pos_detect=0.5, min_score=0.3, margin_confident=0.15, margin_ambiguous=0.05,
             min_cluster_cells=10):
    """Per-cluster call and confidence tier from the score tables. ref: {cluster: (panel type or None, agreement)} or None.
    confident: one type clearly ahead, its positives present, none of its negatives on, reference (if given) agrees.
    probable: a lead exists but is modest, the cluster is very small, or the reference disagrees.
    ambiguous: two types are close, negatives of the leading type are on, or the reference disagrees with an already modest lead.
    unassigned: no type has enough of its positive markers detected."""
    import numpy as np
    import pandas as pd
    rows = []
    for c in score.index:
        elig = [t for t in score.columns if pdet.loc[c, t] >= min_pos_detect and score.loc[c, t] >= min_score]
        base = {"cluster": c, "n_cells": int(sizes[c]) if sizes is not None else None}
        if not elig:
            top = score.loc[c].sort_values(ascending=False)
            rows.append({**base, "label": "unassigned", "tier": "unassigned", "best": top.index[0], "best_score": float(top.iloc[0]), "second": "", "second_score": np.nan,
                         "margin": np.nan, "pos_detect": float(pdet.loc[c, top.index[0]]), "negatives_on": 0, "reference": "", "reference_agreement": np.nan,
                         "reason": "no cell type has at least half of its positive markers detected with a clear signal"})
            continue
        s = score.loc[c, elig].sort_values(ascending=False)
        best = s.index[0]
        second = s.index[1] if len(s) > 1 else ""
        sec_score = float(s.iloc[1]) if len(s) > 1 else 0.0
        margin = float(s.iloc[0]) - sec_score
        reasons = []
        if margin < margin_ambiguous:
            tier = "ambiguous"
            reasons.append(f"{best} and {second} score almost the same (difference {margin:.2f})")
        elif nneg.loc[c, best] > 0:
            tier = "ambiguous"
            reasons.append(f"{int(nneg.loc[c, best])} negative marker(s) of {best} are on")
        elif margin >= margin_confident:
            tier = "confident"
        else:
            tier = "probable"
            reasons.append(f"modest lead over {second} ({margin:.2f})")
        rl, ra = ("", np.nan)
        if ref is not None and c in ref:
            rl, ra = ref[c]
            rl = rl or ""
            if rl and rl != best:
                if tier == "confident":
                    tier = "probable"
                elif tier == "probable":
                    tier = "ambiguous"
                reasons.append(f"reference says {rl}")
        if sizes is not None and sizes[c] < min_cluster_cells and tier == "confident":
            tier = "probable"
            reasons.append(f"only {int(sizes[c])} cells")
        label = best if tier in ("confident", "probable") else "ambiguous"
        rows.append({**base, "label": label, "tier": tier, "best": best, "best_score": float(s.iloc[0]), "second": second, "second_score": sec_score if second else np.nan,
                     "margin": margin, "pos_detect": float(pdet.loc[c, best]), "negatives_on": int(nneg.loc[c, best]), "reference": rl, "reference_agreement": ra,
                     "reason": "; ".join(reasons) if reasons else "clear lead, positives present, no negatives on"})
    return pd.DataFrame(rows)

def ann_legends(info):
    """Plain-language figure legends built from the computed numbers."""
    k, t = info["n_clusters"], info["n_types"]
    return {
        "annotation_scores": (f"Each row is a cluster and each column a cell type from the marker panel ({t} types). The colour shows how strongly the cluster expresses that type's "
                              f"positive marker genes compared with the other clusters (0 = lowest cluster, 1 = highest cluster). A black outline marks the type a cluster was named "
                              f"after. Names are only given when one type is clearly ahead; the label on the right gives the confidence tier (confident, probable, ambiguous or unassigned)."),
        "annotation_umap": (f"Each dot is a cell, coloured by the name of its cluster. Grey means the cluster could not be named with enough confidence: ambiguous (two types fit "
                            f"about equally, or markers conflict) or unassigned (no type fits). {info['n_named']} of {k} clusters are named. UMAP is for viewing: distances between "
                            f"separated groups and group sizes are not quantitative. Very small clusters are drawn with larger outlined dots. A name is a label from the marker panel, not proof of identity."),
        "annotation_markers": (f"Each row is a cluster and each column a marker gene from the panel, grouped by the cell type it belongs to. Dot size is the share of the cluster's cells "
                               f"in which the gene is detected; colour is average expression scaled for each gene (dark = lowest, bright = highest across clusters). A cluster that fits "
                               f"a cell type should light up most genes of that group."),
    }

def ann_figures(adata, info, tab, fig_dir, panel, means, detect, score, ccol="cluster"):
    import os
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rc = {"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
          "axes.spines.top": False, "axes.spines.right": False, "xtick.direction": "out", "ytick.direction": "out",
          "axes.titleweight": "normal", "axes.titlelocation": "left", "legend.frameon": False}
    os.makedirs(fig_dir, exist_ok=True)
    files = []
    names = list(score.index)
    types = list(score.columns)
    k = len(names)
    tier = dict(zip(tab["cluster"], tab["tier"]))
    label = dict(zip(tab["cluster"], tab["label"]))
    ncell = dict(zip(tab["cluster"], tab["n_cells"]))
    with matplotlib.rc_context(rc):
        fig, ax = plt.subplots(figsize=(0.5 * len(types) + 3.2, 0.3 * k + 1.8))
        im = ax.imshow(score.to_numpy(), cmap="viridis", vmin=0, vmax=1, aspect="auto")
        for i, c in enumerate(names):
            if label[c] in types:
                j = types.index(label[c])
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="black", lw=1.6))
        ax.set_xticks(range(len(types)))
        ax.set_xticklabels(types, rotation=60, ha="right")
        ax.set_yticks(range(k))
        ax.set_yticklabels([f"{c} (n={ncell[c]})" for c in names])
        ax.set_ylabel("cluster")
        ax2 = ax.secondary_yaxis("right")
        ax2.set_yticks(range(k))
        ax2.set_yticklabels([tier[c] for c in names])
        ax2.spines["right"].set_visible(False)
        ax.set_title("how well each cluster matches each cell type")
        fig.colorbar(im, ax=ax, label="marker expression (scaled)", shrink=0.6, pad=0.22)
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "annotation_scores.png"), dpi=300)
        plt.close(fig)
        files.append("annotation_scores.png")
        um = adata.obsm["X_umap"]
        lab = adata.obs[ccol].astype(str).to_numpy()
        ct = adata.obs["cell_type"].astype(str).to_numpy()
        named = [t for t in types if t in set(ct)]
        cmap = plt.get_cmap("tab20")
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        for grey, nm in (("#d9d9d9", "unassigned"), ("#9e9e9e", "ambiguous")):
            m = ct == nm
            if m.any():
                ax.scatter(um[m, 0], um[m, 1], s=3, c=grey, linewidths=0, rasterized=True, label=nm)
        for i, t in enumerate(named):
            m = ct == t
            ax.scatter(um[m, 0], um[m, 1], s=3, color=cmap(i % 20), linewidths=0, rasterized=True, label=t)
        span = np.ptp(um, axis=0)
        for c in names:
            sel = lab == c
            m = np.median(um[sel], axis=0)
            if sel.sum() < max(30, 0.01 * len(lab)):
                # very small cluster: enlarge its dots (outlined) and put the number beside them
                col = cmap(named.index(ct[sel][0]) % 20) if ct[sel][0] in named else ("#9e9e9e" if ct[sel][0] == "ambiguous" else "#d9d9d9")
                ax.scatter(um[sel, 0], um[sel, 1], s=14, color=col, edgecolors="black", linewidths=0.4, zorder=3)
                ax.text(m[0], m[1] + 0.05 * span[1], c, fontsize=7, ha="center", va="bottom", weight="bold")
            else:
                ax.text(m[0], m[1], c, fontsize=7, ha="center", va="center", weight="bold", bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.6))
        ax.set_title(f"{info['n_named']} of {k} clusters named")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        arrow = dict(arrowstyle="->", lw=0.8)
        ax.annotate("", xy=(0.14, 0.0), xytext=(0.0, 0.0), xycoords="axes fraction", arrowprops=arrow)
        ax.annotate("", xy=(0.0, 0.14), xytext=(0.0, 0.0), xycoords="axes fraction", arrowprops=arrow)
        ax.text(0.07, -0.03, "UMAP1", transform=ax.transAxes, ha="center", va="top", fontsize=7)
        ax.text(-0.03, 0.07, "UMAP2", transform=ax.transAxes, ha="right", va="center", rotation=90, fontsize=7)
        h, l = ax.get_legend_handles_labels()
        leg = ax.legend(h, l, loc="center left", bbox_to_anchor=(1.0, 0.5), markerscale=3)
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "annotation_umap.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)
        files.append("annotation_umap.png")
        show, groups = [], []
        for r in panel.itertuples():
            g = [x.upper() for x in ann_split(r.positive) if x.upper() in means.columns][:3]
            show += [x for x in g if x not in show]
            groups += [(r.cell_type, x) for x in g if x in show and x not in [y for _, y in groups]]
        if show:
            lo, hi = means[show].min(), means[show].max()
            col = ((means[show] - lo) / (hi - lo).replace(0, np.nan)).fillna(0).to_numpy()
            fr = detect[show].to_numpy()
            fig, ax = plt.subplots(figsize=(0.24 * len(show) + 2.4, 0.28 * k + 1.9))
            gx, cy = np.meshgrid(np.arange(len(show)), np.arange(k))
            sc_ = ax.scatter(gx.ravel(), cy.ravel(), s=4 + 70 * fr.ravel(), c=col.ravel(), cmap="viridis", vmin=0, vmax=1, linewidths=0)
            ax.set_xticks(range(len(show)))
            ax.set_xticklabels(show, rotation=90, fontstyle="italic", fontsize=6.5)
            ax.set_yticks(range(k))
            ax.set_yticklabels([f"{c}: {label[c]}" if label[c] in types else f"{c}: {label[c]}" for c in names])
            ax.set_ylim(k - 0.5, -0.5)
            ax.margins(0.03)
            start = 0
            owner = {x: t for t, x in groups}
            for t in types:
                n = sum(1 for x in show if owner.get(x) == t)
                if n:
                    ax.text(start - 0.3, -0.75, t, rotation=35, ha="left", va="bottom", fontsize=6, clip_on=False)
                    if start:
                        ax.axvline(start - 0.5, color="0.8", lw=0.6)
                    start += n
            ax.set_title("marker genes of each panel cell type", pad=60)
            fig.colorbar(sc_, ax=ax, label="scaled mean expression", shrink=0.6)
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, "annotation_markers.png"), dpi=300)
            plt.close(fig)
            files.append("annotation_markers.png")
    return files

def ann_run(adata, panel, reference=None, ref_label_key="cell_type", label_map=None, cluster_key="cluster", min_detect=0.25, min_pos_detect=0.5,
            min_score=0.3, margin_confident=0.15, margin_ambiguous=0.05, neg_on=0.5, min_cluster_cells=10, min_gene_coverage=0.5,
            min_named_cells=0.7, min_ref_agreement=0.7, suggest_negatives=False, seed=0, fig_dir=None):
    """Name clusters from a marker panel with positive and negative markers, optionally compare with a labelled reference, and assign each
    cluster a confidence tier. Writes obs['cell_type'], obs['annotation_tier'], obs['cell_type_candidates'] and uns['annotation'].
    Returns (adata, records). Thresholds are defaults, not literature values."""
    import os
    import re
    import numpy as np
    import pandas as pd
    if cluster_key not in adata.obs:
        return adata, [ann_record("marker panel", "FAIL", f"obs['{cluster_key}'] not found: run checkpoint 5 first")]
    panel = ann_read_panel(panel)
    recs = []
    auto = []
    if suggest_negatives:
        panel, auto = ann_suggest_negatives(panel)
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    names = sorted(np.unique(labels), key=lambda x: (0, int(x)) if x.isdigit() else (1, x))
    sizes = pd.Series(labels).value_counts()
    genes = [g for r in panel.itertuples() for g in ann_split(r.positive) + ann_split(r.negative)]
    means, detect = ann_gene_stats(adata, genes, labels, names)
    present = set(means.columns)
    cover, weak, usable = {}, [], []
    for r in panel.itertuples():
        pos = [g.upper() for g in ann_split(r.positive)]
        n_present = sum(g in present for g in pos)
        cover[r.cell_type] = n_present / len(pos) if pos else 0.0
        if n_present < 2:
            weak.append(f"{r.cell_type} ({n_present} of {len(pos)} positive genes in the data)")
        if n_present >= 1 and cover[r.cell_type] >= min_gene_coverage:
            usable.append(r.cell_type)
    lost = [t for t in panel["cell_type"] if t not in usable]
    if not usable:
        recs.append(ann_record("marker panel", "FAIL", "no cell type has enough of its positive marker genes in this dataset (gene names or species may not match)"))
        return adata, recs
    recs.append(ann_record("marker panel", "PASS" if not lost else "WARN",
                           f"{len(usable)} of {len(panel)} cell types scored" + ("" if not lost else "; left out for missing genes: " + ", ".join(f"{t} ({100 * cover[t]:.0f}% present)" for t in lost)) +
                           "; sources: " + ("; ".join(sorted(set(re.sub(r"CL:\d+ ?", "", s_) for s_ in panel['source'] if s_))) or "not recorded")))
    recs.append(ann_record("panel strength", "PASS" if not weak else "WARN",
                           "every cell type has at least 2 positive genes in the data" if not weak else "fewer than 2 positive genes: " + ", ".join(weak) + "; add markers to the panel file"))
    no_neg = [r.cell_type for r in panel.itertuples() if not ann_split(r.negative)]
    if auto:
        nd = "suggested from the other types for: " + ", ".join(auto) + " (shared or leaky genes can trigger false conflicts; edit the panel to replace them)"
    elif len(no_neg) == len(panel):
        nd = "none given: only positive evidence is used (add negative markers to the panel to enable the conflict check)"
    elif no_neg:
        nd = "given for most types; missing for: " + ", ".join(no_neg)
    else:
        nd = "negative markers were given for every type"
    recs.append(ann_record("negative markers", "INFO" if (auto or no_neg) else "PASS", nd))
    pn = panel[panel["cell_type"].isin(usable)].reset_index(drop=True)
    sparse = {r.cell_type: [g for g in ann_split(r.positive) if g.upper() in detect.columns and detect[g.upper()].max() < 0.05] for r in pn.itertuples()}
    sparse = {t: g for t, g in sparse.items() if g}
    if sparse:
        recs.append(ann_record("sparse panel genes", "INFO", "set aside because they are detected in under 5% of cells in every cluster: " +
                               "; ".join(f"{t}: {', '.join(g)}" for t, g in sparse.items())))
    thin = []
    for r in pn.itertuples():
        given = ann_split(r.positive)
        left = len([g for g in given if g.upper() in detect.columns and detect[g.upper()].max() >= 0.05])
        if left < max(2, min(3, len(given))):
            thin.append(f"{r.cell_type} ({left} of {len(given)} usable)")
    if thin:
        recs.append(ann_record("enough usable markers", "WARN", "too few detectable positive markers to score these types, so no cluster can be named after them: " + ", ".join(thin) + "; add markers that are detected in this kind of data"))
    score, pdet, nneg, _, neg_genes = ann_score_table(pn, means, detect, min_detect, neg_on)
    ref_map, ref_raw, ref_info = None, {}, None
    if reference is not None:
        rl, rc, ref_info = ann_reference_labels(adata, reference, ref_label_key, seed=seed)
        if rl is None:
            recs.append(ann_record("reference agreement", "WARN", ref_info["error"]))
        else:
            ref_map = {}
            for c in names:
                v = pd.Series(rl[labels == c]).value_counts()
                top, agree = v.index[0], float(v.iloc[0] / v.sum())
                ref_raw[c] = (top, agree, float(np.mean(rc[labels == c])))
                ref_map[c] = (ann_match_label(top, usable, label_map) if agree >= 0.5 else None, agree)
    tab = ann_call(score, pdet, nneg, ref_map, sizes, min_pos_detect, min_score, margin_confident, margin_ambiguous, min_cluster_cells)
    if ref_raw:
        tab["reference_label"] = [ref_raw[c][0] for c in tab["cluster"]]
        tab["reference_agreement"] = [ref_raw[c][1] for c in tab["cluster"]]
    cand = {}
    for c in names:
        e = score.loc[c][(pdet.loc[c] >= min_pos_detect) & (score.loc[c] >= min_score)].sort_values(ascending=False)
        cand[c] = " / ".join(e.index[:2])
    tab["candidates"] = [cand[c] for c in tab["cluster"]]
    m_lab, m_tier = dict(zip(tab["cluster"], tab["label"])), dict(zip(tab["cluster"], tab["tier"]))
    adata.obs["cell_type"] = pd.Categorical([m_lab[c] for c in labels])
    adata.obs["annotation_tier"] = pd.Categorical([m_tier[c] for c in labels], categories=list(TIERS))
    adata.obs["cell_type_candidates"] = pd.Categorical([cand[c] for c in labels])
    adata.uns["annotation"] = tab
    adata.uns["annotation_panel"] = pn
    k = len(names)
    cnt = tab["tier"].value_counts()
    recs.append(ann_record("annotation tiers", "INFO", ", ".join(f"{cnt.get(t, 0)} {t}" for t in TIERS) + f" of {k} clusters"))
    named_cells = int(sizes[tab.loc[tab["tier"].isin(["confident", "probable"]), "cluster"]].sum())
    share = named_cells / adata.n_obs
    left = tab[~tab["tier"].isin(["confident", "probable"])]
    recs.append(ann_record("cells named", "PASS" if share >= min_named_cells else "WARN",
                           f"{100 * share:.0f}% of cells are in named clusters (limit {100 * min_named_cells:.0f}%)" +
                           ("" if left.empty else "; not named: " + ", ".join(f"{r.cluster} ({r.tier}: {r.reason})" for r in left.itertuples())[:400])))
    conf = tab[tab["negatives_on"] > 0]
    recs.append(ann_record("negative marker conflicts", "PASS" if conf.empty else "WARN",
                           "no leading type has negative markers switched on" if conf.empty else "; ".join(
                               f"cluster {r.cluster} ({r.best}): " + ", ".join(neg_genes[r.best][r.cluster]) for r in conf.itertuples())))
    cm = adata.uns.get("cluster_markers")
    if cm is not None and len(cm):
        top = cm.groupby(cm["cluster"].astype(str))["gene"].apply(lambda s: set(g.upper() for g in s)).to_dict()
        posmap = {r.cell_type: set(g.upper() for g in ann_split(r.positive)) for r in pn.itertuples()}
        nm = tab[tab["tier"].isin(["confident", "probable"])]
        unsup = [f"{r.cluster} ({r.label})" for r in nm.itertuples() if not (top.get(str(r.cluster), set()) & posmap[r.label])]
        recs.append(ann_record("data-driven marker support", "PASS" if not unsup else "WARN",
                               "every named cluster has at least one panel gene among its own top markers" if not unsup else
                               "no panel gene among the cluster's own top markers: " + ", ".join(unsup)))
    if ref_map is not None:
        cmp_ = tab[tab["tier"].isin(["confident", "probable", "ambiguous"])].copy()
        mapped = [(r.cluster, r.best, ref_map[r.cluster][0], ref_map[r.cluster][1]) for r in cmp_.itertuples() if ref_map[r.cluster][0] is not None]
        if not mapped:
            recs.append(ann_record("reference agreement", "WARN", "no reference label could be matched to a panel cell type; pass label_map to link them"))
        else:
            wsum = sum(sizes[c] for c, _, _, _ in mapped)
            agree = sum(sizes[c] for c, b, r_, _ in mapped if b == r_)
            dis = [f"{c} (panel {b}, reference {r_})" for c, b, r_, _ in mapped if b != r_]
            frac = agree / wsum
            recs.append(ann_record("reference agreement", "PASS" if frac >= min_ref_agreement else "WARN",
                                   f"{100 * frac:.0f}% of the cells in comparable clusters get the same name from the reference ({len(mapped)} of {k} clusters comparable; "
                                   f"{ref_info['n_shared']} shared genes)" + ("" if not dis else "; disagree: " + ", ".join(dis))))
    lab_counts = tab[tab["tier"].isin(["confident", "probable"])].groupby("label")["cluster"].apply(list)
    multi = {t: v for t, v in lab_counts.items() if len(v) > 1}
    if multi:
        recs.append(ann_record("same type in several clusters", "INFO", "; ".join(f"{t}: clusters {', '.join(v)}" for t, v in multi.items()) + " (subtypes, states or over-splitting; compare their markers)"))
    unseen = [t for t in usable if t not in set(tab["label"])]
    if unseen:
        recs.append(ann_record("panel types not found", "INFO", "no cluster was named: " + ", ".join(unseen)))
    if fig_dir is not None:
        os.makedirs(fig_dir, exist_ok=True)
        tab.to_csv(os.path.join(fig_dir, "annotation_table.csv"), index=False)
        pn.to_csv(os.path.join(fig_dir, "marker_panel_used.csv"), index=False)
        info = {"n_clusters": k, "n_types": len(usable), "n_named": int(tab["tier"].isin(["confident", "probable"]).sum())}
        files = ann_figures(adata, info, tab, fig_dir, pn, means, detect, score, cluster_key) if "X_umap" in adata.obsm else []
        leg = ann_legends(info)
        names_f = {"annotation_scores": "annotation_scores.png", "annotation_umap": "annotation_umap.png", "annotation_markers": "annotation_markers.png"}
        with open(os.path.join(fig_dir, "figure_legends.md"), "w", encoding="utf8") as f:
            f.write("# Figure legends\n\n" + "\n\n".join(f"**{k_.replace('_', ' ')}** ({names_f[k_]}): {v}" for k_, v in leg.items() if names_f[k_] in files) + "\n")
    return adata, recs

def write_ann_report(records, prefix="annotation_report", scope=None):
    """Write <prefix>.json and <prefix>.md; return the overall status (worst of PASS/WARN/FAIL).
    States WHY the verdict was reached, lists every check that ran with what it tests, and lists catalogued checks that did not run."""
    cat = ann_catalog()
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
    lines = [f"# Cell type annotation report: {overall}", "", f"**Why:** {why}", "",
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
