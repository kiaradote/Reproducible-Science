"""Stage functions for the notebook. Each is fn(previous_checkpoint, stage_params, global_params) -> (output, records)."""
import os, urllib.request
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scrna_integrity_kernel as K1
import scrna_qc_kernel as K2
import scrna_features_kernel as K3
import scrna_dimred_kernel as K4
import scrna_cluster_kernel as K5
import scrna_annotate_kernel as K6
import scrna_interpret_kernel as K7
import scrna_compare_kernel as K8
import adapters
import captions
from ckpt import sha256_file

FETCH_CHECKS = ["file size", "file sha256", "gzip stream"]
ASSEMBLE_CHECKS = ["unique cell IDs", "cell ID sets match metadata", "column order vs metadata", "NaN/inf values",
                   "negative values", "detected scale", "scale matches declaration", "published numbers"]


def write_fetch(records, prefix):
    return K1.write_integrity_report(records, prefix, scope=FETCH_CHECKS)


def write_assemble(records, prefix):
    return K1.write_integrity_report(records, prefix, scope=ASSEMBLE_CHECKS)


def write_qc(records, prefix):
    return K2.write_qc_report(records, prefix)


def fetch(prev, p, g):
    raw = Path(g["raw_dir"]); raw.mkdir(parents=True, exist_ok=True)
    recs, files = [], []
    for f in p["files"]:
        dest = raw / f.get("filename", f["name"])
        if not dest.exists():
            print("downloading", f["url"])
            urllib.request.urlretrieve(f["url"], dest)
        for r in K1.check_file(str(dest), f.get("sha256"), f.get("bytes")):
            r["detail"] = f"{f['name']}: {r['detail']}"
            recs.append(r)
        files.append({"name": f["name"], "path": str(dest), "bytes": dest.stat().st_size, "sha256": sha256_file(dest)})
    return {"files": files}, recs


def assemble(prev, p, g):
    paths = {f["name"]: f["path"] for f in prev["files"]}
    adata, ex = getattr(adapters, p["adapter"])(paths, p)
    recs = K1.check_alignment(ex["matrix_ids"], ex["meta_ids"], ex.get("matrix_labels"), ex.get("meta_labels"))
    if ex.get("note"):
        recs.append(K1.integrity_record("alignment note", "INFO", ex["note"]))
    v, detected = K1.check_values(adata.X, declared=p.get("declared_scale"))
    recs += v
    exp = p.get("expected", {})
    if "shape" in exp:
        recs.append(K1.check_expected("shape (cells x genes)", tuple(int(x) for x in adata.shape), tuple(exp["shape"])))
    for col, counts in exp.get("obs_counts", {}).items():
        obs_counts = {str(k): int(v_) for k, v_ in adata.obs[col].value_counts().items()}
        recs.append(K1.check_expected(f"counts of {col}", obs_counts, counts))
    for col, n in exp.get("n_unique", {}).items():
        recs.append(K1.check_expected(f"n unique {col}", int(adata.obs[col].nunique()), n))
    adata.uns["scale"] = detected
    return adata, recs


def linear_row_sums(X, scale, cols=None):
    """Per-cell sum on the linear scale, in row blocks so memory stays small. cols restricts to some genes."""
    Xs = X if cols is None else X[:, cols]
    out = np.zeros(Xs.shape[0])
    for i in range(0, Xs.shape[0], 2000):
        B = Xs[i:i + 2000].tocsr()
        d = B.data.astype(np.float64)
        if scale == "log2tpm":
            d = 2.0 ** d - 1
        elif scale == "log1p_cp10k":
            d = np.expm1(d)
        B = sp.csr_matrix((d, B.indices, B.indptr), shape=B.shape)
        out[i:i + 2000] = np.asarray(B.sum(axis=1)).ravel()
    return out


def qc_plot(obs, p, path):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = [("n_genes", p.get("min_genes"), p.get("max_genes")), ("total", None, None), ("pct_mt", None, p.get("max_pct_mt"))]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.2))
    for ax, (c, lo, hi) in zip(axes, cols):
        v = obs[c].dropna()
        if len(v):
            ax.hist(v, bins=60, color="0.4")
        for t in (lo, hi):
            if t and len(v):
                ax.axvline(t, color="r", lw=1)
        ax.set_xlabel(c); ax.set_ylabel("cells")
        if c == "pct_mt" and not len(v):
            ax.text(0.5, 0.5, "no mitochondrial genes", ha="center", transform=ax.transAxes)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def qc(prev, p, g):
    adata = prev
    scale = adata.uns.get("scale", g.get("scale"))
    recs, mito = K2.check_mito(adata.var_names, g.get("species", "human"))
    X = adata.X.tocsr(); X.eliminate_zeros(); adata.X = X
    adata.obs["n_genes"] = np.diff(X.indptr)
    adata.obs["total"] = linear_row_sums(X, scale)
    if mito:
        mt = linear_row_sums(X, scale, cols=np.where(adata.var_names.isin(mito))[0])
        adata.obs["pct_mt"] = 100 * mt / np.maximum(adata.obs["total"].to_numpy(), 1e-12)
        recs.append(K2.check_mito_signal(adata.obs["pct_mt"], p.get("max_pct_mt", 20)))
    else:
        adata.obs["pct_mt"] = np.nan
    qc_plot(adata.obs, p, str(Path(g["run_dir"]) / "figures" / "qc_before_filtering.png"))
    captions.write_captions(Path(g["run_dir"]) / "figures", captions.cap_qc(adata.n_obs, p, bool(mito)))
    masks = {f"n_genes<{p['min_genes']}": (adata.obs["n_genes"] < p["min_genes"]).to_numpy()}
    if p.get("max_genes"):
        masks[f"n_genes>{p['max_genes']}"] = (adata.obs["n_genes"] > p["max_genes"]).to_numpy()
    if mito:
        masks[f"pct_mt>{p['max_pct_mt']}"] = (adata.obs["pct_mt"] > p["max_pct_mt"]).to_numpy()
    recs += K2.qc_filter_report(adata.obs, masks, sample_col=p.get("sample_col"),
                                max_total_frac=p.get("max_total_frac", 0.2), max_sample_frac=p.get("max_sample_frac", 0.5))
    drop = np.zeros(adata.n_obs, dtype=bool)
    for m in masks.values():
        drop |= m
    adata = adata[~drop].copy()
    gkeep = adata.X.getnnz(axis=0) >= p["min_cells"]
    recs.append(K2.qc_record("genes removed (<min_cells)", "INFO",
                             f"{int((~gkeep).sum())} of {len(gkeep)} genes seen in fewer than {p['min_cells']} cells", "filters"))
    adata = adata[:, gkeep].copy()
    adata.uns["qc_params"] = {k: str(v) for k, v in p.items()}
    return adata, recs


def write_features(records, prefix):
    return K3.write_feat_report(records, prefix)


def features(prev, p, g):
    adata = prev
    out = K3.feat_run(adata, adata.uns["scale"], species=g.get("species", "human"), n_top_genes=p.get("n_top_genes", 2000),
                       target_sum=p.get("target_sum", 10000.0), batch_key=p.get("batch_key"), sample_col=p.get("sample_col"),
                       plot_path=str(Path(g["run_dir"]) / "figures" / "mean_variance.png"))
    ad_out = out[0]
    captions.write_captions(Path(g["run_dir"]) / "figures", captions.cap_features(int(ad_out.var["highly_variable"].sum()), ad_out.n_vars))
    return out


def write_dimred(records, prefix):
    return K4.write_dr_report(records, prefix)


def dimred(prev, p, g):
    fdir = Path(g["run_dir"]) / "figures" / "dimred"
    out = K4.dr_run(prev, batch_key=p.get("batch_key"), correct_batch=p.get("correct_batch", "auto"), n_neighbors=p.get("n_neighbors", 15),
                    max_pcs=p.get("max_pcs", 50), min_pcs=p.get("min_pcs", 5), n_perm=p.get("n_perm", 10), seed=p.get("seed", 0),
                    fig_dir=str(fdir), variance_threshold=p.get("variance_threshold", 0.9), variance_basis=p.get("variance_basis", "signal"))
    captions.write_captions(fdir, captions.cap_dimred(out[0].uns["dimred"]))
    return out


def write_cluster(records, prefix):
    return K5.write_cl_report(records, prefix)


def cluster(prev, p, g):
    fdir = Path(g["run_dir"]) / "figures" / "cluster"
    out = K5.cl_run(prev, resolution=p.get("resolution"), batch_key=p.get("batch_key"), n_neighbors=p.get("n_neighbors", 15),
                    n_boot=p.get("n_boot", 5), seed=p.get("seed", 0), min_stability=p.get("min_stability", 0.8), fig_dir=str(fdir))
    a_ = out[0]
    info = a_.uns["clustering"]
    sizes = a_.obs["cluster"].value_counts().to_numpy()
    rep_text = f"the first {info['n_pcs']} " + ("batch-corrected " if info.get("rep") == "X_pca_harmony" else "") + "principal components"
    captions.write_captions(fdir, captions.cap_cluster(info, sizes, rep_text))
    return out


def write_annotate(records, prefix):
    return K6.write_ann_report(records, prefix)


def annotate(prev, p, g):
    """Name clusters from the marker panel file (p['panel']), optionally compared with a labelled reference (p['reference'], an .h5ad)."""
    import anndata as ad
    panel_path = p.get("panel")
    if not panel_path or not Path(panel_path).exists():
        raise FileNotFoundError(f"marker panel not found: {panel_path}. Create a CSV with columns cell_type, positive, negative, source "
                                "(genes separated by ';'); see the scrna-cell-type-annotation skill for building one from CellGuide.")
    ref = None
    if p.get("reference"):
        if not Path(p["reference"]).exists():
            raise FileNotFoundError(f"reference dataset not found: {p['reference']}. Provide a labelled .h5ad or set 'reference' to None to skip the comparison.")
        ref = ad.read_h5ad(p["reference"])
    fdir = Path(g["run_dir"]) / "figures" / "annotate"
    out = K6.ann_run(prev, K6.ann_read_panel(panel_path), reference=ref, ref_label_key=p.get("reference_label_key", "cell_type"),
                     label_map=p.get("label_map"), suggest_negatives=p.get("suggest_negatives", False), fig_dir=str(fdir))
    a_ = out[0]
    if "annotation" in a_.uns:
        tab = a_.uns["annotation"]
        captions.write_captions(fdir, captions.cap_annotate(len(tab), int(tab["tier"].isin(["confident", "probable"]).sum()),
                                                            len(a_.uns["annotation_panel"]), a_.n_obs))
    return out


LAST_INTERPRET = {}      # findings and caveats of the latest interpret run, read by the report writer


def write_interpret(records, prefix):
    return K7.write_int_report(records, LAST_INTERPRET.get("findings"), LAST_INTERPRET.get("caveats"), prefix)


def interpret(prev, p, g):
    """Interpretation and validation. Parameters: gene_sets (GMT path; with gene_sets_kegg='hsa' the file is fetched from KEGG when missing),
    expectations (CSV path), signature (dict), perturb (bool), max_cells, nmf_k."""
    import pandas as pd
    run = Path(g["run_dir"])
    gs = None
    path = p.get("gene_sets")
    if path:
        path = Path(path)
        if not path.exists() and p.get("gene_sets_kegg"):
            path.parent.mkdir(parents=True, exist_ok=True)
            K7.int_fetch_kegg_gmt(str(path), p["gene_sets_kegg"])
        if path.exists():
            gs = K7.int_read_gmt(str(path))
    ranked = None
    if "comparison_de" in prev.uns:          # checkpoint 7 ran: rank genes by its test statistic for pre-ranked enrichment
        de = prev.uns["comparison_de"]
        ranked = {t: d.dropna(subset=["stat"]).drop_duplicates("gene").set_index("gene")["stat"].astype(float) for t, d in de.groupby("type")}
    exp = pd.read_csv(p["expectations"]) if p.get("expectations") and Path(p["expectations"]).exists() else None
    caveats = K7.int_collect_caveats(str(run / "reports"))
    a_, recs, det = K7.int_run(prev, gene_sets=gs, expectations=exp, ranked=ranked, signature=p.get("signature"), perturb=p.get("perturb", True),
                               max_cells=p.get("max_cells", 5000), caveats=caveats, nmf_k=tuple(p.get("nmf_k", (5, 8, 12))), seed=p.get("seed", 0),
                               fig_dir=str(run / "figures" / "interpret"))
    LAST_INTERPRET.update(findings=det.get("findings"), caveats=caveats)
    if det.get("findings"):
        a_.uns["interpretation_findings"] = pd.DataFrame(det["findings"]).astype(str)
    return a_, recs


LAST_COMPARE = {}        # design, tables and caveats of the latest compare run, read by the report writer


def write_compare(records, prefix):
    return K8.write_cmp_report(records, LAST_COMPARE.get("det") or {}, LAST_COMPARE.get("caveats"), prefix, question=LAST_COMPARE.get("question"))


def compare(prev, p, g):
    """Control vs treatment. Parameters: condition, sample (required; without condition the stage records that no comparison was requested), batch, groups=[control, treatment],
    gene_sets (GMT path; gene_sets_kegg='hsa' fetches it when missing), min_units, min_cells, fdr, min_lfc, method, state, n_swaps, question."""
    run = Path(g["run_dir"])
    if not p.get("condition"):
        LAST_COMPARE.update(det={}, caveats=None, question=None)
        return prev, [K8.cmp_record("groups", "INFO", "no comparison requested for this dataset: set 'condition' and 'sample' in the stage parameters to compare two groups")]
    gs = None
    path = p.get("gene_sets")
    if path:
        path = Path(path)
        if not path.exists() and p.get("gene_sets_kegg"):
            path.parent.mkdir(parents=True, exist_ok=True)
            K7.int_fetch_kegg_gmt(str(path), p["gene_sets_kegg"])
        if path.exists():
            gs = K7.int_read_gmt(str(path))
    caveats = K7.int_collect_caveats(str(run / "reports"), exclude=("compare", "interpret"))
    a_, recs, det = K8.cmp_run(prev, p["condition"], p["sample"], type_col=p.get("type_col", "cell_type"), batch_col=p.get("batch"), groups=p.get("groups"), gene_sets=gs,
                               counts_layer=p.get("counts_layer", "counts"), min_units=p.get("min_units", 3), min_cells=p.get("min_cells", 10), fdr=p.get("fdr", 0.05),
                               min_lfc=p.get("min_lfc", 0.5), method=p.get("method", "auto"), state=p.get("state", True), n_swaps=p.get("n_swaps", 20),
                               fig_dir=str(run / "figures" / "compare"), seed=p.get("seed", 0))
    LAST_COMPARE.update(det=det, caveats=caveats, question=p.get("question"))
    return a_, recs
