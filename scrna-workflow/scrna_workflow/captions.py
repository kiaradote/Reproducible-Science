"""One-sentence captions that say only what each figure shows (axes, colours, counts). Longer plain-language legends explain the methods."""
import json
from pathlib import Path


def write_captions(fig_dir, caps):
    """Merge {filename: caption} into figure_captions.json in fig_dir."""
    d = Path(fig_dir)
    d.mkdir(parents=True, exist_ok=True)
    f = d / "figure_captions.json"
    cur = json.load(open(f, encoding="utf-8")) if f.exists() else {}
    cur.update(caps)
    json.dump(cur, open(f, "w", encoding="utf-8"), indent=1)


def read_captions(fig_dir):
    f = Path(fig_dir) / "figure_captions.json"
    return json.load(open(f, encoding="utf-8")) if f.exists() else {}


def cap_qc(n_cells, p, has_mito):
    t = f"Histograms of genes detected per cell, total signal per cell and mitochondrial percentage for {n_cells:,} cells before filtering; red lines mark the cut-offs"
    t += f" ({p.get('min_genes')} genes minimum" + (f", {p.get('max_genes')} maximum" if p.get("max_genes") else "") + (f", {p.get('max_pct_mt')}% mitochondrial" if has_mito else "") + ")."
    if not has_mito:
        t += " The mitochondrial panel is empty because no mitochondrial percentage could be computed."
    return {"qc_before_filtering.png": t}


def cap_features(n_selected, n_genes):
    return {"mean_variance.png": f"Mean expression (x, log scale) against normalised variability (y) for {n_genes:,} genes; the {n_selected:,} selected variable genes are red and the rest grey."}


def cap_dimred(info):
    n, nmax, nn = info["n_pcs"], info["nmax"], info["n_cells"]
    caps = {
        "pca_variance.png": (f"Left: variation captured by each of {nmax} principal components, with the level expected from shuffled data as a dashed line and the {n} components used marked by a red line. "
                             f"Right: cumulative share of real structure reached by the first components."),
        "pca_covariates.png": f"Link strength between each of the first {n} principal components (columns) and each per-cell measurement (rows); {n - 0 if False else ''}numbers are absolute correlations.".replace("; numbers", "; numbers") ,
        "pca_scatter.png": (f"{nn:,} cells on the first two principal components" + (", before (left) and after (right) batch correction" if info.get("corrected") else "")
                            + f", coloured by {info.get('color_what', 'total signal per cell')}."),
        "umap.png": (f"UMAP of {nn:,} cells" + (", coloured by batch (left) and by " + str(info.get("color_what", "total signal")) + " (right)" if info.get("has_batch") else
                                                f", coloured by {info.get('color_what', 'total signal per cell')}") + "."),
    }
    caps["pca_covariates.png"] = f"Absolute correlation (colour and number) between each of the first {n} principal components (columns) and each per-cell measurement (rows)."
    return caps


def cap_cluster(info, sizes, rep_text):
    k = info["n_clusters"]
    return {
        "cluster_resolution.png": (f"Number of clusters (top), stability under resampling (middle; line = median adjusted Rand index, band = range, dotted line = target) and silhouette (bottom) "
                                   f"at each tested resolution; the red line marks the chosen resolution {info['resolution']}."),
        "cluster_umap.png": f"UMAP of {sum(sizes):,} cells coloured by {k} clusters, with cluster numbers at the cluster centres; clusters were computed on {rep_text}, not on the UMAP.",
        "cluster_sizes.png": f"Left: number of cells in each of {k} clusters (red = under the minimum size). Right: how well each cluster is recovered when 20% of cells are left out (dotted line = threshold).",
        "cluster_markers.png": f"Dot plot of the top marker genes (columns) for each of {k} clusters (rows): dot size = share of the cluster's cells expressing the gene, colour = mean expression scaled per gene.",
    }


def cap_annotate(n_clusters, n_named, n_types, n_cells):
    return {
        "annotation_scores.png": f"Heat map of marker-panel scores for {n_clusters} clusters (rows) against {n_types} cell types (columns); black outlines mark the assigned type and the right-hand label gives the confidence tier.",
        "annotation_umap.png": f"UMAP of {n_cells:,} cells coloured by assigned cell type; grey = ambiguous or unassigned. {n_named} of {n_clusters} clusters are named.",
        "annotation_markers.png": f"Dot plot of panel marker genes (columns, grouped by cell type) across {n_clusters} clusters (rows): dot size = share of cells expressing the gene, colour = mean expression scaled per gene.",
    }
