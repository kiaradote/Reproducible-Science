"""Dataset adapters: turn raw files into an AnnData (cells x genes, sparse) plus what the integrity checks need.
Add one function per dataset. Each returns (adata, extras) where extras has matrix_ids, meta_ids and, when the matrix
file carries its own label row/column, matrix_labels and meta_labels (used to test column order)."""
import re
from itertools import islice
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from loaders import read_genes_by_cells_text, open_text


def sade_feldman(paths, p):
    """GSE120575: TPM matrix (genes x cells, label row) + GEO sample sheet with comment lines before the header."""
    X, ids, genes, labels, layout = read_genes_by_cells_text(paths["matrix"], has_label_row=True,
                                                             chunk_rows=p.get("chunk_rows", 1000))
    with open_text(paths["meta"]) as f:
        head = list(islice(f, 200))
    skip = next(i for i, line in enumerate(head) if line.startswith("Sample name"))
    raw = pd.read_csv(paths["meta"], sep="\t", skiprows=skip, dtype=str, encoding="latin-1")

    def pick(*keys):
        for c in raw.columns:
            if any(k in c.lower() for k in keys):
                return c
        raise KeyError(f"no column containing {keys}; columns are {list(raw.columns)}")

    c_name, c_id = raw.columns[0], pick("title")
    c_pt, c_resp, c_ther = pick("patinet", "patient"), pick("response"), pick("therapy")
    meta = raw[raw[c_name].astype(str).str.startswith("Sample")].rename(
        columns={c_id: "cell_id", c_pt: "patient_timepoint", c_resp: "response", c_ther: "therapy"}).set_index("cell_id")
    obs = meta.reindex(ids)[["patient_timepoint", "response", "therapy"]].copy()
    obs["timepoint"] = obs["patient_timepoint"].str.split("_").str[0]
    obs["patient"] = obs["patient_timepoint"].str.split("_").str[1]
    obs.index = ids
    # The matrix label row marks sorted fractions with a suffix (e.g. Post_P4_T_enriched); the sample sheet does not.
    # Separate it: the sample part is compared with the metadata, the fraction is kept as a covariate.
    suffix = re.compile(r"_(T_enriched|myeloid_enriched)$")
    sample_labels = [suffix.sub("", x) for x in labels]
    obs["sort_fraction"] = [suffix.search(x).group(1) if suffix.search(x) else "unsorted" for x in labels]
    n_suffix = int((obs["sort_fraction"] != "unsorted").sum())
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    adata.uns["loader_layout"] = layout
    return adata, {"matrix_ids": ids, "matrix_labels": sample_labels, "meta_ids": list(meta.index),
                   "meta_labels": meta["patient_timepoint"],
                   "note": (f"{n_suffix} cells carry a sort-fraction suffix in the matrix label row (T_enriched or myeloid_enriched); "
                            "it was removed before comparing with the sample sheet and stored in obs['sort_fraction']")}


def baron_csv(paths, p):
    """GSE84133 per-donor CSV: cells in rows, gene columns, plus barcode and assigned_cluster columns."""
    cols = pd.read_csv(paths["matrix"], nrows=1, index_col=0).columns
    df = pd.read_csv(paths["matrix"], index_col=0, dtype={c: np.int32 for c in cols[2:]})
    genes = list(cols[2:])
    X = sp.csr_matrix(df[genes].to_numpy(dtype=np.float32))
    obs = df[["barcode", "assigned_cluster"]].astype(str).copy()
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    ids = list(df.index)
    return adata, {"matrix_ids": ids, "matrix_labels": obs["assigned_cluster"].to_numpy(), "meta_ids": ids,
                   "meta_labels": obs["assigned_cluster"],
                   "note": "labels and metadata come from the same table, so the order test here is trivially true"}


def kang_mtx(paths, p):
    """GSE96583 batch 2 (Kang 2018): control and IFN-beta stimulated PBMC libraries (10x, genes x cells MatrixMarket), eight donors pooled in each
    library and labelled by demuxlet in a shared metadata table. Keeps singlets only (demuxlet doublets and ambiguous droplets are dropped).
    Barcodes found in both libraries cannot be tied to one metadata row and are dropped. p['max_cells'] optionally subsamples cells at random,
    stratified by donor and condition, to fit small machines (seed p['seed'])."""
    from scipy.io import mmread
    meta = pd.read_csv(paths["meta"], sep="\t", index_col=0)
    genes = pd.read_csv(paths["genes"], sep="\t", header=None)
    parts, bcs, libs = [], [], []
    for lib in ("ctrl", "stim"):
        M = sp.csr_matrix(mmread(paths[f"mtx_{lib}"]).T.tocsr(), dtype=np.float32)       # cells x genes
        bc = pd.read_csv(paths[f"bc_{lib}"], header=None)[0].astype(str).tolist()
        if M.shape[0] != len(bc) or M.shape[1] != len(genes):
            raise ValueError(f"{lib}: matrix {M.shape} does not match {len(bc)} barcodes and {len(genes)} genes")
        parts.append(M); bcs.append(bc); libs += [lib] * len(bc)
    X = sp.vstack(parts).tocsr()
    bc_all = pd.Series(sum(bcs, []))
    lib_all = pd.Series(libs)
    dup = bc_all.duplicated(keep=False).to_numpy()
    in_meta = bc_all.isin(meta.index).to_numpy()
    sing = np.zeros(len(bc_all), bool)
    sing[in_meta] = (meta.loc[bc_all[in_meta], "multiplets"].to_numpy() == "singlet")
    keep = in_meta & ~dup & sing
    n_dup, n_multi = int(dup.sum()), int((in_meta & ~dup & ~sing).sum())
    idx = np.where(keep)[0]
    obs = meta.loc[bc_all[idx]].copy()
    obs["library"] = lib_all[idx].to_numpy()
    obs["sample"] = obs["ind"].astype(str).to_numpy()
    obs["condition"] = obs["stim"].astype(str).to_numpy()
    obs["author_cell_type"] = obs["cell"].astype(str).to_numpy()          # the authors' labels: for evaluation only
    obs = obs[["library", "sample", "condition", "author_cell_type"]]
    obs.index = bc_all[idx].to_numpy()
    note = (f"kept {len(idx)} of {len(bc_all)} droplets: {n_dup} dropped for barcodes present in both libraries, {n_multi} dropped as demuxlet doublets or ambiguous")
    if p.get("max_cells") and len(idx) > p["max_cells"]:
        rng = np.random.default_rng(p.get("seed", 0))
        grp = (obs["sample"] + "_" + obs["condition"]).to_numpy()
        frac = p["max_cells"] / len(idx)
        sel = np.concatenate([rng.choice(np.where(grp == g)[0], size=max(1, int(round(frac * (grp == g).sum()))), replace=False) for g in np.unique(grp)])
        sel.sort()
        idx, obs = idx[sel], obs.iloc[sel]
        note += f"; subsampled to {len(idx)} cells, stratified by donor and condition, seed {p.get('seed', 0)}"
    adata = ad.AnnData(X=X[idx], obs=obs, var=pd.DataFrame(index=genes[1].astype(str).to_numpy()))
    adata.var["ensembl"] = genes[0].astype(str).to_numpy()
    adata.var_names_make_unique()
    ids = list(adata.obs_names)
    return adata, {"matrix_ids": ids, "matrix_labels": adata.obs["library"].to_numpy(), "meta_ids": ids,
                   "meta_labels": pd.Series(meta.loc[ids, "stim"].to_numpy(), index=ids), "note": note}
