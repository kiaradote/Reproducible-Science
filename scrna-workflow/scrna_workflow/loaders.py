"""Memory-safe loaders: large text matrices go straight into a sparse matrix, never a dense table."""
import gzip
import numpy as np
import pandas as pd
import scipy.sparse as sp


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="latin-1", newline="")
    return open(path, "rt", encoding="latin-1", newline="")


def read_genes_by_cells_text(path, sep="\t", has_label_row=False, chunk_rows=1000):
    """Delimited text, genes in rows and cells in columns, one header line of cell names.
    Optional second line of per-cell labels (e.g. sample labels). Returns
    (X cells x genes CSR float32, cell_ids, gene_ids, labels or None, layout dict).
    Layout is detected, not assumed: header with/without a blank gene column, trailing empty fields."""
    with open_text(path) as f:
        header = f.readline().rstrip("\r\n").split(sep)
        line1 = f.readline().rstrip("\r\n").split(sep)
        line2 = f.readline().rstrip("\r\n").split(sep) if has_label_row else None
    first = line2 if has_label_row else line1
    names = header[1:] if header[0] == "" else header
    while names and names[-1] == "":
        names.pop()
    n = len(names)
    if len(first) < n + 1:
        raise ValueError(f"first data row has {len(first)} fields but header names {n} cells")
    extra = [e for e in first[1 + n:] if e.strip()]
    if extra:
        raise ValueError(f"unexpected non-empty values after the last cell column: {extra[:3]}")
    labels = line1[1:1 + n] if has_label_row else None
    layout = {"header_fields": len(header), "data_fields": len(first), "n_cells": n,
              "header_has_blank_gene_column": header[0] == "", "trailing_fields_after_cells": len(first) - 1 - n,
              "label_row": bool(has_label_row)}
    skip = 1 + int(has_label_row)
    dtype = {0: str, **{i: np.float32 for i in range(1, n + 1)}}
    reader = pd.read_csv(path, sep=sep, header=None, skiprows=skip, index_col=0, usecols=range(n + 1), dtype=dtype,
                         chunksize=chunk_rows, encoding="latin-1", engine="c")
    genes, blocks = [], []
    for chunk in reader:
        genes.extend(chunk.index.astype(str))
        blocks.append(sp.csr_matrix(chunk.to_numpy(dtype=np.float32)).T)   # cells x genes_in_chunk (CSC)
        del chunk
    X = sp.hstack(blocks, format="csc")
    del blocks
    X = X.tocsr()
    layout["n_genes"] = len(genes)
    return X, names, genes, labels, layout
