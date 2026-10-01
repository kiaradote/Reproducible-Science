"""Tests for the one-line figure captions."""
import json, os
from pathlib import Path
import pytest
import captions as C

def test_dimred_captions_use_the_computed_numbers():
    info = dict(n_pcs=9, nmax=50, n_cells=1234, corrected=False, has_batch=False, color_what="total signal per cell on a log scale")
    caps = C.cap_dimred(info)
    assert set(caps) == {"pca_variance.png", "pca_covariates.png", "pca_scatter.png", "umap.png"}
    assert "1,234 cells" in caps["umap.png"] and "9" in caps["pca_variance.png"] and "50" in caps["pca_variance.png"]
    both = C.cap_dimred(dict(info, corrected=True, has_batch=True))
    assert "before (left) and after (right) batch correction" in both["pca_scatter.png"] and "batch (left)" in both["umap.png"]

def test_qc_caption_mentions_mito_only_when_available():
    p = dict(min_genes=200, max_genes=None, max_pct_mt=20)
    assert "20% mitochondrial" in C.cap_qc(1000, p, True)["qc_before_filtering.png"]
    assert "empty" in C.cap_qc(1000, p, False)["qc_before_filtering.png"]

def test_cluster_and_annotation_captions_say_what_is_plotted():
    caps = C.cap_cluster(dict(n_clusters=9, resolution=0.5), [10, 20, 30], "the first 12 principal components")
    assert "9 clusters" in caps["cluster_umap.png"] and "0.5" in caps["cluster_resolution.png"] and len(caps) == 4
    a = C.cap_annotate(9, 7, 12, 1937)
    assert "7 of 9 clusters are named" in a["annotation_umap.png"] and len(a) == 3

def test_write_and_merge(tmp_path):
    C.write_captions(tmp_path, {"a.png": "first"})
    C.write_captions(tmp_path, {"b.png": "second"})
    assert C.read_captions(tmp_path) == {"a.png": "first", "b.png": "second"}

@pytest.mark.parametrize("ds", ["baron", "kang"])
def test_every_figure_in_a_real_run_has_a_caption(ds):
    root = Path("runs") / ds / "figures"
    if not root.exists():
        pytest.skip("run folder not present")
    for png in root.rglob("*.png"):
        caps = C.read_captions(png.parent)
        assert png.name in caps and len(caps[png.name]) > 30, f"no caption for {png}"
