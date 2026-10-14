
def manifest_new():
    """Start an empty evidence manifest (list of dicts)."""
    return []


def sha256_file(path):
    """SHA-256 of a local file, read in chunks."""
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest_add(manifest, label, source_url, local_path=None, status="retrieved", parsed=None, note=""):
    """Record one evidence file. status in retrieved/failed/blocked/not_attempted.
    parsed: True/False/None (None = not tried). Size and checksum are computed if local_path exists."""
    import os
    if status not in ("retrieved", "failed", "blocked", "not_attempted"):
        raise ValueError("status must be retrieved/failed/blocked/not_attempted")
    size = None
    sha = None
    if local_path and os.path.exists(local_path):
        size = os.path.getsize(local_path)
        sha = sha256_file(local_path)
    elif status == "retrieved":
        raise ValueError("status=retrieved but local file does not exist: " + str(local_path))
    entry = {"label": label, "source_url": source_url, "local_path": local_path, "status": status,
             "size_bytes": size, "sha256": sha, "parsed": parsed, "note": note}
    manifest.append(entry)
    return entry


def fetch_with_fallback(urls, dest, manifest=None, label=None, timeout=60):
    """Try each URL in order (publisher link first, then alternatives). Standard library only.
    Records every attempt in the manifest; returns local path or None. A 403 from the sandbox proxy
    means: call request_network_access for that domain, do not use mirrors or spoofed headers."""
    import urllib.request
    attempts = []
    for u in urls:
        try:
            with urllib.request.urlopen(u, timeout=timeout) as r:
                data = r.read()
            with open(dest, "wb") as f:
                f.write(data)
            if manifest is not None:
                manifest_add(manifest, label or dest, u, dest, "retrieved", None,
                             "attempts before success: " + str(len(attempts)))
            return dest
        except Exception as e:
            attempts.append(u + " -> " + repr(e)[:120])
    if manifest is not None:
        manifest_add(manifest, label or dest, urls[-1] if urls else "", None, "failed", None,
                     "all sources failed: " + " | ".join(attempts))
    return None


def manifest_save(manifest, path):
    """Write the manifest to CSV."""
    import csv
    cols = ["label", "source_url", "local_path", "status", "size_bytes", "sha256", "parsed", "note"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for row in manifest:
            w.writerow(row)
    return path


def classify_data_files(names):
    """Heuristic file-role guess from names: raw_reads / processed_matrix / sample_metadata / code / unknown.
    A guess only: confirm by opening the file."""
    out = {}
    for n in names:
        s = n.lower()
        if s.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz", ".sra", ".bam", ".cram")):
            out[n] = "raw_reads"
        elif any(k in s for k in ("series_matrix", "metadata", "pheno", "sample_info", "samples", "sdrf", "design")):
            out[n] = "sample_metadata"
        elif any(k in s for k in ("count", "matrix", "expression", "tpm", "fpkm", "rpkm", ".h5", ".mtx", "normalized")):
            out[n] = "processed_matrix"
        elif s.endswith((".py", ".r", ".rmd", ".ipynb", ".sh", ".smk")) or "snakefile" in s:
            out[n] = "code"
        else:
            out[n] = "unknown"
    return out


def reconcile_ids(ids_a, ids_b, label_a="a", label_b="b"):
    """Compare two identifier lists (for example data-file sample IDs vs metadata sample IDs).
    Reports exact matches, matches after case/whitespace/underscore-dash normalisation, and unmatched IDs."""
    def norm(x):
        return str(x).strip().lower().replace("-", "_").replace(" ", "_")
    a = list(ids_a)
    b = list(ids_b)
    sa, sb = set(a), set(b)
    exact = sa & sb
    na = {norm(x): x for x in a if x not in exact}
    nb = {norm(x): x for x in b if x not in exact}
    fuzzy = set(na) & set(nb)
    only_a = [na[k] for k in na if k not in fuzzy]
    only_b = [nb[k] for k in nb if k not in fuzzy]
    return {"n_" + label_a: len(sa), "n_" + label_b: len(sb), "exact": len(exact),
            "matched_after_normalisation": sorted((na[k], nb[k]) for k in fuzzy),
            "only_in_" + label_a: only_a, "only_in_" + label_b: only_b,
            "duplicates_" + label_a: len(a) - len(sa), "duplicates_" + label_b: len(b) - len(sb)}
