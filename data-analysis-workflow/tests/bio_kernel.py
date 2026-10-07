"""Generic helpers for the biological meaning of feature-level results (gene, protein or metabolite lists with effects and statistics): gene-set files, over-representation,
pre-ranked enrichment, redundancy trimming, flags for generic terms, expectation and negative-control checks, and agreement with another study. Prefix bio_.
Needs numpy, scipy and pandas; bio_fetch_kegg_gmt needs internet access. Feature names are compared in upper case."""
import json
import os
import re

HOUSEKEEPING_WORDS = ("ribosome", "proteasome", "spliceosome", "oxidative phosphorylation", "aminoacyl", "rna polymerase", "translation", "protein processing in endoplasmic")
DISEASE_PATTERN = r"\[(hsa|mmu|rno)05\d{3}\]"      # KEGG disease and infection pathways
DISEASE_WORDS = ("infection", "disease", "cancer", "carcinoma", "leukemia", "syndrome", "diabetes", "myocarditis", "virus", "hepatitis", "tuberculosis", "malaria", "asthma", "arthritis", "addiction")

def bio_split(x):
    if x is None or (isinstance(x, float) and x != x):
        return []
    if isinstance(x, (list, tuple, set)):
        return [str(g).strip() for g in x if str(g).strip()]
    return [g.strip() for g in str(x).replace(",", ";").split(";") if g.strip()]

def bio_fetch_kegg_gmt(path, species="hsa"):
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

def bio_read_gmt(path):
    """GMT file -> {set name: set of upper-case gene symbols}. Lines starting with '#' are ignored."""
    out = {}
    for line in open(path, encoding="utf8"):
        if not line.strip() or line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            out[p[0]] = set(g.upper() for g in p[2:] if g)
    return out

def bio_filter_sets(sets, background, min_size=10, max_size=500):
    """Restrict every set to the background and keep sets with min_size to max_size genes left."""
    out = {}
    for k, v in sets.items():
        w = v & background
        if min_size <= len(w) <= max_size:
            out[k] = w
    return out

def bio_bh(p):
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

def bio_ora(genes, sets, background, min_overlap=3):
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
        df["fdr"] = bio_bh(df["p"].to_numpy())
        df = df.sort_values("p").reset_index(drop=True)
    else:
        df["fdr"] = []
    return df

def bio_trim(df, max_jaccard=0.5):
    """Drop a significant term whose overlapping genes largely repeat those of a better term already kept."""
    keep, seen = [], []
    for r in df.itertuples():
        s = set(r.genes.split(";"))
        if all(len(s & t) / max(len(s | t), 1) < max_jaccard for t in seen):
            keep.append(r.Index)
            seen.append(s)
    return df.loc[keep].reset_index(drop=True)

def bio_gsea(ranked, sets, n_perm=1500, min_overlap=10, seed=0, chunk=250):
    """Pre-ranked gene-set enrichment. ranked: Series of scores indexed by upper-case gene (higher = more up). Enrichment score is the running-sum
    statistic weighted by |score|; the null comes from random gene sets of the same size (sizes grouped on a logarithmic grid), drawn n_perm times per size group.
    The p-value is the share of null scores at least as extreme in absolute value; because that cannot go below 1/(n_perm+1), and a correction across hundreds of sets would then
    leave a single truly enriched set above FDR 0.05, the tail beyond what the permutations resolve is extended with an exponential fit to the top 10% of the null scores. Returns term, size, es, nes, p, fdr."""
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
        if p < 20.0 / (len(nul) + 1):                    # few permutations cannot resolve a tiny p: extend the tail with an exponential fit to the top 10% of the null
            q = float(np.quantile(nul, 0.9))
            exc = nul[nul > q] - q
            if len(exc) >= 10 and abs(es) > q and exc.mean() > 0:
                p = float(min(p, max(0.1 * np.exp(-(abs(es) - q) / exc.mean()), 1e-12)))
        rows.append({"term": t, "size": len(idx), "es": float(es), "nes": float(es / max(np.mean(nul), 1e-12)), "p": float(p)})
    df = pd.DataFrame(rows, columns=["term", "size", "es", "nes", "p"])
    if len(df):
        df["fdr"] = bio_bh(df["p"].to_numpy())
        df = df.sort_values("p").reset_index(drop=True)
    else:
        df["fdr"] = []
    return df

def bio_ranked(result, id_col="feature", stat_col="stat"):
    """Series of test statistics indexed by upper-case feature, duplicates and missing values dropped; input for bio_gsea."""
    r = result.dropna(subset=[stat_col]).copy()
    r[id_col] = r[id_col].astype(str).str.upper()
    return r.drop_duplicates(id_col).set_index(id_col)[stat_col].astype(float)

def bio_flag_terms(df, term_col="term"):
    """Add a 'flag' column that says when an enriched term is probably generic: housekeeping machinery or a disease or infection name (membership of such sets often reflects
    common immune, stress or translation genes, not the disease)."""
    out = df.copy()
    flags = []
    for t in out[term_col].astype(str):
        low = t.lower()
        if any(w in low for w in HOUSEKEEPING_WORDS):
            flags.append("housekeeping machinery")
        elif any(w in low for w in DISEASE_WORDS) or re.search(DISEASE_PATTERN, t):
            flags.append("disease or infection name; may reflect generic genes")
        else:
            flags.append("")
    out["flag"] = flags
    return out

def bio_expectations(result, expectations, id_col="feature", effect_col="effect", padj_col="padj", alpha=0.05):
    """Check written-down expectations against a result table. expectations: list of dicts with name, features (list), direction ('up', 'down' or 'none' for a negative control)
    and optional min_fraction (default 0.8: the share of listed features that must move the stated way). Features missing from the result make an expectation 'not testable'
    when fewer than two remain. Returns a list of dicts: expectation, status ('met', 'not met', 'not testable'), found."""
    r = result.copy()
    r[id_col] = r[id_col].astype(str).str.upper()
    r = r.drop_duplicates(id_col).set_index(id_col)
    out = []
    for e in expectations:
        feats = [f.upper() for f in e["features"]]
        have = [f for f in feats if f in r.index]
        if len(have) < 2:
            out.append({"expectation": e["name"], "status": "not testable", "found": f"only {len(have)} of {len(feats)} listed features are in the result"})
            continue
        sub = r.loc[have]
        sig = sub[padj_col] < alpha
        up = (sig & (sub[effect_col] > 0)).mean()
        down = (sig & (sub[effect_col] < 0)).mean()
        need = e.get("min_fraction", 0.8)
        if e["direction"] == "up":
            ok, text = up >= need, f"{up:.0%} of {len(have)} features are significantly higher"
        elif e["direction"] == "down":
            ok, text = down >= need, f"{down:.0%} of {len(have)} features are significantly lower"
        else:
            ok, text = sig.sum() == 0, f"{int(sig.sum())} of {len(have)} control features are significant"
        out.append({"expectation": e["name"], "status": "met" if ok else "not met", "found": text})
    return out

def bio_concordance(a, b, id_col="feature", effect_col="effect", padj_col="padj", alpha=0.05):
    """Agreement between two results (this study and another study, or two conditions). Over the features present in both: how many are significant in each, the overlap with a
    hypergeometric p-value, the Jaccard index, the share of jointly significant features that move the same way, and the Spearman correlation of the effects."""
    import numpy as np
    from scipy.stats import hypergeom, spearmanr
    A = a.copy(); B = b.copy()
    A[id_col] = A[id_col].astype(str).str.upper(); B[id_col] = B[id_col].astype(str).str.upper()
    A = A.drop_duplicates(id_col).set_index(id_col); B = B.drop_duplicates(id_col).set_index(id_col)
    shared = A.index.intersection(B.index)
    if len(shared) < 10:
        return {"n_shared": int(len(shared)), "status": "not testable"}
    A, B = A.loc[shared], B.loc[shared]
    sa, sb = (A[padj_col] < alpha).to_numpy(), (B[padj_col] < alpha).to_numpy()
    both = sa & sb
    p = float(hypergeom.sf(both.sum() - 1, len(shared), sa.sum(), sb.sum())) if sa.sum() and sb.sum() else 1.0
    same = float((np.sign(A[effect_col][both]) == np.sign(B[effect_col][both])).mean()) if both.any() else float("nan")
    rho = float(spearmanr(A[effect_col], B[effect_col]).statistic)
    return {"status": "tested", "n_shared": int(len(shared)), "sig_a": int(sa.sum()), "sig_b": int(sb.sum()), "overlap": int(both.sum()), "overlap_p": p,
            "jaccard": float(both.sum() / max((sa | sb).sum(), 1)), "same_direction": same, "effect_correlation": rho}
