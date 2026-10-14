
def recompute_p(test, stat, df1=None, df2=None, two_sided=True):
    """Recompute a p-value from a reported statistic. test in {'t','z','F','chi2','r'}.
    t: df1=df. F: df1,df2. chi2: df1. r: stat=r, df1=n."""
    from scipy import stats
    if test == "t":
        p = stats.t.sf(abs(stat), df1)
        return 2 * p if two_sided else p
    if test == "z":
        p = stats.norm.sf(abs(stat))
        return 2 * p if two_sided else p
    if test == "F":
        return stats.f.sf(stat, df1, df2)
    if test == "chi2":
        return stats.chi2.sf(stat, df1)
    if test == "r":
        n = df1
        t = stat * ((n - 2) / (1 - stat ** 2)) ** 0.5
        return 2 * stats.t.sf(abs(t), n - 2)
    raise ValueError("unknown test")


def ci_p_consistent(estimate, lo, hi, p_reported, ratio=False):
    """Check a reported 95% CI against a reported p (normal approx). ratio=True for HR/OR/RR (log scale).
    Returns dict with the p implied by the CI and whether it agrees at reported precision (within factor 2)."""
    import math
    from scipy import stats
    if ratio:
        estimate, lo, hi = math.log(estimate), math.log(lo), math.log(hi)
    se = (hi - lo) / (2 * 1.959964)
    z = estimate / se
    p_ci = 2 * stats.norm.sf(abs(z))
    ok = (p_reported / 2.0) <= p_ci <= (p_reported * 2.0)
    return {"p_from_ci": p_ci, "p_reported": p_reported, "consistent": ok}


def grim_check(mean, n, decimals=2, items=1):
    """GRIM: can a mean of n integer-valued observations (x items) round to the reported mean?"""
    import math
    total_n = n * items
    target = round(mean, decimals)
    base = math.floor(mean * total_n)
    for s in (base - 1, base, base + 1, base + 2):
        if round(s / total_n + 1e-12, decimals) == target:
            return True
    return False


def multiplicity_table(pvals, alpha=0.05):
    """Raw, Bonferroni, Holm and Benjamini-Hochberg adjusted p-values; returns DataFrame."""
    import numpy as np
    import pandas as pd
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    holm = np.empty(m)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (m - rank) * p[i])
        holm[i] = min(1.0, run)
    bh = np.empty(m)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, p[i] * m / (rank + 1))
        bh[i] = prev
    df = pd.DataFrame({"p": p, "bonferroni": np.minimum(1, p * m), "holm": holm, "bh_fdr": bh})
    df["sig_raw"] = df.p < alpha
    df["sig_bh"] = df.bh_fdr < alpha
    return df


def p_near_threshold(pvals, lo=0.04, hi=0.05, mid_lo=0.03, mid_hi=0.04):
    """Counts of p-values just below 0.05 vs the adjacent bin (crude caliper check; needs many p-values)."""
    import numpy as np
    p = np.asarray(pvals, dtype=float)
    a = int(((p >= lo) & (p < hi)).sum())
    b = int(((p >= mid_lo) & (p < mid_hi)).sum())
    return {"n_in_0.04_0.05": a, "n_in_0.03_0.04": b, "ratio": (a / b) if b else None}


def smd_ci(m1, sd1, n1, m2, sd2, n2):
    """Cohen d and Hedges g with approximate 95% CI from summary statistics."""
    import math
    sp = math.sqrt(((n1 - 1) * sd1 ** 2 + (n2 - 1) * sd2 ** 2) / (n1 + n2 - 2))
    d = (m1 - m2) / sp
    j = 1 - 3 / (4 * (n1 + n2) - 9)
    se = math.sqrt((n1 + n2) / (n1 * n2) + d ** 2 / (2 * (n1 + n2)))
    return {"d": d, "g": d * j, "ci_lo": d - 1.959964 * se, "ci_hi": d + 1.959964 * se}


def permutation_diff(x, y, n_perm=10000, stat="mean", seed=0):
    """Two-sample permutation test on mean (or median) difference; returns (observed, two-sided p).
    Permute at the level of the independent replicate, not the pseudoreplicate."""
    import numpy as np
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    f = np.mean if stat == "mean" else np.median
    obs = f(x) - f(y)
    pool = np.concatenate([x, y])
    k = len(x)
    cnt = 0
    for _ in range(n_perm):
        rng.shuffle(pool)
        if abs(f(pool[:k]) - f(pool[k:])) >= abs(obs) - 1e-12:
            cnt += 1
    return obs, (cnt + 1) / (n_perm + 1)


def cluster_check(df, cluster_col, value_col, group_col=None):
    """Pseudoreplication diagnostic: ICC(1), design effect, effective n; and, if group_col is given
    (two groups), naive row-level Welch t-test vs test on cluster means."""
    import numpy as np
    from scipy import stats
    d = df[[c for c in [cluster_col, value_col, group_col] if c is not None]].dropna()
    g = d.groupby(cluster_col)[value_col]
    k = g.ngroup().nunique()
    n = len(d)
    ni = g.size().values
    n0 = (n - (ni ** 2).sum() / n) / (k - 1) if k > 1 else np.nan
    gm = d[value_col].mean()
    msb = (ni * (g.mean().values - gm) ** 2).sum() / (k - 1) if k > 1 else np.nan
    msw = ((d[value_col] - g.transform("mean")) ** 2).sum() / (n - k) if n > k else np.nan
    icc = (msb - msw) / (msb + (n0 - 1) * msw) if k > 1 and n > k else np.nan
    deff = 1 + (ni.mean() - 1) * icc if icc == icc else np.nan
    out = {"n_rows": n, "n_clusters": int(k), "icc": icc, "design_effect": deff,
           "effective_n": n / deff if deff == deff and deff > 0 else np.nan}
    if group_col is not None:
        levels = d[group_col].unique()
        if len(levels) == 2:
            a = d[d[group_col] == levels[0]][value_col]
            b = d[d[group_col] == levels[1]][value_col]
            out["naive_p"] = stats.ttest_ind(a, b, equal_var=False).pvalue
            cm = d.groupby([group_col, cluster_col])[value_col].mean().reset_index()
            ca = cm[cm[group_col] == levels[0]][value_col]
            cb = cm[cm[group_col] == levels[1]][value_col]
            out["cluster_level_p"] = stats.ttest_ind(ca, cb, equal_var=False).pvalue
            out["n_clusters_per_group"] = (len(ca), len(cb))
    return out
