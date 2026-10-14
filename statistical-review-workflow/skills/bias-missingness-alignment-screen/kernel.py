
def missingness_report(df, cols=None, by=None):
    """Missing values per column (count, percent) and, if `by` is a group column, percent missing per group
    with a chi-square test of missingness vs group. A significant test shows missingness depends on group
    (not MCAR); a non-significant test does NOT prove MCAR. Reports empty-string and common NA tokens as missing."""
    import pandas as pd
    from scipy import stats
    d = df.copy()
    d = d.replace({"": None, "NA": None, "NaN": None, "nan": None, "NULL": None, "null": None, "N/A": None})
    cols = list(cols) if cols is not None else [c for c in d.columns if c != by]
    rows = []
    for c in cols:
        m = d[c].isna()
        row = {"column": c, "n_missing": int(m.sum()), "pct_missing": 100.0 * m.mean()}
        if by is not None and m.any() and (~m).any():
            tab = pd.crosstab(d[by], m)
            if tab.shape[0] > 1 and tab.shape[1] > 1:
                row["p_missing_vs_" + by] = float(stats.chi2_contingency(tab)[1])
            for g, frac in m.groupby(d[by]).mean().items():
                row["pct_missing[" + str(g) + "]"] = 100.0 * float(frac)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("pct_missing", ascending=False).reset_index(drop=True)


def funnel_compare(eligible, analysed, cols):
    """Selection check: for each numeric column compare the eligible population with the analysed subset.
    Returns n, mean/median in each, and the standardised mean difference (analysed minus eligible, in eligible SD).
    Large differences mean the analysed units are a selected subset; say how the claims could move."""
    import pandas as pd
    out = []
    for c in cols:
        e = pd.to_numeric(eligible[c], errors="coerce").dropna()
        a = pd.to_numeric(analysed[c], errors="coerce").dropna()
        sd = e.std()
        out.append({"column": c, "n_eligible": len(e), "n_analysed": len(a),
                    "pct_retained": 100.0 * len(a) / len(e) if len(e) else None,
                    "mean_eligible": e.mean(), "mean_analysed": a.mean(),
                    "median_eligible": e.median(), "median_analysed": a.median(),
                    "smd": (a.mean() - e.mean()) / sd if sd and sd == sd else None})
    return pd.DataFrame(out)


def p_curve_summary(pvals, alpha=0.05):
    """Crude p-curve and caliper screen on a set of independent reported p-values for the SAME question.
    Right skew (more very small p than p just under alpha) suggests evidential value; a pile-up just under alpha
    is a risk flag only. Needs roughly 10+ independent tests; with fewer it is Not assessable. Cannot demonstrate p-hacking."""
    import numpy as np
    from scipy import stats
    p = np.asarray([x for x in pvals if x is not None], float)
    sig = p[p < alpha]
    n = len(sig)
    if n < 10:
        return {"n_significant": int(n), "status": "Not assessable: fewer than 10 significant p-values"}
    low = int((sig < alpha / 2).sum())
    binom = stats.binomtest(low, n, 0.5, alternative="greater").pvalue
    near = int(((p >= alpha * 0.8) & (p < alpha)).sum())
    above = int(((p >= alpha) & (p < alpha * 1.2)).sum())
    return {"n_significant": int(n), "n_below_alpha_over_2": low, "p_right_skew_binomial": float(binom),
            "n_in_[0.8a,a)": near, "n_in_[a,1.2a)": above,
            "note": "risk flag only; independence of the p-values is assumed and not checked"}


def spec_curve(df, estimate_fn, choices):
    """Forking-paths sensitivity. `choices` maps a choice name to a list of options, e.g.
    {"threshold": [0.01, 0.025], "exclude_outliers": [False, True], "statistic": ["pearson", "spearman"]}.
    `estimate_fn(df, **choice)` must return (estimate, p_value_or_None). Runs every combination, returns a table,
    and the share of specifications whose sign and (if given) significance match the paper's chosen one is up to
    the caller to state against the paper's own specification."""
    import itertools
    import pandas as pd
    names = list(choices)
    rows = []
    for combo in itertools.product(*[choices[n] for n in names]):
        kw = dict(zip(names, combo))
        est, p = estimate_fn(df, **kw)
        r = dict(kw)
        r["estimate"] = est
        r["p"] = p
        rows.append(r)
    return pd.DataFrame(rows)


def alignment_row(claim_id, abstract_wording, results_wording, discussion_wording, figure_or_table,
                  data_supports, tests_run, wording_strength, alignment, note=""):
    """Compare one conclusion across the published record: abstract vs results vs discussion vs figure/table vs data.
    wording_strength: descriptive / associative / causal / predictive / generalised.
    alignment: aligned / overstated / understated / inconsistent / not assessable.
    data_supports: what the recomputed data actually show (or 'not assessed'); tests_run: the tests that produced it.
    'overstated' or 'inconsistent' must name which locations disagree in `note`."""
    strengths = ("descriptive", "associative", "causal", "predictive", "generalised")
    aligns = ("aligned", "overstated", "understated", "inconsistent", "not assessable")
    if wording_strength not in strengths:
        raise ValueError("wording_strength must be one of " + str(strengths))
    if alignment not in aligns:
        raise ValueError("alignment must be one of " + str(aligns))
    if alignment in ("overstated", "inconsistent") and not note:
        raise ValueError("name the locations that disagree in note")
    return {"claim_id": claim_id, "abstract": abstract_wording, "results": results_wording,
            "discussion": discussion_wording, "figure_or_table": figure_or_table, "data_supports": data_supports,
            "tests_run": tests_run, "wording_strength": wording_strength, "alignment": alignment, "note": note}
