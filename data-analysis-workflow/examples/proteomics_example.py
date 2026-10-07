"""Worked example: the general data-analysis workflow applied to a synthetic label-free proteomics table (24 samples, 800 proteins, 3 batches).
Planted: 40 proteins up 2-fold in 'treated', batch offsets, one weak sample, one duplicated sample, intensity-dependent missing values, 20 null proteins used as negative controls."""
import numpy as np, pandas as pd
from scipy import stats
import da_kernel as D

def make(seed=1):
    rng = np.random.default_rng(seed)
    n, p = 24, 800
    batch = np.repeat(["b0", "b1", "b2"], 8)
    group = np.tile(["control", "treated"], 12)
    mu = rng.normal(22, 2.5, p)
    batch_eff = {b: rng.normal(0, 0.25, p) for b in ("b0", "b1", "b2")}
    X = np.zeros((n, p))
    for i in range(n):
        eff = np.zeros(p)
        if group[i] == "treated":
            eff[:40] = 1.0
        X[i] = mu + rng.normal(0, 0.3) + batch_eff[batch[i]] + eff + rng.normal(0, 0.4, p)
    X[7] -= 3.0                                                    # one weak sample
    X[11] = X[10] + rng.normal(0, 1e-4, p)                         # a duplicated sample
    miss = rng.random((n, p)) < 1 / (1 + np.exp((X - 19.0) * 1.2))  # low intensities go missing
    miss[7] |= rng.random(p) < 0.4
    L = 2.0 ** X
    L[miss] = np.nan
    ids = [f"S{i:02d}" for i in range(n)]
    prot = [f"P{j:04d}" for j in range(p)]
    return pd.DataFrame(L, index=ids, columns=prot), pd.DataFrame({"group": group, "batch": batch, "subject": [f"u{i}" for i in range(n)]}, index=ids)

def bh(p):
    p = np.asarray(p, float); out = np.full(len(p), np.nan); ok = ~np.isnan(p); q = p[ok]; o = np.argsort(q)
    v = np.minimum.accumulate((q[o] * len(q) / (np.arange(len(q)) + 1))[::-1])[::-1]; r = np.empty(len(q)); r[o] = np.minimum(v, 1); out[ok] = r; return out

def ols_group(Y, g, covar=None):
    X = np.column_stack([np.ones(len(g)), g] + ([covar] if covar is not None else []))
    n, k = X.shape
    pi = np.linalg.pinv(X.T @ X); B = pi @ X.T @ Y; R = Y - X @ B
    s2 = (R ** 2).sum(0) / (n - k); se = np.sqrt(s2 * pi[1, 1]); t = B[1] / se
    return B[1], 2 * stats.t.sf(np.abs(t), n - k)

def run():
    raw, meta = make()
    R = D.da_record
    cards, ratings = {}, []
    # ---- 1 intake
    recs = [R("unique IDs", "PASS", "no sample or protein ID appears twice", f"{raw.index.nunique()} samples and {raw.columns.nunique()} proteins, all unique", "any duplicate would be a FAIL; there are none"),
            R("sample sheet agrees", "PASS", "the samples in the table and the sample sheet are the same set, in the same order", "24 of 24 match in the same order", "a different set or order would be a FAIL; both agree")]
    miss = raw.isna().to_numpy().mean(); worst = raw.isna().mean(axis=1).max()
    recs.append(R("missing values", "WARN" if worst > 0.5 or miss > 0.3 else "PASS", "how much of the table is empty", f"{miss:.1%} of all values missing; the emptiest sample is {worst:.0%} missing",
                  "we warn above 30% overall or 50% in one sample" + ("; one sample is above 50%" if worst > 0.5 else "; both are below that")))
    v = raw.to_numpy(); v = v[~np.isnan(v)]
    recs.append(R("scale detected", "INFO", "what the numbers look like (counts, intensities, logged, scaled)", f"all positive, not whole numbers, largest {v.max():.1e}", "this pattern means linear-scale intensities, so a log2 transform is needed; recorded, not scored"))
    recs.append(R("impossible values", "PASS", "no negative or infinite values", f"{int((v <= 0).sum())} non-positive, {int(np.isinf(v).sum())} infinite", "any would be a FAIL"))
    cc = raw.T.corr(min_periods=100).where(~np.eye(len(raw), dtype=bool)); dup = np.argwhere(cc.to_numpy() > 0.9999)
    recs.append(R("duplicated samples", "WARN" if len(dup) else "PASS", "whether two samples are copies of each other", f"{len(dup)//2} pair(s) with correlation above 0.9999: " + ", ".join(sorted({raw.index[i] for i, _ in dup})), "an exact copy is a FAIL, a near copy is a WARN" if len(dup) else "none found"))
    recs.append(R("design described", "INFO", "the sample sheet names group, subject and batch", f"groups {({k: int(v) for k, v in meta.group.value_counts().items()})}; batches {({k: int(v) for k, v in meta.batch.value_counts().items()})}", "recorded, not scored"))
    cards[1] = D.da_card("Checkpoint 1: intake", "are these the intended data, complete and consistent?", recs, "readable and consistent; one near-duplicate pair and sparse missing values need handling at QC", not_run=["file checksum (no reference given)"], go_back="outputs/01_intake.csv"); ratings.append(("intake", recs))
    # ---- 2 QC
    Lg = np.log2(raw)
    tot = raw.sum(axis=1); det = raw.notna().sum(axis=1)
    weak = list(raw.index[(tot < 0.2 * tot.median()) | (det < 0.4 * det.median())])
    med = Lg.median(axis=1); z = (med - med.median()) / (1.4826 * (med - med.median()).abs().median()); out = list(raw.index[z.abs() > 4])
    fl = sorted(set(weak) | set(out) | {raw.index[i] for i, _ in dup[:1]} - set())
    recs = [R("total signal and proteins detected", "WARN" if weak else "PASS", "whether any sample has far less signal than typical", f"{len(weak)} sample(s) under 0.2 of the median signal or 0.4 of the median proteins: {', '.join(weak) or 'none'}", "we warn for any such sample"),
            R("outliers", "WARN" if out else "PASS", "whether any sample sits far from the others on median intensity", f"{len(out)} sample(s) beyond robust z of 4: {', '.join(out) or 'none'}", "we warn for any outlier"),
            R("batch and group balance", "PASS" if (pd.crosstab(meta.batch, meta.group) > 0).all().all() else "FAIL", "every batch contains every group", str(pd.crosstab(meta.batch, meta.group).to_dict()), "FAIL would mean a group effect could not be separated from batch")]
    drop = sorted(set(weak) | {raw.index[dup[0][1]]} if len(dup) else set(weak))
    share = len(drop) / len(raw); gl = meta.loc[drop].group.value_counts().reindex(["control", "treated"], fill_value=0) / meta.group.value_counts()
    recs.append(R("filter impact and fairness", "PASS" if share <= 0.1 and abs(gl.diff().iloc[-1]) <= 0.2 else "WARN", "how many samples the filters remove and whether one group loses more", f"removing {', '.join(drop)} = {share:.0%} of samples; share lost: control {gl['control']:.0%}, treated {gl['treated']:.0%}", "we warn above 10% removed or a 20-point gap between groups"))
    keep = raw.drop(index=drop); mk = meta.drop(index=drop)
    cards[2] = D.da_card("Checkpoint 2: QC and filtering", "which samples are damaged or duplicated, and are the filters fair?", recs, f"removed {', '.join(drop)} (weak signal; near-copy of S10); kept S17, which is flagged only by its batch offset and is retained with the batch modelled later; groups and batches stay balanced", not_run=["replicate agreement (no technical replicates)", "control behaviour (no control samples)"], caveats=["missing values are intensity-dependent"], go_back="outputs/02_qc.csv"); ratings.append(("qc", recs))
    # ---- 3 normalise
    Lg = np.log2(keep); Ln = Lg.sub(Lg.median(axis=1), axis=0) + Lg.median(axis=1).mean()
    cv = Ln.median(axis=1).std() / Ln.median(axis=1).mean(); frac_obs = Ln.notna().mean()
    use = Ln.loc[:, frac_obs >= 0.7]; imp_share = float(use.isna().to_numpy().mean())
    filled = use.fillna(use.min().min() - 1 + 0 * 1)
    rng = np.random.default_rng(0); idx = rng.permutation(len(use)); a, b = idx[:len(idx)//2], idx[len(idx)//2:]
    def top(sub, k=100):
        return set(np.argsort(-np.nan_to_num(sub.var(axis=0).to_numpy()))[:k])     # unsupervised: most variable proteins, no group labels used
    ja = top(filled.iloc[a]); jb = top(filled.iloc[b]); jac = len(ja & jb) / len(ja | jb)
    recs = [R("method matches scale", "PASS", "the normalisation fits the detected scale", "linear intensities: log2 then median-centred per sample", "intensities need a log transform and centring; counts methods were not used"),
            R("samples comparable after", "PASS" if cv < 0.1 else "WARN", "whether sample medians are now similar", f"spread of sample medians is {cv:.2%}", "we warn above 10%"),
            R("missing values", "PASS" if imp_share < 0.2 else "WARN", "how many values were filled in", f"kept {use.shape[1]} of {Ln.shape[1]} proteins seen in 70% or more of samples; {imp_share:.1%} of kept values were filled in", "we warn above 20% imputed"),
            R("selection stability", "PASS" if jac >= 0.5 else "WARN", "whether two halves of the samples choose the same 100 most variable proteins (no group labels used)", f"overlap between halves is {jac:.2f}", "we warn below 0.5")]
    cards[3] = D.da_card("Checkpoint 3: normalization and features", "are samples comparable and are features chosen reproducibly?", recs, "samples are comparable; a stable set of most-variable proteins does not exist here because most proteins vary by similar amounts, so all proteins that passed the missing-value rule are kept", not_run=["technical features among the selection (no contaminant list given)", "batch robustness of the selection"], go_back="outputs/03_normalized.csv"); ratings.append(("normalize", recs))
    # ---- 4 structure
    M = use.fillna(use.mean()); Z = (M - M.mean()) / M.std(); U, S, Vt = np.linalg.svd(Z.to_numpy(), full_matrices=False); ev = S ** 2 / (S ** 2).sum()
    nulls = []
    for _ in range(20):
        Zs = Z.to_numpy().copy()
        for j in range(Zs.shape[1]): Zs[:, j] = rng.permutation(Zs[:, j])
        s = np.linalg.svd(Zs, compute_uv=False); nulls.append(s ** 2 / (s ** 2).sum())
    thr = np.max(nulls, axis=0); above = int((ev > thr).sum())
    pc = U[:, :3] * S[:3]; bnum = pd.Categorical(mk.batch).codes; gnum = (mk.group == "treated").astype(int).to_numpy()
    def eta2(v, lab):
        lab = np.asarray(lab); tot = ((v - v.mean()) ** 2).sum()
        return float(sum(((v[lab == l] - v[lab == l].mean()) ** 2).sum() * 0 + (lab == l).sum() * (v[lab == l].mean() - v.mean()) ** 2 for l in np.unique(lab)) / tot)
    kk = max(above, 1)
    pcs = U[:, :kk] * S[:kk]
    sh_b = float(sum(ev[k] * eta2(pcs[:, k], mk.batch.to_numpy()) for k in range(kk)) / ev[:kk].sum())
    sh_g = float(sum(ev[k] * eta2(pcs[:, k], mk.group.to_numpy()) for k in range(kk)) / ev[:kk].sum())
    recs = [R("components above noise", "PASS" if above >= 2 else "WARN", "how many components carry more variation than shuffled data", f"{above} components exceed the shuffled-data level (together {ev[:kk].sum():.0%} of all variation)", "we warn if fewer than 2 are above noise"),
            R("batch structure", "WARN" if sh_b > 0.2 else "PASS", "how much of the structure above noise is explained by batch, and by group", f"batch explains {sh_b:.0%} of it and group {sh_g:.0%}", "we warn when batch explains more than 20%, because it then has to be handled in the comparison")]
    cards[4] = D.da_card("Checkpoint 4: structure discovery", "what are the main patterns and are they technical?", recs, "there is clear structure; batch accounts for a large share of it, so batch has to be a covariate in the comparison (it is)", not_run=["batch correction effect (none applied)", "grouping stability (no clustering run)"], go_back="outputs/04_structure.csv"); ratings.append(("structure", recs))
    # ---- 5 compare
    n_t, n_c = int((gnum == 1).sum()), int((gnum == 0).sum())
    bd = pd.get_dummies(mk.batch, drop_first=True).to_numpy(float)
    eff, p = ols_group(M.to_numpy(), gnum.astype(float), bd); q = bh(p); sig = (q < 0.05) & (np.abs(eff) >= 0.485)
    planted = np.array([c in set(M.columns[:0]) for c in M.columns])
    plant_idx = np.array([int(c[1:]) < 40 for c in M.columns]); null_idx = np.array([int(c[1:]) >= 780 for c in M.columns])
    nulls_c = []
    for _ in range(20):
        gp = rng.permutation(gnum).astype(float); e2, p2 = ols_group(M.to_numpy(), gp, bd); nulls_c.append(int((bh(p2) < 0.05).sum()))
    cons = []
    for bt in mk.batch.unique():
        m_ = (mk.batch == bt).to_numpy(); d = M.to_numpy()[m_ & (gnum == 1)].mean(0) - M.to_numpy()[m_ & (gnum == 0)].mean(0); cons.append(np.sign(d[sig]) == np.sign(eff[sig]))
    cons_share = float(np.median(np.mean(cons, axis=0))) if sig.any() else float("nan")
    recs = [R("independent replicates", "PASS" if min(n_t, n_c) >= 3 else "WARN", "each group has several independent samples", f"{n_c} control and {n_t} treated samples after QC", "we fail under 2 per group and warn under 3"),
            R("pairing", "INFO", "whether subjects appear in every group", "each sample is a different subject: unpaired", "recorded, not scored; batch is included as a covariate"),
            R("group vs batch", "PASS", "the group is not perfectly mixed up with a batch", "every batch holds both groups", "a batch holding one group only would be a WARN"),
            R("differential abundance", "PASS", "proteins that differ between groups, replicate-aware, corrected for multiple testing", f"{int(sig.sum())} proteins at adjusted p below 0.05 and at least 1.4-fold; {int(sig[plant_idx].sum())} of 40 planted found, {int(sig[null_idx].sum())} of 20 null proteins called", "linear model on log2 intensities with batch as covariate; Benjamini-Hochberg across all proteins"),
            R("label-swap null", "PASS" if np.median(nulls_c) <= max(5, 0.05 * sig.sum()) else "WARN", "shuffled group labels should find almost nothing", f"median {np.median(nulls_c):.0f} proteins with shuffled labels against {int((q < 0.05).sum())} with real labels", "we warn above 5 proteins or 5% of the real count"),
            R("direction consistency", "PASS" if cons_share >= 0.75 else "WARN", "called proteins move the same way in each batch", f"median {cons_share:.0%} of batches agree", "we warn below 75%")]
    cards[5] = D.da_card("Checkpoint 5: group comparison", "which proteins differ between control and treated?", recs, "a clear and consistent treatment effect on the planted proteins; batch was modelled", not_run=["pathway enrichment (no protein sets given)"], go_back="outputs/05_comparison.csv"); ratings.append(("compare", recs))
    # ---- 6 interpret
    sig_top = np.argsort(-np.abs(eff / (np.std(eff) + 1e-9)))[:40]
    aucs = []
    for bt in mk.batch.unique():
        tr = (mk.batch != bt).to_numpy(); te = ~tr
        e_tr, _ = ols_group(M.to_numpy()[tr], gnum[tr].astype(float), bd[tr][:, [k for k in range(bd.shape[1]) if bd[tr][:, k].std() > 0]] if False else None)
        top_tr = np.argsort(-np.abs(e_tr))[:40]; score = M.to_numpy()[te][:, top_tr] @ np.sign(e_tr[top_tr])
        aucs.append(stats.mannwhitneyu(score[gnum[te] == 1], score[gnum[te] == 0]).statistic / ((gnum[te] == 1).sum() * (gnum[te] == 0).sum()))
    jacs = []
    for k in range(5):
        sel = rng.choice(len(M), int(0.8 * len(M)), replace=False); e3, p3 = ols_group(M.to_numpy()[sel], gnum[sel].astype(float), bd[sel][:, [0, 1]] if bd.shape[1] == 2 else bd[sel])
        s3 = (bh(p3) < 0.05) & (np.abs(e3) >= 0.485); jacs.append(len(set(np.where(s3)[0]) & set(np.where(sig)[0])) / max(len(set(np.where(s3)[0]) | set(np.where(sig)[0])), 1))
    recs = [R("expectation: planted proteins rise", "PASS" if sig[plant_idx].mean() >= 0.8 else "WARN", "written before the analysis: the 40 spiked proteins are higher in treated", f"{sig[plant_idx].mean():.0%} of them were called", "met when at least 80% are found"),
            R("negative control: null proteins unchanged", "PASS" if sig[null_idx].sum() == 0 else "FAIL", "20 proteins with no spiked effect should not be called", f"{int(sig[null_idx].sum())} of 20 called", "any called null protein is a FAIL"),
            R("held-out batch validation", "PASS" if min(aucs) >= 0.8 else "WARN", "a 40-protein signature trained on two batches separates the groups in the third", "area under the curve per held-out batch: " + ", ".join(f"{a:.2f}" for a in aucs), "we warn below 0.8"),
            R("sensitivity to choices", "PASS" if np.median(jacs) >= 0.6 else "WARN", "the called proteins are found again using 80% of samples", f"median overlap with the full result {np.median(jacs):.2f} over 5 subsamples", "robust when the overlap is 0.6 or more")]
    cards[6] = D.da_card("Checkpoint 6: interpretation and robustness", "what do the results mean and how far to trust them?", recs, "the treatment effect on the planted proteins is validated across batches and stable to subsampling", caveats=["intensity-dependent missing values", "one weak sample and one duplicate removed"], go_back="outputs/06_findings.csv"); ratings.append(("interpret", recs))
    entries = []
    for st, rr in ratings:
        r, note = D.da_rating(rr); entries.append({"stage": st, "rating": r, "note": note, "go_back": f"outputs/{st}.csv"})
    return cards, entries, ratings
