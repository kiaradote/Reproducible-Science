---
name: omics-biological-interpretation
description: "Biological meaning of feature-level results for any omics-type data (bulk or single-cell RNA-seq, proteomics, metabolomics, CRISPR screens): what the changed or marker features mean and how well that is supported. Use after a group comparison or marker finding. Tests the feature list and the full ranking against gene or metabolite sets with the measured features as background (over-representation, up and down separately, and pre-ranked enrichment with enough permutations for many sets), trims redundant terms, flags generic housekeeping or disease-named terms, checks written-down expectations and negative controls, and compares the result with another study (overlap, direction, effect correlation). Plain-English tests with reasons and a PASS/WARN/FAIL verdict."
---

# Biological meaning (checkpoint 6 of the general workflow, usable on its own)

Question answered: what do the changed or marker features mean biologically, how strongly does prior knowledge support that, and does it agree with other studies?
Follow the rules of `data-analysis-workflow` (plain-English card, reasons, quality rating). Input: a result table (feature, effect, test statistic, adjusted p) for each comparison, or marker lists per group, the list of all features measured, and a collection of sets (pathways, categories, gene signatures). Single-cell, bulk and proteomics data all use the same steps.

## Tests, plain definitions and defaults

| Test | What it checks, in plain words | Default rule |
|---|---|---|
| Feature names usable | Features are named the same way as the sets (gene symbols, same species) | WARN if under 50% of measured features appear in any set; FAIL if none do |
| Sets suitable | After keeping only measured features, enough sets of a useful size remain | WARN if fewer than 10 sets keep 10 to 500 members |
| Right background | Enrichment is judged against the features actually measured, not the whole genome or proteome | INFO states the background size; FAIL if the whole genome is used for a targeted assay |
| Over-representation | Features that went up, and those that went down, are tested separately against each set (hypergeometric test, Benjamini-Hochberg across sets) | A term counts at adjusted p below 0.05 with at least 3 shared members |
| Pre-ranked enrichment | The whole ranking by test statistic is tested against each set, so small coordinated shifts are found | Permutation null per set size; use 1,500 or more permutations (a correction over hundreds of sets cannot be passed with 200); WARN if fewer |
| Redundancy | Terms that repeat the same members are merged | Keep the strongest term of any pair sharing more than half their members |
| Generic terms | Terms that mostly reflect housekeeping machinery (ribosome, proteasome) or have disease or infection names | INFO flag; never report these as findings without a specific reason |
| Expectations met | Known biology written before the analysis (positive and negative controls) is reproduced | WARN if any expectation is not met; an untestable one is listed, not counted |
| Agreement with another study | Overlap of significant features, share moving the same way, and correlation of effects | WARN if fewer than 70% of jointly significant features move the same way, or the overlap is not above chance (p over 0.05) |
| Programs (optional) | Co-varying feature programs are reproducible across restarts and tied to groups | Use only when the question asks for programs |
| Regulators (optional) | Regulator activity inferred from a network | Not run unless a network resource is available |

## Where to get sets

- **Human or mouse genes:** KEGG pathways can be downloaded through the KEGG web service (`bio_fetch_kegg_gmt`; check KEGG's terms for non-academic use); Hallmark, Reactome or GO sets come from a GMT file you supply (some sites block scripted downloads).
- **Proteomics:** use gene names of the proteins with the same sets, and measure the background as the proteins quantified.
- **Metabolomics:** use compound-class or pathway sets (KEGG compounds, SMPDB); map names to one identifier system first.
- **CRISPR screens:** rank genes by the screen statistic; essential-gene sets (ribosome, proteasome, spliceosome) should be enriched among depleted genes, and a screen where they are not has likely failed.

## Helpers (kernel.py)

`bio_read_gmt`, `bio_fetch_kegg_gmt`, `bio_filter_sets`, `bio_ora`, `bio_gsea`, `bio_trim`, `bio_bh`, `bio_ranked`, `bio_flag_terms`, `bio_expectations`, `bio_concordance`, `bio_split`. They need numpy, scipy and pandas.

## Pitfalls

- A term name is a hint about shared members, not proof that a pathway is active, and disease-named terms often mean immune or stress genes.
- Mixing up and down features in one list hides opposite effects.
- Using the whole genome as background inflates every enrichment when only part of it was measured.
- Sets built from the same data you test (for example markers from the same experiment) are circular.
- With few features (under about 50) over-representation has little power; say so rather than report no enrichment.

## Tests

`test_bio_kernel.py` (10): multiple-testing values, a planted set found and its near-copy merged, the background changes the result, random features not enriched, pre-ranked enrichment finds one true set among 126 without false alarms, gene-set file reading, generic-term flags, expectations met / not met / not testable, agreement with another study (identical, flipped, unrelated, too few shared), and ranking construction.
