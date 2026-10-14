# statistical-review-workflow

A staged, reproducible workflow for verifying the statistical claims of a paper against its public data. It ships as a set of skills plus an agent profile (`STAT_REVIEWER`). It is domain-agnostic: each paper gets the full sequence, and conclusions are never carried from one paper to another.

## Principles
- A plausible methodological concern is not a demonstrated statistical error. Every concern records its reasoning, a verification procedure and an attempt to falsify it.
- Verdicts: Verified, Contradicted, Supported concern, Not assessable. Verified and Contradicted require a performed reproduction from data that were actually obtained.
- The study design (observational, independent sampling and replication units) is established before any method is criticised.
- Missing data, inaccessible files and unperformed checks are stated explicitly.
- The reviewer runs on the existing environment: no environments are created and no packages are installed (standard library, or ask first).

## Stages and skills (`skills/`)
| Stage | Skill | Purpose |
|---|---|---|
| orchestrator | `statistical-review-workflow` | quality rules, stage order, report template, reporting-guideline list, general helpers (recompute_p, ci_p_consistent, grim_check, multiplicity_table, smd_ci, permutation_diff, cluster_check, p_near_threshold) |
| 1 | `paper-claim-extraction` | six-part claim ledger (claim and location, required evidence, independent calculation, methodological assessment, conclusion, remaining uncertainty) |
| 2 | `data-retrieval-reproducibility` | data availability, publisher supplement links, accession and file-role inspection, label reconciliation, evidence manifest (source, size, sha256, parse status) |
| 3 | `independent-statistical-verification` | recompute, compare (match/minor/major), sensitivity, documented simulations |
| 4 | `statistical-methodology-audit` | design gate, typed concerns with falsification attempts, priority = confidence x impact |
| 4b | `bias-missingness-alignment-screen` | attrition and missing values, selection and bias direction, forking paths and p-hacking risk flags, conclusion alignment across the published record |
| 5 | `adversarial-review-synthesis` | challenge pass, downgrade or withdraw weak concerns, calibrated per-claim verdicts |

The report template (`skills/statistical-review-workflow/references/report_template.md`) ends with a section listing the limits of the review itself and the items the reader should oversee.

## Agent (`agent/STAT_REVIEWER.md`)
System prompt, skill list and connectors (PubMed, omics-archives, literature, bioRxiv) for the specialist agent.

## Example (`examples/bacsik_2023_eLife_86852/`)
A review of Bacsik et al., eLife 2023 (DOI 10.7554/eLife.86852), with the correlation reproduction table, selection funnel and a raw-data diagnostic plot. It was produced with the first version of this workflow, before the data-retrieval gate, the six-part ledger, the bias/alignment screen and the limitations section were added. It illustrates the style of output, not the current report format.

## Installing the skills
Each folder under `skills/` is a skill (`SKILL.md`, optional `kernel.py` helper sidecar, optional `references/`). `kernel.py` files contain only function definitions, so loading a skill defines helpers without running anything.
