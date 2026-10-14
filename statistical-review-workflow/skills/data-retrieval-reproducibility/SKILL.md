---
name: data-retrieval-reproducibility
description: "Stage 2 gate of a statistical claim-verification review: read the paper's Data Availability statement and supplementary links, fetch the actual publisher links, check GEO/SRA accessions and classify files (raw reads, processed matrices, sample metadata, code), reconcile sample identifiers with control/disease labels, and keep an evidence manifest with source, size, checksum and parse status. Ships manifest_add, fetch_with_fallback, manifest_save, classify_data_files, reconcile_ids. Standard library only; installs nothing."
---

# Data retrieval and reproducibility checkpoint

Nothing downstream may claim a reproduction using data that were not obtained. This stage ends with an **evidence manifest** and a list of claims whose evidence is obtained, partial, or missing.

## Steps
1. Read the Data Availability statement and the code availability statement. Copy accessions, repository URLs and supplement references verbatim into the manifest notes.
2. Inspect the publisher's supplementary file index. Use the links the publisher actually provides (open the article page or the publisher's supplementary section); do not guess file names or extensions. Record each file's stated title (for example "Additional file 3: Table S2") next to its URL.
3. For GEO/SRA/ArrayExpress accessions: use the omics-archives connector or the archive's own pages to list the series, samples and supplementary files. Classify each file with `classify_data_files` (a guess; confirm by opening it): raw reads, processed expression matrix, sample metadata, code.
4. Download each needed file in its own cell with `fetch_with_fallback([publisher_url, alternative_url, ...], dest, manifest, label)`. If a source fails, try a reasonable alternative (publisher page link, repository copy, Zenodo/Figshare record, the authors' code repository, GEO supplementary files) before declaring the evidence unavailable. If the sandbox proxy blocks a domain, call request_network_access for it; never use mirrors, caches or a spoofed browser header. A block on a private or reserved address cannot be granted: continue without it.
5. Parse each file with the standard library or already-installed packages. Do NOT create environments or install packages. If a parser is missing (for example an Excel reader), read the file with the standard library (an .xlsx is a zip of XML) or ask the user first. Set `parsed=True/False` in the manifest.
6. Match control/disease (or group) labels to sample identifiers with `reconcile_ids(data_ids, metadata_ids)`. Document every unmatched, duplicated or normalised-only match. If labels cannot be matched, any group comparison is Not assessable.
7. Record every file with `manifest_add(...)`: source URL, size, SHA-256 (computed when the file exists), status (retrieved/failed/blocked/not_attempted), parse status, note. Save with `manifest_save(manifest, "evidence_manifest.csv")` and as an artifact.
8. Write the **evidence-to-claim map**: for each claim in the ledger, state obtained / partial / missing, and for partial or missing, exactly what is lacking.

## Rules
- Never describe a result as reproduced unless the inputs are in the manifest with status retrieved and parsed.
- "Not found at the tried address" is not "unavailable": record what was tried.
- Declare the data level you will reproduce from (raw, processed, summary) and any substitution.
- Subsampled or capped retrievals are "not assessed at full depth".
