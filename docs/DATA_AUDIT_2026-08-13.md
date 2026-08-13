# Data Audit — 2026-08-13

Scope: the 6 datasets acquired/organized during the Aim-1 data-acquisition work of 2026-08-12/13 (`Arutyunyan2023_MFI` + its Visium/multiome siblings, `Greenbaum_NatMed_2024`, `VentoTormo_Nature_2018`, `HumanDevelopmentMultiomicAtlas`, `HPA_trophoblast_proteome`), plus one pre-existing dataset (`2026_human_maternal_fetal_Nature`) whose documentation was found stale during this work. This is **not** a full audit of everything under `DATA/` — just the datasets touched by this project so far.

Method: for every dataset, (1) list the actual files on disk, (2) compare byte-for-byte or via local integrity check against the declared source size, (3) cross-check that `link.md`, `DATA/dataset.index.md`, and this repo's `datasets/<id>/dataset.md` agree with each other and with what's actually on disk, (4) resolve any previously-flagged open questions where feasible.

## Summary table

| Dataset | Files present | Size verified | Docs consistent | Open items |
|---|---|---|---|---|
| `Arutyunyan2023_MFI` (scRNAseq) | ✅ 4/4 h5ad | ✅ byte-exact | ✅ | none |
| `Arutyunyan2023_MFI_Visium` (SpatialTranscriptomics) | ✅ 8/8 tar.gz | ✅ gzip/tar integrity (remote HEAD was flaky at check time, see below) | ✅ | none |
| `Arutyunyan2023_MFI_multiome` (scATACseq) | ✅ skeleton only, as intended | n/a (nothing downloaded, by design) | ✅ | decision to skip ~299GB raw FASTQ still stands, unchanged |
| `Greenbaum_NatMed_2024` (scRNAseq) | ✅ SCP2601 structure intact | not independently re-verified this pass (user-downloaded, structure spot-checked) | ✅ | none found |
| `VentoTormo_Nature_2018` (scRNAseq) | ✅ 1/1 h5ad | not re-verified (single file, no declared remote size to check against) | ✅ | **resolved this pass**: confirmed NOT decidua-only (see below) |
| `HumanDevelopmentMultiomicAtlas` (scRNAseq) | ✅ 7/7 organ `.rds` | ✅ byte-exact (after 1 fix, see below) | ✅ | 5/12 organs deliberately not downloaded (documented, by design) |
| `HPA_trophoblast_proteome` (in `1.Databases/signatures/...`) | ✅ all tier files present | ✅ line counts match script output | ✅ | flagged discrepancy vs. KB-cited numbers, documented, not resolved (see prior dataset.md) |
| `2026_human_maternal_fetal_Nature` (pre-existing) | ✅ data was already there | n/a (pre-existing, not downloaded this project) | ❌ → ✅ **fixed this pass** | `link.md` said "skeleton_only / no files mirrored" — false, now corrected |

## Findings and fixes applied in this audit

### 1. Fixed: stale documentation on `2026_human_maternal_fetal_Nature`

`link.md` claimed "This is currently a clean directory scaffold only. No source files have been mirrored back into the folder yet." This was false — `raw/scPlacenta_host.h5ad` (3.06GB) and `raw/snRNA_raw_counts.h5ad` (4.66GB) have been present since 2026-06-04 (file mtimes), and a first-pass fetal-vs-maternal DE + three-signature analysis already exists in `result/v0.1/`. Also had a blank `Project:` field. Both fixed in `link.md`; `DATA/dataset.index.md` was already corrected for this in an earlier session pass (status `skeleton_only` → `imported_v0.1`).

### 2. Resolved: `VentoTormo_Nature_2018` cell-type coverage

Previously flagged as an open question: does `decidua-v3.h5ad` include trophoblast, or is it decidua (maternal) cells only, despite the paper covering the full interface? Inspected directly via `h5py` (no `anndata` package available, read the HDF5 structure manually):

- `obs['CellType']` — 32 categories, including `EVT` (3,626 cells), `SCT` (1,261), `VCT` (9,479) — trophoblast is well represented (≈14,366 cells)
- `obs['Location']` — `Blood` (11,266) / `Decidua` (40,512) / `Placenta` (18,547)

**Conclusion: not decidua-only.** Safe to use as an independent trophoblast reference. `link.md` and this repo's `datasets/VentoTormo_Nature_2018/dataset.md` updated with the full breakdown.

### 3. Caught and fixed during acquisition (already logged, restated here for completeness): HDMA StomachEsophagus truncated download

Not new to this audit pass (caught and fixed same-session on 2026-08-13 during the download itself, documented in `Worklog.md`), but included here since it's exactly the kind of thing this audit is meant to catch: a `curl | tail` pipeline swallowed a mid-transfer failure, silently truncating `StomachEsophagus_RNA_obj_clustered_final.rds` to ~53% of its expected size while the background task reported success. Caught by comparing on-disk size against Zenodo's declared file size; re-downloaded without the pipe; now byte-exact (19,915,325,339 bytes). Root-cause lesson saved to Claude's cross-session memory (`curl_pipe_swallows_exit_code.md`) so it isn't repeated in future sessions on this or other projects.

### 4. Visium (`Arutyunyan2023_MFI_Visium`) verification method note

At the time of this audit, `ftp.ebi.ac.uk` (the ArrayExpress FTP-fire host) was intermittently failing TLS handshakes (`SSL_ERROR_SYSCALL`) when re-checking `Content-Length` via HEAD requests — unrelated to the original download, which completed fine. Fell back to a network-independent check: `gzip -t` + `tar -tzf` on all 8 tarballs. All 8 passed, each with exactly 79 archive entries (a truncated gzip stream fails this check immediately, so this is a reliable completeness signal even without being able to re-confirm the exact remote byte count).

## Not re-verified in this pass (lower priority / lower risk)

- `Greenbaum_NatMed_2024`: structure was spot-checked (directory tree matches the SCP2601 `file_supplemental_info.tsv` manifest) but individual file sizes weren't re-compared against a remote source in this audit — Broad SCP doesn't expose a simple bulk file-size API the way Zenodo/ArrayExpress do. Low risk: this was a single continuous portal download (not the multi-request pattern that caused the HDMA issue), and internal structure is consistent with the SCP-provided manifest.
- `HPA_trophoblast_proteome`: the underlying `proteinatlas.tsv.zip` → `proteinatlas.tsv` unzip would have failed loudly on corruption (it didn't), and the derived tier files' gene counts match exactly what the filter script printed at run time. Not re-run in this audit pass since nothing changed since generation.
- E-MTAB-12595 (Arutyunyan multiome): still not downloaded, decision unchanged (~299GB raw FASTQ, no processed alternative exists). Out of scope for this audit since there's nothing on disk to verify.

## Net result

All 6 newly-acquired datasets are present, complete, and consistently documented across `link.md` / `DATA/dataset.index.md` / this repo's `datasets/*/dataset.md`. One pre-existing dataset's stale documentation was corrected. One previously-open question (Vento-Tormo cell-type coverage) was resolved with a positive result. One real data-integrity bug was caught and fixed (HDMA StomachEsophagus), with the root cause fixed at the tooling level (documented in cross-session memory) rather than just patched once.

**Data acquisition for Aim 1 is complete and verified.** Next step is unrelated to data: resume the Q1–Q6 framework in `2026-GPT-TWEAKR-Oncofetal.md#定义清楚Placenta的问题` and start building the P1/P2/P3/D signatures.
