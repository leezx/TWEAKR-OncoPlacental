# Data Audit — 2026-08-13

Scope: the 6 datasets acquired/organized during the Aim-1 data-acquisition work of 2026-08-12/13 (`Arutyunyan2023_MFI` + its Visium/multiome siblings, `Greenbaum_NatMed_2024`, `VentoTormo_Nature_2018`, `HumanDevelopmentMultiomicAtlas`, `HPA_trophoblast_proteome`), plus one pre-existing dataset (`2026_human_maternal_fetal_Nature`) whose documentation was found stale during this work. This is **not** a full audit of everything under `DATA/` — just the datasets touched by this project so far.

Method: for every dataset, (1) list the actual files on disk, (2) verify each file against the strongest check actually available for it (see tier legend below — these are not interchangeable, and this report does not claim they are), (3) cross-check that `link.md`, `DATA/dataset.index.md`, and this repo's `datasets/<id>/dataset.md` agree with each other and with what's actually on disk, (4) resolve any previously-flagged open questions where feasible.

## Verification tier legend

Different datasets could only be checked to different strengths, depending on what the source exposes. Listed strongest to weakest — **do not read a lower tier as equivalent to a higher one**:

- **Tier 1 — Byte-exact**: on-disk file size compared against the source's own declared size (Zenodo API `files[].size`, or an equivalent authoritative manifest) and matches exactly.
- **Tier 2 — Archive/structural integrity**: no independently-declared remote size was available to re-check (or the host was unreachable at audit time), so verified via a check that a truncated/corrupted file would fail: `gzip -t`/`tar -tzf` archive integrity, MatrixMarket header declared-nnz vs. actual data-row count, or HDF5 internal dimension cross-consistency (`obs`/`var`/`X` shapes agreeing). This is real evidence of completeness, not a placeholder — but it is evidence the *file is well-formed and internally consistent*, not proof it matches a specific external byte count.
- **Tier 3 — Process/unzip succeeded only**: the acquisition step (download, unzip) completed without error, and any downstream self-generated output is internally consistent with what the generating script reported — but no independent check (remote size or archive integrity beyond a successful unzip) was performed in this audit.

## Summary table

| Dataset | Files present | Verification tier | Docs consistent | Open items |
|---|---|---|---|---|
| `Arutyunyan2023_MFI` (scRNAseq) | ✅ 4/4 h5ad | **Tier 1** — byte-exact vs. portal/Zenodo declared size | ✅ | none |
| `Arutyunyan2023_MFI_Visium` (SpatialTranscriptomics) | ✅ 8/8 tar.gz | **Tier 2** — `gzip -t`+`tar -tzf` pass, 79/79 archive entries each; remote HEAD unreachable at check time (`ftp.ebi.ac.uk` TLS handshake failures), not re-attempted | ✅ | none |
| `Arutyunyan2023_MFI_multiome` (scATACseq) | ✅ skeleton only, as intended | n/a — nothing downloaded, by design | ✅ | decision to skip ~299GB raw FASTQ still stands, unchanged |
| `Greenbaum_NatMed_2024` (scRNAseq) | ✅ SCP2601 structure intact | **Tier 2** — this audit pass: `gzip -t` passed on all 5 `.gz` files; both `matrix.mtx`/`matrix_atac.mtx` (2.6GB/4.9GB, the two largest files) have their MatrixMarket-declared nnz count matched exactly against actual data-row count (92,502,573 and 362,992,363 respectively) | ✅ | none found |
| `VentoTormo_Nature_2018` (scRNAseq) | ✅ 1/1 h5ad | **Tier 2** — this audit pass: HDF5 opened via `h5py`; `obs` index (70,325), `CellType` codes (70,325), and `X` matrix shape ([70325, 31764]) all cross-consistent; `var` index (31,764) matches `X` column count; CSC `indptr` length (31,765) = n_var+1 exactly. No independently-declared remote byte size exists to check against (single portal file, no manifest) | ✅ | **resolved this pass**: confirmed NOT decidua-only (see below) |
| `HumanDevelopmentMultiomicAtlas` (scRNAseq) | ✅ 7/7 organ `.rds` | **Tier 1** — byte-exact vs. Zenodo API declared size (after 1 fix, see below) | ✅ | 5/12 organs deliberately not downloaded (documented, by design) |
| `HPA_trophoblast_proteome` (in `1.Databases/signatures/...`) | ✅ all tier files present | **Tier 2** for the raw source (`unzip -t proteinatlas.tsv.zip` passes clean, this audit pass) / **Tier 3** for the derived tier files (gene counts match exactly what the filter script printed at generation time, but that's internal self-consistency, not independent verification) | ✅ | flagged discrepancy vs. KB-cited numbers, documented, not resolved (see prior dataset.md) |
| `2026_human_maternal_fetal_Nature` (pre-existing) | ✅ data was already there | not re-verified this pass (pre-existing, not acquired by this project; out of scope beyond fixing its docs) | ❌ → ✅ **fixed this pass** | `link.md` said "skeleton_only / no files mirrored" — false, now corrected |

## Findings and fixes applied in this audit

### 1. Fixed: stale documentation on `2026_human_maternal_fetal_Nature`

`link.md` claimed "This is currently a clean directory scaffold only. No source files have been mirrored back into the folder yet." This was false — `raw/scPlacenta_host.h5ad` (3.06GB) and `raw/snRNA_raw_counts.h5ad` (4.66GB) have been present since 2026-06-04 (file mtimes), and a first-pass fetal-vs-maternal DE + three-signature analysis already exists in `result/v0.1/`. Also had a blank `Project:` field. Both fixed in `link.md`; `DATA/dataset.index.md` was already corrected for this in an earlier session pass (status `skeleton_only` → `imported_v0.1`). **Not independently re-verified against a remote source in this audit** — out of scope, this dataset wasn't acquired by this project, only its documentation was touched.

### 2. Resolved: `VentoTormo_Nature_2018` cell-type coverage

Previously flagged as an open question: does `decidua-v3.h5ad` include trophoblast, or is it decidua (maternal) cells only, despite the paper covering the full interface? Inspected directly via `h5py` (no `anndata` package available, read the HDF5 structure manually):

- `obs['CellType']` — 32 categories, including `EVT` (3,626 cells), `SCT` (1,261), `VCT` (9,479) — trophoblast is well represented (≈14,366 cells)
- `obs['Location']` — `Blood` (11,266) / `Decidua` (40,512) / `Placenta` (18,547)

**Conclusion: not decidua-only.** Safe to use as an independent trophoblast reference. This same inspection also produced the Tier 2 structural-integrity evidence in the summary table above (obs/var/X dimensions all cross-consistent). `link.md` and this repo's `datasets/VentoTormo_Nature_2018/dataset.md` updated with the full breakdown.

### 3. Real integrity checks run this pass on `Greenbaum_NatMed_2024` and the HPA raw export

The first version of this audit only spot-checked directory *structure* for `Greenbaum_NatMed_2024` and asserted the HPA zip "would have failed loudly" on corruption without actually testing it — both were caught in review as unverified claims dressed up as verification. Fixed by actually running the checks: `gzip -t` on all 5 `.gz` files in `Greenbaum_NatMed_2024/raw/SCP2601/` (all pass), MatrixMarket declared-nnz vs. actual row count on both matrix files (exact match on both), and `unzip -t` on `proteinatlas.tsv.zip` (clean). Results folded into the summary table and tier legend above.

### 4. Caught and fixed during acquisition (already logged, restated here for completeness): HDMA StomachEsophagus truncated download

Not new to this audit pass (caught and fixed same-session on 2026-08-13 during the download itself, documented in `Worklog.md`), but included here since it's exactly the kind of thing this audit is meant to catch: a `curl | tail` pipeline swallowed a mid-transfer failure, silently truncating `StomachEsophagus_RNA_obj_clustered_final.rds` to ~53% of its expected size while the background task reported success. Caught by comparing on-disk size against Zenodo's declared file size (Tier 1); re-downloaded without the pipe; now byte-exact (19,915,325,339 bytes). Root-cause lesson saved to Claude's cross-session memory (`curl_pipe_swallows_exit_code.md`) so it isn't repeated in future sessions on this or other projects.

## Not independently verified in this pass

- E-MTAB-12595 (Arutyunyan multiome): still not downloaded, decision unchanged (~299GB raw FASTQ, no processed alternative exists). Out of scope for this audit since there's nothing on disk to verify.
- `2026_human_maternal_fetal_Nature`: docs fixed, but the underlying data itself (present since 2026-06-04, not acquired by this project) wasn't re-checked against any remote source in this audit.

## Net result — stated precisely, not blanket

- **Tier 1 (byte-exact)**: `Arutyunyan2023_MFI` (4/4 files), `HumanDevelopmentMultiomicAtlas` (7/7 organs) — matched against an independently-declared remote size.
- **Tier 2 (archive/structural integrity)**: `Arutyunyan2023_MFI_Visium` (8/8), `Greenbaum_NatMed_2024`, `VentoTormo_Nature_2018`, and the HPA raw zip — internally verified as complete and well-formed, but not cross-checked against an independent remote byte count (either none exists, or the host was unreachable at audit time).
- **Tier 3 (weakest)**: the HPA derived tier files — self-consistent with their own generating script, not independently re-verified.
- **Documentation**: all 6 newly-acquired datasets have `link.md` / `DATA/dataset.index.md` / this repo's `datasets/*/dataset.md` cross-checked and consistent. One pre-existing dataset's stale documentation was corrected (data itself not re-verified, out of scope). One previously-open question (Vento-Tormo cell-type coverage) was resolved with a positive result, using the same inspection that produced its Tier 2 evidence.

All 6 datasets have *some* real, executed verification behind them — none are rubber-stamped — but only 2 of the 6 (`Arutyunyan2023_MFI`, `HumanDevelopmentMultiomicAtlas`) can honestly be called byte-exact verified. The other 4 are Tier 2/3, which is the strongest check available for them given what their sources expose, not a gap in effort.

**Data acquisition for Aim 1 is functionally complete.** Next step is unrelated to data: resume the Q1–Q6 framework in `2026-GPT-TWEAKR-Oncofetal.md#定义清楚Placenta的问题` and start building the P1/P2/P3/D signatures.
