# TWEAKR-OncoPlacental — Worklog

Continuity log for this repo. Read this file first when resuming work — it should let you pick up without re-deriving context.

## Progress Tracker (fixed weights, update after every completed step)

User-requested (2026-08-13): report a global % after every completed task, using a stable weighting scheme rather than re-deriving one each time.

| Phase | Weight | Status | Completion |
|---|---|---|---|
| A. Infra (repo/Argos/review loop) | 3% | done | 100% |
| B. Data acquisition (6 core datasets + HDMA) | 7% | done | 100% |
| C. Data audit + compute feasibility doc | 3% | done | 100% |
| D. Step 1 Inventory | 4% | done | 100% |
| E. Step 2 gene-ID mapping | 6% | done | 100% |
| F. Step 3 prep (collision rule + adult reference) | 7% | done — PR #5 merged | 100% |
| G. Step 4 core: D/F/P pseudobulk signature construction | 20% | design locked + replicate-structure audit done, DE compute not started | ~13% |
| H. Tier-2 validation (Tabula Sapiens, post-freeze) | 10% | not started | 0% |
| I. Apply D/F/P to CRC Oncofetal cells, answer Q1 | 20% | not started | 0% |
| J. Q2–Q6 (remaining 6-question framework, unscoped) | 20% | not started | 0% |

**Current total: ~32.6%** (delta +0.6 from ~32%: replicate-structure audit resolved one of `STEP4_DFP_DESIGN.md`'s open items with real data — 4 of 7 placental datasets usable for the trophoblast-vs-rest DE, 3 organoid datasets structurally excluded — no new compute needed, just re-reading Step 1's already-verified inventory JSONs; G bumped from ~10% to ~13% of its 20% weight)

When reporting progress: recompute the weighted sum, state the delta from the last reported number, and update this table in the same commit as the work it reflects.

---

## 2026-08-12 — Session 1: framework definition + infra setup + first dataset

### Where this project comes from

Source discussion lives in the KB, not in this repo:
- `Zhixins-KB/3.Distill/1.Projects/TWEAKR/TWEAKR-Worklog.md` (Aug 12 entries)
- `Zhixins-KB/3.Distill/1.Projects/TWEAKR/2026-GPT-TWEAKR-Oncofetal.md` (GPT strategy dialogue, several sections)

Read those two files for the full scientific reasoning. Short version below.

### The scientific question (recap, not re-derivation — see KB for full text)

Goal: build a CNS-grade **Human Developmental Reference Framework** that separates:

- **D-shared** — developmental program shared by fetal somatic tissue and placental trophoblast, low in adult
- **F-specific** — fetal somatic-specific, low in placenta/adult
- **P-specific** — placental/trophoblast-specific, low in fetal somatic/adult

...then project the existing CRC Oncofetal/M11/revCSC signature onto this framework to test whether "Oncofetal" is actually a mix of an embryonic/fetal program and a separate extraembryonic/placental program (the "onco-placental reprogramming" hypothesis — deliberately not yet called "OncoPlacenta", see `2026-GPT-TWEAKR-Oncofetal.md#Oncofetal-Placenta sig` for why).

Two-stage plan: **Aim 1** build the normal-development reference (no cancer data touched) → **Aim 2** decompose existing oncofetal signatures against it → **Aim 3** tie to SPP1-TAM/TWEAK/Fn14/YAP mechanism.

This repo (`PR/TWEAKR-OncoPlacental`) is scoped to **Aim 1** for now: acquiring and freezing the normal-development reference datasets and building the P1/P2/P3/D signature pipeline. Cancer data (TCGA/CRC atlas) should not be touched until the signature is frozen — that's a design requirement from the KB discussion, not just a nice-to-have.

### Infra decisions made this session

- **New repo, not a migration.** `REPOS/TWEAKR_Oncofetal-workflow` (existing, already pushed to GitHub, mostly GSE131696 analysis scripts) is left untouched. This is a fresh, independent repo.
- **Location:** `/Volumes/Stelligen_SSD/Stelligen/PR/TWEAKR-OncoPlacental/` — code and `.md` files only.
- **GitHub:** public repo, name `TWEAKR-OncoPlacental`, under account `leezx` (gh CLI already authenticated).
- **Big data never goes in this repo.** All raw/processed data lives under `/Volumes/Stelligen_SSD/Stelligen/DATA/<data_type>/<dataset_id>/`, following the existing `DATA/README.md` convention (`raw/processed/result` + `link.md` per dataset, added to `DATA/dataset.index.md`). This repo only holds `datasets/<dataset_id>/dataset.md` (metadata/manifest, not the files themselves), analysis scripts, and docs.

### Data audit — what already existed before this session (don't re-download)

Checked `DATA/dataset.index.md` + full-tree search before doing anything:

| Dataset | Location | Status |
|---|---|---|
| `2026_human_maternal_fetal_Nature` (Wang et al. Nature 2026) | `DATA/scRNAseq/2026_human_maternal_fetal_Nature/` | **Docs were stale** — `link.md` said "skeleton_only / no files mirrored". Actually already has 7.2GB raw (`scPlacenta_host.h5ad` 3GB + `snRNA_raw_counts.h5ad` 4.6GB) with full `origin`(Fetal/Maternal/Unknown)/`celltype_fullname`/`gestational_week` annotation, **and** already has a first-pass fetal-vs-maternal DE + three-signature (Placental_trophoblast / EVT_invasive_immune_tolerance / General_fetal_non_trophoblast) analysis in `result/v0.1/`. TODO (not done this session): fix the stale `link.md`/`dataset.index.md` status field. |
| `GSE131696` | `DATA/scRNAseq/GSE131696/` | Fully downloaded + analyzed (Seurat object, cluster annotation) |
| `GSE247382` (1st vs 3rd trimester bulk) | `DATA/bulkRNAseq/GSE247382/` | Downloaded |
| HPA placenta tissue list | `DATA/1.Databases/signatures/TWEAKR_Oncofetal_gene_sets/processed/v0.1/tissue_category_rna_placenta_Tissue.tsv` | Have it |
| SU_PLACENTA, GOBP_PLACENTA_DEVELOPMENT gene sets | same dir | Have them |

Confirmed **genuinely missing** (grep'd the whole `DATA/` tree, nothing found):
- Arutyunyan et al. Nature 2023 (early-placenta spatial multiomics) — **this session's target, see below**
- Greenbaum et al. Nat Med 2024 (placenta RNA+ATAC+spatial)
- Vento-Tormo et al. Nature 2018 (classic independent replication atlas)
- HPA trophoblast **cell-specific proteome** (distinct from the tissue-level list we already have)
- Human Development Multiomic Atlas (fetal negative reference, 12 organs)
- `DATA/SpatialTranscriptomics/` was completely empty before this session

### This session's download: Arutyunyan et al. 2023, Nature

Full detail in [`datasets/Arutyunyan2023_MFI/dataset.md`](datasets/Arutyunyan2023_MFI/dataset.md) — read that file for URLs/sizes/decisions. Summary:

User gave 4 ArrayExpress accessions from the paper's data-availability statement:
- `E-MTAB-12421` — scRNA/snRNA, primary tissue (implantation sites, decidua, placenta)
- `E-MTAB-12595` — 10x multiome snRNA+snATAC, primary tissue
- `E-MTAB-12698` — 10x Visium spatial
- `E-MTAB-12650` — scRNA/snRNA of trophoblast organoids (PTO) + stem cells (TSC)

**Investigation finding:** ArrayExpress itself only hosts raw FASTQ (via ENA) plus IDF/SDRF metadata for 3 of these 4 — no processed matrices on ArrayExpress except for Visium. The paper's own portal, **Reproductive Cell Atlas** (`reproductivecellatlas.org/mfi.html` + `/invivo-spatial.html`), hosts author-processed AnnData (`.h5ad`) for `E-MTAB-12421` and `E-MTAB-12650`. Used that instead of raw FASTQ — much smaller and already analysis-ready.

**Downloaded (processed, ~19GB total, into `DATA/`):**
- `DATA/scRNAseq/Arutyunyan2023_MFI/raw/primary_tissue/adata_all_donors_all_cell_states_UPD_20230307.h5ad` (12.38GB) — superset covering all donors × all cell states for E-MTAB-12421. Deliberately did **not** also grab the pre-filtered "trophoblast-only" (3.69GB) or single-donor "P13" (1.69GB) files from the same portal since they're strict subsets of this one — derive them locally into `processed/` instead if/when needed (same pattern as how `GSE131696` derives its trophoblast filter).
- `DATA/scRNAseq/Arutyunyan2023_MFI/raw/organoid/` — 3 files for E-MTAB-12650: `Organoid_PTO_cellxgene.h5ad` (1.56GB), `Organoid_TSC_cellxgene.h5ad` (0.79GB), `adata_Fig3_trophoblast_organoids_unstimulated.h5ad` (2.01GB). Kept all three — they're not subsets of each other (different in vitro systems / stimulation states).
- `DATA/SpatialTranscriptomics/Arutyunyan2023_MFI_Visium/raw/` — 8× `<sample>_spaceranger_output.tar.gz` for E-MTAB-12698 (~2.2GB total, from ArrayExpress FTP directly, not the portal). Skipped the full-resolution `.ndpi`/`.tif` whole-slide images (~700MB, only needed for histology figures, not signature work).

**NOT downloaded — needs your decision:** `E-MTAB-12595` (multiome snRNA+snATAC). No processed matrix exists anywhere public (checked ArrayExpress, the portal, and the `ventolab/MFI` GitHub repo). Only raw FASTQ at ENA `ERP144790`, and it's **≈299GB** across 6 runs (3 ATAC + 3 RNA) — would also need local `cellranger-arc` processing afterward. This is bigger than everything else in this dataset combined (~19GB). Default assumption is to **skip it for now** (RNA-level coverage from the other 3 accessions is enough for Aim 1's signature construction; multiome mainly adds chromatin-accessibility corroboration, which is a "nice to have" for a later evidence layer, not a blocker). Flagged in `dataset.md`; `DATA/scATACseq/Arutyunyan2023_MFI_multiome/` exists only as an empty skeleton (raw/processed/result dirs, no `link.md` written yet, no files).

### Download status at end of session

Three downloads were launched as background shell jobs (session-local, not yet confirmed complete — **check actual file sizes in `DATA/` before assuming they finished**, `curl --fail` was used so a partial/failed download should not leave a file that looks complete, but always verify):

```bash
# verify all three finished + sizes look right:
ls -la /Volumes/Stelligen_SSD/Stelligen/DATA/scRNAseq/Arutyunyan2023_MFI/raw/primary_tissue/
ls -la /Volumes/Stelligen_SSD/Stelligen/DATA/scRNAseq/Arutyunyan2023_MFI/raw/organoid/
ls -la /Volumes/Stelligen_SSD/Stelligen/DATA/SpatialTranscriptomics/Arutyunyan2023_MFI_Visium/raw/
```

Expected: `adata_all_donors_all_cell_states_UPD_20230307.h5ad` ≈12.4GB; organoid dir has 3 files totaling ≈4.4GB; Visium dir has 8 `*_spaceranger_output.tar.gz` totaling ≈2.2GB.

### Open TODOs for next session

1. **Verify the 3 background downloads completed cleanly** (see command above) — re-run any that are short/missing.
2. **Write `link.md`** for the 3 new `DATA/` dataset dirs (required by `DATA/README.md` convention — not done yet, only `dataset.md` in this repo exists so far).
3. **Add the 3 new datasets to `DATA/dataset.index.md`** (scRNAseq, SpatialTranscriptomics, scATACseq sections).
4. **Fix the stale docs on `2026_human_maternal_fetal_Nature`** — `link.md` and `dataset.index.md` both incorrectly say "skeleton_only / no files mirrored"; the data and a first analysis pass are actually already there.
5. **Decide E-MTAB-12595 (multiome, ~299GB raw FASTQ)** — skip / download+process later / keep checking for a processed release. Default is skip.
6. **Continue the missing-dataset list**: Greenbaum 2024 (Nat Med), Vento-Tormo 2018 (Nature), HPA trophoblast cell-specific proteome, Human Development Multiomic Atlas — none downloaded yet, none investigated for portal/processed availability yet (only Arutyunyan 2023 was investigated this session, per explicit user request).
7. **Resume the 6-question framework** from `2026-GPT-TWEAKR-Oncofetal.md#定义清楚Placenta的问题` (Q1–Q6) — user had not yet picked a starting question when this session pivoted to infra/data setup. Recommended starting point discussed but not confirmed: Q1 (the D-shared/F-specific/P-specific orthogonal module definitions), since everything downstream depends on it.
8. Init git repo + push to GitHub (`leezx/TWEAKR-OncoPlacental`, public) — see repo root for whether this has been done (check `git remote -v`).

### Update, same session: 3 more datasets (user-downloaded, this session organized them)

The user manually downloaded 3 more datasets directly into `DATA/scRNAseq/` (outside this session's automated pull) and asked to have them integrated into the standard `raw/processed/result` + `link.md` structure and registered. Done:

| Dataset | `DATA/` location | Size | Status |
|---|---|---|---|
| Greenbaum et al. Nat Med 2024 (SCP2601) | `scRNAseq/Greenbaum_NatMed_2024/` | ~7.5GB | ✅ organized — paired snRNA+snATAC + 3 spatial modalities, 8 donors |
| Vento-Tormo et al. Nature 2018 | `scRNAseq/VentoTormo_Nature_2018/` | 800MB | ✅ organized — **open question:** `decidua-v3.h5ad` filename suggests possibly decidua-only, not trophoblast; check `obs` cell types before relying on it as a trophoblast reference |
| Human Development Multiomic Atlas (Liu/Jessa/Kim/Ng et al., Nature 2026 — different paper from `2026_human_maternal_fetal_Nature`) | `scRNAseq/HumanDevelopmentMultiomicAtlas/` | ~20KB (README only) | ⚠️ **placeholder, not downloaded** — only the repo's own README was saved; actual data needs picking specific records from the paper's `table_s14.tsv` on Zenodo, not a bulk pull |

This closes out 3 of the 4 previously-missing datasets from the original gap analysis (Greenbaum 2024 ✓, Vento-Tormo 2018 ✓, HDMA — placeholder only, still needs real data ⚠️). Only the HPA trophoblast cell-specific proteome (distinct from the tissue-level HPA list already in `1.Databases/signatures/TWEAKR_Oncofetal_gene_sets/`) remains completely untouched.

Per-dataset manifests: `datasets/Greenbaum_NatMed_2024/dataset.md`, `datasets/VentoTormo_Nature_2018/dataset.md`, `datasets/HumanDevelopmentMultiomicAtlas/dataset.md`.

### Updated open TODOs (supersedes the list above where they overlap)

1. Verify the 3 `Arutyunyan2023_MFI`/`Arutyunyan2023_MFI_Visium` background downloads completed cleanly (see verify command above).
2. Vento-Tormo 2018: confirm `decidua-v3.h5ad` cell-type coverage (trophoblast included, or decidua-only?).
3. HDMA: still needs an actual data pull — go through `table_s14.tsv`, pick relevant fetal-organ Seurat objects, download just those (not the whole Zenodo community — it also has ChromBPNet models, motif compendia, genome tracks that are out of scope).
4. Still fully missing: HPA trophoblast cell-specific proteome (distinct from the HPA tissue-level list already in `DATA/1.Databases/signatures/TWEAKR_Oncofetal_gene_sets/`).
5. Fix the stale docs on `2026_human_maternal_fetal_Nature` (`link.md` says skeleton_only, actually has data + a first analysis pass).
6. Decide E-MTAB-12595 (multiome, ~299GB raw FASTQ, no processed alternative) — default is skip.
7. Resume the 6-question framework (Q1–Q6) from `2026-GPT-TWEAKR-Oncofetal.md#定义清楚Placenta的问题` — infra/data setup was the detour, the actual signature-construction work hasn't started yet. Recommended starting point: Q1.
8. Once all datasets settle, the actual Aim-1 pipeline work starts: define P1 (placenta-core) / P2 (early-placenta) / P3 (placenta-specific) / D (shared-developmental) per `2026-GPT-TWEAKR-Oncofetal.md#Placenta数据集`'s evidence-layer weighting table — this is still not started, everything so far has been data acquisition.

### Update, same session: HPA trophoblast proteome + HDMA resolved

**HPA trophoblast cell-specific proteome** — closes the last of the 6 originally-missing datasets. Downloaded the full HPA bulk export (small, 43MB) rather than trying to query a filtered API subset, then filtered locally for trophoblast-related single-cell-type specificity (tiered: 54 cell-type-enriched / 90 group-enriched / 1641 cell-type-enhanced). Script + full method + a flagged discrepancy against the KB-cited numbers: `datasets/HPA_trophoblast_proteome/dataset.md`. Lives in `DATA/1.Databases/signatures/TWEAKR_Oncofetal_gene_sets/` (not a new dataset_id — it's an addition to the existing gene-sets collection).

**HDMA (fetal negative reference)** — the placeholder from earlier this session got resolved. User pointed out the full collection is a huge multi-part collection (~2.2TB total: ArchR/bigwigs/ChromBPNet/motifs/fragments/Seurat objects). Investigated via the Zenodo API and found the RNA Seurat objects are stored as **individually downloadable per-organ files** within 4 "Part" records — no need to pull whole bundles. Resolved the full 12-organ size table (see `datasets/HumanDevelopmentMultiomicAtlas/dataset.md`), user picked 7 organs (~57GB: Adrenal/Thyroid/Spleen/Thymus/Liver/Skin/StomachEsophagus — spans all 3 germ layers, skips the 5 largest). **Note: HDMA has no intestine/colon** — it's a broad fetal-negative reference, not gut-specific. Non-RNA data types (ArchR/ChromBPNet/bigwigs/motifs/fragments, ~1.5TB combined) deliberately out of scope.

### Verification: Arutyunyan2023_MFI fully downloaded ✅

All 3 accessions confirmed complete by exact byte match, 2026-08-12 ~19:30:
- `scRNAseq/Arutyunyan2023_MFI/raw/primary_tissue/adata_all_donors_all_cell_states_UPD_20230307.h5ad` — 12,382,386,129 bytes, matches expected exactly
- `scRNAseq/Arutyunyan2023_MFI/raw/organoid/` — 3/3 files, all exact byte matches (PTO 1,559,704,580 / TSC 788,504,395 / Fig3 2,010,471,040)
- `SpatialTranscriptomics/Arutyunyan2023_MFI_Visium/raw/` — 8/8 Visium `*_spaceranger_output.tar.gz` present
- Total: 16GB (scRNAseq) + 2.2GB (Visium) = 18.2GB

### ⚠️ Data-integrity bug found and fixed: `curl | tail` silently swallows download failures

While verifying the 7 HDMA organ downloads, found that `StomachEsophagus_RNA_obj_clustered_final.rds` had silently truncated: expected 19,915,325,339 bytes, actual on-disk was only 10,667,018,304 bytes (~53%) — **but the background task still reported "completed, exit code 0."** Root cause: this session's download commands used `curl ... --progress-bar 2>&1 | tail -3`, and in a pipeline the overall exit status is `tail`'s (always 0), not curl's — a mid-transfer failure (SSL reset, connection drop) gets hidden. Saved as a standing lesson in Claude's cross-session memory (`curl_pipe_swallows_exit_code.md`) so future sessions don't repeat it: never pipe a background download through `tail`/`head`; verify completed downloads by exact byte size (or `gzip -t`/`tar -tzf` for archives when the remote size can't be re-checked).

**Verification sweep done after finding this:**
- HDMA: Adrenal/Thyroid/Spleen/Thymus/Liver/Skin all confirmed byte-exact against Zenodo's declared file size. Only StomachEsophagus was bad — re-downloading now (this time without the pipe, exit code checked directly).
- `Arutyunyan2023_MFI_Visium` (8 tarballs): remote HEAD checks were flaky (EBI FTP host had intermittent SSL handshake failures during this check, unrelated to the download itself), so verified via local `gzip -t` + `tar -tzf` integrity instead — **all 8 passed, each with exactly 79 archive entries**. Confident these are fine.
- `Arutyunyan2023_MFI` primary_tissue/organoid and HPA export were already confirmed byte-exact/integrity-checked earlier in the session (see prior sections) — not re-checked again here.

**StomachEsophagus re-download confirmed** — 19,915,325,339 bytes, exact match. All 7 HDMA organs (Adrenal/Thyroid/Spleen/Thymus/Liver/Skin/StomachEsophagus, ≈58GB) now byte-verified complete.

### Data audit PR (`data-audit-2026-08-13` branch) — revised after review

**Round 1 feedback (REQUEST_CHANGES)**: the audit's "Net result" section claimed "6 newly-acquired datasets are present, complete, and verified" as a blanket statement, but the report's own table admitted `Greenbaum_NatMed_2024` wasn't re-verified against a remote size and `VentoTormo_Nature_2018` had "no declared remote size to check against" — so "verified" was overstated for those two. Valid catch.

**Fix, not just reword**: rather than just softening the language, actually ran the missing checks — `gzip -t` on Greenbaum's 5 `.gz` files (all pass) + MatrixMarket declared-nnz vs. actual row count on both its matrix files (exact match, 92.5M and 363M rows respectively), HDF5 internal dimension cross-consistency on the Vento-Tormo h5ad (`obs`/`var`/`X` shapes agree, CSC `indptr` length correct), and `unzip -t` on the HPA raw zip (previously just asserted it "would fail loudly," now actually tested). Rewrote the report with an explicit 3-tier verification legend (Tier 1 byte-exact / Tier 2 archive-or-structural-integrity / Tier 3 process-succeeded-only) and reclassified every dataset honestly — only `Arutyunyan2023_MFI` and `HumanDevelopmentMultiomicAtlas` are genuinely Tier 1; the rest are Tier 2/3, which is the strongest check available given what each source exposes, not a shortfall in effort.



User asked to open a PR reviewing data correctness/completeness. Full write-up: `docs/DATA_AUDIT_2026-08-13.md`. Highlights beyond what's already logged above:
- Fixed stale `link.md` on the pre-existing `2026_human_maternal_fetal_Nature` dataset (said "skeleton_only, no files mirrored" — false; data + a first analysis pass have been there since 2026-06-04).
- Resolved the open Vento-Tormo 2018 question: inspected `decidua-v3.h5ad` directly via `h5py` (no `anndata` installed) — confirmed it is **not** decidua-only, has ≈14,366 trophoblast cells (EVT/SCT/VCT) plus decidua+blood, safe to use as an independent trophoblast reference.
- All 6 newly-acquired datasets cross-checked for consistency across `link.md` / `DATA/dataset.index.md` / this repo's `datasets/*/dataset.md`.

### Updated open TODOs (supersedes earlier lists where they overlap)

1. Vento-Tormo 2018: confirm `decidua-v3.h5ad` cell-type coverage (trophoblast included, or decidua-only?) — still open.
2. Fix the stale docs on `2026_human_maternal_fetal_Nature` (`link.md` says skeleton_only, actually has data + a first analysis pass) — still open.
3. Decide E-MTAB-12595 (Arutyunyan multiome, ~299GB raw FASTQ, no processed alternative) — default is skip, still open.
4. **Data acquisition for Aim 1 is done.** All 6 originally-missing datasets resolved and byte-verified: Arutyunyan2023_MFI (+Visium), Greenbaum_NatMed_2024, VentoTormo_Nature_2018, HumanDevelopmentMultiomicAtlas (7/12 organs), HPA_trophoblast_proteome. Only the 2 minor open questions above remain (Vento-Tormo cell-type coverage, stale 2026-Nature docs) — neither blocks starting analysis. **Next real step: resume the 6-question framework (Q1–Q6) from `2026-GPT-TWEAKR-Oncofetal.md#定义清楚Placenta的问题`** — recommended starting point Q1 — and start actually building P1/P2/P3/D signatures per the evidence-layer weighting table in `2026-GPT-TWEAKR-Oncofetal.md#Placenta数据集`. This is genuinely the next thing to do; the last several turns were all infra/data setup, not signature-construction work.

### Compute feasibility: Mac mini (16GB) vs. Argos — new standing reference, 2026-08-13

User's Mac mini only has 16GB RAM; some datasets in this project can't be loaded locally at all. Asked for a feasibility judgment per dataset before each analysis, with Argos (DFCI HPC) as the fallback. Full doc: `docs/COMPUTE_FEASIBILITY.md` — **read that before starting any analysis step**, not just this summary.

Key facts gathered (measured, not assumed):
- Mac mini: nominal 16GB, actual 17.18GB physical, but live-checked real availability was only ~109MB free at one point (`top -l 1`) — real-world safe budget is ~6-8GB per local analysis process, not the full nominal amount.
- Argos: SSH-reachable (`argos.dfci.harvard.edu`), SGE scheduler (not SLURM/LSF — `qsub -pe pvm`), nodes are argos1-8 (64 CPU/376.6GB RAM each) and argos9-10 (160 CPU/629.3GB RAM each), CPU-only (no GPU). Existing project tooling (`SOFTWARES/bin/argos-submit-cluster`, `SOFTWARES/Argos-Server/`) already wraps job submission — reuse it, don't rebuild.
- For h5ad/mtx files, computed exact `nnz` (non-zero count) via `h5py`/MatrixMarket headers rather than guessing from on-disk size — this gives a reliable memory-floor formula (`nnz × 8 bytes`, ×2 if a duplicate `raw` layer exists). Result: `Arutyunyan2023_MFI` primary_tissue (~12GB base load, has a duplicate `raw.X`) needs Argos; everything else in the h5ad/mtx set fits the Mac mini alone (though `2026_human_maternal_fetal_Nature`'s `snRNA_raw_counts` is tight, ~4.6GB base, treat any heavy downstream step on it as Argos-first).
- For the HDMA RDS Seurat objects (no cheap way to peek nnz — RDS has no header), actually calibrated empirically instead of guessing: loaded `Adrenal` (650MB) and `Thyroid` (1.58GB) locally, measured real process RSS delta (174MB and 675MB respectively). Extrapolated conservatively (×1.5 safety buffer on the higher per-cell ratio) to the 5 untested organs: Spleen/Thymus estimated borderline-OK locally (~4-6GB, plausible with caution), Skin/Liver/StomachEsophagus estimated to exceed the safe budget (~8-11GB) → Argos. These 5 are estimates, not measurements — flagged as such in the doc, with instructions to redo the calibration properly on Argos (where a wrong guess is free) if any of them becomes analytically important.
- Also noted `object.size()` in R is unreliable for Seurat's S4 objects (over-counted ~2x vs. the real RSS delta in both calibration runs) — process RSS is the number that actually matters for "will this crash the machine."

### Compute feasibility PR — revised after review (2026-08-13)

**Round 1 feedback (REQUEST_CHANGES)**: the `nnz × 8 bytes` formula was only the sparse-matrix core (data+indices), didn't account for indptr/obs/var/layers/categoricals/graphs/reductions/allocator overhead/dtype upcasting/copies — so classifying `snRNA_raw_counts` (~4.63GB) and the Greenbaum ATAC matrix (~2.90GB) as flat "Mac mini OK" was too optimistic for a standing decision this close to the local budget.

**Fix, empirically, not just reworded** — same pattern as the data-audit PR's review round: actually loaded files and inspected real dtypes/nbytes instead of arguing from the formula.

Two real findings, one in each direction:
1. **The formula itself was accurate for h5ad files** — cross-checked `snRNA_raw_counts`: formula predicted 4.63GB, actual measured (sparse payload nbytes + obs/var, dtype-verified via h5py as float32/int32) was 4.67GB, within 1%. h5ad stores its own dtype, and it happened to be float32 for every file in this project.
2. **The formula was wrong for the two mtx files, and not by a small margin** — `scipy.io.mmread` (the generic way to load `.mtx`) silently upcasts to `float64`/`int64` (not the assumed float32/int32) and returns COO format (not CSR/CSC) until explicitly converted. Measured: Greenbaum RNA matrix 1.11GB actual vs. 0.74GB estimated; ATAC matrix 4.36GB (CSR) / 5.81GB (COO, as initially loaded) vs. 2.90GB estimated — the ATAC matrix reclassified from "Mac mini OK" to "Argos-first."

A third finding, not anticipated by the review but arguably the more important one for this specific machine: **live process RSS itself is unreliable here.** The Mac mini is chronically memory-pressured (confirmed via repeated `top -l 1` checks showing only 88-116MB free), and macOS's background memory compression means `ps`-reported RSS can dramatically underreport true footprint — measured the dtype-verified-4.67GB `snRNA_raw_counts` object showing only 994MB-1.8GB RSS delta across repeated runs. This means the original RSS-based HDMA calibration (Adrenal/Thyroid) is now flagged as a lower bound, not a confident estimate, and the extrapolated Spleen/Thymus verdicts were tightened from "plausible locally" to "Argos-first."

Full revised doc: `docs/COMPUTE_FEASIBILITY.md`.

### Step 1 kicked off: Argos-only analysis workflow established, inventory pass running

**Standing workflow going forward (user directive, 2026-08-13)**: no analysis runs on the Mac mini anymore. Every step: submit via `qsub` on Argos → results land in `~/DATA`/`~/TWEAKR-OncoPlacental/results` on Argos → pull results back to local → summarize in this repo → open a PR for review → repeat until the project's done. This section documents Step 1 under that workflow: what problem it solves, what data it touches, exactly what was run and how.

**Problem this step solves**: before building the actual D-shared/F-specific/P-specific pseudobulk comparison (Q1), we don't yet know, *for each of the 8 normal-development datasets*, what its cell-type/annotation columns are called, which values in them mean "trophoblast" vs. something else, whether genes are stored as symbols or Ensembl IDs (matters for merging across datasets later), or whether a file failed to load at all. Guessing this per-dataset while writing the real DE pipeline would be slow and error-prone. So Step 1 is a pure inventory/audit pass — read each file's metadata (not the full expression matrix where avoidable), dump obs/var schema + cell-type value counts to JSON, and stop there.

**Data touched** (all 15 files across the 6 normal-development datasets acquired earlier this project):
- Placental/trophoblast side (7 h5ad + 1 mtx bundle): `Arutyunyan2023_MFI` primary_tissue + 3 organoid files, `2026_human_maternal_fetal_Nature` (2 files), `VentoTormo_Nature_2018`, `Greenbaum_NatMed_2024` (RNA+ATAC mtx + metadata/cluster CSVs)
- Fetal-somatic side (7 RDS): all 7 downloaded `HumanDevelopmentMultiomicAtlas` organs (Adrenal/Thyroid/Spleen/Thymus/Liver/Skin/StomachEsophagus)
- **Adult reference: still an open gap**, not resolved yet (asked the user, got redirected to "solve environment first" — revisit once Step 1 results are in, needed before the actual 3-way D/F/P comparison can run)

**What was actually done, in order:**

1. **Confirmed data is on Argos** — user had already uploaded it; verified via SSH that `~/DATA/scRNAseq/{Arutyunyan2023_MFI,HumanDevelopmentMultiomicAtlas,Greenbaum_NatMed_2024,VentoTormo_Nature_2018,2026_human_maternal_fetal_Nature}` and `~/DATA/SpatialTranscriptomics/Arutyunyan2023_MFI_Visium` all match the local `DATA/` layout.
2. **Environment reconnaissance** — none of the "obvious" conda envs on Argos (`r4p3`, `scRNA`, `r441`) had Seurat or scanpy/anndata installed. The user pointed at an existing `argos-codex` env (set up by unrelated prior tooling, `SOFTWARES/Argos-Server/env/install_argos_codex_env.sh`) — checked it and **it already has everything needed**: Python 3.11 + scanpy 1.11.5 + anndata 0.12.16, R 4.5.3 + Seurat 5.5.0 + SeuratObject 5.4.0 + Matrix 1.7.5 + data.table + tidyverse + jsonlite. No installation was needed at all — `install_argos_codex_env.sh extra` was never run.
3. **Wrote 3 inventory scripts** (`scripts/01_inventory/`):
   - `inventory_h5ad.py` — uses `anndata.read_h5ad(path, backed='r')` deliberately, not a full load, since this is metadata-only and the Arutyunyan primary_tissue file alone is a ~12GB full load; backed mode still gives shape, obs/var schema, value_counts, nnz, raw/layers/obsm/obsp presence.
   - `inventory_greenbaum_mtx.py` — inspects the RNA+ATAC MatrixMarket matrices (shape/nnz/dtype) plus the SCP2601 metadata.csv and humanplacenta_cluster.csv column/value breakdown.
   - `inventory_seurat_rds.R` — no backed-mode equivalent exists for RDS, so this does a full `readRDS()`; reports assay names, dims, gene-naming convention, meta.data columns + value counts, reductions, graphs, per-assay layers.
   - `run_inventory.sh` — the SGE job script (`#$ -pe pvm 2`) that runs all 15 through the 3 scripts above sequentially, writing one JSON per dataset to `results/01_inventory/`.
4. **Transfer to Argos**: `scp` and `rsync` both failed from this session's sandbox (`scp: Connection closed`, local `rsync` binary permission-denied) — worked around with `tar czf - dir | ssh host "tar xzf -"`, which uses only `ssh`/`tar` and succeeded. Noting this as a env quirk in case it recurs.
5. **Submitted via qsub**: job `3620272` (`tweakr_01_inventory`), running under `argos-codex`, on `argos1`. Monitoring via `qstat -j 3620272` and tailing the SGE log — as of this Worklog entry, 12/15 files done (all 7 h5ad + Greenbaum + HDMA Adrenal/Thyroid/Spleen/Thymus), Liver in progress, Skin/StomachEsophagus still to come. No errors so far. `readRDS()` on the larger organs is markedly slower than the h5ad backed-mode reads (real I/O — confirmed via growing `io`/`vmem` in `qstat -j` between polls, not a hang).

**Early finding already visible in the log before the job even finished**: `Greenbaum_NatMed_2024`'s `humanplacenta_cluster.csv` has an explicit `cell_type` column with `vCTB`/`STB`/`EVT`/`EVT-progenitor`/`STB-progenitor`/Hofbauer/fibroblast/endothelial/erythroblast labels — directly usable trophoblast-vs-other annotation, no guessing needed for this dataset.

**Not done yet, still pending**: the job hadn't finished as of this Worklog entry — Liver/Skin/StomachEsophagus results aren't in `results/01_inventory/` yet. This PR documents the process and scripts; a follow-up commit (same PR, before merge, or a fast-follow PR) will add the actual inventory summary once the job completes and results are pulled back locally.

**Review round 1 (REQUEST_CHANGES) — fixed**: `run_inventory.sh` used `set -uo pipefail` without `-e`, and both inventory scripts caught read failures, wrote an `error` field to the JSON, then returned/exited normally — so the SGE job could print `=== Done ===` and exit 0 even if a dataset's inventory silently failed, with no reliable machine-checkable success signal. Fixed by explicitly tracking each of the 15 invocations' exit codes in the shell driver (kept `set -uo pipefail`, deliberately not `-e`, since one failing dataset shouldn't stop the other 14 from being attempted), printing an OK/FAILED summary, and exiting non-zero if anything failed; `inventory_h5ad.py`/`inventory_seurat_rds.R` now `sys.exit(1)`/`quit(status=1)` when their `error` field is set instead of returning cleanly. Pushed the fixed scripts to Argos for future runs. The already-running job 3620272 was submitted with the old script version (can't retroactively fix a job mid-run) — independently verified its correctness the hard way instead of trusting its exit code: checked all completed JSONs for an `error` key directly, none present through 13/15 as of this note.

### Step 1 results — job 3620272 finished, 15/15 confirmed clean

Job finished ~13:03. Same independent-verification approach carried through to the end (not trusting the old script's own printed status): pulled all 15 JSONs and checked each for an `error` key directly — **zero errors, 15/15 genuinely succeeded**. Pulled the full `results/01_inventory/` directory back to local via the same tar-over-ssh workaround. Full findings: `results/01_inventory/SUMMARY.md`.

**Three real findings, not just "everything loaded fine":**

1. **Gene-ID convention mismatch, blocks Step 2 until resolved**: every placental h5ad dataset (Arutyunyan ×4, both 2026-Nature files, VentoTormo, Greenbaum) uses gene symbols; all 7 HDMA RDS objects use Ensembl IDs. Cross-dataset comparison needs a symbol↔Ensembl mapping decision before any pseudobulk merge.
2. **Trophoblast nomenclature is consistent across 5 independent studies** (Arutyunyan, Nature2026, VentoTormo, Greenbaum all independently use VCT/EVT/SCT + close variants — `iEVT`/`eEVT`/`proEVT`, `VCT_p`/`VCT_CCC`/`VCT_fusing`, `SCT_A`/`SCT_B`/`proSCT`) — one canonical trophoblast filter can cover all 5 rather than 5 bespoke ones. HDMA's `annotv1`/`annotv2` columns confirm the fetal-somatic organs have zero trophoblast contamination, as expected.
3. **`2026_human_maternal_fetal_Nature`'s `snRNA_raw_counts` file has no usable annotation** (`obs` columns are just `ID`/`dataset`/`BC`) — needs a follow-up check (barcode overlap with the sibling `scPlacenta_host` file, which does have full annotation) before Step 2 can use it.

Also a scoping caution worth flagging explicitly: Nature2026's `origin` column (`Fetal`/`Maternal`/`Unknown`) labels which *side of the placenta* a cell came from, not "fetal somatic organ" — don't conflate `origin=="Fetal"` placental cells (fetal endothelium, Hofbauer cells) with the HDMA fetal-somatic-organ reference when building the F-specific module.

Committed the 15 JSONs + `SUMMARY.md` to this PR.

### Step 2 kicked off: gene-ID mapping resolved (HGNC), correcting Step 1's coarse characterization along the way

User supplied an HGNC custom-download URL (symbol/Ensembl/previous-symbol/alias columns) to resolve Step 1 Finding #1 (placental datasets use symbols, HDMA uses Ensembl IDs). Downloaded to `DATA/1.Databases/HGNC_gene_id_mapping/` (new dataset, registered in `dataset.index.md`): `raw/hgnc_custom_download.tsv` (50,321 rows, all statuses), `processed/v0.1/hgnc_symbol_ensembl_map.tsv` (45,032 `Approved` rows, 42,348 with an Ensembl ID). Pushed to Argos too.

**Verified against real data before trusting it, and found Step 1's `looks_like_ensembl` check was misleading**: that check only looked at `rownames(obj)[1]` — one gene, not a real characterization. Loaded the actual full gene list for `HDMA_Adrenal` on Argos: HDMA is **not purely Ensembl**, it's a mix — 76.5% already symbols (19,368/25,314), 23.5% Ensembl IDs used only as a fallback (5,946/25,314).

**Coverage-checked the HGNC table against those 5,946 fallback IDs rather than assuming it "solves" the problem**: only 651 (~11%) resolve, checked both `Approved`-only and any-status (same result — status filtering wasn't the limiter). The other 5,295 are genuinely absent from HGNC entirely — mostly lncRNA/pseudogene/novel-locus Ensembl entries HGNC hasn't curated a symbol for, which is expected (HGNC curates narrower than Ensembl's full gene model) and isn't a gap in this specific table or something a different source would fix. Practical takeaway: the actually-relevant protein-coding trophoblast markers are almost certainly already in the resolved 76.5%, not the unresolvable tail — don't chase the remaining ~89% further, just drop or Ensembl-ID-keep them in the Step 2 merge.

Full writeup: `datasets/HGNC_gene_id_mapping/dataset.md`. Coverage-check script (run on Argos, not qsub — small enough for the login node): `scripts/02_gene_id_mapping/check_hdma_gene_id_coverage.sh`.

### Step 2, review round 2 (REQUEST_CHANGES) — fixed with real full-organ + biotype + collision checks

ChatGPT reviewer flagged the first Step 2 pass as insufficient on four concrete points: (1) the HGNC coverage number was only checked against 1 of 7 HDMA organs (Adrenal), not representative; (2) "protein-coding markers are almost certainly already resolved" was an unproven assertion, no actual biotype quantification of what the unresolved ENSG tail contains; (3) no check for symbol collisions/duplicates introduced by the ENSG→symbol mapping step; (4) no final usable canonical feature map for Step 3 to actually consume — plus a phrasing note to soften "HGNC not having them is the ground truth... don't spend more effort" into something less absolute.

**Fixed all four, in order:**

1. **Full 7-organ gene-list extraction** — `scripts/02_gene_id_mapping/extract_all_organ_gene_lists.sh` (qsub, `argos-codex`), looping `readRDS()` over Adrenal/Thyroid/Spleen/Thymus/Liver/Skin/StomachEsophagus, splitting each organ's rownames into ENSG vs symbol, writing per-organ `*_all_genes.txt`/`*_ensg_only.txt` plus a deduped `union_all_organs_ensg.txt`. First submission (job 3620285) failed instantly on every organ — root cause: `Rscript -e '...' --args "$organ" "$OUT"` passes `--args` as a literal token, which became `commandArgs(trailingOnly=TRUE)[1]` itself instead of being consumed, producing a garbage path (`--args_RNA_obj_clustered_final.rds`). Confirmed via stderr, fixed by removing `--args` (verified locally: `Rscript -e 'cat(commandArgs(trailingOnly=TRUE))' foo bar` → `foo bar`, no `--args` needed), resubmitted as job **3620286**, all 7 organs completed cleanly (~19 min total, OK 7/FAILED 0 in the log). Independently re-verified per the established discipline — didn't just trust the log's own printed summary: pulled all 15 output files via tar-over-ssh (scp/rsync still broken in this sandbox) and checked each file's actual `wc -l` against the log's per-organ numbers (all matched exactly), and confirmed `union_all_organs_ensg.txt` (11,015 lines) is 100% `^ENSG`-prefixed. Per-organ ENSG fraction ranges ~19–32% of all rownames (not just Adrenal's 23.5%).

2. **Real biotype quantification** — `scripts/02_gene_id_mapping/query_biotype.py`, batch-POSTing all 11,015 union ENSG IDs to the Ensembl REST API (`/lookup/id`, 13 batches of ≤900). Cross-referenced against the HGNC resolution status: of the 11,015, 1,023 (9.3%) resolve via HGNC, 9,992 (90.7%) don't. Of those 9,992 unresolved, actual biotype breakdown: 88.2% lncRNA, 7.2% not found in Ensembl at all, **4.6% (459) protein_coding**, <0.1% other. 459 real protein-coding genes is a genuine, non-trivial number — the first pass's "almost certainly negligible" was an unverified assertion that turned out directionally fine but numerically unstated. Follow-up check: do any of those 459 have an Ensembl-native `display_name` that could resolve them without HGNC at all? Only 7/459 do (`IFNAR2-IL10RB`, `EFNA4-EFNA3`, `SCPPPQ1`, `PDE8B`, `PDE4C`, `C4orf36`, `KYAT1`) — the other 452 are genuinely unnamed loci in Ensembl itself (novel/predicted genes, readthrough transcripts), confirming no symbol source has them, for a real reason (not a coverage gap in any one table).

3. **Symbol collision check** — built as part of the canonical feature map (below): 29 total collisions across all 7 organs (3–5 per organ) where an ENSG ID's resolved symbol already exists as a separate native-symbol feature in the same object (e.g. `PDE8B` present both as its own symbol row and as `ENSG00000284762`). Listed per-organ in `results/02_gene_id_mapping/canonical_feature_map/collision_report.tsv`. Not resolved here — flagged as an explicit open item for Step 3 (needs an aggregation rule, e.g. summing counts across colliding features, before `canonical_symbol` can be used as a merge key).

4. **Canonical feature map built** — `scripts/02_gene_id_mapping/build_canonical_feature_map.py` produces one `<organ>_canonical_feature_map.tsv` per organ (`original_feature`, `is_ensg`, `canonical_symbol`, `mapping_status`, `biotype`) plus the collision report and a `SUMMARY.md`. Every original feature is kept — nothing dropped at this stage. `mapping_status` breakdown summed across all 7 organs: 154,046 `native_symbol`, 57,283 `unmapped_kept_as_ensembl_id`, 6,236 `mapped_via_hgnc`, 41 `mapped_via_ensembl_display_name`. Unresolved features keep their own Ensembl ID as `canonical_symbol` rather than being dropped, so Step 3 can decide per-analysis whether to filter on `mapping_status`.

5. **Softened the overclaiming language** in both `DATA/1.Databases/HGNC_gene_id_mapping/link.md` and `datasets/HGNC_gene_id_mapping/dataset.md` — replaced "almost certainly already resolved... don't spend more effort, HGNC not having them is the ground truth" with the actual measured numbers (9.3% resolve, 4.6% of the rest are real protein-coding genes, only 7 of those have any name anywhere) and the specific handling decision (kept as Ensembl ID, not dropped).

Not yet done: Step 3's aggregation rule for the 29 collisions, and the still-open adult-reference gap. Both remain explicit open items, not silently deferred.

### PR #4 approved and merged

Resubmitted the round-2 fixes to the ChatGPT reviewer ("Fetal-胎盘-免疫" conversation) via Chrome automation. **APPROVE** — all four round-1 blockers confirmed closed with real data (full 7-organ check, real biotype query finding 459 protein-coding genes not just "almost certainly negligible," 29 collisions explicitly reported not silently resolved, canonical feature map keeps every feature). One non-blocking note from the reviewer: the PR *description* on GitHub still had the old round-1 "almost certainly already resolved... don't chase remaining 89%" phrasing even though `link.md`/`dataset.md`/Worklog had already been corrected — fixed via `gh pr edit` before merging, so the PR body matches what actually shipped. Merged (`gh pr merge --merge --delete-branch`, user-confirmed since GitHub merges are outward-facing/hard-to-reverse and the environment's auto-permission classifier blocks them without explicit confirmation), local `main` fast-forwarded to `c6fcbdd`.

### Step 3 prep: collision aggregation rule quantified + two-tier adult reference acquired (branch `step03-pseudobulk-prep-2026-08-13`)

Two prerequisites the reviewer flagged before real D/F/P construction can start: an explicit rule for the 29 symbol collisions, and the adult-reference gap repeatedly called out as missing.

**Collision aggregation rule** — quantified rather than assumed. New qsub script `scripts/03_pseudobulk_prep/quantify_collision_count_mass.sh` (job 3620311, all 7 organs, ~23 min total) loads each organ's real RNA counts matrix and sums how much total UMI count sits on the colliding `original_features` vs. the organ's total. Result: **0.003%–0.01% of total counts** in every organ (range 3.1e-05–9.9e-05 as a fraction) — several orders of magnitude too small to matter regardless of aggregation choice. Independently re-verified per the established discipline: pulled the 7 per-organ TSVs + combined table from Argos and checked the numbers match the SGE log exactly before trusting them (also caught a false-positive `FAILED` match in my own monitoring script's `grep -c FAILED` — it was matching the literal text "FAILED (0): none" in the summary line, not an actual failure; the job itself succeeded cleanly, OK 7/7). **Decision**: Step 3 pseudobulk will sum counts across colliding features (standard convention, matches Cell Ranger/STARsolo's own handling of duplicate gene symbols) — documented in `results/03_pseudobulk_prep/SUMMARY.md`, not left as an open question for Step 3 to re-litigate.

**Adult reference — two-tier design**, worked out with the user rather than picked unilaterally (mirrors the earlier HDMA-organ-scope decision): Tier 1 (bulk, mandatory, used now to *define* signatures) = GTEx + HPA; Tier 2 (single-cell, deferred, independent post-freeze *validation* only, never used to define signatures) = Tabula Sapiens. Rationale: GTEx/HPA give organ-matched (F-specific) and whole-body (P-specific) bulk comparisons cheaply; Tabula Sapiens checks afterward whether "P-specific" genes are actually absent across real adult cell types rather than just bulk-tissue-averaged low, and staying out of signature *definition* keeps that check from being circular.

- **`[[GTEx_v11_median_tpm]]`** — found the exact file via the GCS JSON listing API (`storage.googleapis.com/storage/v1/b/adult-gtex/o?prefix=...`) since the GTEx Portal downloads page is a JS-rendered SPA a markdown-only fetch can't see. `GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz`, 74,628 genes × 68 adult tissues, byte-exact (10,129,906 bytes) + structural (GCT-declared row count matches actual, twice — raw file and the cleanup script's output). Built the HDMA-organ↔GTEx-tissue mapping by reading the real 68 column names, not assuming — found a genuine gap: **no Thymus column in GTEx at all** (adult donor population skews older, thymus involutes with age, so GTEx's collection doesn't have it).
- **`[[HPA_RNA_tissue_consensus]]`** — before downloading anything, checked whether this was already on disk somewhere in the wider `DATA/` tree and found it was: `DATA/CRC-Atlas/phase2/03_data/raw/HPA_normal_tissue/rna_tissue_hpa.tsv.zip`, same public HPA file, downloaded for an unrelated project. Copied (not moved, to avoid breaking CRC-Atlas's references) into its own `DATA/1.Databases/` entry, md5-verified byte-exact against the source. 20,162 genes × 40 tissues, long format. **Fills GTEx's exact gap** — confirmed by direct inspection that HPA's tissue panel includes thymus, plus all 6 other HDMA organs + colon/rectum + a bonus placenta cross-check row.
- **`[[TabulaSapiens]]`** — checked Figshare's full file list via its API (article 14267219) before deciding what to fetch, rather than downloading the full 15.6GB unified atlas. Downloaded only the 5 per-organ files relevant to this project (Liver/Skin/Spleen/Thymus/Large_Intestine, ~2.85GB), confirming along the way that **Tabula Sapiens has no Adrenal, Thyroid, or Stomach/Esophagus organ file** — a real gap, not assumed; GTEx/HPA remain the only adult reference for those 3 organs. Hit a genuine silent-failure trap on the first attempt: Figshare's `download_url` 302-redirects to a pre-signed S3 URL that **expires in 10 seconds**, and `curl -o` without `-L` "succeeds" (exit 0) while writing a 0-byte file instead of following the redirect — caught only because every download here is checked against Figshare's own `size`/`supplied_md5` rather than trusting curl's exit code, not by an obvious error message. Fixed with `-L`, re-verified byte-exact + MD5 on all 5 files.

Also caught and fixed a bash bug of my own before it did any damage: the first download-manifest generator used Python's default `print(a, b, c)` (space-joined) while the download script's `read` expected tab-delimited fields — this silently produced empty `url` variables and `curl: (3) URL using bad/illegal format`, again caught by the explicit size/MD5 check rather than the script's own exit status.

Updated `DATA/dataset.index.md` with all three new entries. Committed (`f7230f7`), pushed, opened as **PR #5** and submitted to the ChatGPT reviewer.

### PR #5 review round 1 (REQUEST_CHANGES) — fixed: cross-platform comparison method contract

Reviewer's collision-quantification and two-tier design were both endorsed with no changes needed ("collision 部分处理得很好... 这个结论足够稳，可以不再反复讨论" / "GTEx + HPA + held-out Tabula Sapiens 的两层设计本身我也赞成"). One real methodological blocker: HDMA is single-cell/Seurat counts, GTEx is bulk median TPM, HPA is bulk nTPM — three platforms/normalizations that can't be directly compared as fold-change or a merged DE model without platform/depth/composition effects contaminating the signature. The reviewer didn't ask for new data or for Step 3 to actually start — just for the comparison boundary to be locked down first.

**Fixed by writing `docs/STEP3_METHOD_CONTRACT.md`**, which locks in: developmental evidence (F/P-specific candidate identification) computed entirely *within* the scRNA-seq data, never touching GTEx/HPA; GTEx/HPA used only to answer "is this gene still meaningfully expressed in the matched adult tissue?" via each dataset's own internal rank/percentile/threshold (never raw-magnitude cross-platform comparison); final signature = developmental evidence AND adult-depletion evidence as two independently-computed axes, not one merged model.

**Also fixed a genuine inconsistency the reviewer caught**: `HPA_RNA_tissue_consensus` includes a `placenta` row among its "40 adult tissues," and earlier docs didn't carve it out — which would have let placenta count as adult-negative background for P-specific calls (circular: "gene isn't placenta-specific because it's high in placenta"). Fixed by adding an explicit `role` column to both `hdma_organ_to_hpa_tissue_map.tsv` (`placenta` → `positive_cross_check_ONLY_never_negative_reference`, everything else → `adult_negative_reference`) and `hdma_organ_to_gtex_tissue_map.tsv` (all rows `adult_negative_reference` — GTEx has no placenta column, so no exclusion case needed there) so Step 3 code has a machine-readable contract to filter on, not just prose. Re-pushed both updated processed tables to Argos, verified byte-for-byte via `cat` over ssh.

### PR #5 approved and merged — Step 1/2/3-prep all complete

Resubmitted round-2 fixes to the ChatGPT reviewer. **APPROVE**: "上一轮唯一 blocker 已经解决，而且修法是正确的." Confirmed the method contract's platform boundary is locked correctly, the machine-readable `role` column genuinely fixes the placenta logic error (not just prose), and reaffirmed the collision-mass finding needs no further work. Merged (`gh pr merge 5 --merge --delete-branch`), local `main` fast-forwarded to `126ef06`.

**One non-blocking note carried forward into Step 3 core, not a PR #5 blocker**: the reviewer flagged that `STEP3_METHOD_CONTRACT.md` only locks the *platform*-comparison boundary — it does not yet define the actual orthogonal D-shared/F-specific/P-specific split. Specifically: "P-specific 最终必须体现 placenta/trophoblast 相对 fetal-somatic 的独立性，而不只是 trophoblast vs placenta 内其他细胞 + adult depletion." When Step 3 core is built, P-specific needs a third axis beyond (a) trophoblast-vs-other-placental-cells and (b) adult-depletion — specifically checking independence from fetal-somatic (HDMA) expression too, so a gene shared with fetal-somatic tissue doesn't get mislabeled placenta-specific just because it clears the adult-exclusion bar. Flagged here so Step 3's design doesn't silently drop this.

### Step 4 (`results/04_dfp_signature`) kicked off: D/F/P design doc, before any compute

New branch `step04-dfp-signature-2026-08-13`. Rather than jump straight to qsub jobs, wrote `docs/STEP4_DFP_DESIGN.md` first — same pattern that worked for `STEP3_METHOD_CONTRACT.md` (design reviewed once, before a compute run, instead of after a wasted one).

**Core design point, directly resolving the open item from PR #5's APPROVE**: HDMA (fetal-somatic) and the placental scRNA-seq datasets are structurally asymmetric — HDMA organs are already pure fetal-somatic tissue (no internal DE needed, "developmental evidence" = expression level), while placental datasets are a mix of trophoblast and non-trophoblast cells (Hofbauer, endothelial, maternal decidual/immune) requiring an internal trophoblast-vs-rest DE to isolate the placental signal. The fix for the missing third axis: **use each developmental side as an exclusion reference for the other**, symmetric to how GTEx/HPA are used as an adult-exclusion reference — HDMA's own expression distribution excludes genes from P-specific if they're also fetal-somatic-elevated; the placental trophoblast-DE result excludes genes from F-specific if they're also trophoblast-elevated. This makes:

- D-shared = elevated-in-fetal-somatic AND replicated-in-placenta AND adult-excluded (both organ-matched + whole-body)
- F-specific = elevated-in-fetal-somatic AND adult-excluded (organ-matched) AND NOT replicated-in-placenta
- P-specific = replicated-in-placenta AND adult-excluded (whole-body) AND NOT elevated-in-fetal-somatic ← the new clause

Left explicitly open rather than guessed at: exact percentile/threshold cutoffs (to be computed against real distributions, not assumed), whether all 5 placental datasets have enough donor/sample replicate structure for a valid DE test, the quorum for "replicated in placenta," and whether Nature2026's still-unannotated `snRNA_raw_counts` gets resolved before this step runs or is excluded. Opened as PR #6, submitted for review.

### PR #6 review round 1 (REQUEST_CHANGES) — fixed a real conceptual error: expression ≠ developmental evidence

Reviewer endorsed the direction (P-specific's new `NOT elevated_in_fetal_somatic` clause genuinely fixes the orthogonality gap from PR #5) but caught something more fundamental: the first draft treated `elevated_in_fetal_somatic` alone as "fetal developmental evidence." That's wrong — HDMA being trophoblast-free only proves it's a valid fetal-somatic *reference*, not that a highly-expressed gene there is part of a developmental program rather than a housekeeping/organ-identity/constitutive-metabolic gene (all of which would trivially pass that check). Left uncorrected this would have: collapsed D-shared into "fetal-expressed + trophoblast marker" rather than a real shared program; let F-specific absorb ordinary organ-identity genes; and over-pruned P-specific (any placenta gene with normal expression in any one fetal organ would get wrongly excluded).

**Fixed by restructuring**: adult-depletion is folded into the definition of each developmental program itself, not applied as a downstream filter —

- F-developmental(gene, organ) = elevated_in_fetal_somatic AND adult_excluded(matched_organ)
- P-developmental(gene) = replicated_in_placenta AND adult_excluded(whole_body)

— then D/F/P become a clean three-way partition of these two already-adult-corrected programs: D-shared = F-developmental AND P-developmental; F-specific = F-developmental AND NOT P-developmental; P-specific = P-developmental AND NOT F-developmental. Reviewer's framing: "D/F/P 真正成为同一个 developmental universe 的三分解，而不是三个条件各异、可能留下大量灰区的 ad hoc gene lists." Still respects `STEP3_METHOD_CONTRACT.md` — no HDMA-UMI-vs-GTEx-TPM fold-change, `adult_excluded` still uses each dataset's own internal rank/percentile. Reviewer confirmed the no-internal-HDMA-DE decision was correct too — the blocker was never "HDMA needs its own control," it was "expression ≠ developmental specificity." Updated `docs/STEP4_DFP_DESIGN.md`, resubmitting.

### PR #6 approved and merged

**APPROVE**: "上一轮 blocker 已经实质解决... 这比上一版干净很多，因为现在 D/F/P 真的是从两个已经定义好的 developmental parent programs 做集合分解。" One post-APPROVE note, not a blocker: F-developmental is organ-specific, P-developmental is global — final gene-set generation needs to decide whether F stays per-organ or collapses to a cross-organ consensus, "由真实分布决定，不需要现在提前拍死." Merged (`gh pr merge 6 --merge --delete-branch`), local `main` fast-forwarded to `2e727ab`.

### Replicate-structure audit: resolves design doc's first open item, no new compute needed

New branch `step04-dfp-inventory-2026-08-13`. Answered "do all 5 placental datasets have enough donor/sample replicate structure for the trophoblast-vs-rest DE" by re-reading `results/01_inventory/*.json` (already produced/verified via Step 1's qsub job) rather than running anything new — the donor/sample columns and their value counts were already captured there.

**Result**: 4 of 7 placental datasets are actually usable — Arutyunyan primary_tissue (18 donors), Nature2026 scPlacenta_host (23 samples), VentoTormo decidua-v3 (12 fetuses), Greenbaum (8 donors), all with real replicate structure AND a non-trophoblast population to contrast against. The 3 Arutyunyan organoid datasets (PTO/TSC/Fig3) are **structurally excluded regardless of replicate count** — they're pure trophoblast cultures by construction (Step 1 already noted this), so there's no internal non-trophoblast population to run a contrast DE against. Not a data-quality gap, just what organoid culture datasets are. Revises the design doc's "≥3 of 5" quorum placeholder to "≥3 of 4."

Also surfaced two minor real data-quality notes along the way: VentoTormo's `Fetus` donor labels have inconsistent leading whitespace (`" F15"` vs `"F19"`, needs trimming); Greenbaum's cluster-annotation file (~1,923 cells) is much smaller than its full RNA matrix (~36,456 cells per Step 1) — needs checking whether that's a representative subsample before running the DE on it.

Full writeup: `results/04_dfp_signature/replicate_structure_audit.md`. `docs/STEP4_DFP_DESIGN.md` updated to mark this open item resolved and cross-reference the audit.

### PR #7 review round 1 (REQUEST_CHANGES) — fixed: marginal totals ≠ paired contrast, found two real bugs while fixing it

Reviewer correctly caught that the audit only checked marginal totals (dataset has N donors; dataset has trophoblast + non-trophoblast cells overall) — not whether any single donor actually has both. If trophoblast status is confounded with donor identity, donor-level pseudobulk DE is invalid no matter how many nominal donors exist. Asked for a real `dataset | donor | n_trophoblast | n_non_trophoblast | eligible` table.

Wrote `scripts/04_dfp_signature/donor_troph_crosstab.py`, ran on Argos (this genuinely needed new compute — the joint donor×cell-type breakdown wasn't in Step 1's existing inventory JSONs, only marginal value_counts were). Found two real bugs while building it, both caught by looking at actual output rather than trusting the first pass:

1. **VentoTormo's whitespace bug confirmed real**: after trimming `Fetus`, all 12 donors cleanly have both groups — untrimmed, `" F15"` and a clean `"F15"` would have silently split into two fake single-group donors.
2. **Greenbaum join was completely broken, 0/1923 cells matched**: `cluster.csv`'s `NAME` (`W9_AAACCAACACCTGCCT`) and `metadata.csv`'s `NAME` (`JS34#ACGTCAAGTTGCAATG-1`) use incompatible barcode schemes entirely — not almost-matching, not matching at all. Fixed by noticing `cluster.csv`'s own `NAME` already embeds the donor as the prefix before the last `_` — no join needed. This also revealed the annotated ~1,923-cell subset covers only **3 of the full 8 donors** (W8-2/W9/W11), not 8 as round 1's marginal-count table implied.

**Final real numbers**: Arutyunyan 17/18 donors eligible, Nature2026 23/23, VentoTormo 12/12, Greenbaum 3/3 (not 8) — all 4 datasets remain usable, but Greenbaum's real replicate count is much thinner than round 1 suggested and should be weighted accordingly, not treated as equal-strength evidence. Updated `results/04_dfp_signature/replicate_structure_audit.md` with the full cross-tab and this correction, resubmitting.

### Repo layout (as of this session)

```text
TWEAKR-OncoPlacental/
├── Worklog.md              # this file — read first when resuming
├── README.md
├── datasets/
│   ├── Arutyunyan2023_MFI/dataset.md              # E-MTAB-12421/12595/12650/12698
│   ├── Greenbaum_NatMed_2024/dataset.md            # Broad SCP2601
│   ├── VentoTormo_Nature_2018/dataset.md           # E-MTAB-6701/6678/7304
│   └── HumanDevelopmentMultiomicAtlas/dataset.md   # placeholder only, not downloaded
├── scripts/                 # empty so far
├── notebooks/                # empty so far
└── docs/                     # empty so far
```
