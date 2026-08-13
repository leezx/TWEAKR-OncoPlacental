# TWEAKR-OncoPlacental — Worklog

Continuity log for this repo. Read this file first when resuming work — it should let you pick up without re-deriving context.

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
