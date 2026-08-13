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

### Updated open TODOs (supersedes earlier lists where they overlap)

1. **Confirm the StomachEsophagus re-download finished and matches 19,915,325,339 bytes exactly** — was in progress as of this update.
2. Vento-Tormo 2018: confirm `decidua-v3.h5ad` cell-type coverage (trophoblast included, or decidua-only?) — still open.
3. Fix the stale docs on `2026_human_maternal_fetal_Nature` (`link.md` says skeleton_only, actually has data + a first analysis pass) — still open.
4. Decide E-MTAB-12595 (Arutyunyan multiome, ~299GB raw FASTQ, no processed alternative) — default is skip, still open.
5. **Data acquisition for Aim 1 is essentially done** (6/6 originally-missing datasets resolved, modulo the StomachEsophagus re-download and the 2 minor open questions above). Next real step: resume the 6-question framework (Q1–Q6) from `2026-GPT-TWEAKR-Oncofetal.md#定义清楚Placenta的问题` — recommended starting point Q1 — and start actually building P1/P2/P3/D signatures per the evidence-layer weighting table in `2026-GPT-TWEAKR-Oncofetal.md#Placenta数据集`.

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
