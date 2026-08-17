# Step 7: CLiM/CLuM external data acquisition — design

New project phase, initiated 2026-08-16 after Step 6 (this project's
original scope) reached 100% completion (`Worklog.md`). Sources a new
external dataset: **"Spatial multi-omics landscape of colorectal cancer
macro- and micrometastases"** (Liu, Jadhav, Pan, et al.; Kopetz, Wang,
Maru, senior/corresponding; *Cancer Cell*, 2026;
DOI [10.1016/j.ccell.2026.06.009](https://doi.org/10.1016/j.ccell.2026.06.009);
PMC [PMC13390483](https://pmc.ncbi.nlm.nih.gov/articles/PMC13390483/)).
The paper profiles 49 tumors from 19 patients (primary CRC + liver
metastases [CLiM] + lung metastases [CLuM]) across 5 modalities, and
separately reuses several existing public cohorts for validation.

## Scope (user-confirmed 2026-08-16)

**In scope**: all public scRNA-seq and bulk RNA-seq/microarray cohorts
this paper cites (both its own newly-generated data where publicly
accessible, and the external cohorts it reuses).

**Explicitly out of scope**: all spatial transcriptomics (Visium ST,
Visium HD), LCM-WGS, and PhenoCycler-Fusion protein imaging — for two
independent reasons: (1) user-confirmed descope this round; (2) the
majority of this paper's own newly-generated data across *every*
modality — including 9 of its 12 Visium HD samples, its 4 in-house
scRNA-seq samples, all LCM-WGS, and possibly its 117 MDACC bulk
RNA-seq samples — sits in a **Zenodo record requiring login/authentication**
(`10.5281/zenodo.19043057`, 257.2GB, license CC-BY-4.0 but access-gated),
not a plain public download. Per this session's standing rule (no
account creation, no credential handling), this repo is not accessible
without the user's own Zenodo access. Not attempted here — flagged as a
distinct future item if the user later provides access.

## Real inventory findings (this round — every accession checked
directly, not assumed from the paper's abstract-level description)

| Dataset | Paper's stated scope | Verified structure | Notes/discrepancies found |
|---|---|---|---|
| **GSE231559** | scRNA-seq, 9 CLiM + 6 primary CRC | 26 samples total, `GSE231559_RAW.tar` (692.1MB), MTX/TSV (10x-style) | Clean match to this project's existing loader pattern (bare MTX/barcodes/features, same as HTAN/CRLM). Raw counts available. |
| **GSE225857** | scRNA-seq, 4 CLiM + 4 primary CRC | 8 samples listed at series level, `GSE225857_RAW.tar` (607.0MB) mixing MTX/TSV **with JPG/PNG/JSON** (spatial-image-adjacent files) | **Real finding**: this series appears to co-package spatial-transcriptomics-related files alongside the scRNA-seq matrices (its own title mentions "single-cell and spatial transcriptome analysis"). Must isolate the scRNA-seq-only files during download/inventory — not assumed to be scRNA-seq-only just because the paper cites it that way. Raw sequencing unavailable (patient privacy); only processed matrices provided (per GEO's own notice), deposited separately in CNGB for raw. |
| **GSE285991** (SuperSeries) | scRNA-seq, 10 liver metastases | **Not scRNA-seq at series level.** Two subseries: `GSE285989` (CUT&Tag/ATAC-seq, mouse+human) and `GSE285990` (RNA-seq, mouse+human). `GSE285990`'s 10 human samples (`P01_LM`–`P10_LM`) are labeled only as "RNA-seq" — GEO's own record does not confirm these are single-cell rather than bulk/pseudobulk, and explicitly states **raw sequence data for the human samples was withheld from GEO for patient-privacy reasons** (only processed expression matrices deposited). | **Real, unresolved discrepancy**: the paper describes this as a scRNA-seq source; GEO's own metadata does not confirm single-cell resolution for the human arm and may in fact be bulk. **Must be resolved by directly inspecting the downloaded processed-matrix file structure** (a true per-cell matrix will have far more "samples"/barcodes than 10 govern-level entries; a bulk matrix will have exactly ~10 columns) before this dataset is used for anything beyond inventory — flagged as a blocking verification item, not assumed either way. |
| **GSE131418** | bulk microarray, 170 liver-met samples | SuperSeries-scale, **1,135 samples total** (333+545 primary tumors, 184+73 metastases, discovery+validation cohorts combined), CEL files (5.3GB) + processed TXT/CSV | The paper's "170" is a **subset** of this 1,135-sample series (the liver-metastasis-site subset specifically) — needs to be identified via the series' own clinical/site metadata (sample characteristics fields) during inventory, not assumed to be the whole series. |
| **GSE17538** | bulk microarray, primary CRC, n=232 | SuperSeries, 244 samples total (human **and mouse** mixed, 2 platforms: `GPL570` human, `GPL1261` mouse), CEL files (1.8GB) | Paper's 232 is the human-only subset of the 244; mouse samples must be excluded during inventory via the platform/organism field. |
| **GSE21510** | bulk microarray, primary CRC, n=146 | 148 samples total, `GPL570`, CEL files (735.4MB) | Close match (148 vs. paper's 146) — likely 2 excluded for QC in the paper; not investigated further, reported as-is. |
| **TCGA-CRC** | bulk RNA-seq, n=610 | Not a GEO accession — public, open-access via the [GDC Data Portal](https://portal.gdc.cancer.gov/) (`TCGA-COAD` + `TCGA-READ` projects combined). Gene-level expression counts/FPKM are open-tier (no dbGaP application needed); only raw BAMs/germline variants are controlled-access. Typical combined COADREAD RNA-seq sample count in public reporting is ~570–630 depending on exact QC filter — paper's 610 is within this range, not independently re-derived here. | Access method differs from every other cohort in this table (GDC API/`gdc-client`, not a GEO `_RAW.tar`) — noted explicitly so the download script doesn't assume one uniform mechanism. |

## What "preliminary analysis" means for this PR (scoped explicitly)

Matching this project's own Step 1 precedent (`docs/` — the original
data-inventory phase that preceded any scoring compute by several PRs):
**this PR is inventory only** — real download, byte-verified, structural
characterization (sample/cell counts, gene-ID format, raw-vs-processed
status, cell-type annotations if any, resolution of the GSE225857
spatial-file and GSE285990 bulk-vs-single-cell discrepancies found
above). **Not in scope for this PR**: running the existing D/F/P/revCSC
scoring pipeline against any of these cohorts — that is a natural,
directly-analogous follow-on (mirroring exactly how `HTAN_CRC_
progressive_plasticity`/`CRLM_NMP_ATLAS` were inventoried in the
original Step 6 design round before any scoring ran), but is
deliberately deferred to a separate PR once this inventory is reviewed,
to keep this PR's scope checkable in one pass rather than conflating
"did the download work" with "does the scoring pipeline apply cleanly."

## Download mechanics (per this project's standing discipline)

- Runs on Argos, `argos-qsub1`, network download (not qsub-wrapped —
  I/O-bound, not compute-bound, same precedent as this project's
  original HDMA/Arutyunyan acquisition phase).
- Every download verified by **exact byte size** against the source's
  declared size (GEO's own reported file sizes above; GDC's own
  manifest-reported sizes for TCGA) — never piped through `tail`/`head`
  (the standing `curl_pipe_swallows_exit_code` lesson from this
  project's own history: a pipe silently swallows `curl`'s real exit
  code, so a truncated download can report success).
- Archive integrity additionally checked via `tar -tzf`/`gzip -t` where
  applicable, matching the same precedent.
- New directory structure: `results/07_clim_external_data/` (metadata,
  inventory reports) + a `DATA/` staging area on Argos (not committed to
  git — matching this project's existing convention that raw data lives
  under `/home/zz950/DATA/`, only small inventory/summary tables get
  committed).

## What this design does not do

Does not download anything from the Zenodo-restricted record (§Scope
above). Does not attempt to resolve the GSE225857/GSE285990
discrepancies here — they are locked as **required verification items**
during the actual download/inventory compute, not assumed away. Does
not run any scoring/analysis beyond structural inventory. Does not touch
this project's already-frozen D/F/P/revCSC gene sets or scoring
machinery (would be reused as-is in a later follow-on PR, not
re-derived).
