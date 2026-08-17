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
| **GSE231559** | scRNA-seq, 9 CLiM + 6 primary CRC | 26 samples total, `GSE231559_RAW.tar` (692.1MB declared), MTX/TSV (10x-style) | Clean match to this project's existing loader pattern (bare MTX/barcodes/features, same as HTAN/CRLM). Raw counts available. **Round-1 review correction**: technical loadability is not the same as matching the paper's cited 9+6 subset — reproducing that exact 15-sample cohort from the 26-sample series (via patient-ID/tissue-site metadata) is a required inventory acceptance criterion, not assumed automatic. |
| **GSE225857** | scRNA-seq, 4 CLiM + 4 primary CRC | **Round-1 review correction (real, sharper finding than originally reported)**: at GSM level, this series contains only **2 genuine scRNA-seq samples** — `GSM7058754` ("immune cells") and `GSM7058755` ("nonimmune cells"), each an aggregated pool across all patients (41,892 CD45− + 196,473 CD45+ cells total), **not 8 per-patient samples**. The remaining 6 GSMs (`GSM7058756`–`GSM7058761`) are separate spatial-transcriptomics samples (4 primary + 2 CLiM, FFPE) that also carry `barcodes.tsv`/`features.tsv`/`matrix.mtx`-shaped files, not just images — so file-extension/MTX-shape alone cannot distinguish scRNA-seq from spatial in this series. | Downloader **must select by exact GSM accession** (`GSM7058754`/`GSM7058755` only), never by file extension or "looks 10x-shaped." Recovering per-patient (4 CLiM + 4 primary) breakdown from these 2 pooled immune/non-immune files, if possible at all, requires the cells' own per-cell patient-of-origin metadata. **Round-2 review correction**: this metadata is not "inside the matrix" — each scRNA-seq GSM deposits it as a **separate accompanying file** (confirmed directly: `GSM7058754` carries `GSM7058754_immune_counts.txt.gz` (213.9MB, the counts matrix) *and* `GSM7058754_immune_meta.txt.gz` (9.6MB, separate per-cell metadata); `GSM7058755` follows the same two-file pattern) — a required inventory verification item, not assumed to exist. Raw sequencing unavailable (patient privacy) for the scRNA-seq GSMs; only processed matrices provided. |
| **GSE285991** (SuperSeries) | scRNA-seq, 10 liver metastases | Two subseries: `GSE285989` (CUT&Tag/ATAC-seq) and `GSE285990` (RNA-seq). **Round-1 review correction**: the first draft called `GSE285990`'s 10 human samples (`P01_LM`–`P10_LM`) "genuinely unresolved" single-cell-vs-bulk — **false, confirmed directly by fetching the individual GSM record** (`GSM8714595`/P01_LM): `Library source: transcriptomic single cell`, supplementary files are genuine per-cell `barcodes.tsv.gz`/`features.tsv.gz`/`matrix.mtx.gz` (82.3MB matrix for this one sample alone — far too large to be a 10-column bulk matrix). **These are confirmed real single-cell data**, not an open question. **Round-2 review correction**: the same GSM record's platform/chemistry attribution is internally inconsistent within GEO's own metadata — its `Description` field says `BD Rhapsody` while its `Extraction protocol` field says libraries were "prepared using the Chromium Single Cell 3' Reagent Kit v2." The single-cell-vs-bulk determination does not depend on resolving this (it rests on `Library source: transcriptomic single cell` plus the genuine per-cell barcode/feature/matrix files, independent of which chemistry produced them); the specific chemistry/platform is left unresolved rather than asserted as Chromium. | Corrected acceptance criterion: verify matrix integrity, real per-cell dimensions (thousands of barcodes, not ~10), and the 10-sample→per-patient mapping — a QA step, not a "decide if this is single-cell at all" step (that question is already answered). Human raw sequence still withheld from GEO for privacy (processed matrices only, as originally noted — this part was correct). |
| **GSE131418** | bulk microarray, 170 liver-met samples | SuperSeries-scale, **1,135 samples total** (333+545 primary tumors, 184+73 metastases, discovery+validation cohorts combined), CEL files (5.3GB declared) + processed TXT/CSV | The paper's "170" is a **subset** of this 1,135-sample series. **Round-1 review correction (strengthened acceptance criterion)**: identifying "liver metastasis" by site keyword alone is not sufficient — the inventory must reproduce the paper's specific 170-sample cohort together with the treatment/status metadata the paper actually analyzes, using the series' own clinical annotation fields, not a keyword match. |
| **GSE17538** | bulk microarray, primary CRC, n=232 | **Round-1 review correction (real subsetting error, not just wording)**: the first draft assumed 232 = 244-total minus mouse (238) — **wrong, confirmed directly**: `GSE17538` is a SuperSeries of 4 subseries — `GSE17536` (177 human CRC, Moffitt) + `GSE17537` (55 human CRC, Vanderbilt, confirmed exact sample count directly) = **232 exactly**, matching the paper's cited n; `GSE19072` (human colonic adenomas) and `GSE19073` (6 mouse) are excluded entirely, not merely filtered by organism. | Corrected reconstruction: download and use `GSE17536` + `GSE17537` specifically (not the SuperSeries' full sample set, and not a naive organism filter on the combined 244). |
| **GSE21510** | bulk microarray, primary CRC, n=146 | 148 samples total, `GPL570`, CEL files (735.4MB declared) | **Round-1 review correction**: the first draft's "likely 2 excluded for QC, not investigated further" was itself an unverified assumption — exactly the kind of thing this project's discipline requires checking, not guessing. Corrected: stated as an **unresolved 148→146 discrepancy** — the exact 2 excluded samples and the exclusion basis must be identified from the series' own clinical/QC metadata during inventory, not assumed. Downloading all 148 is fine; declaring the paper's cohort reconstructed is not, until this is resolved. |
| **TCGA-CRC** | bulk RNA-seq, n=610 | Not a GEO accession — public, open-access via the [GDC Data Portal](https://portal.gdc.cancer.gov/) (`TCGA-COAD` + `TCGA-READ` projects combined). Gene-level expression counts/FPKM are open-tier (no dbGaP application needed); only raw BAMs/germline variants are controlled-access. | Access method differs from every other cohort in this table (GDC API/`gdc-client`, not a GEO `_RAW.tar`) — noted explicitly so the download script doesn't assume one uniform mechanism. **Round-1 review correction**: "610 is within the typical ~570–630 public range" is not itself a reconstruction of the paper's cohort — inventory must establish which specific sample/aliquot/sample-type filter (primary tumor only vs. all sample types, one aliquot per patient, etc.) reproduces the paper's n=610 as closely as achievable, not stop at "plausible magnitude." |

## What "preliminary analysis" means for this PR (scoped explicitly)

Matching this project's own Step 1 precedent (`docs/` — the original
data-inventory phase that preceded any scoring compute by several PRs):
**this PR is inventory only** — real download, byte-verified, structural
characterization (sample/cell counts, gene-ID format, raw-vs-processed
status, cell-type annotations if any). **Round-2 review correction**:
the GSE225857 scRNA-vs-spatial file-type question and the GSE285990
single-cell-vs-bulk question are **already resolved** (both fixed round
1, above) — not open items this PR still needs to "resolve." What
remains genuinely open for this PR's inventory pass is narrower:
GSE225857's per-patient (4 CLiM + 4 primary) reconstruction from the
2 pooled immune/non-immune matrices' accompanying per-cell metadata
files, and GSE285990's matrix-integrity/dimension/sample-mapping QA
(confirming real per-cell structure, not deciding single-cell-vs-bulk).
**Not in scope for this PR**: running the existing D/F/P/revCSC
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
- **Round-1 review correction, byte-verification mechanics**: the first
  draft proposed verifying against GEO's human-readable declared sizes
  above (e.g. "692.1 Mb") — these are rounded UI display values, **not**
  an exact byte oracle, and cannot be used for byte-exact verification.
  Corrected: verification uses the exact integer byte count from
  machine-readable source metadata at download time — the HTTP
  `Content-Length` header (or FTP `SIZE` response) for GEO downloads,
  and the GDC manifest's own reported file size + checksum for TCGA —
  never the GEO webpage's rounded MB/GB display. Never piped through
  `tail`/`head` regardless (the standing `curl_pipe_swallows_exit_code`
  lesson from this project's own history: a pipe silently swallows
  `curl`'s real exit code, so a truncated download can report success).
- **Round-1 review correction, archive-integrity mechanics**: the first
  draft proposed `tar -tzf` for every archive — wrong for GEO's ordinary
  `_RAW.tar` files, which are **plain (non-gzip-compressed) tar archives
  containing individually-gzipped members** (per GEO's own documented
  convention: `tar -xf GSExxxx_RAW.tar` then decompress the `.gz`
  members inside). Corrected: `tar -tf` (not `-tzf`) to verify the outer
  archive's own integrity/listing, then `gzip -t` on each extracted
  `.gz` member — `tar -tzf` reserved only for genuinely gzip-compressed
  tarballs (`.tar.gz`/`.tgz`), which none of these `_RAW.tar` files are.
- New directory structure: `results/07_clim_external_data/` (metadata,
  inventory reports) + a `DATA/` staging area on Argos (not committed to
  git — matching this project's existing convention that raw data lives
  under `/home/zz950/DATA/`, only small inventory/summary tables get
  committed).

## What this design does not do

Does not download anything from the Zenodo-restricted record (§Scope
above). Does not attempt to reconstruct the paper's exact per-cohort
sample counts (GSE231559's 9+6, GSE225857's per-patient CLiM/primary
split, GSE131418's 170, GSE21510's 146, TCGA's 610) here — locked as
**required inventory acceptance criteria** during the actual
download/inventory compute, not assumed satisfied by "the series
downloads and loads." (`GSE285990`'s single-cell-vs-bulk question,
originally listed here as unresolved, is no longer open — confirmed
genuine scRNA-seq directly against its own GSM record, round 1.) Does
not run any scoring/analysis beyond structural inventory. Does not touch
this project's already-frozen D/F/P/revCSC gene sets or scoring
machinery (would be reused as-is in a later follow-on PR, not
re-derived).

## Review history

- **Round 1 (REQUEST_CHANGES — 2 real factual errors plus several
  acceptance-criteria gaps, all independently re-verified against live
  GEO records before fixing, no download performed yet so nothing to
  re-run)**: (1) `GSE285990`'s human samples were called "genuinely
  unresolved" single-cell-vs-bulk — false; fetched the individual GSM
  record (`GSM8714595`/P01_LM) directly, confirmed `Library source:
  transcriptomic single cell`, genuine per-cell
  barcodes/features/matrix files — corrected to state this is confirmed
  scRNA-seq, not an open question (round 2 further found this same GSM
  record's chemistry attribution internally inconsistent — see below —
  so "Chromium kit" is removed from this entry rather than repeated).
  (2) `GSE17538`'s "232 = 244 total
  minus mouse" reconstruction was wrong — confirmed directly that
  `GSE17538` is itself a 4-subseries SuperSeries, and the paper's 232 is
  exactly `GSE17536` (177) + `GSE17537` (55, confirmed exact count) —
  corrected the reconstruction rule. (3) `GSE225857`'s scRNA-seq/spatial
  isolation rule was underspecified — confirmed at GSM level that only
  2 GSMs (`GSM7058754`/`GSM7058755`, pooled immune/non-immune across all
  patients, not 8 per-patient samples) are scRNA-seq, and that the other
  6 GSMs are spatial samples that also carry MTX-shaped files — fixed to
  require exact-GSM-accession selection, not file-shape heuristics.
  (4)/(5)/(6) `GSE21510`'s "likely 2 excluded for QC" was itself an
  unverified assumption — corrected to an explicit unresolved
  discrepancy requiring identification during inventory; added explicit
  paper-cohort-reconstruction acceptance criteria for `GSE231559`,
  `GSE131418`, and TCGA-CRC (previously implied, not stated as a
  requirement). (7) Byte-verification mechanics corrected from GEO's
  rounded MB/GB display values to exact `Content-Length`/manifest byte
  counts. (8) Archive-integrity mechanics corrected from `tar -tzf`
  (wrong — GEO's `_RAW.tar` files are plain tar containing gzipped
  members) to `tar -tf` + per-member `gzip -t`.
- **Round 2 (REQUEST_CHANGES — narrow, 2 residual wording/provenance
  issues, no methodological or scope blocker; both independently
  re-verified against live GEO records before fixing)**: (1) the
  "preliminary analysis" scoping paragraph still described the
  GSE225857 and GSE285990 discrepancies as open items to "resolve" in
  this PR, contradicting the rest of the document (both were already
  fixed in round 1) — corrected to state plainly what's already
  resolved vs. what's still genuinely open (GSE225857's per-patient
  reconstruction; GSE285990's matrix QA). (2) Two overly specific
  metadata claims were corrected after direct re-fetch: `GSE225857`'s
  per-cell patient-of-origin metadata was said to live "inside the
  matrix" — false; confirmed directly that each scRNA-seq GSM deposits
  it as a **separate file** (`GSM7058754_immune_counts.txt.gz` +
  `GSM7058754_immune_meta.txt.gz`, and the analogous pair for
  `GSM7058755`) — corrected. `GSE285990`'s P01 sample was said to be
  "prepared with the Chromium Single Cell 3' Reagent Kit v2" as
  settled fact — confirmed directly that the same GSM record's own
  `Description` field says `BD Rhapsody` while its `Extraction
  protocol` field says Chromium, i.e. GEO's own metadata is internally
  inconsistent on chemistry/platform — corrected to leave chemistry
  unresolved rather than assert Chromium (the single-cell-vs-bulk
  determination doesn't depend on it).

Submitting for round-3 review before download, same discipline as every
prior step this session.
