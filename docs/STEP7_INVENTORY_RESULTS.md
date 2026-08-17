# Step 7: CLiM/CLuM external data inventory — real compute results

Executes `docs/STEP7_CLIM_DATA_ACQUISITION_DESIGN.md` (PR #35, APPROVE
after 3 review rounds) for real: downloads every public scRNA-seq and
bulk RNA-seq/microarray cohort locked in that design, byte-verifies each
download, checks archive integrity, and structurally characterizes each
dataset (sample/cell counts, gene-ID format, raw-vs-processed status,
paper-cohort-reconstruction status). **Inventory only** — no D/F/P/revCSC
scoring against any of these cohorts here, per the design's explicit
scope boundary; that is a separate, deliberately deferred follow-on PR.

All compute ran on Argos (`argos-qsub1`, plain network I/O for downloads,
`argos-codex` conda env for the inventory characterization script — not
qsub-wrapped, same precedent as this project's original HDMA/Arutyunyan
acquisition phase). Every download byte-verified against the exact
integer `Content-Length` (GEO) or GDC manifest `file_size`+`md5sum`
(TCGA) — never GEO's rounded MB/GB webpage display, never piped through
`tail`/`head` (standing `curl_pipe_swallows_exit_code` lesson).

**Went through 3 review rounds** (`PR #36`) — round 1: 4 blockers plus 1
completeness gap; round 2 and round 3: 2 successive residual correctness
issues in GSE225857's patient-cohort reconstruction, the second of which
this project could not fully resolve either way (see round 3). All
independently re-verified against real committed data/live sources
before fixing. See "Round-1 review fixes," "Round-2 review fix," and
"Round-3 review fix" below.

## Scripts

- `scripts/07_clim_external_data/download_geo.sh` — downloads all 6 GEO
  cohorts (GSE231559, GSE225857, GSE285990, GSE131418, GSE17536,
  GSE17537, GSE21510), exact-Content-Length byte verification per file.
- `scripts/07_clim_external_data/download_tcga.py` — TCGA-COAD+READ via
  the GDC API (distinct access method from every GEO cohort — not a GEO
  `_RAW.tar`), size+md5 verification against the GDC manifest.
- `scripts/07_clim_external_data/inventory_clim_data.py` — structural
  characterization + archive integrity + paper-cohort-reconstruction
  attempt, per dataset.

## Download results: 0 failures across all 7 cohorts

| Dataset | Bytes downloaded | Verification |
|---|---|---|
| GSE231559 | 725,678,080 (RAW.tar) | Exact `Content-Length` match |
| GSE225857 | 4 files, 326,770,005 total | Exact `Content-Length` match, each |
| GSE285990 | 30 files (10 GSMs × 3), 615,512,611 total | Exact `Content-Length` match, each |
| GSE131418 | 5,661,460,480 (RAW.tar) + 4 clinical/annotation files | Exact `Content-Length` match, each |
| GSE17536 | 1,600,727,040 (RAW.tar) | Exact `Content-Length` match |
| GSE17537 | 295,290,880 (RAW.tar) | Exact `Content-Length` match |
| GSE21510 | 771,153,920 (RAW.tar) | Exact `Content-Length` match |
| TCGA-CRC | 701 files via GDC | Exact `file_size` + `md5sum` match, each — **independently re-verified exhaustively, round-1 review fix, see below** |

Archive integrity (`tar -tf` on the outer plain-tar archive + `gzip -t`
on **every** extracted `.gz` member — exhaustive, round-1 review fix,
see below): **all clean**, 0 failures, for every `_RAW.tar` (GSE231559:
78/78 members; GSE131418: 1,138/1,138; GSE17536: 177/177; GSE17537:
55/55; GSE21510: 148/148).

## Round-1 review fixes (`PR #36`, all independently re-verified before
fixing)

**Blocker 1 — GSE225857's locked per-patient reconstruction acceptance
criterion was never actually executed.** PR #35 explicitly left the 4
CLiM + 4 primary reconstruction, from the two pooled scRNA metadata
files, as required work for this inventory pass. The first draft of
`inventory_gse225857()` only read metadata column headers, regex-matched
for `patient|sample|donor|origin`, and reported
`per_patient_reconstruction_possible=true` without ever reading a single
value — the results doc even said reconstruction was "not attempted."
Confirmed real by re-reading the committed function directly. **Fixed**:
now parses the actual `patients`/`patients_organ` values from both
GSMs' metadata files. **Real, honest result**: 7 unique patients have
liver-cancer (CLiM) tissue data (`LCL`/`LCT` prefixes), 6 unique
patients have colon-cancer (primary CRC) tissue data (`CCL`/`CCT`
prefixes) — **neither matches the paper's cited 4+4 exactly**. These are
pooled dissociated-cell fractions (immune vs. non-immune) across all
profiled patients, not the paper's discrete per-tumor-block sample
structure, so a clean 4+4 subset isn't recoverable by tissue-prefix
counting alone. Reported as PARTIAL/UNRESOLVED, not forced to match.
**Superseded by the round-2 fix below** — the 7/6 split here turned out
to rest on an unsupported classification assumption, corrected in round
2.

**Blocker 2 — GSE131418 declared unresolved without opening the clinical
metadata files the PR itself downloads.** `download_geo.sh` downloads
`GSE131418_GEO_submission_Recurrence_meta_data_V2_Updated.xls.gz` and
`...stage4_meta_data_V2_Updated.xls.gz`, specifically for cohort
reconstruction — but `inventory_gse131418()` never opened either file,
only the series matrix. Confirmed real by re-reading the committed
function. **Fixed**: both legacy-binary `.xls` files (OLE2/Composite
Document Format, read via `xlrd`) are now actually opened. **Real,
honest result**: both files are exclusively `mcc.prim.*` samples
(Recurrence file: 134/134 primary; Stage4/survival file: 40/40 primary)
— **neither contains a single metastasis sample**. They are primary-
tumor clinical-outcome/survival annotation files, not annotations of the
liver-metastasis samples themselves, so they structurally cannot resolve
which subset of the 197 GEO-identified liver-met samples matches the
paper's cited 170. The 170-sample reconstruction **remains UNRESOLVED**
— but now for a directly-verified reason, not a premature dismissal.

**Blocker 3 — archive integrity was weaker than the locked design.**
PR #35 locked `tar -tf` + `gzip -t` on every extracted member; the first
draft's `gzip -t` only ran on a 20-file spot-check (via a slow
per-member `tar -xOf <big-tar> <member> | gzip -t`, which re-scans the
whole archive each call — infeasible to run exhaustively that way on a
1,138-member/5.66GB archive), while the PR body and this doc still
called every archive "clean" and "per the design's corrected mechanics"
— overstating what a 20/1,138 spot-check actually established. Confirmed
real by re-reading the committed function. **Fixed**: each archive is
now extracted once (a single sequential pass), then `gzip -t` runs on
**every** extracted `.gz` member directly off local disk (fast per-file,
no re-scanning) — genuinely exhaustive. Re-run confirms 0 failures
across all 1,596 `.gz` members checked (78+1,138+177+55+148).

**Blocker 4 — TCGA's "verified" claim didn't actually check md5.**
`download_tcga.py`'s resume-skip logic and `inventory_tcga()`'s
verification count both only compared on-disk file size against the
manifest, never recomputed md5 — so a stale/same-length-but-corrupted
file could have silently passed, despite the PR's claim of "701/701
verified (size+md5)." Confirmed real by re-reading both committed
functions. **Fixed**: the resume-skip predicate now requires size AND
md5 match; `inventory_tcga()` independently re-hashes **every** one of
the 701 files against the manifest's own md5sum (not sampled). Re-run
confirms: 701/701 verified, 0 missing, 0 size mismatches, 0 md5
mismatches.

**Completeness gap — bulk/TCGA structural characterization was
shallower than the design's declared dimensions.** The design commits to
characterizing "sample/cell counts, gene-ID format, raw-vs-processed
status" per dataset; the scRNA loaders did this, but
`inventory_cel_series()` only counted CEL archive members and
`inventory_tcga()` never opened a single expression file. **Fixed**:
`inventory_cel_series()` now extracts a real Affymetrix array-type
string directly from a sample CEL file's own binary header (confirmed
`HG-U133_Plus_2` for all 3 CEL series — matches GPL570, independently
confirming the series-matrix `!Sample_platform_id` field rather than
just trusting it) and records CEL's raw-probe-level status explicitly.
`inventory_tcga()` now opens a sample `augmented_star_gene_counts.tsv`
file directly and confirms real Ensembl gene IDs (60,660/60,664 rows are
genes, versioned Ensembl format, the other 4 are STAR's own
`N_unmapped`/`N_multimapping`/`N_noFeature`/`N_ambiguous` summary rows)
and that both raw (`unstranded`, integer STAR counts) and processed
(`tpm_unstranded`/`fpkm_unstranded`/`fpkm_uq_unstranded`) representations
coexist in the same file.

## Round-2 review fix (`PR #36`, independently re-verified before fixing)

**Blocker — GSE225857's round-1 patient-cohort reconstruction rested on
an unsupported classification assumption.** The round-1 fix unioned
`LCL`(immune)`|LCT`(non-immune) into a single "liver-cancer patient"
count and `CCL|CCT` into a single "colon-cancer patient" count, treating
the immune-cell-sorting-fraction prefixes (`LCL`/`CCL`) as equally
tumor-defining as the non-immune/tumor-containing fraction (`LCT`/`CCT`)
— explicitly justified as "the same L/C convention confirmed for
GSE231559." That justification doesn't hold: GSE231559 and GSE225857 are
unrelated datasets with independent naming schemes, and the resulting
7-liver-cancer-patient count directly conflicts with GSE225857's own
series summary, which states "6 CRC patients with liver metastasis were
enrolled" / "27 samples of 6 CRC patients." Confirmed real by directly
re-fetching GSE225857's series summary text.

**Fixed** with a properly-grounded reconstruction: the non-immune
fraction's `LCT` (liver-tumor) and `CCT` (colon-tumor) patient sets are
IDENTICAL — `{s0107, s0115, s0813, s0920, s1231}`, n=5 — a genuinely
paired primary+liver-met tumor-tissue cohort, not an assumption borrowed
from an unrelated series. This also explains the 7-vs-6 discrepancy
precisely: `s1125` has immune-fraction (CD45+) presence in *both* liver
and colon but no tumor-fraction data in either — adding it to the 5
paired patients gives 6 patients with disease representation in both
organs, matching GEO's cited "6 CRC patients" exactly. `s0816` has
liver+blood representation only, with zero colon representation of any
kind (immune or tumor fraction) — the 7th total unique patient across
all organ prefixes, structurally distinct from the other 6, consistent
with being a non-CRC inclusion (flagged as a hypothesis, not confirmed
from public data alone). The corrected result (5 paired tumor-tissue
patients) is closer to the paper's cited 4+4 than the round-1 result (7)
but still doesn't match exactly — reported honestly as unresolved, not
forced to match. **This section's `s1125`/`s0816` directional hypothesis
was itself corrected in round 3, below** — kept here as an accurate
historical record of what round 2 actually fixed and asserted, not
retroactively edited.

## Round-3 review fix (`PR #36`, independently investigated before
fixing — the finding this time could not be confirmed either way)

**Blocker — the round-2 fix's directional hypothesis about which patient
is "likely non-CRC" was itself an unsupported claim.** Round 2 asserted
`s1125` (immune-fraction presence in both organs) "completes" GEO's 6
cited CRC patients and that `s0816` (liver+blood only) is "consistent
with being a non-CRC inclusion." The reviewer's round-3 finding disputed
this directly, citing "the original authors' own public analysis code"
as constructing the tumor-cell input from 6 patient IDs that include
`s0816` but *exclude* `s1125` — the opposite assignment.

**This project attempted to independently verify the reviewer's specific
counter-claim and could not.** Searched for any publicly accessible
GitHub/Zenodo/code repository tied to GSE225857's source publication
(Wang F, Long J et al., *Science Advances* 2023, PMID 37327339) via
GEO's own contact/publication metadata and targeted web search for the
specific variable/ID strings cited (`tumor_merge`, `immune_merge`,
`s0816_1`) — found no such repository reachable without authentication.
Per this session's standing discipline (never accept a specific factual
claim, including a reviewer's, without independent confirmation), this
project does not adopt the reviewer's counter-claim as fact either.

**Fixed** by removing the directional assertion entirely rather than
flipping it: both `s1125` and `s0816` are now reported neutrally, with
their objective, directly-computed organ-prefix membership stated
plainly and *no* CRC/non-CRC label assigned to either. Which one (if
either) corresponds to GEO's 6th cited "CRC patient with liver
metastasis" is explicitly left unresolved — this project's own public
GEO metadata and public-search access are insufficient to settle it
either way, and the honest record is to say so rather than assert a
direction this project cannot itself verify.

## Per-dataset structural characterization + cohort-reconstruction results

### GSE231559 — scRNA-seq, 26 total samples

Clean 10x-style MTX loader pattern (bare/versioned Ensembl gene IDs,
100% integer-valued sampled nonzeros in every sample — genuine raw
counts). **Paper's cited 9 CLiM + 6 primary CRC subset reconstructs
EXACTLY**: 9 liver-tumor samples (`L1T, L4T, L6T, L8T1, L8T2, L9T, L10T,
L11T, L12T`) = CLiM exactly; 6 colon-tumor samples (`C1T-C6T`) = primary
CRC exactly, classified by joining the RAW.tar's extracted GSM
accessions against the GEO series matrix's own `Sample_title` field (not
the tar's internal library-ID filenames, which use an unrelated naming
scheme). The remaining 11 samples (8 liver-normal, 3 colon-normal) are
paired-normal reference tissue, not part of the cited 9+6 tumor cohort.
Full per-sample table: `results/07_clim_external_data/GSE231559_inventory.tsv`.

### GSE225857 — scRNA-seq, 2 genuine scRNA GSMs

`GSM7058754` (immune, 196,473 cells) and `GSM7058755` (non-immune,
41,892 cells) — matches the design's locked exact-GSM-accession
selection (the other 6 GSMs in this series are spatial, correctly not
downloaded). Both GSMs' accompanying `*_meta.txt.gz` files carry real
per-cell metadata including `patients`/`sampletag`/`patients_organ`
columns, confirming the design's round-2 fix (patient-of-origin metadata
lives in a **separate file**, not "inside the matrix"). **Real
per-patient reconstruction result** (round-2/round-3 review fixes, see
above): the non-immune fraction's liver-tumor (`LCT`) and colon-tumor
(`CCT`) patient sets are IDENTICAL — 5 paired patients
(`s0107, s0115, s0813, s0920, s1231`) — with `s1125` (immune-fraction
presence in both organs, no tumor-fraction data) and `s0816`
(liver+blood only, zero colon representation) accounting for the
remaining 2 of GSE225857's 7 total unique patients. One of these two,
added to the 5 paired patients, would give 6 — matching GEO's own
series-summary count of "6 CRC patients" — but **which of the two is
not resolvable from public data** (round-3 fix: neither is asserted to
be the CRC/non-CRC one). Close to but not an exact match for the paper's
cited 4+4 either way, reported as PARTIAL/UNRESOLVED. Full table:
`results/07_clim_external_data/GSE225857_inventory.tsv`.

### GSE285990 — scRNA-seq, 10/10 human liver-metastasis samples confirmed

All 10 `P01_LM`-`P10_LM` GSMs (`GSM8714595`-`GSM8714604`) load cleanly:
consistent axes, 100% bare/versioned-Ensembl gene IDs, genuine integer
raw counts, cell counts ranging 9,203-17,088 per sample (real, not
degenerate). **Matches the paper's cited 10-sample cohort exactly.**
This round's inventory additionally confirmed — directly, via each
GSM's own record, not assumed — that this GEO series also contains a
**separate Mus musculus Kupffer-cell/FOLFOX mechanistic sub-study**
(`GSM8714605`-`GSM8714609`: `NT_LM`/`RS_LM`/`RL_LM`/`WT_LM`/`WT_DTR_LM`,
plus the series-level `KCs_gene_exp.txt.gz`/`WT_DTR_LM_*` files) and 2
**Mus musculus Stereo-seq spatial** samples (`GSM8714610`-`8714611`) —
both correctly excluded from download (wrong species for the human
cohort; the Stereo-seq pair also wrong modality). Full table:
`results/07_clim_external_data/GSE285990_inventory.tsv`.

### GSE17536 + GSE17537 — bulk microarray, primary CRC

177 + 55 = **232 CEL files exactly**, matching the paper's cited n=232
(subseries-count reconstruction already confirmed during the design-
review round). This round confirms, directly from a sample CEL file's
own binary header (not assumed from GEO's `!Sample_platform_id` field
alone): both series use the **Affymetrix HG-U133_Plus_2** array (GPL570)
— real raw probe-level intensity data, not summarized/normalized
expression values. Archive integrity: exhaustive `tar -tf` + `gzip -t`
clean on both `_RAW.tar` files (177/177, 55/55).

### GSE131418 — bulk microarray, 1,135 total samples; paper cites 170 liver-met

Real, honest attempt to reconstruct the paper's cited 170-sample
liver-metastasis subset, **not forced to match by construction**:

- Two sub-cohorts distinguishable by sample-title prefix:
  `consortium` (618 samples) and `mcc` (517 samples), 618+517=1,135 ✓
- `site of metastasis: LIVER` — **197 total** (141 MCC + 56 Consortium)
- `site of metastasis: LIVER` restricted to `treatment status: PRE` only
  — 53
- **Round-1 review fix**: the 2 clinical XLS files (`Recurrence`: 134
  rows; `Stage4`/survival: 40 rows) were actually opened this round —
  both are exclusively primary-tumor samples, contain zero metastasis
  samples, and therefore cannot help resolve the liver-met subset
  question at all.
- **None of these groupings reproduce 170 exactly.**

**Reported as UNRESOLVED, now with a directly-verified reason** — the
paper's exact 170-sample subset likely requires its own supplementary
sample list not present in any of the publicly downloadable GEO
material for this series. Downloading all 1,135 samples (done,
byte-verified, archive-integrity-clean, exhaustive gzip-t 1,138/1,138)
is sufficient per the design ("downloading all X is fine; declaring the
cohort reconstructed is not"); the exact cohort match remains an open
item. Full per-sample metadata:
`results/07_clim_external_data/GSE131418_sample_metadata.tsv`.

### GSE21510 — bulk microarray, 148 total samples; paper cites 146

Real reconstruction attempt against GEO's own series-matrix
characteristics:

- Tissue-prep breakdown: `cancer, LCM`=104, `normal, homogenized`=25,
  `cancer, homogenized`=19 — sums to 148 exactly, no natural
  146-sized subgroup among these 3 categories or any pairwise
  combination.
- 107 unique patients total (41 with >1 sample); 104 unique patients
  have a `cancer, LCM` sample (no duplicates within that category).
- Array type confirmed directly from a sample CEL header:
  **Affymetrix HG-U133_Plus_2** (GPL570).

**Reported as UNRESOLVED** — the exact 2 excluded samples and the
exclusion basis (QC failure, sample swap, etc.) cannot be identified
from GEO's public series-matrix metadata alone; this would require the
paper's own supplementary methods/QC exclusion list. Downloading all
148 (done, byte-verified, exhaustive archive-integrity-clean 148/148) is
sufficient per the design; the exact 146-sample cohort remains an open
item. Full per-sample metadata:
`results/07_clim_external_data/GSE21510_sample_metadata.tsv`.

### TCGA-CRC — bulk RNA-seq via GDC, 701 open-access files; paper cites 610

TCGA-COAD (524 files) + TCGA-READ (177 files) = 701, gene-level STAR-
Counts `augmented_star_gene_counts.tsv`, all open-access (no dbGaP
application needed). **Round-1 review fix**: every one of the 701 files
is now independently re-hashed against the GDC manifest's own md5sum
(not sampled, not size-only) — 701/701 verified, 0 missing, 0 size
mismatches, 0 md5 mismatches. Structural characterization confirmed
directly from a sample file: 60,660 of 60,664 data rows are real genes
(versioned Ensembl IDs, e.g. `ENSG00000000003.15`), the other 4 are
STAR's own summary rows (`N_unmapped` etc.); both raw (`unstranded`,
integer STAR counts) and processed (`tpm_unstranded`/`fpkm_unstranded`/
`fpkm_uq_unstranded`) representations coexist in the same file.

Sample-type breakdown: 647 Primary Tumor, 51 Solid Tissue Normal, 2
Recurrent Tumor, 1 Metastatic. Restricting to Primary-Tumor files gives
647 files across **624 unique cases** (13 cases have >1 Primary-Tumor
aliquot — real technical replicates, not paper-cohort members counted
twice).

**Reported as PARTIAL, not exact** — 624 unique primary-tumor cases is
close to but does not exactly match the paper's cited 610; the specific
filter needed (which aliquot to keep per multi-aliquot case, any
additional QC/data-completeness exclusion the paper applied) is not
reconstructable from GDC file metadata alone. This is the closest
achievable reconstruction from public metadata, reported honestly as
such, not declared exact.

## What this PR does not do

Runs no D/F/P/revCSC scoring against any of these cohorts — deliberately
deferred to a separate follow-on PR, per the locked design. Does not
resolve the 3 explicitly-unresolved cohort-reconstruction items
(GSE131418's 170, GSE21510's 146, TCGA's exact 610, GSE225857's exact
4+4) — these remain open, honestly documented, not silently assumed or
forced to match. Does not touch the Zenodo-restricted record
(`10.5281/zenodo.19043057`) or any spatial/LCM-WGS/protein-imaging data
— out of scope per the user-confirmed descope and the design's scope
boundary.
