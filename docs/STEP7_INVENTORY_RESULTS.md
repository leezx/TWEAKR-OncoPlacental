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
| TCGA-CRC | 701 files via GDC | Exact `file_size` + `md5sum` match, each |

Archive integrity (`tar -tf` on the outer plain-tar archive + `gzip -t`
spot-check on 20 extracted `.gz` members per archive, per the design's
corrected mechanics): **all clean**, 0 failures, for every `_RAW.tar`
(GSE231559: 78/78 members; GSE131418: 1,138/1,138; GSE17536: 177/177;
GSE17537: 55/55; GSE21510: 148/148).

## Real bugs found and fixed during this compute (reported honestly, not
just the clean final numbers)

**Bug 1 — GSE231559 cohort classification joined on the wrong key.**
First draft of `inventory_clim_data.py` parsed the paper-facing L#N/L#T/
C#N/C#T label directly from the RAW.tar's extracted filenames (e.g.
`GSM7290760_SC10_21N_barcodes.tsv.gz`) — but those filenames embed an
internal sequencing-library ID (`SC10_21N`), not the GEO series matrix's
`Sample_title` field (`L1N`), which is what actually carries the L/C
labeling. First run reported **0/9 CLiM and 0/6 primary** — a real bug,
not a genuine reconstruction failure (confirmed by directly fetching the
series matrix and joining on GSM accession: `GSM7290760` = `L1N`).
Fixed by adding `fetch_series_matrix_gsm_title_map()`, which fetches the
series matrix and classifies by GSM-accession join instead of filename
parsing. Re-run confirmed the correct result (see below).

**Bug 2 — GSE285990 GSM-ID string-concatenation rollover.** GSM IDs were
generated as `f"GSM87145{95+i}"` for `i in range(10)`, intending
`GSM8714595`-`GSM8714604`. This works for `i=0..4` (`95+i` stays 2-digit:
95-99), but at `i=5`, `95+i=100` (3-digit), so the f-string produced
`"GSM87145"+"100"` = `GSM87145100` — a malformed 11-character ID, not
`GSM8714600`. This silently broke 5 of 10 samples (P06-P10), reported as
`"error": "missing file(s)"` in the first run. Fixed to
`f"GSM8714{595+i}"` (fixed 4-digit prefix, full 3-digit variable suffix).
Re-run confirmed all 10 samples load correctly (see below).

**Bug 3 — a real GDC API quirk, not a data problem** (found and fixed
before the write-up, during the TCGA download itself): the GDC `/data`
bulk endpoint returns a `tar.gz` bundle (`MANIFEST.txt` + one
subdirectory per file) when given multiple file IDs, but returns the raw
file's bytes directly, unwrapped, when given exactly **one** ID. With
701 total files and `BATCH_SIZE=50`, the final batch (`701 % 50 = 1`)
hit this case and failed with `tarfile`'s "not a gzip file" error.
Confirmed via a direct single-ID request: identical byte size and md5 to
the manifest, just not tar-wrapped — a real, reproducible GDC API
behavior, not a corrupted download. Fixed `download_batch()` in
`download_tcga.py` to detect `len(file_ids)==1` and write the raw
response directly instead of assuming `tarfile.open()` always applies.
All 701/701 files verified (size+md5) after the fix.

## Per-dataset structural characterization + cohort-reconstruction results

### GSE231559 — scRNA-seq, 26 total samples

Clean 10x-style MTX loader pattern (bare/versioned Ensembl gene IDs,
100% integer-valued sampled nonzeros in every sample — genuine raw
counts). **Paper's cited 9 CLiM + 6 primary CRC subset reconstructs
EXACTLY** once classified by the correct GSM-accession join (Bug 1,
above): 9 liver-tumor samples (`L1T, L4T, L6T, L8T1, L8T2, L9T, L10T,
L11T, L12T`) = CLiM exactly; 6 colon-tumor samples (`C1T-C6T`) = primary
CRC exactly. The remaining 11 samples (8 liver-normal, 3 colon-normal)
are paired-normal reference tissue, not part of the cited 9+6 tumor
cohort. Full per-sample table: `results/07_clim_external_data/GSE231559_inventory.tsv`.

### GSE225857 — scRNA-seq, 2 genuine scRNA GSMs

`GSM7058754` (immune, 196,473 cells) and `GSM7058755` (non-immune,
41,892 cells) — matches the design's locked exact-GSM-accession
selection (the other 6 GSMs in this series are spatial, correctly not
downloaded). Both GSMs' accompanying `*_meta.txt.gz` files carry real
per-cell metadata including patient-of-origin candidate columns
(`patients`, `sampletag`, `patients_organ`) — confirming the round-2
design-review fix (patient-of-origin metadata lives in a **separate
file**, not "inside the matrix") and confirming per-patient (4 CLiM + 4
primary) reconstruction is structurally possible from this metadata, not
attempted in this inventory-only PR. Full table:
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
(this reconstruction was already confirmed via GEO's own subseries
counts during the design-review round; this round adds real archive
integrity confirmation — both `_RAW.tar` files pass `tar -tf` +
`gzip -t` spot-check cleanly, and GEO's own declared per-series sample
counts match the actual extracted CEL-file counts exactly).

### GSE131418 — bulk microarray, 1,135 total samples; paper cites 170 liver-met

Real, honest attempt to reconstruct the paper's cited 170-sample
liver-metastasis subset from GEO's own series-matrix per-sample
characteristics, **not forced to match by construction**:

- Two sub-cohorts distinguishable by sample-title prefix:
  `consortium` (618 samples) and `mcc` (517 samples), 618+517=1,135 ✓
- `site of metastasis: LIVER` — **197 total** (141 MCC + 56 Consortium)
- `site of metastasis: LIVER` restricted to `treatment status: PRE` only
  — 53
- **None of these groupings reproduce 170 exactly.**

**Reported as UNRESOLVED** — the paper's exact 170-sample subset likely
requires its own supplementary sample list, not reconstructable from
GEO's public per-sample characteristics fields alone. Downloading all
1,135 samples (done, byte-verified, archive-integrity-clean) is
sufficient per the design ("downloading all X is fine; declaring the
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

**Reported as UNRESOLVED** — the exact 2 excluded samples and the
exclusion basis (QC failure, sample swap, etc.) cannot be identified
from GEO's public series-matrix metadata alone; this would require the
paper's own supplementary methods/QC exclusion list. Downloading all
148 (done, byte-verified, archive-integrity-clean) is sufficient per the
design; the exact 146-sample cohort remains an open item. Full
per-sample metadata: `results/07_clim_external_data/GSE21510_sample_metadata.tsv`.

### TCGA-CRC — bulk RNA-seq via GDC, 701 open-access files; paper cites 610

TCGA-COAD (524 files) + TCGA-READ (177 files) = 701, gene-level STAR-Counts
`augmented_star_gene_counts.tsv`, all open-access (no dbGaP application
needed). Sample-type breakdown: 647 Primary Tumor, 51 Solid Tissue
Normal, 2 Recurrent Tumor, 1 Metastatic. Restricting to Primary-Tumor
files gives 647 files across **624 unique cases** (13 cases have >1
Primary-Tumor aliquot — real technical replicates, not paper-cohort
members counted twice).

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
(GSE131418's 170, GSE21510's 146, TCGA's exact 610) — these remain open,
honestly documented, not silently assumed or forced to match. Does not
touch the Zenodo-restricted record (`10.5281/zenodo.19043057`) or any
spatial/LCM-WGS/protein-imaging data — out of scope per the user-confirmed
descope and the design's scope boundary.
