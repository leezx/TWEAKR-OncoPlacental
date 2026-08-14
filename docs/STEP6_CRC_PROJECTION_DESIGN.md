# Step 6 (Phase I): Projecting the Frozen D/F/P Signature onto Real CRC Oncofetal Data

## Purpose

This is the step that finally answers the project's original Q1 (`docs/PROJECT_SUMMARY.md:9-13`):

> Is CRC's "Oncofetal" malignant epithelial cell state actually one program, or a
> mixture of two distinct normal-development programs — a fetal-somatic
> (embryonic organ development) program and a separate placental/trophoblast
> program — plus a possible third, genuinely shared program between the two?

Inputs (all frozen, closed, reviewer-approved, and Tier-2-validated — Steps 1-5):
`results/04_dfp_signature/dfp_gene_sets/D_shared_FINAL.txt` (6 genes),
`F_specific_FINAL.txt` (2,504 genes), `P_specific_FINAL.txt` (78 genes).

This step does **not** re-open or re-calibrate any Step 4 cutoff. It is the first
time cancer data is touched in this project — by explicit design (`Worklog.md:103`:
"Cancer data (TCGA/CRC atlas) should not be touched until the signature is
frozen").

## Real inventory findings (this round, job 3620558, `results/06_crc_projection/inventory/crc_dataset_inventory.json`)

No CRC dataset lives inside TWEAKR-OncoPlacental yet. Four real candidates were
found staged on Argos (mirrored from a sibling, unrelated project's data lake —
none of this was previously registered inside this repo) and inventoried directly
(same discipline as Steps 1/5 — never assumed structure or raw-counts
availability from documentation elsewhere):

| Dataset | Cells | Genes | Raw counts | Cell-type annotation | Patients/donors |
|---|---|---|---|---|---|
| **CRC_single_cell_atlas_2025** (meta-atlas) | 665,473 | 28,476 (Ensembl) | `layers['counts']`, 100% integer (confirmed) | `atlas_cell_type_middle`: `Cancer cell` (509,421) / `CRLM` (156,052); also coarse/fine tiers | 54 constituent studies, 9 platforms, rich clinical/molecular metadata (CMS subtype, MSI, KRAS/BRAF/TP53 status, tumor stage) |
| **HTAN_CRC_progressive_plasticity** | 47,107 (epithelial-only export) | 25,344 (Ensembl) | `.raw.X`, 100% integer (confirmed; `X` itself is normalized) | `cell_type`: `malignant cell` (26,551) vs. named normal epithelial subtypes (early colonocyte, secretory, ISC, goblet, BEST4+, tuft, enteroendocrine); `Tumor Status`/`Sample Type`: Primary/Metastasis/Non-Tumor | 29 patients |
| **CRLM_NMP_ATLAS** | 75,104 | 30,257 (Ensembl) | `.raw.X`, 100% integer (confirmed) | `cell_type`: `malignant cell` only 4,051 of 75,104 — dataset is TME/immune-focused (T cell, neutrophil, NK, macrophage dominate) | 6 donors, liver-metastasis-focused (`timepoint`: Before/After NMP) |
| **GSE178318** | 140,281 | 33,694 (Ensembl+symbol pairs in `genes.tsv`) | Raw 10x mtx, 100% integer (confirmed) | **None** — raw matrix only, no cell-type calls; barcode suffix encodes sample (e.g. `_COL07_CRC`), but no metadata table inventoried yet | 9 patients (15 GSM samples per prior registry note, not independently re-verified this round) |

**No dataset carries a pre-existing "Oncofetal" (or "revCSC"/"M11") label.** A
keyword scan of every low-cardinality obs column for
`oncofetal|fetal|placent|trophoblast|revcsc|m11` returned only coincidental
substring matches inside random-looking auto-generated join IDs and one dataset
name (`...LM112`) — not a real annotation. This confirms the Explore-agent
finding from before this design: "Oncofetal" must be **operationally defined**
via the D/F/P projection itself, not read off an existing column. (The
KB-strategy doc's M11/revCSC meta-program, computed separately on this same
meta-atlas via NMF, was not found as an `obsm`/`uns` entry in `adata_nmf.h5ad` —
it likely lives in a separate results file in that directory tree; worth locating
later as an independent cross-check, but not required to start this step.)

## Proposed dataset plan

**Primary: `CRC_single_cell_atlas_2025` meta-atlas.** Largest sample size (509K
malignant "Cancer cell" calls), broadest external validity (54 independent
studies, 9 sequencing platforms), richest clinical/molecular annotation
(CMS subtype, MSI status, driver-gene mutation status — useful for later
subgroup checks), and raw counts already confirmed. Existing `atlas_cell_type_*`
annotation is used as-is (not re-derived) — same principle as Step 5's
"Tabula Sapiens is held-out validation, not re-tuned" discipline: this project
does not re-cluster or re-annotate someone else's atlas, it uses the
already-published cell-type calls.

**Secondary / independent replication: `HTAN_CRC_progressive_plasticity`.**
Smaller but very cleanly annotated — importantly, it separates malignant cells
from *named, specific* normal epithelial subtypes (not just "normal" as one
bucket) and from non-tumor tissue, and explicitly separates Primary vs.
Metastasis. Different underlying studies than the meta-atlas (avoids
double-dipping on the same source data when checking whether a finding
replicates).

**Tertiary / liver-metastasis context: `CRLM_NMP_ATLAS`.** Smaller malignant
population (4,051 cells) and TME-focused, but relevant given the project's
liver-metastasis angle and its `Before/After_NMP` timepoint structure — kept as
a secondary check, not a primary analysis population, given the small malignant
n.

**Deferred: `GSE178318`.** No cell-type annotation exists yet in the files
inventoried this round — would need either a separate metadata/annotation file
(not yet located) or de novo clustering, which this project has consistently
avoided doing itself (using published annotations only). Not blocking — flagged
for a later round if the primary/secondary datasets leave an open question this
one could resolve.

## Open items before compute (analogous to Step 2's gene-ID mapping precedent)

1. **Gene-ID mapping**: all 4 candidates use Ensembl gene IDs; the frozen D/F/P
   gene lists are `canonical_symbol`-based (from HDMA, Step 2). Needs an
   Ensembl→symbol map. Candidate reuse: Step 2's HGNC-based resolution table
   (`datasets/HGNC_gene_id_mapping/`) — needs a coverage check against these
   4 datasets' specific Ensembl ID sets before assuming it's sufficient (same
   discipline as Step 2's original 9.3%/4.6% coverage audit, not assumed).
2. **Scoring method for single-cell gene-set projection**: Step 4's original
   percentile-based scoring was designed for bulk/pseudobulk data
   (organ-level, not cell-level). For single-cell, propose `scanpy.tl.score_genes`
   (Seurat `AddModuleScore`-equivalent: mean expression of the gene set minus
   mean expression of an expression-level-matched control gene set) — directly
   analogous in spirit to the expression-matched permutation-null design just
   built and reviewer-approved in Step 5, and it is the field-standard method
   for this exact "does this cell express this curated gene program" question.
   Proposed for D-shared, F-specific, and P-specific scored separately per cell.
3. **Aggregation / donor-awareness**: per the standing donor/sample-aware
   discipline (Steps 4/5), per-cell scores will be aggregated to
   per-patient/per-sample pseudobulk-style summaries (median score, fraction of
   malignant cells scoring above a to-be-calibrated threshold) before any
   cross-patient claim is made — a single patient's cells should not be able to
   drive a "CRC malignant cells are F-dominant" conclusion. The meta-atlas's 54
   constituent studies and 9 platforms also need to be tracked as a covariate
   (same lesson as Step 5's cross-platform sequencing-depth confound) — a
   finding that only shows up in one study/platform is a weaker claim than one
   replicated across studies.
4. **Primary hypothesis test**: within-patient malignant-vs-matched-normal-epithelial
   contrast (available directly in HTAN's annotation; the meta-atlas would need
   normal/polyp samples cross-referenced via `medical_condition`/`sample_type`) —
   do malignant cells show elevated D/F/P scores relative to their own patient's
   normal epithelial cells, controlling for patient identity? This is the direct
   test of "oncofetal reactivation" and avoids simple between-patient confounds.

## What this step is not

Not a re-derivation of "Oncofetal" via new clustering/NMF on this project's own
compute — existing published cell-type calls are used as-is. Not a re-calibration
of Step 4's frozen D/F/P cutoffs. Not yet a spatial-transcriptomics analysis
(no spatial CRC dataset has been inventoried in this round — single-cell only
for this first pass).

Submitting this design for review before running any real gene-ID mapping or
scoring compute.
