# Step 6 (Phase I): Projecting the Frozen D/F/P Signature onto Real CRC Oncofetal Data

## Purpose

**PR #16 review round 2 (REQUEST_CHANGES) caught a definitional-circularity blocker**
in everything below this point, before it was ever run: the original plan defined
"Oncofetal cells" *via* the D/F/P projection itself ("Oncofetal must be
operationally defined via the D/F/P projection"), then proposed asking whether
those same D/F/P-selected cells are F/P/D-dominant — a circular measurement,
since the yardstick used to find "Oncofetal" cells was itself F/P/D. Fixed by
adding an explicit **Layer 1 / Layer 2 structure** (see the new section below,
"Independent Oncofetal anchor (Layer 1)") — Layer 1 locates/defines Oncofetal
status using a gene signature with zero overlap in construction with D/F/P;
only Layer 2 projects the frozen D/F/P signature onto Layer-1-identified cells.

**PR #16 review round 3 (REQUEST_CHANGES, then withdrawn once verified)**
asked for concrete proof that the M11/revCSC anchor is actually independent —
construction-independent, not just "found in a different directory" — and for
the exact provenance of the 297,307-cell subset it lives on. Both were checked
directly (reading the actual NMF pipeline script, checking file mtimes,
querying the atlas for the subset's cell-type/study/platform composition) and
are folded into the "Independent Oncofetal anchor" section below; reviewer
confirmed this closes the circularity concern and also refined the primary
question and analysis hierarchy (see "Revised analysis structure" below).

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

**No dataset carries a pre-existing "Oncofetal" (or "revCSC"/"M11") *column*.** A
keyword scan of every low-cardinality obs column for
`oncofetal|fetal|placent|trophoblast|revcsc|m11` returned only coincidental
substring matches inside random-looking auto-generated join IDs and one dataset
name (`...LM112`) — not a real annotation.

## Independent Oncofetal anchor (Layer 1) — fixes the circularity blocker

**Round 1 of this design proposed defining "Oncofetal" via the D/F/P projection
itself, then asking whether D/F/P-selected cells are F/P/D-dominant — reviewer
correctly identified this as circular** (the yardstick used to find "Oncofetal"
cells would be the same yardstick used to characterize them). Fixed by locating
real, pre-existing, independently-computed evidence in the
`CRC_single_cell_atlas_2025` meta-atlas's own directory tree (same machine, not
previously checked in this repo) — **not re-derived by this project, used as-is**,
same discipline as reusing published cell-type calls:

- `NMF/metamodule_fnmf/MM_alt_clean_byscore_wardD2/deliver.mm_top_genes.csv` —
  21 NMF meta-programs (M01-M21) independently decomposed from the meta-atlas's
  malignant-cell expression variance, each with a ranked gene list (Ensembl IDs,
  50 genes/MP in the "top50" version, larger top100/top200 versions also exist).
  This decomposition used **no input from this project's D/F/P developmental
  atlas** — it is pure NMF on CRC malignant-cell expression.
- `NMF/dual_workflow_annotation/CRC_21MP_dual_workflow_annotation.csv` —
  independent functional annotation of each MP by two different published
  classification schemes (a TNBC-derived and a 3CA-derived taxonomy). **M11's
  representative genes are `KRT19, KRT8, KRT18, ANXA2, S100A11, S100A10`**,
  annotated by both schemes as "Basal / partial EMT" / "EMT-like /
  mesenchymal-like."
- `NMF/csc_signature_jaccard/` + `NMF/revCSC_overlap_tables/` — a **prior,
  independent Jaccard-overlap validation** of all 21 MPs against a published
  mouse-intestinal-regeneration-derived CRC cancer-stem-cell subtype signature
  set (`Lgr5CSC`/`proCSC`/`revCSC`, human-ortholog-mapped, 31/42 revCSC genes
  successfully mapped). **M11 is the best- or near-best-matching MP to `revCSC`
  across every robustness check run** (unique-gene and repeated-gene top-gene
  list versions, at cutoffs 50/100/200 — `revCSC_overlap_summary.md`): e.g. at
  top50 (unique-gene version), M11 has the highest Jaccard (0.085,
  overlap genes `ANXA1;KRT18;SFN;TMSB4X;TNFRSF12A`) of all 21 MPs. `TNFRSF12A`
  is the TWEAK receptor (TWEAKR) — directly relevant to this project's broader
  SPP1-TWEAK-TWEAKR mechanistic interest, independently of the D/F/P question.

**This is a legitimate Layer-1 anchor**: M11's gene identity and its match to
the published revCSC signature were both established using zero information
from this project's D/F/P developmental-atlas construction (HDMA/GTEx/HPA/
Arutyunyan/Nature2026/etc. never enter this NMF or Jaccard-overlap computation).

**Independence verified concretely, on three separate axes, per PR #16 review
round 3 (REQUEST_CHANGES → withdrawn once verified)** — the reviewer required
direct evidence rather than a plausibility claim, so each axis was checked
against the actual files/scripts, not asserted:

1. **Construction independence**, confirmed by reading `metamodule_fnmf.s2c.crc.R`
   directly: M11 is a purely expression-only, unsupervised discovery. Per-sample
   fastNMF is run independently on each of the meta-atlas's constituent samples
   (input: the expression matrix only); the resulting per-sample NMF factors are
   compared pairwise by Jaccard similarity; hierarchical clustering (`hclust`) on
   that similarity matrix, cut (`cutree`) at a chosen module count, produces
   M01-M21. **No gene set, no fetal/placenta/trophoblast reference, and no
   revCSC signature enters this construction step at any point.**
2. **Annotation independence**: the revCSC Jaccard-overlap comparison
   (`NMF/csc_signature_jaccard/`, `NMF/revCSC_overlap_tables/`) was run *after*
   M01-M21 were already fixed by the clustering above — it asks "which of these
   21 already-frozen modules looks most like the published revCSC signature,"
   not "which genes should define a module called M11." revCSC is a post-hoc
   biological annotation of a pre-existing module, not a construction input.
3. **Chronological independence**, confirmed by file mtimes: M11's top-gene
   table (`deliver.mm_top_genes.csv`) is dated 2026-07-15; the revCSC-overlap
   validation (`revCSC_overlap_summary.md`) is dated 2026-07-22; this project's
   own `F_specific_FINAL.txt` (the frozen D/F/P signature) is dated 2026-08-14 —
   roughly a month later. D/F/P could not have influenced M11's discovery.

**Net effect on the project's causal/analytical framing** (this is the version
to carry into any write-up): CRC expression data → unsupervised discovery of
M11 → independent finding that M11 resembles the published revCSC signature →
(a month later, independently) construction of the D/F/P normal-development
framework → *only now* asking what developmental ancestry explains M11. This is
a stronger design than "pick a published Oncofetal signature and score cancer
cells with it" — the cancer data grew its own M11 state first; D/F/P is used
after the fact to interpret it, not to find it.

**One real gap found, not glossed over**: the readily-available per-cell module
score table (`NMF/viz_signature_MM_alt_clean_byscore_wardD2/addmodulescore.df.tsv`,
297,307 rows) covers 13 of 21 MPs (M01, M04, M05, M08, M10, M12, M13, M15, M16,
M17, M18, M19, M20) — **M11 is not among them** (so are M02/M03/M06/M07/M09/M14/M21;
likely dropped in a redundancy-collapse step unrelated to Oncofetal relevance,
not investigated further since M11 is the one that matters here). Rather than
searching further for a possibly-nonexistent precomputed M11 score, Layer 1
compute will score M11 directly from its own already-independently-derived gene
list (top50/top100/top200 versions available) using the same `scanpy.tl.score_genes`
+ expression-matched-null-calibration mechanics being built for D/F/P (Layer 2) —
**sharing the scoring *method* across layers is fine and intentional (consistency,
comparability); what must stay independent is the *gene content* defining
"Oncofetal," which it does: M11's gene list has no construction-time dependency
on D/F/P.**

**297,307-cell subset provenance, checked directly (not assumed)**: this is
**not** a random slice of the 665,473-cell atlas. Every one of the 297,307
cells maps to the h5ad's `obs_names` by exact barcode string (100% overlap,
verified directly: `n_ids in modulescore file = 297307`, `overlap = 297307`).
It is a **malignant-only NMF discovery subset**: 198,126 `Cancer cell` +
99,181 `CRLM` (no normal/immune/polyp cells), COAD (152,712) + READ (30,300)
only, drawn from 14 of the meta-atlas's 54 constituent studies and 3 of 9
platforms (10x 5p, 10x 3p, BD Rhapsody). **Per reviewer instruction, the reason
for this narrower-than-full-atlas coverage is stated honestly as unresolved,
not guessed**: *the pre-existing NMF analysis was available for 297,307
malignant cells from 14 constituent studies and three platforms; the reason
for the restriction relative to the full atlas remains to be established.*
This does not block using M11 as an anchor — it does mean M11 status is only
defined for this 297,307-cell subset, not claimed for the full 665,473-cell
atlas.

## Revised analysis structure: two distinct cell populations, three analyses

Per reviewer guidance, **two populations must be kept explicit and not
conflated**, though they validate each other:

- **M11 analysis cohort** = the 297,307-cell malignant-only NMF discovery
  subset. Answers: what is the developmental identity of the independently
  discovered M11/revCSC meta-program?
- **Full-atlas D/F/P cohort** = all malignant cells passing gene/scoring QC
  across the full 665,473-cell atlas (not restricted to the NMF subset).
  Answers: what F/P/D developmental landscape exists across the CRC malignant
  compartment overall — independent of whether a cell has an M11 score?

**Primary question, reframed** (per review — stronger and less presuppositional
than the original "are CRC malignant cells F-dominant or P-dominant"): *What is
the developmental identity of the independently discovered CRC M11/revCSC
meta-program?* This does not presuppose the answer is P-high (placental) —
if M11 turns out to be strongly F-specific with P-specific defining a separate,
orthogonal malignant state, that is an equally (or more) valuable result: it
would mean "OncoFetal" and "OncoPlacental" are not the same thing renamed, but
two genuinely separable axes of developmental reactivation, plus whatever
D-shared core they hold in common.

1. **Primary analysis: M11/revCSC ↔ D/F/P continuous decomposition**, on the
   297,307-cell M11 analysis cohort. Correlate the (null-calibrated) M11 score
   against the (null-calibrated) D-shared / F-specific (global + per-organ) /
   P-specific scores **continuously — not binarized into M11-high/M11-low up
   front** — with patient- and study-level validation (no single
   patient/study driving the correlation). This is the direct, least-assumption
   answer to "what developmental ancestry explains M11."
2. **Secondary analysis: M11-high cells' developmental composition.** Using a
   calibrated M11-high threshold (not picked post-hoc to flatter a result),
   check whether M11-high malignant cells are preferentially P-specific, a
   particular F-lineage (especially GI/intestinal, given CRC's tissue of
   origin), D-shared, or in fact split into multiple separable developmental
   substates rather than one program.
3. **Tertiary analysis: full-atlas, M11-independent D/F/P landscape.** Score
   D/F/P across the full-atlas D/F/P cohort without reference to M11 at all.
   This can surface a state the "Oncofetal" framework itself might miss — e.g.
   a P-high/M11-low malignant population that exists outside what M11 captures,
   which would itself be a genuinely new finding, not just "M11 has placental
   genes in it."

Only after these three analyses run can the project's original Q1 be answered
with a non-circular measurement. If Layer 1's M11 anchor turns out to have some
unresolvable problem once compute starts (e.g., too few cells clear a sensible
Oncofetal-high threshold), the fallback is to rename the deliverable "CRC
developmental-program projection" rather than claim to have resolved Oncofetal
composition — per the reviewer's explicit fallback instruction.

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

**Secondary: `HTAN_CRC_progressive_plasticity`.**
Smaller but very cleanly annotated — importantly, it separates malignant cells
from *named, specific* normal epithelial subtypes (not just "normal" as one
bucket) and from non-tumor tissue, and explicitly separates Primary vs.
Metastasis. **Not yet labeled "independent replication" of the meta-atlas** —
per review, the meta-atlas's 54 constituent studies were not checked for
whether they include HTAN's own cohort, and CRC meta-atlases commonly
integrate other public studies wholesale. A **study-provenance overlap
audit** (cross-referencing HTAN's patient/sample identifiers against the
meta-atlas's `dataset`/`study_id`/`NCBI_BioProject_accession`/etc. obs
columns) is required before any claim of independent replication; until that
audit runs, HTAN is used as a secondary **technical/external validation**
dataset, not asserted to be a non-overlapping source.

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
2. **Scoring method for single-cell gene-set projection — revised after review
   (blocker fix, no compute run yet).** The original proposal was to compute
   `scanpy.tl.score_genes` (Seurat `AddModuleScore`-equivalent) separately for
   D-shared/F-specific/P-specific and compare the raw scores directly (e.g.
   "F score > P score" implying "more fetal-somatic than placental"). Reviewer
   correctly caught that this is invalid: D-shared (6 genes), F-specific
   (2,504 genes) and P-specific (78 genes) differ ~32x in size between F and P
   alone, and have fundamentally different compositional structure — F-specific
   is a union across 7 fetal organs, internally lineage-heterogeneous;
   P-specific is a single, strict whole-body-adult-depleted program.
   `score_genes` answers "is this cell enriched for this one gene set relative
   to its own matched control," not "how does this gene set's enrichment
   compare in magnitude to a different, differently-sized, differently-structured
   gene set" — raw scores across signatures share no common scale, so a naive
   comparison could reflect signature architecture, not biology.

   **Revised scoring contract**: for each dataset, independently build an
   expression-matched null distribution *per signature* (D-shared, F-specific,
   P-specific, each per-organ F module below, and Layer 1's M11 signature) —
   directly analogous to the matched-permutation-null design just approved in
   Step 5: sample many random gene panels of the same size from the same
   expression-detectability strata as the real signature, compute the same
   `score_genes`-style statistic for each, and convert every cell/sample's
   *observed* raw score into a common-scale value before any cross-signature
   comparison. **Empirical percentile within each signature's own null is the
   primary common scale** (per review: D-shared has only 6 genes, so its null
   distribution is more discrete/sparse and a z-score would be more sensitive
   to null-distribution shape than a percentile is; percentile's interpretation
   stays stable across signatures of very different size). Z-score is computed
   as a secondary sensitivity check, not the primary comparison metric. Only
   these null-calibrated, common-scale values (not raw `score_genes` output)
   may be used to say a cell/sample looks more F-like vs. P-like vs. D-like vs.
   M11(Oncofetal)-like.

   **Lineage-resolved F modules, not one monolithic F-specific score
   (blocker fix)**: Step 4 already established F-developmental is strongly
   organ-specific (`F_developmental_<Organ>.txt`, 7 separate per-organ sets
   underlying the union that became `F_specific_FINAL.txt`) — many genes
   belong to only one fetal organ. Scoring only the pooled 2,504-gene union
   risks diluting a real, narrow signal (e.g. if CRC's fetal-like reactivation
   is specifically intestinal/GI-developmental) against the other six organs'
   unrelated genes. Compute will score **both** the global F-specific set
   (kept for set-logic/completeness) **and** each of the 7 per-organ F modules
   separately (each restricted to genes also in `F_specific_FINAL.txt`, i.e.
   organ module ∩ F-specific, to stay within the frozen, P-deduplicated
   signature) — biological interpretation of "is CRC's Oncofetal state
   fetal-somatic" should center on the lineage-resolved (per-organ, especially
   GI/intestinal) modules, with the global F score used only as a coarse,
   secondary summary.
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
   do malignant cells show elevated null-calibrated D/F/P enrichment (per item 2,
   not raw scores) relative to their own patient's normal epithelial cells,
   controlling for patient identity? This is the direct test of "oncofetal
   reactivation" and avoids simple between-patient confounds.
5. **Study-provenance overlap audit (added after review)**: before calling
   `HTAN_CRC_progressive_plasticity` an independent-replication check against
   the `CRC_single_cell_atlas_2025` meta-atlas, cross-reference HTAN's
   patient/sample identifiers against the meta-atlas's per-cell `dataset`/
   `study_id`/`NCBI_BioProject_accession`/`SRA_sample_accession`/etc. obs
   columns (already inventoried this round, see table above) to check whether
   HTAN's own cohort is one of the meta-atlas's 54 constituent studies. If
   overlap is found, HTAN is downgraded from "replication" language to
   "technical/external validation" (already the working label pending this
   audit, per the dataset-plan section above) and the meta-atlas's 54-study
   count should be reported net of any HTAN overlap where relevant.

## What this step is not

Not a re-derivation of "Oncofetal" via new clustering/NMF on this project's own
compute — both the malignant-cell-type calls *and* the M11 Oncofetal anchor
(Layer 1) are existing, independently-computed, published/prior work, used
as-is, not re-derived. Not a re-calibration of Step 4's frozen D/F/P cutoffs.
Not yet a spatial-transcriptomics analysis (no spatial CRC dataset has been
inventoried in this round — single-cell only for this first pass).

Submitting this design for review before running any real gene-ID mapping or
scoring compute.
