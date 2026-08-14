# Step 5 design: Tier-2 validation of the frozen D/F/P signature against real adult single cells

Per the PR #13 reviewer's explicit next step: independent validation
using Tabula Sapiens (adult single-cell atlas), held out of Step 4's
signature *definition* specifically so this check is not circular. Design
submitted before any qsub compute, same discipline as every prior step.

## What this validates, and why it's a real independent check

Step 3/4's `adult_excluded` used GTEx/HPA **bulk** tissue data — a
tissue-level average. A gene could look "low" in bulk Liver while
actually being highly expressed in one specific, rare cell type within
Liver that gets diluted out by averaging across the whole tissue's cell
mixture. Tabula Sapiens is single-cell, annotated by real cell type
(`cell_ontology_class`, standardized Cell Ontology terms) — this is the
first check of whether the frozen D-shared/F-specific/P-specific genes
are genuinely absent across real adult **cell types**, not just low on
average across a tissue's mixed-cell bulk signal.

## What's actually on disk (checked directly, `results/05_tier2_validation/inventory/`)

5 organs downloaded: Liver, Skin, Spleen, Thymus, Large_Intestine.
Confirmed via direct inventory (not assumed):

- **Gene IDs are native symbols already** (e.g. `DDX11L1`, `WASH7P`) —
  matches HDMA's `canonical_symbol` convention directly, no remapping
  needed.
- **Raw counts confirmed available**: `X` is normalized (only ~1-2% of
  nonzero values integer-valued), but `.raw.X` and `layers['raw_counts']`
  are both **100% integer-valued** in every organ — genuine raw counts,
  checked directly per organ, not assumed from one file.
- **Cell-type annotation**: `cell_ontology_class` (standardized, e.g.
  `hepatocyte`, `cd4-positive helper t cell`) is the primary column to
  use for cross-organ consistency; `free_annotation` is a secondary,
  dataset-specific label; `compartment` (epithelial/immune/endothelial/
  stromal) is a coarser grouping.
- **Real cell-type diversity per organ**: Liver 13 types (hepatocyte,
  macrophage, endothelial, etc.), Skin 25, Spleen and Thymus large
  immune-cell repertoires (memory B cell, various T cell subsets,
  thymocyte stages), Large_Intestine includes real epithelial types
  (enterocyte, goblet cell, Paneth cell) relevant to the eventual CRC
  application.
- **Confirmed gap** (already documented in `datasets/TabulaSapiens/dataset.md`,
  reconfirmed here): no Adrenal, Thyroid, or Stomach/Esophagus organ file
  exists in Tabula Sapiens. **F-specific genes from those 3 HDMA organs
  cannot get Tier-2 single-cell validation from this source** — stated
  explicitly, not glossed over. GTEx/HPA bulk remains the only adult
  reference for those 3 organs.

## Method

- **Pseudobulk per (organ, cell_ontology_class)**: sum raw counts
  (`layers['raw_counts']`) across all cells sharing a `cell_ontology_class`
  label within each organ, CPM-normalize within organ (same
  never-cross-platform-raw-magnitude discipline as every prior step —
  comparing within Tabula Sapiens' own cell-type distribution, not against
  HDMA/GTEx/HPA magnitudes). Cell-type groups below a minimum cell-count
  floor (proposed: ≥20 cells, to avoid noisy pseudobulk from a handful of
  cells) excluded, reported not silently dropped.
- **Detection criterion**: same not-detected-floor (CPM<1) design used
  throughout Step 4, for consistency.
- **Organ-matched check (F-specific genes)**: for each of the 4
  HDMA-organ-matching Tabula Sapiens files (Liver, Skin, Spleen, Thymus),
  check each organ's own frozen `F_developmental_<Organ>.txt` gene list
  against that organ's real per-cell-type detection — report what
  fraction of (gene, cell type) pairs show unexpected detection, and flag
  any specific gene/cell-type combination that looks surprisingly high
  (candidate false positive worth re-examining, not silently averaged
  away in a summary statistic).
- **Whole-body-style check (D-shared and P-specific genes, both organ-
  agnostic by construction)**: check across **all 5 organs' cell types
  combined** (not just organ-matched) — since P-developmental's own
  `adult_excluded` was whole-body, not organ-restricted, its Tier-2 check
  should be too. Large_Intestine included here specifically as the
  CRC-adjacent adult reference.
- **Reporting, not pre-committing to a pass/fail threshold**: this is a
  validation pass, not a re-calibration of Step 4's already-frozen
  cutoffs. Report real detection rates per gene set, flag individual
  genes with concerning single-cell-level detection for manual review,
  but do not silently drop genes from the frozen signature based on this
  alone — any removal would be a separate, reviewed decision.

## What this does NOT do

- Does not touch or re-open Step 4's frozen gene sets — this is
  validation, run after the freeze, as designed from the start (PR #5's
  two-tier reference design).
- Does not attempt Adrenal/Thyroid/Stomach single-cell validation (no
  data source) — flagged as a known, permanent gap for this project's
  current data scope, not something to silently work around.
- Does not yet touch CRC Oncofetal data — that's the step after this one.
