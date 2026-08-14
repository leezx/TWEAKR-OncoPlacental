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

## Method (revised per round-1 review: donor-aware, not organ-pooled)

The first draft pooled raw counts directly by `(organ, cell_ontology_class)`,
discarding the donor dimension entirely. Caught by review as a real
inconsistency with the donor/sample-aware discipline used throughout
Step 4 (HDMA per-sample pseudobulk, P-developmental's per-donor pairing):
donor counts per organ are real but small (checked directly, not
assumed: Liver 2, Skin 2, Spleen 3, Thymus 2, Large_Intestine 2), so
pooling across donors lets a single donor drive an apparent "detected"
hit with no way to tell whether it's real replicated biology or one
donor's idiosyncrasy — cell count alone (e.g. 20+ cells) doesn't
guarantee independent replication if all 20 cells come from one donor.

- **Primary unit: `(organ, donor, cell_ontology_class)` pseudobulk.** Sum
  raw counts (`layers['raw_counts']`) within each donor×cell-type group.
  CPM computed **per donor-celltype pseudobulk's own library size**
  independently (not "normalized within organ" — organ is the comparison
  *scope*, each donor-celltype library is its own normalization
  denominator; this doc's first draft conflated the two).
- **Minimum cell-count floor applied at the donor×cell-type level**
  (proposed: ≥20 cells), not after organ-level pooling — a cell type with
  20+ cells that all come from a single donor must not look the same as
  20+ cells genuinely split across 2 donors.
- **Aggregate to cell-type level for reporting**, keeping donor structure
  visible rather than collapsing it into one pooled number:
  - how many eligible donors contribute this cell type;
  - in how many of those donors the gene clears CPM≥1 (donor-level
    detection fraction);
  - median/range CPM across contributing donors.
  - **Cell types backed by only 1 donor are explicitly labeled
    `single-donor evidence`** and never presented as equal-confidence to
    multi-donor evidence (most Tabula Sapiens organs here have only 2
    donors total, so many cell types will land in this category — stated
    honestly, not hidden in a summary average).
- **Detection criterion**: same not-detected-floor (CPM<1) design used
  throughout Step 4, for consistency.
- **Organ-matched check (F-specific genes)**: for each of the 4
  HDMA-organ-matching Tabula Sapiens files (Liver, Skin, Spleen, Thymus),
  check each organ's own frozen `F_developmental_<Organ>.txt` gene list
  against that organ's real per-donor-celltype detection — report what
  fraction of (gene, cell type) pairs show unexpected detection (with
  donor-level detection fraction, not a flattened count), and flag any
  specific gene/cell-type combination that looks surprisingly high
  (candidate false positive worth re-examining, not silently averaged
  away in a summary statistic).
- **Whole-body-style check (D-shared and P-specific genes, both organ-
  agnostic by construction)**: evidence unit is **gene × cell type ×
  donor** across all 5 organs — not 5 organs' cell types treated as one
  flat pile of equal-weighted pseudobulks. Large_Intestine included as
  the CRC-adjacent adult reference but must not get disproportionate
  statistical weight relative to the other 4 organs just because it's
  more relevant to the eventual application.
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
