# Step 3 prep — collision aggregation decision + adult reference acquisition

## 1. Symbol collision aggregation rule — decided, quantified not assumed

PR #4's `canonical_feature_map` found 29 symbol collisions across the 7 HDMA organs (cases where an unresolved ENSG ID's HGNC/Ensembl-native symbol already exists as a separate native-symbol feature in the same object — e.g. `PDE8B` present both as its own symbol row and as `ENSG00000284762`). The reviewer explicitly deferred the aggregation-rule decision to Step 3 rather than requiring it in PR #4.

**Question**: does it matter whether Step 3 sums counts across colliding features vs. picks one arbitrarily? Quantified via `scripts/03_pseudobulk_prep/quantify_collision_count_mass.sh` (qsub job 3620311, all 7 organs) rather than assumed:

| organ | total_counts | collision_counts | n_collision_features | fraction |
|---|---|---|---|---|
| Adrenal | 26,762,748 | 833 | 6 | 3.11e-05 |
| Thyroid | 16,436,601 | 1,625 | 6 | 9.89e-05 |
| Spleen | 57,991,267 | 4,294 | 10 | 7.40e-05 |
| Thymus | 145,949,688 | 7,655 | 8 | 5.24e-05 |
| Liver | 158,772,055 | 11,380 | 8 | 7.17e-05 |
| Skin | 132,856,857 | 9,227 | 10 | 6.95e-05 |
| StomachEsophagus | 240,435,683 | 8,903 | 10 | 3.70e-05 |

Colliding features account for **0.003%–0.01% of an organ's total RNA counts** — several orders of magnitude below anything that could materially shift a pseudobulk profile regardless of aggregation choice.

**Decision**: Step 3's pseudobulk construction will **sum counts across colliding original_features** before mapping to `canonical_symbol` (the standard convention, matching how Cell Ranger/STARsolo handle duplicate gene symbols at the reference-building stage). Given the negligible count mass, this is a safe default — not a consequential modeling choice — and doesn't need further debate or a more elaborate rule (e.g. picking the higher-expressed duplicate, or annotation-based selection).

Full per-organ output: `Adrenal_collision_mass.tsv` … `StomachEsophagus_collision_mass.tsv`, combined in `all_organs_collision_mass.tsv`.

## 2. Adult normal reference — Tier-1 (bulk) acquired, Tier-2 (single-cell) partially acquired

Per project discussion (2026-08-13), the adult reference the ChatGPT reviewer flagged as a Step-3 prerequisite is designed as two tiers, not one dataset:

- **Tier 1 — mandatory, used to define D/F/P signatures now**: `[[GTEx_v11_median_tpm]]` (74,628 genes × 68 adult tissues) + `[[HPA_RNA_tissue_consensus]]` (20,162 genes × 40 adult tissues). Two different comparisons, not one:
  - **F-specific**: organ-matched fetal (HDMA) vs. adult (GTEx/HPA) — same organ, developmental persistence question.
  - **P-specific**: trophoblast vs. **all** adult tissues (whole-body), not just the 7 HDMA-matched organs — avoids wrongly calling a gene "placenta-specific" just because it's absent from the 7 organs happening to have fetal data, when it might be high in some unrelated adult tissue (e.g. testis).
  - **CRC sanity check**: normal adult colon/rectum columns, since the eventual target is CRC malignant epithelial cells.
  - GTEx has **no Thymus tissue** (donor population skews older, thymus involutes) — HPA fills that specific gap; HPA also serves as a general cross-check on the other 6 organs.
- **Tier 2 — deferred, independent validation only, NOT used to define signatures**: `[[TabulaSapiens]]` (adult single-cell, organ-matched subset: Liver/Skin/Spleen/Thymus/Large_Intestine, ~2.85GB, not the full 15.6GB atlas). Held out of signature *definition* entirely so "P-specific genes are absent across real adult cell types" can be reported as independent post-freeze validation rather than a circular check. Downloaded now (pre-fetched) but not to be touched by Step 3's signature-defining code.
  - Confirmed gap: Tabula Sapiens has no Adrenal, Thyroid, or Stomach/Esophagus organ file — GTEx/HPA remain the only adult reference for those 3 HDMA organs; no cell-type-level validation available for them.

Dataset details, verification, and provenance: `DATA/1.Databases/GTEx_v11_median_tpm/link.md`, `DATA/1.Databases/HPA_RNA_tissue_consensus/link.md`, `DATA/1.Databases/TabulaSapiens/link.md`.

## Open items carried into Step 3 proper

- No GTEx/Tabula-Sapiens adult reference for HDMA's Thymus other than HPA (single source, no cross-validation for that one organ's F-specific comparison).
- No Tabula Sapiens cell-type validation available for Adrenal/Thyroid/StomachEsophagus.
- Tier-2 (Tabula Sapiens) validation logic itself is not yet built — deferred until after a first D/F/P signature version exists to validate.
