# Step 6 gut re-anchor: primary scoring compute design

Implements the **primary analysis** from the already-approved Step 6
design (`docs/STEP6_CRC_PROJECTION_DESIGN.md`, section "Revised analysis
structure", item 1): the continuous revCSC ↔ D/F/P decomposition, scored
on `CRC_single_cell_atlas_2025` (the primary dataset, 665,473 malignant
cells), using the gut-specific Layer 2 contract locked in PR #25
(`docs/STEP6_GUT_REANCHOR_DELTA.md`). The design's method itself
(null-calibrated percentile scoring, overlap-exclusion, donor/study-aware
aggregation) is unchanged and already approved — this doc locks the
concrete implementation parameters needed to make a 665K-cell compute job
tractable, submitted for review before any qsub job runs.

## Scope of this compute (explicitly bounded)

- **Dataset**: `CRC_single_cell_atlas_2025` only
  (`/home/zz950/DATA/scRNAseq/meta_study/CRC_single_cell_atlas_2025/adata_nmf.h5ad`).
  Secondary/tertiary datasets (`HTAN_CRC_progressive_plasticity`,
  `CRLM_NMP_ATLAS`) and the secondary/tertiary analyses (revCSC-high
  cells' developmental composition + M11 concordance; full-atlas
  revCSC-independent landscape) are **out of scope**, separate follow-on
  work.
- **Cells**: all 665,473 obs — confirmed directly (not assumed) that
  `atlas_cell_type_middle` has exactly 2 values (`Cancer cell` 509,421 +
  `CRLM` 156,052 = 665,473 = `n_obs` exactly), so this file is *already*
  malignant-cell-only; no further cell-type filtering needed.
- **Signatures scored** (7 total, all null-calibrated independently):
  revCSC primary (27 genes) and 3 overlap-adjusted variants (see below);
  `D_Gut-shared` (8), `F_Gut-specific` (2,192), `F_Colon-specific` (1,451,
  primary regional), `F_SI-specific` (1,452, secondary regional),
  `P_Gut-specific` (76). `F_Gut-core` and the individual regional D/P
  sets are **not** scored here — per PR #25's own finding, they are
  descriptive/tertiary, not part of the primary Layer 2 axes.

## Gene-ID mapping (reuses PR #25's contract exactly, no new lookup)

`CRC_single_cell_atlas_2025`'s `var_names` are bare Ensembl gene IDs
(confirmed: `ENSG00000186092`, etc., no version suffix) — the **same
bare-Ensembl-ID space** as `var_id_map.tsv`'s `gene_id` column and
`revCSC_human_FINAL.tsv`'s `human_ensembl_id` column, both already used
in PR #25. Every gut D/F/P gene set is mapped `var_name → gene_id` via
`var_id_map.tsv` (unchanged from PR #25); revCSC is used directly via its
`human_ensembl_id`. Each signature's final testable gene list is
`{signature Ensembl IDs} ∩ {atlas var_names}`, with `n_input` /
`n_mapped_to_ensembl` / `n_present_in_atlas` reported for every signature
— same auditable-denominator discipline as every prior gene-ID contract
this project has used.

## revCSC overlap-exclusion variants (per PR #25's locked contract)

Four revCSC scoring variants, used only where the target signature
requires exclusion (per `docs/STEP6_GUT_REANCHOR_DELTA.md`):

| Variant | Genes | Used against |
|---|---|---|
| `revCSC_full27` | 27 | `D_Gut-shared`, `P_Gut-specific` (zero overlap, no exclusion needed) — also the standing sensitivity check reported alongside every comparison |
| `revCSC_minus_CLU` (26) | 27 − `CLU` | `F_Colon-specific` |
| `revCSC_minus_ASS1` (26) | 27 − `ASS1` | `F_SI-specific` |
| `revCSC_minus_CLU_ASS1` (25) | 27 − `CLU` − `ASS1` | `F_Gut-specific` |

## Scoring method — null-calibrated percentile (per the approved design, parameters locked here)

For each signature (real gene set) and each of `N_PERM` matched-null
draws: bin all atlas genes into `N_BINS` expression-detectability strata
(fraction of cells with `counts > 0`, computed once from `layers['counts']`,
the same real, measured covariate class already used for every
permutation-null design this project has built — Step 5, Step 4a/PR#22,
PR#24). Each null draw samples one random gene per real-signature gene
from the same stratum (without replacement within a draw), then computes
`scanpy.tl.score_genes` for both the real signature and each null draw on
the full 665,473-cell object. Every cell's **empirical percentile** among
its own `N_PERM` null draws is the primary common-scale value (per the
approved design's own stated reasoning: percentile is stable across
signatures of very different sizes — `D_Gut-shared` at 8 genes vs.
`F_Gut-specific` at 2,192 — where a z-score would be more sensitive to
null-distribution shape for the smallest sets). Raw `score_genes` z-score
is retained as a secondary sensitivity value, not the primary metric.

**`N_PERM = 100`** (reduced from Step 5/Step 4a's `N_PERM = 500`
precedent) — a real, stated tradeoff: this compute scores 7 signatures ×
(1 real + 100 null draws) × 665,473 cells × up to 2,192 genes, an order
of magnitude larger than any prior permutation job this project has run.
100 draws still gives an empirical-percentile floor of ~1% per cell,
sufficient for a continuous correlation analysis (the primary use here)
rather than a strict per-cell significance call. If a later analysis
needs finer per-cell resolution, `N_PERM` can be raised in a follow-up
run — flagged explicitly rather than silently under-resolved.
**`N_BINS = 20`** (matches the Step 4a/PR#22 precedent exactly).
**Fixed seed `20260815`**, stated for reproducibility, same convention
as every prior permutation script this project has built.

## Primary analysis (this compute's deliverable)

Continuous correlation: each cell's null-calibrated revCSC percentile
(using the correctly-overlap-excluded variant per comparison) against
each cell's null-calibrated D/F/P percentile — Pearson + Spearman,
computed **within-study** first (36 `study_id` values) and pooled only
after confirming no single study/patient drives the pooled correlation
(per the approved design's explicit donor/study-aware requirement and
Step 5's cross-platform-confound lesson). `platform` (9 values) and
`patient_id` are carried as covariates in the reported output, not
collapsed away. Reported as real, honest numbers — including if the
correlation is weak, null, or inconsistent across studies; this compute
does not gate or reshape the frozen D/F/P or revCSC sets regardless of
outcome.

## What this does NOT do

- Does not touch the secondary (revCSC-high developmental composition +
  M11 concordance) or tertiary (full-atlas revCSC-independent landscape)
  analyses — separate, later work.
- Does not run on `HTAN_CRC_progressive_plasticity` or `CRLM_NMP_ATLAS`
  — separate, later work.
- Does not re-derive or re-calibrate any frozen gut D/F/P or revCSC gene
  set — purely applies them.
- Does not use `F_Gut-core` or the individual regional D/P sets as
  scoring inputs — per PR #25, they are descriptive/tertiary only.

Submitting for review before any qsub job runs, same discipline as every
prior step.
