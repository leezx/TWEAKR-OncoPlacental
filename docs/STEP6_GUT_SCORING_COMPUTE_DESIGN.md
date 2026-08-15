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
- **Signatures scored** (13 distinct gene sets total, all null-calibrated
  independently — corrected count, see "revCSC overlap-exclusion
  variants" below for the full accounting): 8 revCSC panels (27-gene
  primary + 28-gene extended/sensitivity, each in a full and up to 3
  overlap-excluded forms) and 5 D/F/P panels: `D_Gut-shared` (8),
  `F_Gut-specific` (2,192), `F_Colon-specific` (1,451, primary regional),
  `F_SI-specific` (1,452, secondary regional), `P_Gut-specific` (76).
  `F_Gut-core` and the individual regional D/P sets are **not** scored
  here — per PR #25's own finding, they are descriptive/tertiary, not
  part of the primary Layer 2 axes.

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

## revCSC scoring inventory: 8 panels (corrected — restores the frozen extended sensitivity variant)

**Round 1 correction**: the previous version of this doc scored only the
27-gene primary revCSC and dropped the 28-gene extended/sensitivity
variant (27 `one2one` + `Ly6a`). That was a real gap — the approved Step
6 design (`docs/STEP6_CRC_PROJECTION_DESIGN.md`, PR #17 round 2) locks
the extended variant as "reported alongside every primary result, not
optionally," and PR #25's delta doc restates this carries over unchanged
("27-gene primary/1-gene extended-only sensitivity addition"). Restored
here. PR #25's own overlap audit already checked the extended-only gene
(`Ly6a`→`LY6A`) against every gut D/F/P set and found **zero** overlap
everywhere (`overlap_extended_only_n = 0` in
`results/06_crc_projection/revcsc_gut_overlap_audit/revcsc_gut_dfp_overlap_audit.tsv`)
— so the extended variant needs exactly the same overlap-exclusion
pattern as the primary variant (dropping only `CLU`/`ASS1` where the
primary does), never an `Ly6a`/`LY6A` exclusion.

| Variant | Genes | Used against |
|---|---|---|
| `revCSC_primary27_full` | 27 | `D_Gut-shared`, `P_Gut-specific` (zero overlap, no exclusion needed) — also the standing sensitivity check reported alongside every comparison |
| `revCSC_primary27_minus_CLU` | 26 | `F_Colon-specific` |
| `revCSC_primary27_minus_ASS1` | 26 | `F_SI-specific` |
| `revCSC_primary27_minus_CLU_ASS1` | 25 | `F_Gut-specific` |
| `revCSC_extended28_full` | 28 | `D_Gut-shared`, `P_Gut-specific` — reported alongside every primary-variant result on these two axes, per the locked extended-sensitivity contract |
| `revCSC_extended28_minus_CLU` | 27 | `F_Colon-specific` — alongside `revCSC_primary27_minus_CLU` |
| `revCSC_extended28_minus_ASS1` | 27 | `F_SI-specific` — alongside `revCSC_primary27_minus_ASS1` |
| `revCSC_extended28_minus_CLU_ASS1` | 26 | `F_Gut-specific` — alongside `revCSC_primary27_minus_CLU_ASS1` |

**Total scored gene sets this compute: 13** (8 revCSC panels above + 5
D/F/P panels: `D_Gut-shared`, `F_Gut-specific`, `F_Colon-specific`,
`F_SI-specific`, `P_Gut-specific`). The compute-budget estimate below is
recalculated against this corrected count, not the previous "7
signatures."

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
null-distribution shape for the smallest sets).

**Round 1 correction**: the previous version of this doc mislabeled the
secondary sensitivity value as "raw `score_genes` z-score" — this is
wrong, `scanpy.tl.score_genes` does not output a z-score. The approved
contract's secondary value is a **null-calibrated z-score**: for each
cell, `(observed_score - mean(null_draw_scores)) / std(null_draw_scores)`,
i.e. the observed real-signature score standardized against that same
cell's `N_PERM`-draw permutation-null distribution — not a transform of
the real score alone. Corrected here; both the empirical percentile
(primary) and this null-calibrated z-score (secondary sensitivity value)
are computed from the same `N_PERM` null draws per signature, no
additional compute cost.

**`N_PERM = 100`** (reduced from Step 5/Step 4a's `N_PERM = 500`
precedent) — a real, stated tradeoff: this compute scores 13 signatures ×
(1 real + 100 null draws) × 665,473 cells × up to 2,192 genes (per-panel
gene counts vary 8–2,192; see the 13-panel inventory above), an order of
magnitude larger than any prior permutation job this project has run.
100 draws still gives an empirical-percentile floor of ~1% per cell,
sufficient for a continuous correlation analysis (the primary use here)
rather than a strict per-cell significance call. If a later analysis
needs finer per-cell resolution, `N_PERM` can be raised in a follow-up
run — flagged explicitly rather than silently under-resolved.

**Round 1 addition — convergence check (required before trusting
`N_PERM=100` at full scale)**: on a fixed representative subset (stated
here: 20,000 cells, stratified-sampled proportionally by `study_id` to
preserve the 36-study composition, drawn with the same fixed seed
`20260815`), run the full null-calibration pipeline at both `N_PERM=100`
and `N_PERM=500` for all 13 signatures. Report, per signature: the
correlation between the two runs' per-cell empirical percentiles
(expected ~1.0 if 100 draws is adequate), and the resulting
Pearson/Spearman D/F/P↔revCSC correlations computed at each `N_PERM`
side-by-side. If any signature's two `N_PERM` settings produce
materially different pooled correlations (threshold: Pearson r differs
by >0.02, stated here as a concrete, checkable number rather than left
open), that signature is re-scored at `N_PERM=500` for the full
665,473-cell run; otherwise `N_PERM=100` is used throughout as planned.
This check runs first, cheaply, before the full-scale job — its own
result gates whether any signature needs the larger draw count.

**`N_BINS = 20`** (matches the Step 4a/PR#22 precedent exactly).
**Fixed seed `20260815`**, stated for reproducibility, same convention
as every prior permutation script this project has built.

## Primary analysis (this compute's deliverable)

Continuous correlation: each cell's null-calibrated revCSC percentile
(using the correctly-overlap-excluded variant per comparison, both
primary and extended reported alongside each other per the 8-panel
inventory above) against each cell's null-calibrated D/F/P percentile —
Pearson + Spearman, computed **within-study** first (36 `study_id`
values) and pooled only after the donor/study-aware validation below
confirms no single study/patient drives the pooled correlation (per the
approved design's explicit donor/study-aware requirement and Step 5's
cross-platform-confound lesson).

**Round 1 correction**: the previous version of this doc said pooling
happens "after confirming no single study/patient drives" the result but
did not specify how that confirmation is made, and treated carrying
`platform`/`patient_id` as reported covariates as if that were itself
donor-aware inference — it is not; hundreds to thousands of cells from
one donor remain pseudoreplicates that can dominate a cell-level
correlation regardless of what covariates are reported alongside it.
Corrected to three concrete, executable outputs, computed for every
revCSC↔D/F/P comparison in the 13-panel inventory:

**Round 2 correction — donor unit must be composite, not bare
`patient_id`**: the atlas is a 36-study meta-atlas; bare `patient_id`
values (e.g. `P1`, `Patient01`) are set by each constituent study
independently and are not documented as globally unique across studies.
Using bare `patient_id` as the donor unit risks silently merging
unrelated patients from different studies who happen to share an ID
string into one "patient" — which would corrupt exactly the per-patient
and leave-one-patient-out safeguards this design adds. Fixed by defining
the donor unit as the composite `donor_key = (study_id, patient_id)`
everywhere below (per-patient table, equal-donor study summaries, and
leave-one-donor-out sensitivity) — this is safe regardless of whether
`patient_id` happens to already be globally unique, so no separate
uniqueness audit is needed as a prerequisite.

1. **Per-donor correlation table**: Pearson + Spearman computed
   separately within each `donor_key` (`study_id`, `patient_id`
   composite), reported with that donor's `n_cells`. Donors whose
   correlation is mathematically non-estimable
   (e.g. `n_cells` too small, or zero variance in either variable within
   that donor) are reported as `NOT_ESTIMABLE` with the reason stated
   — never silently dropped, and never filtered based on the resulting
   correlation value.
2. **Equal-donor-weighted study summaries**: within each `study_id`,
   the study-level correlation is the unweighted mean of its
   per-donor correlations (from #1) — since `donor_key` is
   study-scoped, this is unambiguous — not a cell-pooled correlation
   that implicitly weights each donor by their cell count. The
   cell-pooled within-study correlation is also reported alongside it,
   explicitly labeled as cell-weighted, so the two are never conflated.
3. **Leave-one-out pooled sensitivity**: the full-cohort pooled
   Pearson/Spearman recomputed once per left-out `donor_key` and once
   per left-out `study_id` (leave-one-donor-out, leave-one-study-out),
   reporting the resulting correlation's range and which single
   donor/study (if any) shifts it most. A result is only reported as
   "robust to no single donor/study" if this range stays qualitatively
   stable (same sign, no threshold crossing against the pooled value);
   if it does not, that instability is reported as the finding, not
   smoothed over.

`platform` (9 values) is retained as a reported covariate (not folded
into the leave-one-out procedure itself, since it is a technical rather
than biological grouping) so any platform-correlated pattern remains
visible in the output. Reported as real, honest numbers — including if
the correlation is weak, null, inconsistent across studies, or sensitive
to a single donor; this compute does not gate or reshape the frozen
D/F/P or revCSC sets regardless of outcome.

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

## Review history

- **Round 1 (REQUEST_CHANGES)**: reviewer confirmed `N_PERM=100` and
  empirical percentile as the primary common-scale metric are both
  defensible and consistent with the approved Step 6 design — no change
  to either. Four real fixes required and applied: (1) restored the
  frozen 28-gene extended/sensitivity revCSC variant, dropped from the
  original draft's scoring inventory; (2) corrected the scored-gene-set
  count from "7 signatures" to the accurate 13, and recalculated the
  compute-budget language accordingly; (3) corrected the secondary
  sensitivity metric's definition from "raw `score_genes` z-score" (not
  a real thing `score_genes` produces) to the approved contract's actual
  null-calibrated z-score; (4) replaced the vague "confirming no single
  study/patient drives the pooled correlation" with three concrete,
  executable outputs (per-patient correlation table, equal-patient-
  weighted study summaries, leave-one-patient/study-out sensitivity) plus
  an `N_PERM=100` vs. `N_PERM=500` convergence check on a fixed
  representative subset.
- **Round 2 (REQUEST_CHANGES, one blocker)**: reviewer re-cloned the
  repo at head `0718765` and confirmed all four round-1 fixes present.
  One remaining real issue: the round-1 donor-aware validation used bare
  `patient_id`, which is not documented as globally unique across this
  36-study meta-atlas — a real cross-study ID-collision risk that could
  silently merge unrelated patients into one "donor," corrupting the
  very safeguard being added. Fixed by defining the donor unit as the
  composite `donor_key = (study_id, patient_id)` throughout (per-donor
  table, equal-donor study summaries, leave-one-donor-out sensitivity) —
  safe regardless of whether `patient_id` happens to already be globally
  unique, so no separate uniqueness audit is required. `N_PERM=100`, the
  13-panel inventory, and empirical percentile as the primary metric
  were all re-confirmed with no further changes requested.
