# Step 6 gut re-anchor: primary scoring compute — results

Real compute, per the approved design (`docs/STEP6_GUT_SCORING_COMPUTE_DESIGN.md`,
PR #26, 2 review rounds). Projects the 13-panel gut D/F/P + revCSC
inventory onto all 665,473 malignant cells of `CRC_single_cell_atlas_2025`,
using null-calibrated empirical percentile scoring, and reports the
primary continuous revCSC↔D/F/P correlation with the locked
donor/study-aware validation. All numbers below are real qsub output,
pulled back and verified byte-exact (md5) against the Argos-side files
before being trusted or written up — same discipline as every prior step.

## Job provenance

| Stage | Job ID | Scope | Outcome |
|---|---|---|---|
| Gene-set build | 3620942 | 13-panel Ensembl-ID inventory | Clean; 100% resolution, 0 unmapped, exact expected panel sizes |
| Convergence check | 3620943 | 20,000-cell stratified subset, N_PERM=100 vs 500, real `scanpy.tl.score_genes` | Clean; all 13 panels' percentile correlation ≥0.986; 2/10 comparison pairs flagged for N_PERM=500 |
| Full run (1st attempt, CSR) | 3620971 | 665,473 cells | Killed intentionally after 1 panel — see "Compute-cost fixes found mid-run" below |
| Full run (2nd attempt, CSC + `score_genes_fast`) | 3620977 | 665,473 cells, all 13 panels + primary analysis | **Completed cleanly**, ~50 min total |

Full-run per-cell score files (`crc_gut_scoring_all_panels.parquet`, 55MB;
`crc_gut_scoring_cell_metadata.parquet`, 5.8MB) are kept on Argos only
(`results/06_crc_projection/gut_scoring_full/` — outside this repo's
data-file size norms, per `.gitignore`'s "no data files in this repo"
policy), reproducible from the frozen scripts + these committed
parameters. md5 for provenance: `8c942c06923ff8b8c417d67efff88110`
(`crc_gut_scoring_all_panels.parquet`), `5c7250403e7229a4f33fec131ae5434d`
(`crc_gut_scoring_cell_metadata.parquet`). The real, primary deliverable —
the donor/study-aware correlation analysis — is committed in full
(`results/06_crc_projection/gut_scoring_primary_analysis/`, 55 files, 1.2MB).

## Convergence check: N_PERM=100 confirmed adequate

Per-panel correlation between N_PERM=100 and N_PERM=500's per-cell
empirical percentiles, all 13 panels, 20,000-cell subset:

| Panel | Pearson r (100 vs 500) | Spearman ρ |
|---|---|---|
| revCSC_primary27_full | 0.9936 | 0.9867 |
| revCSC_primary27_minus_CLU | 0.9925 | 0.9865 |
| revCSC_primary27_minus_ASS1 | 0.9934 | 0.9813 |
| revCSC_primary27_minus_CLU_ASS1 | 0.9937 | 0.9861 |
| revCSC_extended28_full | 0.9910 | 0.9893 |
| revCSC_extended28_minus_CLU | 0.9930 | 0.9841 |
| revCSC_extended28_minus_ASS1 | 0.9930 | 0.9826 |
| revCSC_extended28_minus_CLU_ASS1 | 0.9928 | 0.9807 |
| D_Gut-shared | 0.9864 | 0.9843 |
| F_Gut-specific | 0.9980 | 0.9792 |
| F_Colon-specific | 0.9976 | 0.9862 |
| F_SI-specific | 0.9964 | 0.9896 |
| P_Gut-specific | 0.9933 | 0.9889 |

All well above any reasonable stability bar — confirms `N_PERM=100`'s
per-cell empirical percentile is not meaningfully different from what
500 draws would give, for the continuous-correlation use case this
compute targets (per the design's own stated limits: not intended for
fine per-cell significance calls).

Per the locked 0.02-Pearson-r gate on the 10 revCSC↔D/F/P comparison
pairs: **2/10 pairs exceeded the threshold** —
`revCSC_primary27_minus_CLU` vs `F_Colon-specific` (Δ=0.0265) and
`revCSC_extended28_minus_CLU` vs `F_Colon-specific` (Δ=0.0235). Per the
contract, this flagged 3 panels for `N_PERM=500` in the full run:
`F_Colon-specific`, `revCSC_extended28_minus_CLU`,
`revCSC_primary27_minus_CLU` (union of the flagged panels across both
pairs) — all other panels used `N_PERM=100` as planned.

## Compute-cost fixes found mid-run (not assumed, profiled directly)

Two real feasibility problems were found and fixed while this compute
was actually running, not anticipated in the design:

1. **Naive `scanpy.tl.score_genes` is infeasible at 665K-cell scale.**
   A timing probe (not a guess) showed per-call cost scaling with cell
   count (1.6s/9.6s/21.5s at 20k/100k/300k cells), extrapolating to
   ~45-55s/call at 665,473 — ~18h+ for the full 13-panel × 101-draw job.
   Traced by reading `scanpy` 1.11's `_score_genes_bins` source directly
   on Argos: it recomputes a full-genome (28,476-gene) average-expression
   binning from scratch on *every* call, regardless of the scored
   gene-list's size. Fixed with `score_genes_fast`
   (`crc_gut_scoring_core.py`): precomputes that binning once per scored
   population, reused across all calls. **Numerically validated before
   use** (`validate_score_genes_fast.py`): the selected control-gene
   *set* is byte-identical to real `scanpy.tl.score_genes`'s internal
   selection in every checked case (8/1442/2173-gene panels × real
   signature + 3 null draws each), and resulting scores `allclose` to
   1e-4 relative tolerance (the tiny residual is pure float64
   summation-order noise from a different code path, not an algorithm
   mismatch).
2. **CSR column-slicing is a second, independent bottleneck.** Even with
   `score_genes_fast`, the full run was observed running at ~10-14s/call
   *regardless of panel size* — far slower than a first (incomplete)
   timing probe run only against the largest panel had suggested for
   small panels. Traced directly: `layers['counts']` (and `X` rebuilt
   from it) is CSR (row-oriented), and `score_genes_fast`'s repeated
   column-slicing (`x[:, gene_pos]`, `x[:, ctrl_pos]`, every one of
   ~1,300+ calls) is a well-known scipy inefficiency on CSR regardless of
   column count. A targeted probe confirmed: converting `X` to CSC
   (compressed sparse column) once, before the repeated calls, cuts cost
   from 14.4s/call to 0.9s/call (27-gene panel) and ~11.9s/call to
   3.2s/call (2,173-gene panel) at full scale. One-time conversion cost
   ~80s. The first full-run attempt (job 3620971, CSR) was killed after
   confirming this and resubmitted (job 3620977, CSC) — the one panel
   already checkpointed was safely reused (CSR→CSC conversion changes no
   computed values, only access performance).

Combined effect: **the full run completed in ~50 minutes**, down from an
~18h+ naive projection.

**Minor bug found and fixed in the same pass**: `DataFrame.attrs` does
not survive a `to_parquet`/`read_parquet` round-trip, so the checkpoint-
reuse path was silently reporting `n_testable=-1` for any panel loaded
from a checkpoint instead of freshly scored — this affected exactly one
row (`revCSC_primary27_full`, reused from the killed CSR run) in
`n_testable_genes_per_panel.tsv`. Fixed in `crc_gut_scoring_core.py` by
writing `n_testable` to a small sidecar file instead of relying on
`.attrs` (applies to future runs); this run's one wrong value was
corrected using the number the original qsub log directly reported
(`n_testable=27`, confirmed real, not fabricated) — the underlying score
values themselves were never affected, only this one provenance field.

## Final scoring inventory (13 panels, all cells)

| Panel | n_testable genes | N_PERM used |
|---|---|---|
| revCSC_primary27_full | 27 | 100 |
| revCSC_primary27_minus_CLU | 26 | 500 (gated) |
| revCSC_primary27_minus_ASS1 | 26 | 100 |
| revCSC_primary27_minus_CLU_ASS1 | 25 | 100 |
| revCSC_extended28_full | 28 | 100 |
| revCSC_extended28_minus_CLU | 27 | 500 (gated) |
| revCSC_extended28_minus_ASS1 | 27 | 100 |
| revCSC_extended28_minus_CLU_ASS1 | 26 | 100 |
| D_Gut-shared | 8 | 100 |
| F_Gut-specific | 2,173 | 100 |
| F_Colon-specific | 1,442 | 500 (gated) |
| F_SI-specific | 1,437 | 100 |
| P_Gut-specific | 76 | 100 |

All 13 panels reach 100% Ensembl-ID resolution against the atlas
(per `results/06_crc_projection/gut_scoring/gut_scoring_gene_id_mapping_summary.tsv`),
consistent with PR #25's gene-ID contract.

## Primary analysis: results (real, honest — not gated on outcome)

10 revCSC↔D/F/P comparison pairs, all 665,473 cells, 529 donors
(composite `donor_key = study_id||patient_id`, 475 estimable with
`n_cells ≥ 10`; the 54 non-estimable donors are **all** small-`n_cells`
cases — 1 to 9 cells each — not zero-variance edge cases, a real, expected
long tail across a 36-study meta-atlas, never dropped, always reported
transparently):

| Comparison | Pooled Pearson r | Pooled Spearman ρ | Robust to leave-one-donor/study-out? |
|---|---|---|---|
| revCSC_primary27_full ↔ D_Gut-shared | 0.012 | 0.017 | **Yes** |
| revCSC_primary27_full ↔ P_Gut-specific | -0.031 | -0.024 | **Yes** |
| revCSC_primary27_minus_CLU ↔ F_Colon-specific | 0.012 | -0.017 | **No** |
| revCSC_primary27_minus_ASS1 ↔ F_SI-specific | 0.100 | 0.013 | **No** |
| revCSC_primary27_minus_CLU_ASS1 ↔ F_Gut-specific | 0.119 | 0.004 | **No** |
| revCSC_extended28_full ↔ D_Gut-shared | 0.022 | 0.026 | **Yes** |
| revCSC_extended28_full ↔ P_Gut-specific | -0.036 | -0.027 | **Yes** |
| revCSC_extended28_minus_CLU ↔ F_Colon-specific | 0.046 | 0.019 | **No** |
| revCSC_extended28_minus_ASS1 ↔ F_SI-specific | 0.151 | 0.036 | **Yes** |
| revCSC_extended28_minus_CLU_ASS1 ↔ F_Gut-specific | 0.176 | 0.045 | **Yes** |

**Reported as real, honest numbers — not smoothed over**:

- **All 10 correlations are weak in absolute magnitude** (|r| ≤ 0.18).
  Whatever revCSC↔D/F/P relationship exists at the single-cell level
  across the full, heterogeneous 665,473-cell atlas is not a strong,
  dominant signal — this compute does not find a large "Oncofetal = pure
  F-program" or "= pure P-program" effect at this scale/resolution.
- **D and P comparisons are the weakest and are the ones that ARE
  robust** to single-donor/study removal (|r| ≤ 0.036 throughout, no
  sign flip on leave-one-out). This is consistent with the negligible
  gene-overlap finding from PR #25's audit (D/P axes share zero genes
  with revCSC) — a weak, stable, near-null association is the expected,
  clean result here.
- **F comparisons are directionally consistent (all positive) and
  somewhat larger** (r = 0.01 to 0.18), strongest for
  `F_Gut-specific`/`F_SI-specific` with the extended revCSC variant
  (r = 0.15-0.18) and weakest for `F_Colon-specific` (r = 0.01-0.05)
  despite `F_Colon-specific` being the *primary* regional F axis per
  Step 4a's locked hierarchy — a real, not obviously expected, pattern.
- **4 of the 6 F-comparison pairs are NOT robust** to
  leave-one-donor-or-study-out: all 3 pairings using the primary
  (27-gene) revCSC (`F_Colon-specific`, `F_SI-specific`, `F_Gut-specific`)
  **plus** the `extended28`↔`F_Colon-specific` pairing — not "primary
  revCSC" alone, corrected from an earlier wrong count/attribution in
  this doc. Only `extended28`'s pairings with `F_SI-specific` and
  `F_Gut-specific` are robust among the F comparisons; all 4 D/P pairs
  (both revCSC variants) are robust. For `revCSC_primary27_minus_CLU` ↔
  `F_Colon-specific` specifically,
  excluding the single study `Terekhanova_2023_Nature` shifts the pooled
  Pearson r by -0.052 — the largest single-study effect of any pair,
  more than double the pair's own pooled r (0.012). This means the
  already-weak `F_Colon-specific` correlation is itself substantially
  driven by one study's contribution, not a consistent cross-cohort
  signal. This instability is the finding for these 4 pairs, not
  something to interpret as "revCSC correlates with F" — per the
  design's own stated contract, non-robust pairs are reported as
  unstable, not as evidence of association.
- Equal-donor-weighted vs. cell-weighted within-study summaries and the
  full per-donor tables are committed in full
  (`results/06_crc_projection/gut_scoring_primary_analysis/`) for anyone
  wanting to inspect a specific study or donor's contribution directly.

**What this does NOT show**: a strong, robust, single-axis "Oncofetal is
dominated by the fetal-gut program" or "...by the placental program"
result. The honest read of this primary compute is a weak, largely
donor/study-heterogeneous association, with the D/P axes' near-null,
stable signal being the cleanest result and the F axes' larger-but-
less-stable signal needing the secondary (revCSC-high subset) and
tertiary (full-atlas revCSC-independent) analyses — explicitly out of
scope for this PR — before any stronger claim is defensible.

## What this does NOT do

Same explicit scope boundary as the approved design: no secondary
(revCSC-high developmental composition + M11 concordance) or tertiary
(full-atlas revCSC-independent) analyses; no other CRC datasets
(`HTAN_CRC_progressive_plasticity`, `CRLM_NMP_ATLAS`); no re-derivation
of any frozen gut D/F/P or revCSC gene set.

Submitting for compute review before merge, same discipline as every
prior step.
