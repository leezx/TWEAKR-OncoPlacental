# Step 6 gut re-anchor: primary scoring compute — results

Real compute, per the approved design (`docs/STEP6_GUT_SCORING_COMPUTE_DESIGN.md`,
PR #26, 2 review rounds). Projects the 13-panel gut D/F/P + revCSC
inventory onto all 665,473 malignant cells of `CRC_single_cell_atlas_2025`,
using null-calibrated empirical percentile scoring, and reports the
primary continuous revCSC↔D/F/P correlation with the locked
donor/study-aware validation. All numbers below are real qsub output,
pulled back and verified byte-exact (md5) against the Argos-side files
before being trusted or written up — same discipline as every prior step.

**This is the reproducible re-run** (round 2 of PR #27 review) — the
first run (job 3620943/3620977) used a non-deterministic per-panel RNG
seed (see "Seed-determinism fix" below) and its outputs were discarded
entirely rather than partially reused, so every number in this document
comes from the second, deterministically-seeded run.

## Job provenance

| Stage | Job ID | Scope | Outcome |
|---|---|---|---|
| Gene-set build | 3620942 | 13-panel Ensembl-ID inventory | Clean; 100% resolution, 0 unmapped, exact expected panel sizes |
| Convergence check (round 1, non-deterministic seed) | 3620943 | 20,000-cell subset | Discarded — seed bug (see below) |
| Full run (round 1, CSR, non-deterministic seed) | 3620971 | 665,473 cells | Killed intentionally after 1 panel — CSC fix found mid-run |
| Full run (round 1, CSC, non-deterministic seed) | 3620977 | 665,473 cells | Completed, but discarded once the seed bug was found in review |
| **Convergence check (round 2, deterministic seed)** | **3621023** | 20,000-cell subset, N_PERM=100 vs 500, real `scanpy.tl.score_genes` | **Clean; used for this doc** |
| **Full run (round 2, deterministic seed)** | **3621066** | 665,473 cells, all 13 panels + primary analysis | **Completed cleanly, ~40 min; used for this doc** |

Full-run per-cell score files (`crc_gut_scoring_all_panels.parquet`,
`crc_gut_scoring_cell_metadata.parquet`) are kept on Argos only
(`results/06_crc_projection/gut_scoring_full/` — outside this repo's
data-file size norms, per `.gitignore`'s "no data files in this repo"
policy), reproducible from the frozen scripts + these committed
parameters + the fixed seed. The real, primary deliverable — the
donor/study-aware correlation analysis — is committed in full
(`results/06_crc_projection/gut_scoring_primary_analysis/`, 55 files).

## Seed-determinism fix (PR #27 review round 1 — real blocker, not cosmetic)

Per-panel RNG seeds were originally derived as
`hash((seed, panel_name)) % (2**32)`. Python's built-in `hash()` on
`str`/`tuple` is process-randomized (`PYTHONHASHSEED`) by default, so
this was **not actually reproducible** from the nominal fixed seed
`20260815` across separate qsub process invocations — each individual
run's null draws were still statistically valid (a real, legitimate
permutation test), but re-running "the same" job would silently draw a
different set of null genes every time, violating the design's explicit
"fixed seed, stated for reproducibility" contract. The reviewer
correctly flagged this as a real blocker (not a style nit): the round-1
convergence gate's flagged pairs were close enough to the 0.02 threshold
that a re-run wasn't certified to reproduce the same 3-panel gate.

Fixed with a deterministic hash (`hashlib.sha256` over a plain string,
unaffected by `PYTHONHASHSEED`) — verified directly to return the
identical integer across separate Python processes before trusting it.
**Both the round-1 convergence check and full run were discarded
entirely** (not patched/reused) and re-run from scratch with the fix;
all stale Argos-side checkpoints were cleared first so nothing from the
non-deterministic run could leak into the reproducible one.

**The re-run's gate is genuinely different from the discarded round-1
gate** — direct confirmation that the seed bug was a real, not
theoretical, problem: round 1 (buggy) flagged `F_Colon-specific`,
`revCSC_extended28_minus_CLU`, `revCSC_primary27_minus_CLU`; round 2
(fixed) flags `F_Colon-specific`, `F_Gut-specific`,
`revCSC_extended28_minus_CLU`, `revCSC_primary27_minus_CLU_ASS1` — 3 of
4 panels differ. The **qualitative primary-analysis conclusion did not
change** (every pair's robust/non-robust flag is identical between the
two runs — see below), but the exact correlation values shifted by up to
~0.02-0.03, consistent with the convergence check's own finding that
these particular pairs sit near the edge of `N_PERM=100`'s resolution.

## Convergence check: N_PERM=100 confirmed adequate (reproducible run)

Per-panel correlation between N_PERM=100 and N_PERM=500's per-cell
empirical percentiles, all 13 panels, 20,000-cell subset:

| Panel | Pearson r (100 vs 500) | Spearman ρ |
|---|---|---|
| revCSC_primary27_full | 0.9916 | 0.9888 |
| revCSC_primary27_minus_CLU | 0.9941 | 0.9879 |
| revCSC_primary27_minus_ASS1 | 0.9935 | 0.9804 |
| revCSC_primary27_minus_CLU_ASS1 | 0.9929 | 0.9867 |
| revCSC_extended28_full | 0.9929 | 0.9875 |
| revCSC_extended28_minus_CLU | 0.9923 | 0.9766 |
| revCSC_extended28_minus_ASS1 | 0.9941 | 0.9847 |
| revCSC_extended28_minus_CLU_ASS1 | 0.9932 | 0.9831 |
| D_Gut-shared | 0.9891 | 0.9872 |
| F_Gut-specific | 0.9985 | 0.9809 |
| F_Colon-specific | 0.9976 | 0.9876 |
| F_SI-specific | 0.9981 | 0.9888 |
| P_Gut-specific | 0.9939 | 0.9916 |

All well above any reasonable stability bar — confirms `N_PERM=100`'s
per-cell empirical percentile is not meaningfully different from what
500 draws would give, for the continuous-correlation use case this
compute targets (per the design's own stated limits: not intended for
fine per-cell significance calls).

Per the locked 0.02-Pearson-r gate on the 10 revCSC↔D/F/P comparison
pairs:

| revCSC panel | D/F/P panel | r (N=100) | r (N=500) | \|Δr\| | Exceeds 0.02? |
|---|---|---|---|---|---|
| revCSC_primary27_full | D_Gut-shared | -0.0037 | -0.0068 | 0.0031 | No |
| revCSC_primary27_full | P_Gut-specific | -0.0394 | -0.0446 | 0.0053 | No |
| revCSC_primary27_minus_CLU | F_Colon-specific | 0.0701 | 0.0712 | 0.0010 | No |
| revCSC_primary27_minus_ASS1 | F_SI-specific | 0.1570 | 0.1475 | 0.0094 | No |
| revCSC_primary27_minus_CLU_ASS1 | F_Gut-specific | 0.1070 | 0.1367 | 0.0297 | **Yes** |
| revCSC_extended28_full | D_Gut-shared | 0.0032 | 0.0008 | 0.0024 | No |
| revCSC_extended28_full | P_Gut-specific | -0.0442 | -0.0484 | 0.0041 | No |
| revCSC_extended28_minus_CLU | F_Colon-specific | 0.1217 | 0.1012 | 0.0205 | **Yes** |
| revCSC_extended28_minus_ASS1 | F_SI-specific | 0.1900 | 0.1923 | 0.0023 | No |
| revCSC_extended28_minus_CLU_ASS1 | F_Gut-specific | 0.1839 | 0.1745 | 0.0094 | No |

**2/10 pairs exceeded the threshold**, flagging 4 panels for `N_PERM=500`
in the full run: `F_Colon-specific`, `F_Gut-specific`,
`revCSC_extended28_minus_CLU`, `revCSC_primary27_minus_CLU_ASS1` (union
of the flagged panels across both pairs) — all other panels used
`N_PERM=100` as planned.

## Compute-cost fixes found mid-run (not assumed, profiled directly)

Two real feasibility problems were found and fixed while the round-1
compute was running (both carried forward into this reproducible re-run
unchanged — they are performance fixes, not correctness fixes, and were
never in question):

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
   ~80s.

Combined effect: both reproducible-seed runs (convergence check + full
run) completed in a few hours total, down from an ~18h+ naive projection
for the full run alone.

## Final scoring inventory (13 panels, all cells, reproducible run)

| Panel | n_testable genes | N_PERM used |
|---|---|---|
| revCSC_primary27_full | 27 | 100 |
| revCSC_primary27_minus_CLU | 26 | 100 |
| revCSC_primary27_minus_ASS1 | 26 | 100 |
| revCSC_primary27_minus_CLU_ASS1 | 25 | 500 (gated) |
| revCSC_extended28_full | 28 | 100 |
| revCSC_extended28_minus_CLU | 27 | 500 (gated) |
| revCSC_extended28_minus_ASS1 | 27 | 100 |
| revCSC_extended28_minus_CLU_ASS1 | 26 | 100 |
| D_Gut-shared | 8 | 100 |
| F_Gut-specific | 2,173 | 500 (gated) |
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
| revCSC_primary27_full ↔ D_Gut-shared | 0.008 | 0.012 | **Yes** |
| revCSC_primary27_full ↔ P_Gut-specific | -0.025 | -0.015 | **Yes** |
| revCSC_primary27_minus_CLU ↔ F_Colon-specific | 0.005 | -0.029 | **No** |
| revCSC_primary27_minus_ASS1 ↔ F_SI-specific | 0.092 | 0.006 | **No** |
| revCSC_primary27_minus_CLU_ASS1 ↔ F_Gut-specific | 0.101 | -0.013 | **No** |
| revCSC_extended28_full ↔ D_Gut-shared | 0.021 | 0.022 | **Yes** |
| revCSC_extended28_full ↔ P_Gut-specific | -0.030 | -0.022 | **Yes** |
| revCSC_extended28_minus_CLU ↔ F_Colon-specific | 0.048 | 0.018 | **No** |
| revCSC_extended28_minus_ASS1 ↔ F_SI-specific | 0.190 | 0.111 | **Yes** |
| revCSC_extended28_minus_CLU_ASS1 ↔ F_Gut-specific | 0.172 | 0.030 | **No** |

**Round 2 correction (found by reviewer, real methodological gap, not
cosmetic)**: `robust_to_no_single_donor_or_study` was originally computed
from **Pearson sign stability only**, even though Spearman ρ is reported
throughout as the co-equal rank-based metric. Direct check against the
committed per-donor/per-study leave-one-out files found a real case where
this mattered: `revCSC_extended28_minus_CLU_ASS1` ↔ `F_Gut-specific`
stayed Pearson-sign-stable under leave-one-out but its Spearman ρ flips
sign when the same influential study (`Terekhanova_2023_Nature`) is
excluded — a pair reported "robust" that is not actually rank-stable.
Fixed: `robust` now requires same-sign stability for **both** Pearson
and Spearman, both leave-one-donor-out and leave-one-study-out
(`crc_gut_scoring_primary_analysis.py`). Re-run from the already-scored,
already-verified per-cell data (no re-scoring needed — only the
robustness *definition* changed, not any correlation value). This is the
only pair whose flag changed; **9 of 10 pairs' flags are identical**
between the Pearson-only and the corrected Pearson-and-Spearman
criterion — the seed-bug re-run's qualitative stability (see above) still
holds under the corrected definition too.

**Reported as real, honest numbers — not smoothed over**:

- **All 10 correlations are weak in absolute magnitude** (|r| ≤ 0.19).
  Whatever revCSC↔D/F/P relationship exists at the single-cell level
  across the full, heterogeneous 665,473-cell atlas is not a strong,
  dominant signal — this compute does not find a large "Oncofetal = pure
  F-program" or "= pure P-program" effect at this scale/resolution.
- **D and P comparisons are the weakest and are the ones that ARE
  robust** to single-donor/study removal (|r| ≤ 0.030 throughout, no
  sign flip on leave-one-out). This is consistent with the negligible
  gene-overlap finding from PR #25's audit (D/P axes share zero genes
  with revCSC) — a weak, stable, near-null association is the expected,
  clean result here.
- **F comparisons are directionally consistent (all positive) and
  somewhat larger** (r = 0.005 to 0.19), strongest for
  `F_SI-specific`/`F_Gut-specific` with the extended revCSC variant
  (r = 0.17-0.19) and weakest for `F_Colon-specific` (r = 0.005-0.05)
  despite `F_Colon-specific` being the *primary* regional F axis per
  Step 4a's locked hierarchy — a real, not obviously expected, pattern.
- **5 of the 6 F-comparison pairs are NOT robust** to
  leave-one-donor-or-study-out (Pearson-and-Spearman criterion, see round
  2 correction above): all 3 pairings using the primary (27-gene) revCSC
  (`F_Colon-specific`, `F_SI-specific`, `F_Gut-specific`) **plus** both
  `extended28`↔`F_Colon-specific` and `extended28`↔`F_Gut-specific`. Only
  `extended28`'s pairing with `F_SI-specific` is robust among the F
  comparisons; all 4 D/P pairs (both revCSC variants) are robust. For
  `revCSC_primary27_minus_CLU` ↔ `F_Colon-specific` specifically,
  excluding the single study `Terekhanova_2023_Nature` shifts the pooled
  Pearson r by -0.054 — the largest single-study effect of any pair,
  more than 11× the pair's own pooled r (0.005). This means the
  already-weak `F_Colon-specific` correlation is itself substantially
  driven by one study's contribution, not a consistent cross-cohort
  signal. This instability is the finding for these 5 pairs, not
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

## Review history

- **Round 1 (REQUEST_CHANGES)**: two real issues found. (1) Wording bug —
  the results prose miscounted the F-comparison denominator (said "5",
  actually 6) and misattributed a non-robust pairing to the wrong revCSC
  variant; corrected, verified directly against the raw overview table.
  (2) Seed-determinism bug (real blocker) — per-panel RNG seeds used
  Python's process-randomized `hash()` instead of a deterministic hash,
  so the nominal fixed seed wasn't actually reproducible; fixed, and
  both the convergence check and full run were discarded and re-run from
  scratch rather than patched. This document reflects the reproducible
  re-run throughout.
- **Round 2 (fresh full review after round-1 fixes)**: reviewer re-fetched
  the PR at head and, while independently re-deriving numbers from the
  committed TSVs, caught two more real issues. (1) A prose/data
  mismatch — the results text cited a -0.064 single-study Pearson-r
  shift, but the committed overview table said -0.054; traced to a
  column-indexing mistake made while manually re-deriving the number for
  prose (read `pooled_spearman_rho_excl` instead of
  `delta_pearson_r_vs_full`) — the committed table's -0.054 was correct
  all along, only the prose was wrong; fixed, along with the resulting
  "13×" ratio claim (corrected to ~11×). (2) A real methodological gap —
  `robust_to_no_single_donor_or_study` was computed from Pearson sign
  stability only; direct verification against the per-donor/per-study
  leave-one-out files found `revCSC_extended28_minus_CLU_ASS1` ↔
  `F_Gut-specific` is Pearson-stable but Spearman-sign-flips when
  `Terekhanova_2023_Nature` is excluded. Fixed the robustness definition
  to require both metrics; re-ran the primary-analysis step only (the
  underlying per-cell scores were untouched, so no re-scoring needed) —
  exactly 1 of 10 pairs' flag changed (that one, True→False), confirming
  the fix was narrow and correctly scoped.

Submitting for compute review before merge, same discipline as every
prior step.
