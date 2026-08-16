# Step 6 secondary analysis: revCSC-high developmental composition + M11 concordance — design

Implements the **secondary analysis** from the already-approved Step 6 design
(`docs/STEP6_CRC_PROJECTION_DESIGN.md`, "Revised analysis structure", item 2),
explicitly deferred out of scope in PR #27 (primary analysis, merged). Reuses
the primary compute's already-scored, already-verified per-cell percentiles
(`crc_gut_scoring_all_panels.parquet`, job 3621066) wherever possible —
this analysis does **not** re-score revCSC or gut D/F/P; it only scores the
one new gene set (M11) and analyzes already-computed data.

**Round 1 review caught 5 real, structural problems in the first draft**,
all found by the reviewer independently re-deriving numbers from the
committed data rather than accepting the prose, and all confirmed directly
before fixing (not just taken on the reviewer's word):

## 1. revCSC-high threshold — completely redesigned (real bug, not wording)

**The first draft's "top decile = percentile ≥90" claim was empirically
wrong, not just imprecisely worded.** Checked directly on the full
665,473-cell cohort: the null-calibrated empirical percentile's own
distribution is heavily right-skewed — **median percentile is already
~90.2**, so "percentile ≥90" selects **~51% of all cells**, not 10%. This
makes sense once stated plainly: the null-calibrated percentile answers "is
this cell's revCSC expression higher than its own matched-null expectation,"
not "how does this cell rank against other cells" — most CRC malignant
cells apparently score above their own null expectation for this gene set,
a real, honest, and separately interesting observation in its own right
(reported explicitly, not just fixed and dropped).

**Fixed with a genuinely different construction: cross-cell rank of the
null-calibrated z-score**, not a fixed cutoff on the null-calibrated
percentile. `revCSC_primary27_minus_CLU_ASS1_zscore` (see §4 for why this
variant, not `_full`, is the cohort-defining one) is ranked across all
665,473 cells (`pandas.rank(pct=True)`, continuous, no ties from the
percentile's own ~1,000-value discretization ceiling) — **top decile by
this construction is exactly 10.0% of cells by definition** (verified:
66,548/665,473), **pre-specified before this secondary analysis and
outcome-independent** (the ranking procedure, not a specific score value,
is fixed before compute; not literal preregistration in the strict sense,
since PR #27's per-cell scores already existed and were inspected before
this construction was chosen — round 2 review correctly flagged the
"genuinely pre-registered" wording as overclaiming). Z-score (not
percentile) is used specifically for this cross-cell ranking step because
it is continuous across its full range (median percentile is capped near a
ceiling many cells share, given only 100-500 null draws; z-score has no
such ceiling) — the percentile remains the reported, primary per-cell
common-scale value for every other purpose in this project, this is the
one place a different (still null-calibrated) metric is used, for a
stated, checkable reason.

**Two pre-specified sensitivity cutoffs, same construction**: top 5% and
top quintile by the same cross-cell z-score rank — reported alongside the
primary top-decile result, never substituted for it.

**Round 2 fix — `revCSC_extended28` sensitivity cohort corrected for
overlap-safety**: the first round-2 draft still specified
`revCSC_extended28_full`'s equivalent (rank of `revCSC_extended28_full_zscore`)
as the extended-cohort sensitivity. This was internally inconsistent with
this design's own overlap-safe cohort rule (§4): `revCSC_extended28_full`
contains `CLU`/`ASS1`, so an F-axis composition/enrichment test run against
an `extended28_full`-defined cohort would reintroduce exactly the
mechanical overlap the primary cohort was fixed to exclude. **Fixed**: the
extended-cohort sensitivity uses `revCSC_extended28_minus_CLU_ASS1`
throughout (zero overlap with every gut F axis, same as the primary
cohort's `_minus_CLU_ASS1` variant) — `revCSC_extended28_full` is no longer
used for any F-facing analysis in this design. **Disclosed limitation**:
PR #27's own output only scored `revCSC_extended28_minus_CLU_ASS1` at
`N_PERM=100` (unlike the primary `revCSC_primary27_minus_CLU_ASS1`, which
happens to be at `N_PERM=500`). Given the real z-rank Jaccard(100,500)
finding below (§ convergence check) — 0.77–0.88, not the ~90% previously
claimed — this extended sensitivity cohort's hard-threshold membership is
disclosed as lower-precision than the primary cohort's; it is reported and
interpreted strictly as a sensitivity check on the *threshold construction
itself* (does the qualitative composition/enrichment conclusion survive an
independently-defined revCSC variant), not as an equal-precision
replication. Re-scoring `extended28_minus_CLU_ASS1` at `N_PERM=500` on the
full atlas was considered and rejected as disproportionate new compute for
a secondary sensitivity check on an already-answered question (Q1's
primary result already rests on the `N_PERM=500`-gated panels).

## 2. Developmental composition of revCSC-high cells — axis structure fixed

**The first draft treated `F_Colon-specific` and `F_SI-specific` as
categories parallel to and mutually exclusive with `D_Gut-shared`/
`P_Gut-specific` in a 4-way argmax.** This is structurally wrong against
PR #25's own locked contract: there is exactly **one** global D axis, **one**
global F axis (`F_Gut-specific`, coarse), **one** global P axis — only F is
region-resolved, and `F_Colon-specific`/`F_SI-specific` are *regional
interpretations of that one F axis*, not independent parallel categories.
The two regional F sets also share roughly 711 genes between them (an
explicit reviewer-provided estimate, not yet independently re-verified this
round but directionally consistent with both being built from the same
`F_Colon-developmental ∪ F_SI-developmental`-derived construction, PR #25) —
running a 4-way argmax across two highly-overlapping-in-content sets plus D
and P would mechanically inflate the appearance of "multiple separable
substates" by splitting what is substantially one shared program into two
argmax-eligible categories.

**Fixed (round 1) with a two-step assignment, matching PR #25's actual axis
hierarchy** — this corrected which axes compete, but round 2 review
correctly identified a second, independent problem still present: an
**unconditional** 3-way argmax forces every cell into exactly one of
D/F/P even when all three axes show weak evidence, and forces a single
label onto a cell that is genuinely high on two axes at once. Argmax alone
cannot distinguish a real discrete substate from coactivation or from the
absence of any clear developmental program — it was never valid evidence
for "multiple separable substates," only for "which axis is relatively
larger, if any of them mean anything for this cell."

**Round 2 fix — axis-supported status computed first, argmax demoted to a
descriptive summary**:
- **Step 0 (axis-supported status, the actual substate-structure result)**:
  for each revCSC-high cell, per axis in `{D_Gut-shared, F_Gut-specific,
  P_Gut-specific}`, "supported" = null-calibrated percentile ≥90 (the same
  pre-specified null-calibrated threshold used for every other hard cutoff
  in this project — evidence versus that cell's own matched null, not a
  cross-cell "top decile" claim). Each cell is then classified into exactly
  one of 8 mutually exclusive, jointly exhaustive categories: **none / D
  only / F only / P only / D+F / D+P / F+P / D+F+P**. This is reported as
  the primary composition result and directly exposes coactivation
  (multi-axis-supported cells) and absence of a clear program
  (none-supported cells) instead of hiding both inside a forced single
  label.
- **Step A (coarse 3-way argmax, retained as a descriptive summary only)**:
  per revCSC-high cell, the axis with the highest percentile among
  `{D_Gut-shared, F_Gut-specific, P_Gut-specific}`. Reported alongside Step
  0, restricted to (or at minimum cross-tabulated against) axis-supported
  status — **explicitly not used, on its own, as evidence that revCSC-high
  cells form "multiple separable substates."** A formal separability claim
  would require an additional pre-specified clustering/mixture criterion,
  out of scope here; absent that, this design's claim is limited to
  "axis-defined composition/substructure."
- **Step A′ (predominant axis among supported cells only)**: restricted to
  cells with exactly one or more supported axes (i.e. excluding
  "none"-status cells, for whom "predominant axis" is not a meaningful
  question), the axis with the highest percentile — but only reported as
  "predominant" if it leads the second-highest axis by a pre-specified
  margin of **Δ5 percentile points** (same margin concept already used for
  the Colon/SI regional call below); cells not meeting the margin are
  reported as "no clear predominance among supported axes" rather than
  forced to a single label.
- **Step B (regional refinement, only within Step 0/A's F-supported or
  F-predominant cells)**: for those cells, which of `F_Colon-specific` /
  `F_SI-specific` has the higher percentile — reported as "of the
  F-assigned cells, X% skew Colon-specific, Y% skew SI-specific, Z% show no
  clear regional skew (percentiles within Δ5 of each other)," explicitly
  framed as a regional breakdown *within* the F category, described
  descriptively ("Colon-biased"/"SI-biased"), not itself claimed as
  evidence of two separable substates (the two regional F panels still
  share ~711 genes).

**Donor/study-aware aggregation (unchanged from the first draft, still
required)**: every step above computed pooled and as an unweighted mean
across donors (composite `donor_key`) and studies.

## 3. M11 scoring (new gene set, only new compute this analysis needs)

Unchanged from the first draft: 297,307-cell M11 subset identified by exact
barcode match (re-verified directly this round: 297,307/297,307, see the
convergence-probe result below); M11 gene list from `deliver.mm_top_genes.csv`'s
`M11` column (top50, primary) + `unique_top_versions/deliver.mm_top_genes.unique{100,200}.csv`
(secondary/sensitivity) — both already Ensembl-ID-native, no symbol mapping
needed (confirmed directly, §4). **Noted directly while building these
files** (`build_m11_gene_sets.py`): the `unique100`/`unique200` files are
*deduplicated-across-modules* versions, so their M11 columns contain fewer
than 100/200 genes after removing genes already claimed by
higher-priority modules — real counts are **39** (top100) and **62**
(top200), not 100/200; the panel names are kept as-is (`M11_top100`,
`M11_top200`) since that is the source files' own naming, but the true
gene counts are reported everywhere a count is shown. Scoring method:
`score_genes_fast`, `N_BINS=20`, fixed seed `20260815`.

**Round 1 correction**: `N_PERM` raised to **500** (not 100) for this M11
scoring specifically, since a hard-threshold (M11-high) use case is more
sensitive to draw count than the continuous-correlation use case
`N_PERM=100` was originally validated for. At 297,307 cells (vs. 665,473
for the primary compute), `N_PERM=500` for one gene set is a small,
affordable cost (well under the full primary compute's total budget), not
worth trading precision for.

**Round 2 correction (real bug, not wording)**: the round-1 draft's "~90%
Jaccard between `N_PERM=100`/`500`" claim was computed using the *old*
percentile-threshold cohort construction (percentile≥90), not the
cross-cell z-score-rank construction this design actually uses after the
round-1 threshold fix — those are different cohort definitions, and the
number did not actually validate the method in use. **Recomputed directly**
on the existing 20,000-cell convergence-check parquet files
(`results/06_crc_projection/gut_scoring_convergence_check/scores_nperm{100,500}.parquet`,
already pulled back, no new compute), using the real construction — cross-cell
rank of `revCSC_primary27_minus_CLU_ASS1_zscore` — at top5%/10%/20%:

| cutoff | N_PERM=100 n | N_PERM=500 n | intersection | union | Jaccard |
|---|---|---|---|---|---|
| top 5% | 1,001 | 1,001 | 874 | 1,128 | **0.7748** |
| top 10% | 2,001 | 2,001 | 1,803 | 2,199 | **0.8199** |
| top 20% | 4,001 | 4,001 | 3,743 | 4,259 | **0.8788** |

**This replaces the incorrect ~90% figure** — the real agreement is lower,
which if anything strengthens (does not weaken) the case for `N_PERM=500`
over `100`. But as the reviewer correctly pointed out, 100-vs-500 agreement
only shows that 100 is *insufficient*; it does not show that 500 itself has
*converged*. A separate 500-vs-1000 probe is required and was run — see the
dedicated convergence-probe result below.

**M11 N_PERM=500 vs. 1000 convergence probe (new compute, run before the
297,307-cell M11 job)**: `scripts/06_crc_projection/m11_nperm_convergence_probe.py`
scores the 45-gene `M11_top50_minus_revCSC_overlap` panel at `N_PERM=500`
and `N_PERM=1000` on a fixed 20,000-cell stratified-by-`study_id` sample of
the 297,307-cell M11 population itself (not the full-atlas convergence
subset — M11 is only defined on its own barcode-matched subset), then
reports the same z-rank cohort-membership Jaccard at top5/10/20% between
the two draw counts.

**Result (job 3621108, byte-exact pulled back and md5-verified, 9 files)**:
barcode match to the M11 subset confirmed exact (297,307/297,307).
`N_PERM=500` vs. `1000` z-rank membership agreement:

| cutoff | N_PERM=500 n | N_PERM=1000 n | intersection | union | Jaccard |
|---|---|---|---|---|---|
| top 5% | 1,001 | 1,001 | 972 | 1,030 | **0.9437** |
| top 10% | 2,001 | 2,001 | 1,967 | 2,035 | **0.9666** |
| top 20% | 4,001 | 4,001 | 3,964 | 4,038 | **0.9817** |

Continuous percentile Pearson r (500 vs. 1000) = **0.9997**. This is
substantially higher agreement than the revCSC 100-vs-500 z-rank numbers
above (0.77–0.88) — **`N_PERM=500` is certified as adequate for the
45-gene `M11_top50_minus_revCSC_overlap` panel** on the real 297,307-cell
run; the hard-threshold-precision claim this design relies on is now
directly validated for the actual construction in use, not inferred from a
different cohort definition.

## 4. M11 × revCSC gene overlap — real blocker found and fixed

**The first draft only audited M11 against gut D/F/P and asserted M11 could
be scored "directly" (no exclusion) against revCSC for the concordance
test, without ever checking M11 against revCSC itself.** This was a real
gap: M11 was *originally identified* via a Jaccard-overlap match to
revCSC (per the approved Step 6 design's own "Independent Oncofetal anchor"
section) — checked directly, M11's top50/top100 gene lists share **5
genes** with `revCSC_primary27_full` (and all its `_minus_CLU`/`_minus_ASS1`
variants, since none of the 5 are `CLU`/`ASS1`): `ANXA1`, `KRT18`, `SFN`,
`TMSB4X`, `TNFRSF12A` (top200 adds a 6th, `PMEPA1`) — nearly 1 in 5 of
revCSC's 27 primary genes. Using unmodified M11 to test "does M11 correlate
with / enrich alongside revCSC" would be answering a question partly
mechanical-by-construction (5 shared genes contributing to both scores),
not purely biological.

**Fixed**: `M11_minus_revCSC_overlap` (45 genes at top50, dropping the 5
shared genes) is the primary variant for **every M11↔revCSC comparison in
this analysis** (both the continuous correlation and the enrichment test,
§5); full (non-excluded) M11 is retained only as a sensitivity check,
reported alongside, never as the primary result. The M11×gut-D/F/P overlap
audit from the first draft (real, F-axis overlap found, zero D/P overlap)
is kept as a reported transparency check — this design's actual compute
never correlates M11 against D/F/P directly, so it does not gate an
exclusion contract here, but the numbers are kept in the doc since they are
real, already-computed, and relevant context for anyone extending this
analysis later.

**Round 2 fix — overlap audit committed as a script/artifact, not just
design prose** (per review: "I would also commit the overlap audit as a
script/result artifact"): `scripts/06_crc_projection/build_m11_gene_sets.py`
builds the M11 top50/100/200 Ensembl-ID gene sets directly from their
already-Ensembl-native NMF source files (`deliver.mm_top_genes.csv` /
`unique_top_versions/deliver.mm_top_genes.unique{100,200}.csv` — confirmed
directly, no symbol mapping needed, unlike the D/F/P panels), both the full
and `_minus_revCSC_overlap` variant of each, and writes the M11×revCSC
overlap audit to `results/06_crc_projection/m11_revcsc_overlap_audit/m11_revcsc_overlap_audit.tsv`
as a committed, re-derivable artifact.

**Also real, also fixed**: the first draft defined "revCSC-high" using
unmodified `revCSC_primary27_full`, then asked whether that cohort enriches
for F-axis percentiles — but `F_Colon-specific`/`F_SI-specific`/
`F_Gut-specific` all share `CLU`/`ASS1` with `revCSC_primary27_full` (the
exact reason PR #27's primary analysis used the `_minus_CLU`/`_minus_ASS1`/
`_minus_CLU_ASS1` variants for every F comparison). Selecting cells by
`_full` and then testing F-axis enrichment would reintroduce precisely the
mechanical overlap PR #27 was built to exclude. **Fixed**: the
cohort-defining revCSC score throughout this entire analysis (§1's
threshold, §2's composition, §5's enrichment test) is
`revCSC_primary27_minus_CLU_ASS1` — zero overlap with every gut D/F/P axis
— a single, uniformly-safe cohort definition rather than tracking several
per-comparison variants. **Round 2 factual correction**: the design's first
round-2 draft additionally claimed this variant has zero overlap "and with
M11" — **this is wrong and has been corrected.** `revCSC_primary27_minus_CLU_ASS1`
only removes `CLU`/`ASS1`; it still contains all 5 of the M11-overlap genes
(`ANXA1`, `KRT18`, `SFN`, `TMSB4X`, `TNFRSF12A` — none of the 5 is `CLU` or
`ASS1`). The M11↔revCSC comparison (§5) is gene-disjoint for a different
reason, entirely on the M11 side: `M11_minus_revCSC_overlap` removes those
same 5 genes from M11's gene list (§4 below). The actual primary
comparison (`M11_minus_revCSC_overlap` vs. `revCSC_primary27_minus_CLU_ASS1`)
is therefore still correctly gene-disjoint — only the stated rationale for
*why* was wrong, which could have misled a future implementer into
thinking the revCSC side alone was already M11-safe. `revCSC_primary27_full`
is retained only
as a stated sensitivity check on the threshold/cohort definition itself,
not used for any F-axis or M11 comparison. This variant is also already
scored at `N_PERM=500` in PR #27's full run (it was one of the 4 gated
panels) — the higher-precision option was already available at zero extra
compute cost.

## 5. M11 concordance test (round 2: enrichment cutoffs matched + MH robustness discipline added)

Within the 297,307-cell subset only:

1. **Continuous correlation**: `M11_minus_revCSC_overlap` percentile (full
   M11 as sensitivity) vs. `revCSC_primary27_minus_CLU_ASS1` percentile
   (already computed at `N_PERM=500`, no re-scoring) — Pearson + Spearman,
   same within-study-then-pooled and donor/study-aware leave-one-out
   validation as the primary analysis (reusing
   `crc_gut_scoring_primary_analysis.py`'s machinery directly).
2. **Enrichment test — round 2 fix: matched cutoff pairs, not a
   fixed-M11/varying-revCSC design.** The round-1 draft varied the revCSC
   cutoff across its 3 pre-specified thresholds while holding M11-high
   fixed at top decile — a different (weaker) robustness question than
   intended. **Fixed**: revCSC-high and M11-high are defined by the *same*
   cross-cell z-score-rank construction (§1's fix, applied identically on
   each side — not a fixed percentile cutoff), evaluated as **matched
   pairs**: **10%×10% primary**, **5%×5%** and **20%×20%** sensitivities —
   both sides move together, never one fixed while the other varies.
   M11-high is defined on `M11_minus_revCSC_overlap`'s null-calibrated
   z-score, ranked across the 297,307-cell subset.
3. **Effect estimate**: a 2×2 odds ratio computed within each donor,
   pooled via Mantel-Haenszel (donor-stratified common OR).
4. **Round 2 fix — MH robustness discipline, matching PR #27's standard
   (hard blocker, unresolved after round 1)**: MH stratification controls
   donor-specific baseline rate differences (avoiding a pooled Simpson's-
   paradox problem), but it does **not** turn the many cells within one
   donor into independent observations, and does not by itself prevent a
   single large donor from dominating the common OR — both the `statsmodels`
   implementation and the standard CMH definition (conditional association
   within strata) are explicit about this, correctly flagged by the
   reviewer. **The MH common OR is kept as the primary effect-size
   estimate, but is no longer described as having solved donor-level
   pseudoreplication**, and is now required to carry the same robustness
   discipline PR #27's continuous correlations already use:
   - the full **per-donor 2×2 table and per-donor OR**, reported (not just
     the pooled MH statistic), with donors whose 2×2 table is
     non-estimable (a zero cell) explicitly listed and excluded from the MH
     pooling with the exclusion count stated, not silently dropped;
   - **leave-one-donor-out** and **leave-one-study-out common OR**,
     mirroring `crc_gut_scoring_primary_analysis.py`'s existing
     leave-one-out machinery, to check no single donor/study drives the
     pooled result;
   - **no conclusion based solely on the cell-level MH p-value/CI** — a
     donor-cluster bootstrap CI (resampling donors, not cells, with
     replacement, recomputing the common OR each resample) is added as the
     donor-cluster-aware inferential-uncertainty companion to the point
     estimate, replacing reliance on the standard MH CI alone.

## Scope boundary (explicit, unchanged)

- Still `CRC_single_cell_atlas_2025` only — no other dataset.
- Does not re-score or re-derive revCSC or any gut D/F/P set.
- Does not touch the tertiary analysis (full-atlas revCSC-independent D/F/P
  landscape) — separate, later work.
- M11's top50/100/200 gene lists and the 297,307-cell subset are used
  exactly as already independently derived — no re-clustering, no
  re-running the NMF pipeline.

## Review history

- **Round 1 (REQUEST_CHANGES, 5 real issues, all confirmed directly before
  fixing)**: (1) "top decile = percentile≥90" was empirically wrong
  (actually ~51% of cells) — replaced with cross-cell z-score-rank
  construction, verified exact 10.0% at the primary cutoff. (2) The
  predominant-axis assignment double-counted `F_Colon-specific`/
  `F_SI-specific` as parallel to D/P despite PR #25's one-global-F-axis
  contract and their ~711-gene mutual overlap — fixed to a two-step
  coarse-then-regional assignment. (3) M11's real 5-gene overlap with
  revCSC itself (the genes that originally identified M11 as revCSC-like)
  was never checked — fixed with an `M11_minus_revCSC_overlap` variant as
  primary for all M11↔revCSC comparisons. (4) The revCSC-high cohort used
  the unmodified `_full` variant, reintroducing the exact CLU/ASS1
  mechanical overlap with F axes that PR #27 excluded — fixed to use
  `revCSC_primary27_minus_CLU_ASS1` uniformly as the cohort-defining score.
  (5) Hard-threshold cohort membership stability at `N_PERM=100` was never
  checked (PR #27's convergence validation was for continuous correlation
  only) — checked directly on the already-computed convergence-check data:
  ~90% Jaccard agreement between `N_PERM=100`/`500` cohorts at a
  percentile-threshold construction; M11 scoring raised to `N_PERM=500`
  accordingly (affordable at 297,307 vs. 665,473 cells), and the
  cohort-defining revCSC variant already happens to be at `N_PERM=500` in
  PR #27's own output.

- **Round 2 (REQUEST_CHANGES, "substantially closer" — 3 of 5 round-1 fixes
  confirmed structurally sound, 2 substantive blockers + 3 smaller
  inconsistencies found, all confirmed directly before fixing)**:
  (1) **Argmax alone cannot show "substates"** — round 1's two-step
  coarse-then-regional fix corrected *which* axes compete but left an
  unconditional 3-way argmax forcing every cell to one label regardless of
  evidence strength, which cannot distinguish discrete substates from
  coactivation/no-program — fixed by computing axis-supported status
  (none/D/F/P/D+F/D+P/F+P/D+F+P via a pre-specified percentile≥90
  per-axis threshold) as the primary composition result, demoting argmax to
  a descriptive summary explicitly disclaimed as substate evidence, and
  adding a Δ5-margin "predominant axis among supported cells only" step.
  (2) **MH donor-stratification does not solve cell-level pseudoreplication**
  — still described that way after round 1 — fixed by keeping MH as the
  point estimate only, adding per-donor 2×2/OR reporting (non-estimable
  donors exposed), leave-one-donor/study-out common OR, and a donor-cluster
  bootstrap CI, matching PR #27's existing robustness discipline.
  (3) **The reported ~90% N_PERM=100-vs-500 Jaccard used the wrong (old
  percentile-threshold) cohort construction** — recomputed directly on the
  existing convergence-check parquet using the actual z-rank construction:
  real numbers are 0.7748/0.8199/0.8788 at top5/10/20%, lower than
  previously claimed; a further new probe (500-vs-1000, on the actual
  45-gene M11 panel and its own 20k-cell population subset, job 3621108) was
  run since 100-vs-500 agreement alone cannot certify 500's own convergence
  — see §3 for the result. (4) **`revCSC_extended28_full` reintroduced
  CLU/ASS1** when used as the extended sensitivity cohort for an F-facing
  test, contradicting the primary cohort's own overlap-exclusion rule —
  fixed to `revCSC_extended28_minus_CLU_ASS1`, with its `N_PERM=100`
  precision limitation (vs. the primary cohort's 500) disclosed rather than
  masked. (5) **A factual error**: `revCSC_primary27_minus_CLU_ASS1` was
  incorrectly described as zero-overlap with M11 (it still contains all 5
  M11-overlap genes; the CLU/ASS1 exclusion is unrelated to the M11 overlap)
  — the actual comparison is still correctly gene-disjoint, but for the
  right reason: `M11_minus_revCSC_overlap` removes those 5 genes on the M11
  side — corrected. Two smaller terminology/design locks also applied:
  "genuinely pre-registered" reworded to "pre-specified before the
  secondary analysis" throughout (PR #27's scores already existed and were
  inspected before this construction was chosen, so it is not literal
  preregistration); and the enrichment test's revCSC/M11 cutoffs changed
  from fixed-M11×varying-revCSC to matched pairs (10%×10% primary,
  5%×5%/20%×20% sensitivities). The M11×revCSC overlap audit was also
  committed as a script/artifact (`build_m11_gene_sets.py` +
  `m11_revcsc_overlap_audit.tsv`) rather than living only in design prose.

- **Round 3 (REQUEST_CHANGES, "no remaining statistical/design blocker"
  after this fix — 1 real issue, confirmed directly before fixing)**: all 5
  round-2 fixes confirmed substantively closed. One remaining blocker in
  `build_m11_gene_sets.py`: the M11×revCSC overlap audit was **circular** —
  it hardcoded the expected 5/6 overlap symbols, converted only those to
  Ensembl IDs, and intersected M11 against that pre-specified set, so an
  unexpected 7th overlap gene could never have been discovered and would
  have silently remained in `*_minus_revCSC_overlap`. **Confirmed directly**
  by reading the script — the reviewer's description was exactly right.
  **Fixed**: `found_overlap` is now derived as
  `genes_set & revcsc_ens_set` (the *full* 27-gene `revCSC_primary27`
  Ensembl set), independently, for every M11 cutoff; the previously-found
  5/6-gene list is retained only as an `assert` check against this
  independently-derived result (fails loudly if a future M11 gene-list
  revision introduces an unexpected overlap), never as the source of the
  exclusion. Re-ran on Argos: the independently-derived overlap is
  identical to the previously (circularly) reported one — 5 genes at
  top50/top100, 6 at top200 — confirming the earlier numbers were correct,
  just for the wrong (circular) reason; all output files byte-identical
  (md5-verified) to the round-2 versions, only the script changed.

Submitting for round-4 review before the 297,307-cell M11 job runs, same
discipline as every prior step.
