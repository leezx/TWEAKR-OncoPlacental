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
66,548/665,473), genuinely pre-registered and outcome-independent (the
ranking procedure, not a specific score value, is fixed before compute).
Z-score (not percentile) is used specifically for this cross-cell ranking
step because it is continuous across its full range (median percentile is
capped near a ceiling many cells share, given only 100-500 null draws;
z-score has no such ceiling) — the percentile remains the reported,
primary per-cell common-scale value for every other purpose in this
project, this is the one place a different (still null-calibrated) metric
is used, for a stated, checkable reason.

**Two pre-registered sensitivity cutoffs, same construction**: top 5% and
top quintile by the same cross-cell z-score rank — reported alongside the
primary top-decile result, never substituted for it.

**`revCSC_extended28_full`'s equivalent (rank of
`revCSC_extended28_full_zscore`, no overlap-exclusion needed for the D/P
axes — see §4) is also computed and reported alongside**, not treated as a
separate primary threshold.

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

**Fixed with a two-step assignment, matching PR #25's actual axis
hierarchy**:
- **Step A (coarse, 3-way, matches the locked global-axis structure)**: per
  revCSC-high cell, the axis with the highest percentile among
  `{D_Gut-shared, F_Gut-specific, P_Gut-specific}` — one global axis each,
  no double-counting. This is the primary composition result, directly
  answering "predominantly D-shared, F (any lineage), or P-specific."
- **Step B (regional refinement, only within Step A's F-assigned cells)**:
  for cells assigned to F in Step A, which of `F_Colon-specific` /
  `F_SI-specific` has the higher percentile — reported as "of the
  F-assigned cells, X% skew Colon-specific, Y% skew SI-specific, Z% show no
  clear regional skew (percentiles within Δ5 of each other)," explicitly
  framed as a regional breakdown *within* the F category, not a 4th
  independent substate.

**Donor/study-aware aggregation (unchanged from the first draft, still
required)**: both steps computed pooled and as an unweighted mean across
donors (composite `donor_key`) and studies.

## 3. M11 scoring (new gene set, only new compute this analysis needs)

Unchanged from the first draft: 297,307-cell M11 subset identified by exact
barcode match (already verified in the approved design); M11 gene list from
`deliver.mm_top_genes.csv`'s `M11` column (top50, primary) +
`unique_top_versions/...unique100/200.csv` (secondary/sensitivity). Scoring
method: `score_genes_fast`, `N_BINS=20`, fixed seed `20260815`.

**Round 1 correction**: `N_PERM` raised to **500** (not 100) for this M11
scoring specifically — the membership-stability check below (§4) found
~90% Jaccard cohort agreement between `N_PERM=100` and `500` for a
percentile-threshold cohort on the 20,000-cell convergence-check subset,
meaning a hard-threshold (M11-high) use case is more sensitive to draw
count than the continuous-correlation use case `N_PERM=100` was originally
validated for. At 297,307 cells (vs. 665,473 for the primary compute),
`N_PERM=500` for one gene set is a small, affordable cost (well under the
full primary compute's total budget), not worth trading precision for.

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
*and* with M11 (none of `CLU`/`ASS1` are among the 5 M11-overlap genes
either) — a single, uniformly-safe cohort definition rather than tracking
several per-comparison variants. `revCSC_primary27_full` is retained only
as a stated sensitivity check on the threshold/cohort definition itself,
not used for any F-axis or M11 comparison. This variant is also already
scored at `N_PERM=500` in PR #27's full run (it was one of the 4 gated
panels) — the higher-precision option was already available at zero extra
compute cost.

## 5. M11 concordance test (revised cohort/variant references only; method unchanged)

Within the 297,307-cell subset only:

1. **Continuous correlation**: `M11_minus_revCSC_overlap` percentile (full
   M11 as sensitivity) vs. `revCSC_primary27_minus_CLU_ASS1` percentile
   (already computed at `N_PERM=500`, no re-scoring) — Pearson + Spearman,
   same within-study-then-pooled and donor/study-aware leave-one-out
   validation as the primary analysis (reusing
   `crc_gut_scoring_primary_analysis.py`'s machinery directly).
2. **Enrichment test**: using the same 3 revCSC-high thresholds (§1,
   cross-cell z-score rank of `revCSC_primary27_minus_CLU_ASS1`) restricted
   to this subset, test whether revCSC-high cells are enriched for M11-high
   status — M11-high defined by the *same* cross-cell z-score-rank
   construction as §1 (top decile of `M11_minus_revCSC_overlap`'s
   null-calibrated z-score, ranked across the 297,307-cell subset — not a
   fixed percentile cutoff, for the identical reason §1 required this fix)
   — via a 2×2 odds ratio, computed within each donor and pooled via
   Mantel-Haenszel (donor-stratified), not a single pooled contingency
   table.

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

Submitting for review before any qsub job runs, same discipline as every
prior step.
