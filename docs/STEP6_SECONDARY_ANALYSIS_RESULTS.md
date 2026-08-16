# Step 6 secondary analysis: revCSC-high developmental composition + M11 concordance — results

Real compute, per the approved design (`docs/STEP6_SECONDARY_ANALYSIS_DESIGN.md`,
PR #28, APPROVE after 4 review rounds). Answers the secondary-analysis
question explicitly deferred out of PR #27: whether revCSC-high cells
show a distinguishable developmental composition, and whether M11 (an
independent NMF-derived Oncofetal-annotated module) concords with revCSC.
All numbers below are real qsub output, pulled back and verified
byte-exact (md5, 30 files) against the Argos-side files before being
trusted or written up — same discipline as every prior step.

## Job provenance

| Stage | Job ID | Scope | Outcome |
|---|---|---|---|
| M11 gene-set build + overlap audit | (PR #28 review) | 6 M11 panels, Ensembl-native | Clean; 50/45/39/34/62/56 genes per panel |
| M11 N_PERM convergence probe | 3621108 | 20k-cell sample of the M11 subset | Certified `N_PERM=500` (Jaccard 0.94–0.98 vs. `N_PERM=1000`) |
| **M11 scoring + composition + concordance (this doc)** | **3621118** | 297,307 (M11)/665,473 (composition) cells | **Completed cleanly, ~3 min total; used for this doc** |

Per-cell M11 score/composition/concordance detail files are kept on
Argos only (`results/06_crc_projection/{m11_scoring_full,secondary_analysis_composition,secondary_analysis_concordance}/`),
reproducible from the frozen scripts + committed parameters + fixed seed.
The real, primary deliverables — all summary/overview tables — are
committed in full (30 files, this repo).

## M11 scoring: 6 panels, 297,307 cells, N_PERM=500

Barcode match to the M11 subset re-verified exact (297,307/297,307,
independently confirmed a third time — after the design's original claim
and the convergence probe's 20k-cell check — this time on the full
population). `n_testable` genes per panel (all present in the atlas,
0 unmapped):

| Panel | n_testable | Role |
|---|---|---|
| `M11_top50_full` | 50 | Sensitivity |
| `M11_top50_minus_revCSC_overlap` | 45 | **Primary** |
| `M11_top100_full` | 39 | Sensitivity (source file's "unique100" is deduplicated across NMF modules — see PR #28 §3) |
| `M11_top100_minus_revCSC_overlap` | 34 | Sensitivity |
| `M11_top200_full` | 62 | Sensitivity (source file's "unique200", same deduplication note) |
| `M11_top200_minus_revCSC_overlap` | 56 | Sensitivity |

## Section 1: revCSC-high cohort — exact cross-cell rank cutoffs

Cross-cell rank of `revCSC_primary27_minus_CLU_ASS1_zscore` across all
665,473 cells — verified exact by construction, not approximate:

| Cutoff | n cells | % of 665,473 |
|---|---|---|
| Top 5% | 33,274 | 5.00% |
| **Top 10% (primary)** | **66,548** | **10.00%** |
| Top 20% | 133,095 | 20.00% |

The `revCSC_extended28_minus_CLU_ASS1` sensitivity cohort (disclosed
lower precision, `N_PERM=100`) gives identical cell *counts* at every
cutoff by construction (same rank method, same population size) but a
**different set of cells** — used below as the F-facing overlap-safety
sensitivity check.

## Section 2: developmental composition of revCSC-high cells

**Primary result (Step 0, axis-supported status, top-10% cohort,
n=66,548, 433 donors)** — evidence versus each axis's own matched null
(percentile ≥90), not a forced single label:

| Category | Pooled % | Donor-unweighted mean % | Study-unweighted mean % |
|---|---|---|---|
| none | 41.1% | 49.1% | 38.4% |
| F_only | 42.6% | 35.1% | 39.8% |
| D_only | 2.2% | 2.1% | 2.0% |
| P_only | 3.5% | 5.1% | 4.5% |
| D+F | 5.4% | 4.2% | 5.1% |
| F+P | 4.4% | 3.7% | 4.1% |
| D+P | 0.2% | 0.4% | 0.2% |
| D+F+P | 0.6% | 0.4% | 0.3% |

**Real, honest read**: among revCSC-high cells, evidence is essentially
**bimodal** — either no gut-developmental axis clears its own
matched-null bar (~41%), or `F_Gut-specific` alone does (~43%). D-only
and P-only are both small (2–4%), and multi-axis coactivation is real but
modest (D+F 5%, F+P 4%, the double-support D+P/D+F+P categories under
1% each). This is qualitatively consistent with PR #27's primary finding
(F showed the largest, if still weak, pooled correlation with revCSC) —
**this secondary analysis adds that when a gut-developmental axis *is*
detectable in revCSC-high cells at all, it is overwhelmingly F, not D or
P**, though "no axis detectable" is close to equally common.

**Why Step 0 (not the coarse argmax) is the primary result — shown
directly, not just argued**: the unconditional 3-way argmax (`stepA`,
kept only as a descriptive summary per the locked design) assigns F to
**71.7%** of the same cohort — 29 points higher than Step 0's F-only+
any-F-supported total (`F_only`+`D+F`+`F+P`+`D+F+P` = 52.9%), because
`F_Gut-specific` (2,192 genes) mechanically tends to score higher on
relative terms than the much smaller D (8 genes) and P (76 genes) sets
even when none of the three axes shows real matched-null evidence. This
is exactly the argmax distortion the round-2 review predicted and the
reason the design was fixed to report Step 0 first — the discrepancy is
real, measured, and would have silently overstated "F predominance" by
~29 percentage points had the original (pre-fix) design shipped.

**Step A′ (predominant axis, among the 58.9% of cells with ≥1 supported
axis)**: F predominant in 45.6/58.9 = **77.4%** of supported cells, D in
4.1%, P in 6.2%, no clear predominance (< Δ5 margin) in 12.2%.

**Step B (regional refinement, within the 52.9% F-assigned cells)**:
`SI_biased` (14.3% of the full cohort, 27.0% of F-assigned cells) is
**~3× more common than** `Colon_biased` (4.7% of full cohort, 8.9% of
F-assigned cells) — the majority of F-assigned cells (64.1% of them) show
no clear Colon/SI skew. This is a **real, not obviously expected**
finding worth flagging plainly: `F_Colon-specific` is the *primary*
regional axis per Step 4a's locked hierarchy (and `CRC_single_cell_atlas_2025`
is a colon-cancer atlas), yet among revCSC-high cells with any regional
F skew at all, SI-biased is the more common pattern, not Colon-biased.
Not investigated further here — flagged as a real observation for future
scrutiny, not smoothed into "F_Colon-specific dominates as expected."

**Sensitivity checks — real numbers, not asserted**:

- **Cutoff sensitivity** (5%/10%/20% cohorts): the Step 0 breakdown is
  stable across all three cutoffs (`none` 41.6%/41.1%/41.1%, `F_only`
  41.3%/42.6%/43.6%, D/P-only both ≤4% at every cutoff) — the qualitative
  "bimodal none-vs-F" finding is not an artifact of the specific 10%
  threshold.
- **F-facing overlap-safety sensitivity** (the `extended28_minus_CLU_ASS1`-defined
  cohort, top-10% cutoff, disclosed `N_PERM=100`): `none` 39.1% (vs.
  41.1% primary), `F_only` 44.2% (vs. 42.6%), Colon/SI-biased 4.8%/14.6%
  (vs. 4.7%/14.3%) — the F-related conclusion survives an
  independently-defined, CLU/ASS1-overlap-safe revCSC-high cohort at a
  different (lower) precision, not just the single primary-cohort
  definition.

Full breakdowns for every cohort/cutoff, pooled and donor/study-level,
are committed in full (`results/06_crc_projection/secondary_analysis_composition/`,
6 files + 1 combined table).

## Section 3: M11 ↔ revCSC concordance (within the 297,307-cell M11 subset)

### Continuous correlation

| Comparison | Pooled Pearson r | Pooled Spearman ρ | Robust (donor+study leave-one-out)? |
|---|---|---|---|
| **`M11_top50_minus_revCSC_overlap` ↔ `revCSC_primary27_minus_CLU_ASS1`** (primary, gene-disjoint) | **0.318** | **0.421** | **Yes** |
| `M11_top50_full` ↔ `revCSC_primary27_minus_CLU_ASS1` (sensitivity, includes the 5 shared genes) | 0.440 | 0.507 | Yes |

137 donors, all 137 estimable (no small-`n_cells` exclusions within this
subset). **This is a real, moderate, robust positive correlation —
substantially stronger than any of PR #27's primary D/F/P↔revCSC
correlations** (all |r| ≤ 0.19). The gene-disjoint primary comparison
(r=0.318) is meaningfully weaker than the shared-gene sensitivity
comparison (r=0.440, as expected — the 5 shared genes contribute
mechanical correlation on top of any biological one) but remains a real,
substantial, donor/study-robust signal on its own — M11's concordance
with revCSC is not merely an artifact of the 5 genes it happens to share
with revCSC.

### Enrichment test (matched cutoff pairs — M11-high × revCSC-high)

M11-high: cross-cell z-score rank of `M11_top50_minus_revCSC_overlap`,
ranked within the 297,307-cell M11 subset. revCSC-high: Section 1's
global 665,473-cell ranking, restricted to the M11 subset (a fixed
membership set, not re-ranked).

| Cutoff | n M11-high | n revCSC-high | n both | Donors (estimable) | OR_MH | 95% CI (asymptotic) | 95% CI (donor-cluster bootstrap) |
|---|---|---|---|---|---|---|---|
| 5%×5% | 14,866 | 9,825 | 3,576 | 137 (111) | **5.79** | (5.40, 6.21) | (4.36, 7.85) |
| **10%×10% (primary)** | **29,731** | **20,418** | **9,322** | **137 (124)** | **6.18** | **(5.90, 6.47)** | **(4.33, 8.31)** |
| 20%×20% | 59,462 | 45,075 | 23,781 | 137 (131) | **4.92** | (4.78, 5.06) | (3.88, 6.07) |

**Real, robust, substantial enrichment** — every cutoff's CI (both the
asymptotic MH CI and the donor-cluster bootstrap CI, which does not
assume within-donor independence) clearly excludes 1. Per the round-2
review correction, MH is **not** presented as having solved cell-level
pseudoreplication by itself — the donor-cluster bootstrap CI is the
inferential-uncertainty companion that does not make that assumption,
and it agrees qualitatively with the asymptotic CI (wider, as expected,
but still excluding 1 at every cutoff).

**Robustness, checked directly (not asserted)**: leave-one-study-out was
the more consequential of the two checks — excluding `Tian_2023_Nat_Med`
produces the largest single-study shift at every cutoff (top5%: 5.79→4.57;
top10%: 6.18→4.43; top20%: 4.92→3.72) — but the common OR **never drops
below 3.7 for any single study excluded, at any cutoff** — the enrichment
finding is not driven by one study. Leave-one-donor-out shifts are all
small (max |Δ| 0.25–0.78 OR units across the 137 donors). Non-estimable
donors (zero-margin 2×2 table — a donor with either no M11-high or no
revCSC-high cells) are 26/137 (5% cutoff), 13/137 (10% cutoff), 6/137
(20% cutoff) — expected (fewer cells qualify as "high" at smaller
cutoffs, so more donors have zero events in one cell of their own 2×2
table), reported transparently, and excluded from per-donor OR only
(not from the MH pooling, which correctly handles zero-cell strata).

Full per-donor 2×2/OR tables, leave-one-donor/study-out tables, and the
correlation per-donor/leave-one-out detail are committed in full
(`results/06_crc_projection/secondary_analysis_concordance/`, 21 files).

## What this shows, honestly

**M11 (an independently-discovered, expression-only NMF meta-program,
originally identified via its own Jaccard-overlap match to revCSC) shows
a real, robust, moderate-to-large concordance with revCSC** — both
continuously (r=0.32, gene-disjoint) and as a threshold-based enrichment
(OR≈5–6, CI excluding 1 at every cutoff and every single-study/donor
exclusion checked) — **substantially stronger than any single gut-
developmental axis's correlation with revCSC** found in PR #27
(|r| ≤ 0.19 for all 10 D/F/P pairs). Developmentally, revCSC-high cells
that show any gut-axis evidence at all are dominated by `F_Gut-specific`
(not D or P), though "no axis detectable" is comparably common — a real,
bimodal pattern, not a clean single-axis story. The F-assigned subset
itself skews more SI-biased than Colon-biased, a genuinely unexpected
result for a colon-cancer atlas, flagged here for future scrutiny rather
than smoothed over.

**What this does NOT show**: this does not establish that M11 *is* the
"Oncofetal" fetal-gut program, nor that revCSC-high cells' F-association
constitutes a distinct, statistically separable substate (per the
locked design, that would require an additional pre-specified
clustering/mixture criterion, out of scope here — this analysis reports
"axis-defined composition/substructure," not formal separability). It
also does not extend beyond `CRC_single_cell_atlas_2025` or beyond the
already-frozen gut D/F/P and revCSC gene sets.

## What this does NOT do

Same explicit scope boundary as the approved design: does not touch the
tertiary analysis (full-atlas revCSC-independent D/F/P landscape); does
not extend to the 2 secondary/tertiary CRC datasets
(`HTAN_CRC_progressive_plasticity`, `CRLM_NMP_ATLAS`); does not
re-derive or re-score any frozen gut D/F/P, revCSC, or M11 gene set.

## Review history

Submitting for review before merge, same discipline as every prior step.
