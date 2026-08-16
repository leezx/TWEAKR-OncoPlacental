# Step 6 secondary analysis: revCSC-high developmental composition + M11 concordance — results

Real compute, per the approved design (`docs/STEP6_SECONDARY_ANALYSIS_DESIGN.md`,
PR #28, APPROVE after 4 review rounds). Answers the secondary-analysis
question explicitly deferred out of PR #27: whether revCSC-high cells
show a distinguishable developmental composition, and whether M11 (an
independent NMF-derived Oncofetal-annotated module) concords with revCSC.
All numbers below are real qsub output, pulled back and verified
byte-exact (md5) against the Argos-side files before being trusted or
written up — same discipline as every prior step.

**Round 1 review caught 2 real implementation deviations from the
approved design plus 4 real write-up errors** (arithmetic, an
unsupported causal claim, a mislabeled study, and an overstated
cross-population comparison), all independently verified against real
data before fixing — see "Review history."

## Job provenance

| Stage | Job ID | Scope | Outcome |
|---|---|---|---|
| M11 gene-set build + overlap audit | (PR #28 review) | 6 M11 panels, Ensembl-native | Clean; 50/45/39/34/62/56 genes per panel |
| M11 N_PERM convergence probe | 3621108 | 20k-cell sample of the M11 subset | Certified `N_PERM=500` (Jaccard 0.94–0.98 vs. `N_PERM=1000`) |
| M11 scoring + composition + concordance (round 1 draft) | 3621118 | 297,307 (M11)/665,473 (composition) cells | Completed cleanly, ~3 min; superseded by the rerun below (composition/enrichment logic fixed, M11 scores unchanged) |
| **Composition + concordance rerun (this doc)** | **3621126** | 665,473/297,307 cells, M11 scores reused from job 3621118 (not rescored — round-1 review confirmed no re-score needed) | **Completed cleanly, <1 min; used for this doc** |

Per-cell M11 score/composition/concordance detail files are kept on
Argos only (`results/06_crc_projection/{m11_scoring_full,secondary_analysis_composition,secondary_analysis_concordance}/`),
reproducible from the frozen scripts + committed parameters + fixed seed.
The real, primary deliverables — all summary/overview/cross-tab tables —
are committed in full to this repo (33 files pulled back and md5-verified
this round; 3 large per-cell parquet files excluded per the size-norm
precedent PR #27 set, kept Argos-only).

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

**Why Step 0 (not the coarse argmax) is the primary result — quantified
directly via a Step0×StepA cross-tab, not just argued**: the round-1
draft of this document compared argmax's 71.7% F share against Step 0's
F-supported total (52.9%) and reported the gap as "~29 points," but
71.7% − 52.9% = 18.8 points, not 29 — a real arithmetic error, caught in
round-1 review and fixed. The cross-tab (`composition_primary27_minus_CLU_ASS1_top10pct__step0_x_stepA_crosstab.tsv`)
gives the precise, correct breakdown instead of an arithmetic difference
of two marginals: of the cohort's 47,747 argmax-F cells (71.7%), **14,208
(21.4 percentage points of the cohort, ~30% of all argmax-F assignments)
come from cells with `step0_status = none`** — i.e. cells where *none* of
D/F/P clears its own matched-null bar, and F is simply the
least-unconvincing of three unconvincing scores. The remaining 33,539
argmax-F cells (50.4% of the cohort) do come from genuinely F-supported
categories (`F_only`+`D+F`+`F+P`+`D+F+P`). **What this does NOT establish**:
the round-1 draft additionally claimed this happens "because
`F_Gut-specific` has 2,192 genes" (mechanically outscoring the much
smaller D/P sets) — round-1 review correctly flagged this as an
unsupported causal claim; the compute null-calibrates each signature
against its own matched-null gene sets specifically to control for
gene-set size, so this result does not by itself demonstrate that size is
the mechanism. The defensible, demonstrated claim is narrower: an
**unconditional argmax forces a relative label even when no axis clears
its own matched-null evidence threshold** — the mechanism is not
separately tested here.

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
plus the Step0×StepA cross-tab per cohort, are committed in full
(`results/06_crc_projection/secondary_analysis_composition/`, 5 composition
tables + 5 cross-tab tables + 2 combined tables).

## Section 3: M11 ↔ revCSC concordance (within the 297,307-cell M11 subset)

### Continuous correlation

| Comparison | Pooled Pearson r | Pooled Spearman ρ | Robust (donor+study leave-one-out)? |
|---|---|---|---|
| **`M11_top50_minus_revCSC_overlap` ↔ `revCSC_primary27_minus_CLU_ASS1`** (primary, gene-disjoint) | **0.318** | **0.421** | **Yes** |
| `M11_top50_full` ↔ `revCSC_primary27_minus_CLU_ASS1` (sensitivity, includes the 5 shared genes) | 0.440 | 0.507 | Yes |

137 donors, all 137 estimable (no small-`n_cells` exclusions within this
subset). This is a real, moderate, robust positive correlation — the
gene-disjoint primary comparison (r=0.318) is meaningfully weaker than
the shared-gene sensitivity comparison (r=0.440, as expected — the 5
shared genes contribute mechanical correlation on top of any biological
one) but remains real and donor/study-robust on its own.

**Round-1 fix — same-population comparison, not a cross-population one**:
the round-1 draft called this "substantially stronger than any of PR #27's
primary D/F/P↔revCSC correlations (all |r| ≤ 0.19)" — correctly flagged
in review as comparing r=0.318 (computed within the 297,307-cell M11
subset) against PR #27's values (computed on the full 665,473-cell
atlas), a cross-population, not same-population, comparison. **Fixed by
directly computing the D/F/P↔revCSC correlations restricted to the exact
same 297,307-cell M11 subset** (no new scoring needed — the D/F/P
percentiles were already scored for all 665,473 cells in PR #27; this
just restricts them to the M11 subset and reuses the same correlation
machinery):

| Comparison (M11 subset only, n=297,307) | Pearson r | Spearman ρ |
|---|---|---|
| `revCSC_primary27_minus_CLU_ASS1` ↔ `D_Gut-shared` | -0.002 | -0.001 |
| `revCSC_primary27_minus_CLU_ASS1` ↔ `F_Gut-specific` | 0.012 | -0.182 |
| `revCSC_primary27_minus_CLU_ASS1` ↔ `F_Colon-specific` | -0.074 | -0.151 |
| `revCSC_primary27_minus_CLU_ASS1` ↔ `F_SI-specific` | 0.005 | -0.140 |
| `revCSC_primary27_minus_CLU_ASS1` ↔ `P_Gut-specific` | -0.057 | -0.031 |

**This is a genuine same-population comparison, and it makes the finding
stronger, not weaker**: within the exact same 297,307-cell M11 subset,
every D/F/P↔revCSC Pearson r is essentially null (|r| ≤ 0.074 — even
smaller in magnitude than PR #27's full-atlas values), while M11's
gene-disjoint correlation with revCSC in the identical population is
r=0.318. **M11 correlates with revCSC far more than any gut-developmental
axis does in the same 297,307-cell population** — a cleaner, more
defensible claim than the round-1 draft's cross-population one.

### Enrichment test (matched cutoff pairs — M11-high × revCSC-high)

M11-high: cross-cell z-score rank of `M11_top50_minus_revCSC_overlap`,
ranked within the 297,307-cell M11 subset. revCSC-high: Section 1's
global 665,473-cell ranking, restricted to the M11 subset (a fixed
membership set, not re-ranked).

**Round-1 fix — MH pooling now excludes non-estimable donors, per the
locked design**: the round-1 draft computed each donor's "estimable"
flag (zero-cell 2×2 table) for the per-donor reporting table but then
passed *all* donor strata (estimable or not) into the pooled MH
estimate, its CIs, and the leave-one-out/bootstrap sensitivities —
directly contradicting `docs/STEP6_SECONDARY_ANALYSIS_DESIGN.md`'s
explicit instruction that non-estimable donors be "excluded from the MH
pooling." Confirmed as a real deviation by re-reading both the design
text and the code. **Checked directly before fixing**: whether
unconditionally retaining zero-cell strata would actually be more
statistically standard is a genuinely separate question (it would — MH's
formula handles zero-cell strata without bias, and roughly half of the
"non-estimable" donors per cutoff have an informative, not fully
degenerate, margin), but a PR explicitly executing an approved design
must implement that design, not a different estimator — so this is fixed
to comply literally: MH pooling, its CIs, and leave-one-out/bootstrap
are now restricted to estimable donors only, with the exclusion count
reported. The shift this causes is small (never changes which side of 1
any CI falls on):

| Cutoff | n M11-high | n revCSC-high | n both | Donors (estimable) | OR_MH | 95% CI (asymptotic) | 95% CI (donor-cluster bootstrap) |
|---|---|---|---|---|---|---|---|
| 5%×5% | 14,866 | 9,825 | 3,576 | 137 (111, 26 excluded) | **5.51** | (5.13, 5.91) | (4.25, 7.57) |
| **10%×10% (primary)** | **29,731** | **20,418** | **9,322** | **137 (124, 13 excluded)** | **6.09** | **(5.82, 6.38)** | **(4.37, 8.27)** |
| 20%×20% | 59,462 | 45,075 | 23,781 | 137 (131, 6 excluded) | **4.92** | (4.78, 5.06) | (3.78, 6.02) |

**Real, robust, substantial enrichment** — every cutoff's CI (both the
asymptotic MH CI and the donor-cluster bootstrap CI, which does not
assume within-donor independence) clearly excludes 1. Per the round-2
review correction, MH is **not** presented as having solved cell-level
pseudoreplication by itself — the donor-cluster bootstrap CI is the
inferential-uncertainty companion that does not make that assumption,
and it agrees qualitatively with the asymptotic CI (wider, as expected,
but still excluding 1 at every cutoff). **Note on scope**: leave-one-donor/study-out
sensitivity below is checked on the point-estimate (OR) scale only — CIs
were not recomputed at every excluded donor/study, so "excludes 1" is
stated for the full-cohort CIs above, not claimed for every individual
leave-one-out scenario.

**Round-1 fix — leave-one-study-out study attribution corrected**: the
round-1 draft stated "excluding `Tian_2023_Nat_Med` produces the largest
single-study shift at every cutoff," which is wrong at 2 of 3 cutoffs —
confirmed directly against the committed per-study tables. The largest
**absolute** shift is study-specific: `Joanito_2022_Nat_Genet` at 5%
(OR 5.51→9.04 excluded, Δ≈+3.53), `Tian_2023_Nat_Med` at 10% (OR
6.09→4.33, Δ≈−1.76), `Chen_2024_Cancer_Cell` at 20% (OR 4.92→6.58 excluded,
Δ≈+1.66). What **is** true at every cutoff: `Tian_2023_Nat_Med` gives the
largest *downward* shift (the lowest post-exclusion OR) at all three —
4.26 (5%), 4.33 (10%), 3.72 (20%). **The common OR never drops below 3.7
for any single study excluded, at any cutoff** — the enrichment finding
is not driven by one study, under either framing. Leave-one-donor-out
shifts are all small (max |Δ| 0.25–0.76 OR units across the 137
estimable donors). Non-estimable donors (zero-margin 2×2 table) are
26/137 (5% cutoff), 13/137 (10% cutoff), 6/137 (20% cutoff) — expected
(fewer cells qualify as "high" at smaller cutoffs, so more donors have
zero events in one cell of their own 2×2 table), reported transparently
in the per-donor table, and now excluded from MH pooling per the locked
design (see above).

Full per-donor 2×2/OR tables (all 137 donors, estimable flag included),
leave-one-donor/study-out tables (estimable donors only, per the fix
above), the correlation per-donor/leave-one-out detail, and the
same-population D/F/P sensitivity table are committed in full
(`results/06_crc_projection/secondary_analysis_concordance/`).

## What this shows, honestly

**M11 (an independently-discovered, expression-only NMF meta-program,
originally identified via its own Jaccard-overlap match to revCSC) shows
a real, robust, moderate-to-large concordance with revCSC** — both
continuously (r=0.32, gene-disjoint, same-population comparison against
D/F/P shows every D/F/P↔revCSC r in the identical 297,307-cell subset is
essentially null, |r|≤0.074) and as a threshold-based enrichment
(OR≈5–6, every cutoff's CI — asymptotic and donor-cluster bootstrap —
excluding 1, common OR never dropping below 3.7 under any single-study
exclusion checked). Developmentally, revCSC-high cells that show any
gut-axis evidence at all are dominated by `F_Gut-specific` (not D or P),
though "no axis detectable" is comparably common (41.1% vs. 42.6%) — a
real, bimodal pattern, not a clean single-axis story. Precisely
quantified via a Step0×StepA cross-tab: 21.4 of the unconditional
argmax's 71.7% F share (about 30% of all argmax-F assignments) comes
from cells with *no* supported developmental axis at all — direct,
measured evidence for why Step 0 (not argmax) is the design's primary
composition result. The F-assigned subset itself skews more SI-biased
than Colon-biased, a genuinely unexpected result for a colon-cancer
atlas, flagged here for future scrutiny rather than smoothed over.

**What this does NOT show**: this does not establish that M11 *is* the
"Oncofetal" fetal-gut program, nor that revCSC-high cells' F-association
constitutes a distinct, statistically separable substate (per the
locked design, that would require an additional pre-specified
clustering/mixture criterion, out of scope here — this analysis reports
"axis-defined composition/substructure," not formal separability). It
also does not demonstrate *why* the unconditional argmax favors F
(gene-set size is a plausible but untested mechanism — see Section 2).
It also does not extend beyond `CRC_single_cell_atlas_2025` or beyond the
already-frozen gut D/F/P and revCSC gene sets.

## What this does NOT do

Same explicit scope boundary as the approved design: does not touch the
tertiary analysis (full-atlas revCSC-independent D/F/P landscape); does
not extend to the 2 secondary/tertiary CRC datasets
(`HTAN_CRC_progressive_plasticity`, `CRLM_NMP_ATLAS`); does not
re-derive or re-score any frozen gut D/F/P, revCSC, or M11 gene set.

## Review history

- **Round 1 (REQUEST_CHANGES, 2 real implementation deviations from the
  locked design + 4 real write-up errors, all independently verified
  against real data before fixing)**: (1) **Blocker**: `stepA` coarse
  argmax was reported as a separate marginal, not cross-tabulated against
  Step 0's axis-supported status as the locked design explicitly
  requires — fixed by adding an explicit Step0×StepA cross-tab
  (`step0_x_stepA_crosstab.tsv` per cohort), which also let the argmax
  "distortion" be quantified precisely (21.4 points from `none`-status
  cells) instead of by an arithmetic difference of two marginals. (2)
  **Blocker**: MH enrichment pooling passed all donor strata (including
  non-estimable, zero-cell ones) into the pooled OR/CI/leave-one-out/
  bootstrap, contradicting the locked design's explicit "excluded from
  the MH pooling" instruction — confirmed as a real deviation (checked
  directly whether unconditional retention would actually be more
  statistically standard: it would, but the design must be executed as
  approved) — fixed to filter to estimable donors only for all pooled
  quantities, exclusion count reported; the numerical shift is small
  (OR_MH 5.79→5.51 at 5%, 6.18→6.09 at 10%, 4.92→4.92 at 20%), no
  qualitative change. (3) The composition write-up's "~29 point"
  argmax-inflation claim was arithmetically wrong (71.7%−52.9%=18.8, not
  29) — fixed with the precise cross-tab-derived number (21.4 points from
  `none`-status cells) instead. (4) The write-up's causal claim
  ("because `F_Gut-specific` has 2,192 genes") was not demonstrated by
  this compute (the scoring method null-calibrates against matched-null
  gene sets specifically to control for size) — removed, narrowed to the
  defensible observational claim. (5) The leave-one-study-out prose
  misattributed the largest shift to `Tian_2023_Nat_Med` at every cutoff;
  the largest *absolute* shift is actually study-specific
  (`Joanito_2022_Nat_Genet` at 5%, `Tian_2023_Nat_Med` at 10%,
  `Chen_2024_Cancer_Cell` at 20%) — `Tian_2023_Nat_Med` does give the
  largest *downward* shift at all three, which is what survives in the
  fixed prose. (6) "M11 substantially stronger than any D/F/P axis"
  compared M11's within-M11-subset r against PR #27's full-atlas r — a
  cross-population comparison — fixed by directly computing the same
  D/F/P↔revCSC correlations restricted to the identical M11 subset (no
  new scoring needed); the same-population comparison turns out even
  more favorable to the original claim (D/F/P↔revCSC all |r|≤0.074 within
  the M11 subset, vs. M11's r=0.318). Two smaller wording fixes: narrowed
  the "CI excludes 1 under every exclusion" claim to the full-cohort CIs
  specifically (per-exclusion CIs were not recomputed); corrected the
  committed/verified file-count wording. M11 scoring itself (job 3621118)
  was not rerun — round 1 review confirmed this was unaffected by any of
  the 6 issues; only composition/concordance (job 3621126) reran.

Submitting for round-2 review before merge, same discipline as every
prior step.
