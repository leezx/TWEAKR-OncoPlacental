# Step 6 tertiary analysis: full-atlas, revCSC-independent D/F/P landscape — design

Implements the **tertiary analysis** item from the already-approved Step 6
design (`docs/STEP6_CRC_PROJECTION_DESIGN.md`, "Revised analysis structure",
item 3): *"Score D/F/P across the full-atlas cohort without reference to
revCSC at all. This can surface a state the 'Oncofetal' framework itself
might miss — e.g. a P-high population that exists outside the operational
revCSC-high cohort used in PR #28/#29."*

**Also directly answers the axis-composition component of Q3 of the
project's original 6-question framework**
(`2026-GPT-TWEAKR-Oncofetal.md#定义清楚Placenta的问题`, read in full for the
first time this session, user-confirmed 2026-08-16 scope): *does the
malignant cell population — not just the revCSC-defined subset — show
axis-defined F(Fetal) × P(Placental) [× D(Shared)] developmental
composition/substructure, rather than being one uniform program?* Step 6
secondary (PR #28/#29) already answered this **within the revCSC-high
subset**; this analysis extends it to **every malignant cell in the atlas,
independent of revCSC status**, which is the genuinely new, non-circular
check — a population could show real F/P/D quadrant occupancy that
revCSC's own 27-gene definition never selects for at all.

**Round-1 review correction (applied from the start, not discovered after
compute)**: this design does **not** claim to establish that these
quadrants are *separable states* in the formal statistical sense (distinct
modes/clusters). Per PR #28's own locked boundary — restated here rather
than re-litigated — axis-supported status and a threshold-based quadrant
table can show *axis-defined composition/substructure* (**multi-axis
support vs. single-axis support vs. none** — round-2 correction: not
"dominance," which is specifically what Step A′'s Δ5-margin predominance
call establishes, not what the 8-category support table alone
establishes), but **cannot by themselves distinguish that from a smooth
continuous F×P gradient** thresholded at the same 90th-percentile bar; a
genuine separability claim would require an additional pre-specified
clustering/mixture criterion, out of scope here. Every result in this
design is described as "quadrant occupancy" / "axis-defined composition"
throughout, never as "separable states" or "dominance" outside Step A′.

**No new scoring required.** All D/F/P and revCSC per-cell percentiles for
all 665,473 `CRC_single_cell_atlas_2025` malignant cells already exist,
computed and byte-verified in PR #27's full run
(`results/06_crc_projection/gut_scoring_full/crc_gut_scoring_all_panels.parquet`).
This is a new *analysis* of already-computed data, not a new compute job —
same "reuse before rescoring" discipline as Step 6 secondary.

## 1. Population: all 665,473 cells, unconditional on revCSC

Unlike Step 6 secondary (which restricted to the revCSC-high cohort), this
analysis runs on **every malignant cell that was scored** — no revCSC-based
filtering of any kind at the cohort-definition stage. revCSC's own score is
used only *afterward*, as a cross-tabulated variable to connect this
analysis back to Q1/Q2 (see §4) — never to define which cells are included.

## 2. Axis-supported status (primary result) — same construction as Step 6 secondary, same reasons

Per PR #28/#29's locked, review-tested methodology (not re-litigated here):

- **Axis-supported** = null-calibrated percentile ≥90 (evidence versus each
  cell's own matched null), independently for each of
  `{D_Gut-shared, F_Gut-specific, P_Gut-specific}` — the same threshold
  already used throughout Step 6, not a new pick.
- **8 mutually exclusive, jointly exhaustive categories**: none / D_only /
  F_only / P_only / D+F / D+P / F+P / D+F+P. This is the primary result —
  it directly answers Q3's axis-composition question (**multi-axis
  support vs. single-axis support vs. none**, not "dominance" — see the
  round-2 correction above), unlike a forced single-label argmax.
- **Coarse 3-way argmax retained only as a cross-tabulated descriptive
  summary** (Step0×StepA cross-tab, per PR #29 round-1's fix) — never used
  alone as evidence of "multiple separable states." PR #29 already showed
  directly, on real data, that argmax alone overstates a single axis's
  share by double-digit percentage points relative to axis-supported
  status; the same distortion risk applies here and is guarded against the
  same way from the start, not discovered again in review.
- **Predominant axis among supported cells only**, Δ5 percentile-point
  margin rule (identical to PR #28/#29's Step A′) — "no clear
  predominance" reported explicitly, not forced.
- **Regional refinement** (`F_Colon-specific` vs. `F_SI-specific`) within
  F-supported cells, same Δ5 margin rule (identical to Step B).

## 3. The explicit Q3 quadrant framing — F × P 2D space, D reported alongside

Q3's original framing is a 2D space (Fetal-specific score × Placenta-specific
score) with 4 quadrants: Fetal-high/Placenta-low, Fetal-high/Placenta-high,
Fetal-low/Placenta-high, double-low. This project's D/F/P framework has a
third axis (D-shared) that Q3's original 2D sketch doesn't have a slot for,
so the quadrant table is reported **both ways**, not by silently dropping D:

- **Primary table**: the full 8-category axis-supported breakdown from §2
  (D/F/P jointly), which is strictly more informative than a 2D quadrant
  and does not force D-shared evidence to be ignored.
- **Q3-quadrant table (explicit, for direct comparability to the original
  question)**: collapsed to F×P only (`F_Gut-specific` axis-supported ×
  `P_Gut-specific` axis-supported, ignoring D), the 4 cells being
  Fetal-high/Placenta-low, Fetal-high/Placenta-high, Fetal-low/Placenta-high,
  double-low — using the exact same percentile≥90 threshold, not a
  different cutoff picked for this table. This is reported as a
  **derived, secondary view of the same primary data**, not an
  independent analysis with its own threshold-picking risk.

## 4. Cross-tabulation against revCSC status (connects back to Q1/Q2, does not gate) — two distinct tables, not one

**Round-1 review correction**: the first draft conflated two different
quantities — "not in the top-10% cross-cell-rank cohort" and "revCSC does
not have evidence for this cell" — which PR #28 already established are
not equivalent (the null-calibrated percentile's own distribution is
heavily right-skewed, median ≈90.2 across this atlas, so a cell can fail
the cross-cell top-decile cut while still showing strong revCSC evidence
versus its own matched null). Fixed by reporting **two separate,
correctly-labeled tables**, using the same matched-null support semantics
(percentile≥90) already used for D/F/P in §2, applied identically to
revCSC:

- **Table 4a — axis-supported category × operational revCSC-high cohort**
  (`revCSC_primary27_minus_CLU_ASS1_zscore`, cross-cell rank, top
  10%/5%/20%, identical construction to Step 6 secondary §1). Interpreted
  strictly as: *"does a P-supported (or F+P, D+F+P) population exist
  outside the operational revCSC-high cohort PR #28/#29 used for the
  secondary analysis?"* — a statement about a specific, already-defined
  cohort, not about revCSC evidence in general.
- **Table 4b — axis-supported category × revCSC matched-null support**
  (`revCSC_primary27_minus_CLU_ASS1_percentile` ≥90, the same
  null-calibrated-evidence semantics as every D/F/P axis in §2, not a
  cross-cell rank). **Round-2 wording correction**: percentile <90 means
  the cell does not meet the pre-specified revCSC-support threshold — it
  is not literally "absence of evidence," which would overclaim a sharp
  biological boundary that a threshold call does not establish. Interpreted
  as: *"does P-supported biology exist in cells that are revCSC-not-supported
  at the same percentile≥90 threshold used for D/F/P?"* — the correct
  construction for the design's actual motivating question, requiring no
  new scoring since `revCSC_primary27_minus_CLU_ASS1_percentile` is
  already computed for all 665,473 cells.

Neither table gates which cells are included in §2/§3's population-level
results — both are reported cross-tabs computed after the fact, exactly
as in the first draft; only the labeling and the addition of Table 4b are
new.

## 5. Donor/study-aware aggregation

Every table in §2–§4 **except the Step0×StepA cross-tab** (see the
round-2 correction immediately below — that one table is pooled-only by
construction, same as PR #29's precedent) is reported pooled, plus
unweighted mean across donors (composite `donor_key`) and across studies
— identical machinery to `secondary_analysis_composition.py`'s
`categorical_donor_study_summary`, reused directly (not reimplemented) on
the full 665,473-cell population instead of the revCSC-high subset.

**Round-2 wording correction**: `step0_x_stepA_crosstab()` (PR #29) is
**explicitly pooled-only by construction** — it does not, and is not
required to, give donor/study aggregation, matching PR #29's own
precedent for this exact descriptive/cross-tabulated table (argmax is
never more than a descriptive summary, so its cross-tab against Step 0
is reported at the pooled level only). The
Q3 quadrant table (§3) and both revCSC joint tables (§4a/§4b) are **not**
implemented via that pooled-only crosstab helper; instead, each joint
category (e.g. `"P_only|revCSC_high_top10pct"`) is built as a single
composite string label and passed through `categorical_donor_study_summary`
directly, so every joint table gets the same pooled + donor-unweighted +
study-unweighted treatment as §2's tables — no new statistical procedure
invented, just consistent reuse of the one function that already does
donor/study aggregation correctly.

## 6. Scope boundary (explicit)

- Still `CRC_single_cell_atlas_2025` only — the 2 additional CRC datasets
  (`HTAN_CRC_progressive_plasticity`, `CRLM_NMP_ATLAS`) are separate,
  already-scoped future work, not folded into this analysis.
- Does not re-score or re-derive any frozen gut D/F/P or revCSC gene set,
  and does not touch M11 (M11 is only defined on its own 297,307-cell
  subset, not the full atlas — out of scope here by construction).
- Reuses PR #27's already-verified per-cell percentiles and PR #28/#29's
  already-reviewed statistical methodology (axis-supported status, Δ5
  margin, donor/study-aware aggregation) without modification — this
  design's only new contribution is applying that methodology to the
  full, revCSC-unconditional population and adding the explicit Q3
  quadrant view + revCSC cross-tab.

## Compute plan

One lightweight script (no qsub required for the analysis itself — reads
the already-scored parquet already on Argos; run via qsub for consistency
with standing discipline anyway, since it touches the full-atlas file):
loads `crc_gut_scoring_all_panels.parquet` + `crc_gut_scoring_cell_metadata.parquet`
(no new scoring), computes §2's 8-category status + Step0×StepA cross-tab,
§3's Q3-quadrant table, §4a/§4b's two revCSC cross-tabs, and §5's
donor/study aggregation for every table **except the intentionally
pooled-only Step0×StepA cross-tab** (via composite joint-labels through
`categorical_donor_study_summary` for the Q3/4a/4b tables) — writes
summary/overview tables only (no new per-cell files), same "small
deliverable" size norm as Step 6 secondary's composition output.

## Review history

- **Round 1 (REQUEST_CHANGES, 2 real conceptual blockers, 1 implementation
  watch item — all confirmed directly against the design doc's own
  wording before fixing)**: (1) the first draft's Q3 framing ("split into
  separable Fetal×Placental quadrant states") overclaimed formal
  statistical separability that a threshold-based quadrant table cannot
  establish on its own (a smooth continuous F×P gradient thresholded at
  the same bar produces the identical table) — fixed by reframing every
  claim in this design as "axis-defined composition/substructure" /
  "quadrant occupancy," explicitly stating the separability boundary
  rather than re-discovering it after compute, matching PR #28's already
  locked precedent. (2) §4's revCSC cross-tab conflated "not in the
  top-10% cross-cell-rank cohort" with "revCSC has no evidence for this
  cell" — two different quantities per PR #28's own established
  right-skew finding (median revCSC percentile ≈90.2) — fixed by
  splitting into two correctly-labeled tables: 4a (operational
  revCSC-high cohort membership, PR #28/#29's exact construction) and 4b
  (revCSC's own matched-null percentile≥90 support, the same semantics
  already used for D/F/P), the latter being the construction that
  actually answers the design's motivating question. (3) Implementation
  watch: the new Q3/revCSC joint tables must use
  `categorical_donor_study_summary` (donor/study-aware) via composite
  joint-category labels, not the pooled-only `step0_x_stepA_crosstab`
  helper — made explicit in §5 rather than left to be caught in a later
  compute review.

- **Round 2 (REQUEST_CHANGES, "narrowly" — both round-1 blockers confirmed
  resolved, 3 small wording/contract cleanups, no new methodology
  needed)**: (1) leftover "coactivation vs. discrete dominance vs. no
  program" phrasing in the round-1 correction and §2 incorrectly implied
  the 8-category support table establishes *dominance*, which is
  specifically Step A′'s job (Δ5-margin predominance) — fixed to "multi-
  axis support vs. single-axis support vs. none" throughout. (2) Table
  4b's "lack matched-null evidence" / "without revCSC evidence" wording
  overclaimed a sharp absence-of-evidence boundary that a percentile<90
  threshold call does not establish — fixed to "revCSC-not-supported at
  percentile≥90." (3) §5's "every table in §2–§4" claim contradicted the
  intentionally pooled-only `step0_x_stepA_crosstab` (PR #29's own
  precedent for that specific descriptive table) — fixed by explicitly
  exempting it rather than leaving an internal contradiction.

- **Round 3 (REQUEST_CHANGES, one residual internal-contradiction leftover
  only)**: the "Compute plan" section still said §5 gives donor/study
  aggregation "for every table ... not the pooled-only crosstab helper,"
  contradicting §5's own just-fixed explicit exemption of the Step0×StepA
  cross-tab — a one-line miss when applying round 2's fix, not a new
  design issue. Fixed to state the exemption consistently in both places.

Submitting for round-4 review before any qsub job runs, same discipline
as every prior step.
