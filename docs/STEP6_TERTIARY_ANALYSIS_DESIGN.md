# Step 6 tertiary analysis: full-atlas, revCSC-independent D/F/P landscape — design

Implements the **tertiary analysis** item from the already-approved Step 6
design (`docs/STEP6_CRC_PROJECTION_DESIGN.md`, "Revised analysis structure",
item 3): *"Score D/F/P across the full-atlas cohort without reference to
revCSC at all. This can surface a state the 'Oncofetal' framework itself
might miss — e.g. a P-high/revCSC-low malignant population that exists
outside what revCSC captures."*

**Also directly answers Q3 of the project's original 6-question framework**
(`2026-GPT-TWEAKR-Oncofetal.md#定义清楚Placenta的问题`, read in full for the
first time this session, user-confirmed 2026-08-16 scope): *does the
malignant cell population — not just the revCSC-defined subset — split into
separable Fetal(F) × Placental(P) [× Shared(D)] developmental quadrant
states, rather than being one program?* Step 6 secondary (PR #28/#29)
already answered this question **within the revCSC-high subset**; this
analysis extends it to **every malignant cell in the atlas, independent of
revCSC status**, which is the genuinely new, non-circular check — a
population could show a real F/P/D quadrant structure that revCSC's own
27-gene definition never captures at all.

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
  it directly answers Q3's actual question (coactivation vs. discrete
  dominance vs. no program), unlike a forced single-label argmax.
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

## 4. Cross-tabulation against revCSC status (connects back to Q1/Q2, does not gate)

For every cell already classified in §2/§3, additionally report its
`revCSC_primary27_minus_CLU_ASS1` cross-cell-rank status (same construction
as Step 6 secondary §1: top 10%/5%/20% by cross-cell z-score rank, computed
identically here since the same score column and same full-atlas population
are being used). This produces a joint table: axis-supported category ×
revCSC-high/not-revCSC-high — directly testable for the design's stated
target finding: **does a P-high (or F+P, D+F+P) population exist that is
NOT revCSC-high?** (I.e., a "P-high/revCSC-low" cell state existing outside
what revCSC's own 27-gene signature captures — the design's explicit
motivating question, checked directly rather than assumed.)

This cross-tab is reported, not used to define any cohort — revCSC status
never gates which cells are included in §2/§3's population-level results.

## 5. Donor/study-aware aggregation

Every table in §2–§4 reported pooled, plus unweighted mean across donors
(composite `donor_key`) and across studies — identical machinery to
`secondary_analysis_composition.py`'s `categorical_donor_study_summary`,
reused directly (not reimplemented) on the full 665,473-cell population
instead of the revCSC-high subset.

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
§3's Q3-quadrant table, §4's revCSC cross-tab, §5's donor/study
aggregation for each — writes summary/overview tables only (no new
per-cell files), same "small deliverable" size norm as Step 6 secondary's
composition output.

Submitting for review before any qsub job runs, same discipline as every
prior step.
