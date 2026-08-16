# Step 6 dataset extension: real results — HTAN_CRC_progressive_plasticity + CRLM_NMP_ATLAS

Real compute, per the approved design (`docs/STEP6_DATASET_EXTENSION_DESIGN.md`,
APPROVE after 3 review rounds). Closes the last remaining item of the
original 3-dataset plan (`docs/STEP6_CRC_PROJECTION_DESIGN.md`) and,
with it, the last substantive item toward this project's "100%" scope.
All numbers below are real qsub output, pulled back and verified
byte-exact (md5, 62 files across two jobs) before being trusted or
written up, and every headline number was independently hand-recomputed
from the committed per-patient/per-donor tables before being reported —
same discipline as every prior step.

## Job provenance

| Stage | Job ID | Scope | Outcome |
|---|---|---|---|
| **Scoring** (loader, coverage check, marker check, provenance audit, 3 scoring passes) | **3621230** | HTAN 26,551 malignant cells (pass 1) + 47,107 malignant+normal cells (pass 2); CRLM 4,051 malignant cells | Completed cleanly, <7 min; 13 files, byte-exact |
| **Analysis** (correlation + patient-matched paired contrast) | **3621236** | Reads pass-1/2/3 scores, no new scoring | Completed cleanly, <1 min; 49 files, byte-exact |

No per-cell files committed beyond the scoring passes' own outputs
(consistent with the "small deliverable" norm from every prior Step 6
round). All output tables committed in full
(`results/06_crc_projection/dataset_extension/`,
`results/06_crc_projection/dataset_extension_analysis/`).

## Pre-compute gates (all required by the design, all run)

**Gene-ID axis contract**: both datasets confirmed to have `raw.var_names`
byte-identical to the working `var_names` (bare Ensembl IDs, no version
suffixes) — the version-suffix canonicalization branch is a real,
executed assertion that resolved to a no-op, not a skipped check.

**Coverage check** (`coverage_check.tsv`, all 13 panels × 2 datasets):
coverage is high for both datasets on every panel except one real
exception — CRLM's `P_Gut-specific` panel resolves only 40/76 genes
(52.6%), the lowest coverage figure of any panel/dataset pair (next
lowest: CRLM `D_Gut-shared` 7/8 = 87.5%). All other panels clear ≥89%
coverage on both datasets. This is reported as-is; CRLM's already-
exploratory framing (see below) absorbs this, not smoothed over.

**Canonical-marker sentinel check** (`canonical_marker_sentinel_check.tsv`):
EPCAM/KRT8/KRT18/KRT19 are elevated (mean log1p expression ~1.4–3.5)
across every epithelial cell type in both datasets (malignant and every
named normal subtype), and `PTPRC` (CD45, immune marker) is uniformly
near-zero (0.003–0.12) across all groups — confirms both datasets are
genuinely epithelial-compartment data with no immune contamination, and
that the loader's raw-counts/normalization split is biologically sane,
not just structurally correct.

**HTAN study-provenance overlap audit** (`htan_provenance_overlap_audit.tsv`,
`..._conclusion.txt`): **NO_OVERLAP_FOUND**, checked programmatically
against every identifier column (`dataset`, `sample_id`, `patient_id`,
`donor_id`, and 5 accession columns) across all 36 of
`CRC_single_cell_atlas_2025`'s constituent studies — not just the one
whose name suggested HTAN provenance. That study, `HTAPP_HTAN`, turned
out on inspection to be a **different HTAN sub-cohort**: its 19,711
cells carry `HTA1_`-prefixed identifiers sourced from Pelka et al.
2021's Synapse project (`syn24181445`), while
`HTAN_CRC_progressive_plasticity`'s 29 patients carry exclusively
`HTA8_`-prefixed `donor_id` values — zero string-level overlap across
every identifier column checked. Per the design's locked claim-gate
contract, this HTAN dataset can now be described as external validation
(not merely a same-provenance internal check), though it is still not
formally "independent replication" of any specific claim, only
confirmed non-overlapping cohort membership.

## Section 1: primary-extension correlation (revCSC↔D/F/P, malignant cells only)

**HTAN** (26,551 malignant cells, 29 patients, 28 estimable): all 10
comparison pairs remain weak, |r|≤0.10 — consistent with PR #27's
primary-atlas finding that no clean single-axis "Oncofetal = F" or "= P"
relationship exists. 3 of 5 primary-revCSC pairs are robust to
leave-one-donor-out (both D/P pairs, plus the broad `F_Gut-specific`
pair); the two regional-F pairs (`F_Colon-specific`, `F_SI-specific`)
are not robust — the same qualitative pattern (D/P more stable than
regional F) PR #27 found on the full 665,473-cell atlas, now
independently observed on a genuinely non-overlapping population.

**A real, honestly-reported nuance**: correlation *signs* are not
consistent between HTAN and the primary atlas for the F/P pairs, even
though magnitudes are weak in both. `revCSC_primary27_full`↔`P_Gut-specific`
is r=−0.025 in the primary atlas (`primary_analysis_overview.tsv`) but
r=+0.102 in HTAN; `revCSC_primary27_minus_CLU_ASS1`↔`F_Gut-specific` is
r=+0.101 in the primary atlas but r=−0.067 in HTAN. This is not a
contradiction of the core "weak" finding (both magnitudes are ≤0.10),
but it does mean the *sign* of these already-weak F/P correlations
should not be treated as a stable, generalizable feature — reported here
plainly rather than cherry-picking the consistent pairs.

**CRLM** (4,051 malignant cells, 6 donors, exploratory only per the
locked design — no robustness classification computed): correlations
are somewhat larger in magnitude (up to r=−0.21 Pearson /
ρ=−0.32 Spearman for `revCSC_extended28_minus_CLU`↔`F_Colon-specific`),
but with only 6 donors this is explicitly not interpreted as evidence of
a stronger relationship — per the design, no LODO-based robust/non-robust
call is made for CRLM at all (the `robust_*` columns are blank in the
committed overview, by design, not omitted from the table).

## Section 2: HTAN patient-matched malignant-vs-normal-epithelial contrast

**This is the most direct, and most striking, real finding of the
dataset extension.** Using the joint-calibrated scoring pass (malignant
+ normal epithelium scored together, so percentiles are on the same
scale) and the locked patient-level paired contract (one malignant
summary + one normal-epithelial summary per patient, patients — not
cells — as the statistical unit, 26 of 29 patients had both cell types
present and were paired):

| Panel | n patients paired | Median malignant − normal Δ (percentile pts) | Wilcoxon p | Paired-t p |
|---|---|---|---|---|
| `F_Gut-specific` | 26 | **+82.4** | 3.3×10⁻⁵ | 6.1×10⁻⁸ |
| `F_Colon-specific` | 26 | **+72.4** | 6.0×10⁻⁸ | 9.3×10⁻¹⁰ |
| `F_SI-specific` | 26 | **+80.0** | 3.3×10⁻⁵ | 1.2×10⁻⁷ |
| `revCSC_primary27_minus_CLU_ASS1` | 26 | +11.4 | 6.6×10⁻⁴ | 1.8×10⁻³ |
| `D_Gut-shared` | 26 | +0.65 | 0.451 | 0.384 |
| `P_Gut-specific` | 26 | −2.0 | 0.648 | 0.999 |

**All three F-developmental axes show a massive, highly significant
shift** — malignant cells sit ~72–82 percentile points higher than that
same patient's own normal epithelium, both by the conservative
rank-based Wilcoxon test and the parametric paired t-test. `revCSC`
itself shows a real but much smaller shift (+11.4 points, still
significant). **`D_Gut-shared` and `P_Gut-specific` show no significant
difference at all** — patient-matched malignant and normal epithelial
cells are statistically indistinguishable on these two axes (p=0.45–0.99
across both tests).

This is the single most direct evidence of "oncofetal reactivation" in
this entire project: within the same patient, controlling for
patient-level confounding, malignant cells are dramatically more
fetal-somatic-gut-like than their own normal epithelium — but they are
**not** more developmentally-shared or more placental/trophoblast-like.
The finding is F-specific, not D-shared or P-specific, and it holds
consistently across both regional F axes (Colon and SI), not just the
coarse F_Gut-specific union.

**Per-named-subtype sensitivity check** (`htan_patient_matched_contrast_by_normal_subtype.tsv`,
descriptive only per the locked design, not a second primary analysis):
the F-panel malignant-vs-normal shift replicates against nearly every
individual normal epithelial subtype (BEST4+, goblet, colonocyte, early
colonocyte, secretory, tuft — all p<0.001, deltas 58–97 points), with
one partial exception: `intestinal crypt stem cell of colon` shows a
much smaller, non-significant F-panel delta (e.g. F_Gut-specific
Δ=6.5, p=0.37) — plausibly because stem cells already carry some
fetal-like developmental character themselves, making them a less sharp
contrast to malignant cells than fully differentiated subtypes. This is
flagged as an interesting pattern, not investigated further here. One
small-n subtype result stands out for D/P (enteroendocrine cells,
n=18: D_Gut-shared Δ=22.75, p=7.6×10⁻⁶; P_Gut-specific Δ=8.1, p=2.0×10⁻⁴)
— reported for completeness but not treated as overturning the primary
all-subtypes-pooled null result for D/P, given the much smaller n and
its status as a sensitivity check, not the primary contrast.

## What this shows, honestly

**The primary-extension correlation results generalize PR #27's core
finding** (weak, not-single-axis revCSC↔D/F/P relationships) to a
dataset with zero confirmed provenance overlap with the primary atlas —
a real, if modest, piece of independent corroboration for that finding's
generality. The sign-instability of the already-weak F/P correlations
between datasets is reported honestly rather than smoothed over.

**The patient-matched contrast is a genuinely new, and genuinely
striking, result this project has not had before**: a *within-patient*,
confound-controlled test showing malignant cells are massively more
F-developmental (fetal-somatic-gut) than their own normal epithelium,
while showing no such shift on D-shared or P-specific axes. This
sharpens, rather than contradicts, the project's core developmental
question — the reactivation this project set out to characterize looks,
by this most direct test available, F-specific rather than D-shared or
P-specific.

**What this does NOT show**: it does not establish that this
malignant-vs-normal F-shift is *causally* driven by any specific
mechanism (SPP1+ TAM/TWEAK/Fn14/YAP or otherwise — Q5, explicitly out of
scope). It is a strong test but not a fully confound-free one — HTAN's
malignant and normal-epithelial samples per patient are not guaranteed
to be anatomically/treatment-matched beyond sharing the same patient
(per the design's own stated caveat). It does not extend to CRLM (no
comparable normal population exists there) or establish anything about
formal biological separability of an F-high malignant state (that would
need a pre-specified clustering/mixture model, out of scope here, same
boundary maintained throughout Step 6).

## What this does NOT do

Does not re-derive or re-score any frozen D/F/P or revCSC gene set
(reused as-is, confirmed via the coverage check rather than assumed).
Does not re-run the secondary (revCSC-high composition) or tertiary
(full-atlas composition) analysis structure on either dataset (locked
scope decision, design section 6). Does not extend to `GSE178318` (no
cell-type annotation exists, already deferred). Does not touch Q5/Q6 of
the original 6-question framework (explicitly out of scope per the
2026-08-16 user confirmation, `docs/Q1_Q2_Q4_CROSS_REFERENCE.md`).

## Project-completion note

This closes the last remaining substantive item toward this project's
user-confirmed "100%" scope (`docs/Q1_Q2_Q4_CROSS_REFERENCE.md`,
"100% = close what's reachable with current data"). Q5/Q6 remain
explicitly out of scope, documented as a distinct future aim.
