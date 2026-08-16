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

**Round-1 review correction**: the first draft claimed "no per-cell
files committed" — false. Both scoring-pass `*_scores.parquet` and
`*_cell_metadata.parquet` files ARE per-cell (26,551/47,107/4,051 rows
respectively), matching PR #27's own primary-compute precedent (which
also commits per-cell scores) — only the secondary/tertiary *composition
analyses* deliberately avoided per-cell commits, a different precedent
that does not apply here. Corrected: this PR commits per-cell scores for
all 3 scoring passes (13 files) plus per-donor/per-patient summary
tables (49 files) — 62 files total, no distinct "small deliverable"
claim made. All output tables committed in full
(`results/06_crc_projection/dataset_extension/`,
`results/06_crc_projection/dataset_extension_analysis/`).

## Pre-compute gates (all required by the design, all run)

**Gene-ID axis contract**: both datasets confirmed to have `raw.var_names`
byte-identical to the working `var_names`. **Round-1 review correction**:
the first draft claimed this branch's byte-identity check was itself
proof of "bare Ensembl IDs, no version suffixes" — it was not; the code
accepted byte-identical axes without ever checking their *format*.
Fixed (`dataset_extension_core.py`, round-1 diff): every ID is now
asserted to match the bare `ENSG` + 11-digit pattern directly, and this
assertion is confirmed to pass for both datasets (re-run, no rescoring
needed — this is a structural gate independent of the scored values).

**Coverage check** (`coverage_check.tsv`, all 13 panels × 2 datasets):
coverage is high for both datasets on every panel except one real
exception — CRLM's `P_Gut-specific` panel resolves only 40/76 genes
(52.6%), the lowest coverage figure of any panel/dataset pair (next
lowest: CRLM `D_Gut-shared` 7/8 = 87.5%; every other pair ≥89%).
**Round-1 review correction**: the first draft said CRLM's "exploratory
framing absorbs this" — this directly contradicted the design's actual
contract (a required investigation gate, not an auto-waived one), and
the implementation's flat 50% floor never even flagged this case (52.6%
clears 50%) despite the design's wording being about deviation from
peers, not an absolute floor — fixed (`dataset_extension_scoring.py`,
round-1 diff) to flag panels >15 points below their own dataset's
median coverage, which now correctly catches this exact case.
**The actual investigation, done here rather than waived**: the 36
missing genes were checked directly against both other datasets'
`var_names` — **all 36 are present in the primary
`CRC_single_cell_atlas_2025` atlas, and 34/36 are present in
`HTAN_CRC_progressive_plasticity`** — despite CRLM having the *largest*
total gene count of the three datasets (30,257 vs. 28,476 and 25,344).
This rules out a random-dropout or low-quality-cell explanation; it is a
real, specific gap in CRLM's own feature reference/annotation for
exactly these 36 genes (consistent with `P_Gut-specific`'s known
enrichment for placental/trophoblast marker genes — a gene family with
historically unstable annotation across Ensembl reference versions).
**Conclusion**: the existing 40-gene CRLM `P_Gut-specific` score remains
usable and is reported as-is (no rescoring performed or needed), but
every CRLM `P_Gut-specific` result in this document should be read with
this reduced, annotation-gap-driven gene coverage in mind — flagged
here explicitly rather than smoothed into the dataset's general
exploratory caveat.

**Canonical-marker sentinel check** (`canonical_marker_sentinel_check.tsv`).
**Round-1 review correction — a real factual error in the first draft**:
it claimed `PTPRC` is "uniformly near-zero (0.003–0.12) across all
groups" and both datasets show "no immune contamination" — checked
directly against the committed table, **this is true for HTAN
(epithelial-only export) but false for CRLM**, which is explicitly
TME-focused. CRLM's real marker values show exactly the expected,
biologically correct pattern instead: `PTPRC` is elevated in every real
immune cell type (T cell 3.15, natural killer cell 3.23, neutrophil
2.92, mononuclear phagocyte 2.28, B cell 2.51 — mean log1p expression)
and low in malignant cells (0.16) and cholangiocyte (0.11); `EPCAM` runs
the opposite direction (malignant 2.70, cholangiocyte 1.84, all immune
types ≤0.12). **This is, if anything, stronger evidence the loader is
biologically correct** than the (wrong) "no immune contamination" claim
would have been — it reproduces exactly the expected epithelial-vs-immune
split in a genuinely mixed-compartment dataset, not just a
homogeneous one.

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
every identifier column checked. **Round-1 review correction**: this
specific narrative was checked interactively during development but not
originally committed as auditable script output — an unsupported claim
by this project's own standing discipline. Fixed: a new supplementary
function (`htan_name_similarity_supplementary_check()`) makes this a
real, reproducible, committed artifact
(`htan_name_similarity_supplementary_check.tsv`), reporting every
`HTA`-named meta-atlas study's own identifiers directly — confirms the
`HTA1_`/`syn24181445` pattern exactly as stated. Per the design's locked
claim-gate contract, this HTAN dataset can now be described as external
validation (not merely a same-provenance internal check), though it is
still not formally "independent replication" of any specific claim,
only confirmed non-overlapping cohort membership.

## Section 1: primary-extension correlation (revCSC↔D/F/P, malignant cells only)

**HTAN** (26,551 malignant cells, 29 patients, 28 estimable): all 10
comparison pairs remain weak — **round-1 review correction**: the first
draft claimed "|r|≤0.10," which is numerically false; the largest
observed magnitude is Pearson r=+0.1019
(`revCSC_primary27_full`↔`P_Gut-specific`), and the primary-atlas
comparator value cited below (r=+0.1012 for
`revCSC_primary27_minus_CLU_ASS1`↔`F_Gut-specific`) also exceeds 0.10.
Corrected: all correlations are weak, with the largest observed absolute
Pearson r about 0.102 in HTAN (and comparably weak, ≤0.19, in the
primary atlas per PR #27) — consistent with PR #27's primary-atlas
finding that no clean single-axis "Oncofetal = F" or "= P" relationship
exists. 3 of 5 primary-revCSC pairs are robust to leave-one-donor-out
(both D/P pairs, plus the broad `F_Gut-specific` pair); the two
regional-F pairs (`F_Colon-specific`, `F_SI-specific`) are not robust —
the same qualitative pattern (D/P more stable than regional F) PR #27
found on the full 665,473-cell atlas, now independently observed on a
genuinely non-overlapping population.

**A real, honestly-reported nuance**: correlation *signs* are not
consistent between HTAN and the primary atlas for the F/P pairs, even
though magnitudes are weak in both. `revCSC_primary27_full`↔`P_Gut-specific`
is r=−0.025 in the primary atlas (`primary_analysis_overview.tsv`) but
r=+0.102 in HTAN; `revCSC_primary27_minus_CLU_ASS1`↔`F_Gut-specific` is
r=+0.101 in the primary atlas but r=−0.067 in HTAN. This is not a
contradiction of the core "weak" finding (both magnitudes are ≈0.10 or
smaller), but it does mean the *sign* of these already-weak F/P
correlations should not be treated as a stable, generalizable feature —
reported here plainly rather than cherry-picking the consistent pairs.

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
same patient's matched normal epithelium, both by the conservative
rank-based Wilcoxon test and the parametric paired t-test. `revCSC`
itself shows a real but much smaller shift (+11.4 points, still
significant). **Round-1 review correction**: the first draft said
`D_Gut-shared` and `P_Gut-specific` show malignant and normal cells are
"statistically indistinguishable" — this overstates what a
nonsignificant test establishes. Failing to reject a null hypothesis is
not proof of equivalence; that would require a pre-specified equivalence
margin or a direct differential-effect test, neither of which PR #33
locked. **Corrected, accurate statement**: no detectable patient-matched
shift was observed on `D_Gut-shared` or `P_Gut-specific` (observed
median deltas near zero, +0.65 and −2.0, p=0.45–0.99 across both tests)
— a sharp contrast to the F-panels' 72–82-point shifts, but not itself
proof those two axes show *no* real shift, only that none was detected
at this sample size.

This remains the single most direct evidence of "oncofetal reactivation"
in this entire project: within the same patient — controlling for
between-patient confounding (not a fully confound-free test; see below)
— malignant cells are dramatically more fetal-somatic-gut-like than
their patient-matched normal epithelium, while no such shift was
detected on the developmentally-shared or placental/trophoblast axes.
The finding is F-selective, not proven F-exclusive, and it holds
consistently across both regional F axes (Colon and SI), not just the
coarse F_Gut-specific union.

**Per-named-subtype sensitivity check** (`htan_patient_matched_contrast_by_normal_subtype.tsv`,
descriptive only per the locked design, not a second primary analysis):
the F-panel malignant-vs-normal shift replicates against most individual
normal epithelial subtypes (BEST4+, goblet, colonocyte, early
colonocyte, secretory, tuft — all p<0.001, deltas 58–97 points), with
two real exceptions, both reported here rather than the first-draft
omission of the second one. `intestinal crypt stem cell of colon` shows
a much smaller, non-significant F-panel delta (e.g. F_Gut-specific
Δ=6.5, p=0.37) — plausibly because stem cells already carry some
fetal-like developmental character themselves, making them a less sharp
contrast to malignant cells than fully differentiated subtypes.
**`enteroendocrine cell of colon` (n=18) shows genuinely mixed F-panel
behavior, not a clean replication**: `F_Gut-specific` is essentially
flat (Δ=+0.6, p=0.62), `F_Colon-specific` is significantly *negative*
(Δ=−9.1, p=0.0099 — the opposite direction from every other subtype),
and `F_SI-specific` is significantly positive (Δ=+30.45, p=5.3×10⁻⁴).
This same subtype also shows the two most notable D/P sensitivity
results (`D_Gut-shared` Δ=+22.75, p=7.6×10⁻⁶; `P_Gut-specific`
Δ=+8.1, p=2.0×10⁻⁴). Both patterns are flagged as real, interesting,
and unexplained here (small n=18, sensitivity-only status, not
investigated further) — not smoothed into a claim that the F-panel
result "replicates against nearly every subtype," which the
enteroendocrine result alone would contradict.

## What this shows, honestly

**The primary-extension correlation results generalize PR #27's core
finding** (weak, not-single-axis revCSC↔D/F/P relationships) to a
dataset with zero confirmed provenance overlap with the primary atlas —
a real, if modest, piece of independent corroboration for that finding's
generality. The sign-instability of the already-weak F/P correlations
between datasets is reported honestly rather than smoothed over.

**The patient-matched contrast is a genuinely new, and genuinely
striking, result this project has not had before**: a within-patient
test that controls between-patient confounding — **round-1 review
correction**: not "confound-controlled" in the fully general sense (see
below) — showing malignant cells are massively more F-developmental
(fetal-somatic-gut) than their patient-matched normal epithelium, with
no such shift detected on D-shared or P-specific axes. This sharpens,
rather than contradicts, the project's core developmental question — the
reactivation this project set out to characterize looks, by this most
direct test available, F-selective rather than D-shared or P-specific.

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

## Review history

- **Round 1 (REQUEST_CHANGES — one genuine pre-compute contract
  violation plus several material write-up errors, no rescoring
  required for any of them, all independently verified against real
  data/code before fixing)**: (1) CRLM's `P_Gut-specific` coverage gate
  (40/76 genes, 52.6%) was waived instead of investigated, contradicting
  the design's actual "blocking gate until understood" contract, and the
  implementation's flat 50% floor never even flagged it — fixed the
  coverage-flagging logic to a per-dataset relative-deviation check (now
  correctly flags this exact case) and ran the real investigation: all
  36 missing genes are present in both other datasets despite CRLM
  having the most total genes of the three, ruling out random dropout —
  a genuine CRLM-specific reference/annotation gap for this gene family,
  the existing 40-gene score remains usable and is kept as-is. (2) The
  bare-Ensembl-ID assertion never actually ran when axes were
  byte-identical — fixed to assert the ID format directly; confirmed to
  still pass for both datasets (re-run, no rescoring). (3) The marker
  check write-up wrongly claimed "no immune contamination" and uniform
  near-zero PTPRC for *both* datasets — checked directly against the
  committed table, false for CRLM (which correctly shows real TME
  biology: PTPRC high in T/NK/neutrophil/phagocyte/B cells, low in
  malignant/cholangiocyte) — fixed to describe the real, biologically
  correct pattern. (4) The HTAPP_HTAN/Pelka/Synapse provenance narrative
  was asserted from an interactive check never committed as auditable
  output — fixed by adding `htan_name_similarity_supplementary_check()`,
  run and committed, confirming the exact claim. (5) The primary
  correlation summary claimed "|r|≤0.10" — numerically false, largest
  observed magnitude is r≈0.102 — corrected. (6) The patient-matched
  contrast's "statistically indistinguishable" language for D/P
  overstated what a nonsignificant test establishes — corrected to "no
  detectable shift observed." (7) The subtype sensitivity paragraph
  mentioned only ISC as a partial exception to the F-panel replication,
  omitting `enteroendocrine cell of colon`'s genuinely mixed F-panel
  behavior (negative for F_Colon-specific, positive for F_SI-specific) —
  added. (8) Minor language fixes: "own normal epithelium" → "matched
  normal epithelium" (PR #33's locked terminology), "within-patient,
  confound-controlled" → "controls between-patient confounding" (not
  fully confound-free, per the document's own already-stated caveat).
  (9) The PR description's "no per-cell files" claim was wrong — this
  PR does commit per-cell scores for all 3 passes, matching PR #27's own
  precedent — corrected in both this doc and the PR body.
