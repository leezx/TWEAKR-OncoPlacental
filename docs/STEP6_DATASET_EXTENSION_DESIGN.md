# Step 6 dataset extension: `HTAN_CRC_progressive_plasticity` + `CRLM_NMP_ATLAS`

Executes the last remaining item from the original 3-dataset plan
(`docs/STEP6_CRC_PROJECTION_DESIGN.md`, "Proposed dataset plan"): the
primary/secondary/tertiary analysis structure has so far run only on
`CRC_single_cell_atlas_2025` (Steps 6 primary/secondary/tertiary, PRs
#27/#29/#31, all merged, all real results). This design locks the
compute contract for extending the same projection to the two smaller
datasets already inventoried in that same design doc but never scored:
`HTAN_CRC_progressive_plasticity` and `CRLM_NMP_ATLAS`. No new dataset
discovery — both were already found, inventoried, and structurally
characterized (job 3620558, `results/06_crc_projection/inventory/`) in
the original design round; this doc only locks how they get scored.

## What's already locked and directly reused (not re-decided here)

Per this project's standing discipline (reuse already-reviewed machinery,
don't reinvent), every one of the following is carried over unchanged
from Step 6 primary/secondary — cited, not re-litigated:

- **Frozen gene panels**: all 13 panels (8 revCSC variants + `D_Gut-shared`/
  `F_Gut-specific`/`F_Colon-specific`/`F_SI-specific`/`P_Gut-specific`),
  already stored as bare-Ensembl-ID `.ensembl.txt` files
  (`results/06_crc_projection/gut_scoring/`). **No new gene-set
  construction or gene-ID mapping step is needed** — the panels were
  already built directly against `CRC_single_cell_atlas_2025`'s Ensembl
  `var_names` (Step 6's Ensembl-ID-primary switch, PR #25), and both
  `HTAN_CRC_progressive_plasticity` and `CRLM_NMP_ATLAS` also use
  Ensembl gene IDs (confirmed at inventory time: 25,344 / 30,257 genes
  respectively). This is a real simplification relative to the original
  design doc's Step 2 gene-ID-mapping concern, which was written before
  the Ensembl-ID-primary switch existed.
- **Null-calibration scoring method**: `score_genes_fast` +
  expression-detectability-matched null draws, empirical percentile as
  primary common-scale metric, z-score secondary
  (`crc_gut_scoring_core.py`, PR #26/#27, unchanged).
- **revCSC↔D/F/P comparison-pair contract**: the same 10 pairs
  (`COMPARISON_PAIRS` in `crc_gut_scoring_core.py`), same overlap-exclusion
  logic (`CLU`/`ASS1`).
- **Donor/study-aware aggregation discipline**: per-patient pseudobulk
  summaries before any cross-patient claim, matching the standing rule
  used in every Step 6 analysis so far.

## What's genuinely new and needs a decision here

### 1. Raw-counts location differs from the primary atlas — a real technical difference, not assumed away

`crc_gut_scoring_core.py`'s `load_atlas()` rebuilds `adata.X` from
`layers['counts']` — that is `CRC_single_cell_atlas_2025`'s specific,
confirmed-raw layer. Per the original inventory (job 3620558), the two
new datasets instead carry confirmed raw counts in **`.raw.X`**, with
`adata.X` itself already normalized. Reusing `load_atlas()` unmodified
would silently re-normalize already-normalized data — a real bug, not a
style choice.

**Round-1 review correction**: the first draft's fix (copy `raw.X` into
`X`, normalize, done) was incomplete — `crc_gut_scoring_core.py`'s
`compute_detectability()` independently reads `adata.layers["counts"]`
to build the expression-detectability-matched null strata that the
*entire* null-calibration scoring machinery depends on. A loader that
only touches `X` leaves `layers["counts"]` empty and silently breaks
that machinery, not just the normalization step. **Corrected fix**: the
new loader (`load_extension_dataset(path)`) must standardize each
extension dataset into the *same internal contract* the primary atlas
loader produces, not just visually similar output — specifically: (a)
read `adata.raw.X` together with `adata.raw.var_names`/`adata.raw.var`
(the actual axis the count matrix is indexed against, not assumed to
align positionally with `adata.var_names` — dimension-matching alone
does not prove ordering or gene identity, must be checked explicitly);
(b) assert the raw feature axis is what's expected before proceeding;
(c) write that raw count matrix into `adata.layers["counts"]` (matching
the primary atlas's contract exactly, so `compute_detectability()` and
every downstream function work unmodified); (d) only then rebuild `X`
via the identical `normalize_total(target_sum=1e4)` + `log1p` recipe.
Canonical-marker verification (elevated epithelial/malignant markers in
the annotated compartment) remains a required *biological* sentinel
check, but is not a substitute for the structural axis-alignment
assertions in (a)/(b) — a loader could pass a marker sanity check while
still having a subtle gene-axis misalignment elsewhere in the panel.

### 2. Gene-ID version-suffix check (verify, don't assume)

The inventory table records both datasets' gene IDs as "Ensembl" but
does not record whether `var_names` carry version suffixes (e.g.
`ENSG00000141510.4` vs. bare `ENSG00000141510`) — the frozen panels are
bare-ID. **Round-1 review correction**: the check must run against the
feature axis that actually indexes the count matrix used for scoring —
per §1's corrected loader contract, that is `adata.raw.var_names`, not
necessarily `adata.var_names` (these are not guaranteed identical; if
they differ, checking the wrong one silently validates the wrong gene
axis). **Required first compute step, before any scoring**: confirm
`raw.var_names` format directly (not assumed from the inventory doc's
shorthand) — either assert byte-exact equality between `raw.var_names`
and `var_names` (simplest, if true), or, if versioned, strip suffixes
with the same logic already used elsewhere in this project *and assert
the stripped IDs contain no duplicate collisions* (two versioned IDs
stripping to the same bare Ensembl ID would silently corrupt the
gene-to-column mapping). Report per-panel `n_testable = |panel ∩
raw.var_names|` for all 13 panels on both datasets before proceeding —
if coverage is unexpectedly low for any panel, that blocks compute until
understood, not silently accepted.

### 3. HTAN's malignant-vs-matched-normal-epithelial contrast — a genuinely new analysis this dataset uniquely enables

Unlike `CRC_single_cell_atlas_2025` (malignant cells only,
`atlas_cell_type_middle`), `HTAN_CRC_progressive_plasticity` carries
**named, specific normal epithelial subtypes** (early colonocyte,
secretory, ISC, goblet, BEST4+, tuft, enteroendocrine) from the *same
patients* as the malignant cells, per the original inventory. This is
the direct within-patient malignant-vs-own-normal-epithelial contrast
flagged as an open item in the original design doc (item 4, never
executed).

**Round-1 review correction**: the first draft's "compare each patient's
malignant cells' percentile distribution against that same patient's
normal epithelial cells, using a paired test" was ambiguous about the
statistical unit — read literally, it could be (mis)implemented as a
cell-level paired test across thousands of cells per patient, which
would be pseudoreplication (treating within-patient cells as independent
observations). **Corrected, locked contract**: (a) score the *entire*
HTAN epithelial dataset (malignant + normal, jointly) using one common
null-calibration construction per panel — the null strata depend on the
population being scored, so malignant and normal cells must **not** be
independently recalibrated against separate null distributions and then
have their resulting percentile scales treated as directly comparable;
(b) for each of the 5 D/F/P panels plus the primary revCSC panel, and
for every patient with both malignant and normal-epithelial cells
present, collapse to **one malignant summary value and one
normal-epithelial summary value per patient per panel** (e.g. median
percentile); (c) the paired (within-patient) test runs across these
per-patient summary pairs — **patients, not cells, are the independent
statistical unit**. All-subtypes-pooled is the primary normal-epithelial
comparator; per-named-subtype summaries (early colonocyte, secretory,
ISC, goblet, BEST4+, tuft, enteroendocrine) are descriptive/sensitivity
only, reported where per-subtype n permits, not a second primary
analysis. **Also corrected**: "matched-normal"/"own-normal" language is
weakened to **patient-matched normal epithelium** — HTAN's annotation
includes primary tumor, metastasis, and non-tumor specimens per patient,
so same-patient matching removes patient-level confounding but does not
by itself rule out anatomical-site, treatment, or epithelial-subtype
composition effects between a patient's malignant and normal samples.
This remains the most direct oncofetal-reactivation test in this
project, but it is a strong test, not a confound-free one — stated
honestly rather than oversold. It only exists for `HTAN_CRC_
progressive_plasticity` — `CRC_single_cell_atlas_2025` has no matched
normal population, and `CRLM_NMP_ATLAS` is TME-focused with too few
malignant cells for a similar contrast (see below).

### 4. Study-provenance overlap audit for HTAN (required before "replication" language — flagged since the original design, never run)

The original design doc explicitly requires this before calling HTAN
"independent replication" of the meta-atlas: cross-reference HTAN's
patient/sample identifiers against `CRC_single_cell_atlas_2025`'s
per-cell `dataset`/`study_id`/`NCBI_BioProject_accession`/`SRA_sample_
accession` obs columns. **Locked here as a required pre-writeup step,
not optional**: run the audit before either dataset's results are
described as "replication" anywhere in this project's docs. **Round-1
review correction (strengthened, not just restated)**: this is a
**claim gate, not a compute gate** — the scores themselves remain
useful and get computed/reported regardless of the audit's outcome; only
the *language* used to describe them is gated. Until the audit runs, the
write-up avoids not just "replication" but also "external validation"
(both imply the audit already ran) and instead calls this the "**HTAN
dataset-extension analysis**" throughout — if the audit finds no
overlap, the stronger "external validation"/"replication" language may
be promoted afterward, not assumed up front.

### 5. `CRLM_NMP_ATLAS` is explicitly scoped as exploratory only

Per the original inventory: only 4,051 of 75,104 cells are malignant
(the rest is TME — T cell, neutrophil, NK, macrophage), across 6 donors.
**Locked scope**: run the same revCSC↔D/F/P null-calibrated scoring on
its 4,051 malignant cells — **round-1 review correction**: the first
draft said "correlation + composition, same construction as the
primary/secondary analyses" here, directly contradicting §6's "no
default composition re-run" — §6 is the correct, intended contract.
Corrected: **primary continuous correlation only** (the same
primary-analysis-equivalent construction used for HTAN), explicitly
**not** the secondary revCSC-high composition or tertiary full-atlas
composition machinery. Every result is reported explicitly labeled
exploratory/small-n — no donor-cluster bootstrap CI, no "replication" or
"confirms" language, matching this project's standing fallback-naming
discipline ("if too few cells clear a sensible threshold, rename the
deliverable rather than claim resolved," carried over from the original
design's explicit fallback instruction). **Round-1 review correction on
the rationale itself**: the first draft framed this as driven by "only
4,051 malignant cells" — that overstates the weakness of the cell count,
which is actually adequate for estimating cell-level score
distributions. The real limiting factor is **6 donors** — the effective
sample size for any donor-level claim — compounded by the
liver-metastasis/TME-focused sampling and `Before/After_NMP` structure,
exactly as the original inventory itself already describes. Corrected
statement: CRLM is exploratory because it has only 6 biological
replicates in a specialized metastatic cohort, not because it has "too
few" malignant cells. **Round-1 review correction, leave-one-donor-out**:
the design previously said both "no leave-one-donor-out robustness
claim" (here) and that LODO correlation would still be reported (compute
plan) — not actually contradictory, but under-specified. Clarified: LODO
is computed and shown as a **descriptive sensitivity diagnostic only** —
it does not generate a robust/non-robust classification the way it does
for the two much-larger datasets. The `Before/After_NMP` timepoint
structure is recorded but not analyzed as a covariate here — a genuine
gap, flagged for future scrutiny only if the small-n result turns out
interesting enough to warrant it.

### 6. No new secondary/tertiary composition analysis by default

Given both datasets are far smaller than `CRC_single_cell_atlas_2025`
(47,107 and 75,104 cells vs. 665,473), this design scopes compute to a
**primary-analysis-equivalent** (continuous revCSC↔D/F/P correlation,
donor/study-aware) plus HTAN's unique malignant-vs-normal contrast
(§3) plus CRLM's exploratory version (§5) — not a full re-run of the
secondary (revCSC-high composition) or tertiary (full-atlas,
revCSC-independent) analysis structure on each dataset. Rationale: the
secondary/tertiary analyses' main scientific payoff (bimodal composition
generalizing population-wide, SI:Colon skew) was already established at
full statistical power on the 665,473-cell primary atlas; re-running the
identical composition machinery on two much smaller populations mainly
adds replication value, not new structure, and can be scoped as a
follow-up if the primary-equivalent correlation here surfaces something
that specifically warrants it — avoiding scope creep for its own sake,
per this project's standing efficiency discipline.

## Compute plan

1. **Coverage check** (§2): confirm gene-ID format, report per-panel
   `n_testable` for both datasets, all 13 panels — blocking gate before
   any scoring.
2. **Loader verification** (§1): confirm `.raw.X`-based normalization
   reproduces expected canonical-marker behavior (e.g. epithelial/
   malignant markers elevated in the annotated malignant/epithelial
   compartment) before trusting scored output.
3. **Study-provenance overlap audit** (§4): run before any HTAN write-up
   uses "replication" or "external validation" language.
4. **Primary-analysis-equivalent scoring**: null-calibrated percentile +
   z-score, all 13 panels, both datasets. **Round-1 review correction**:
   the first draft claimed `N_PERM=100` was already certified for "this
   panel set" and needed no new convergence probe — false on two counts,
   confirmed directly against `docs/STEP6_GUT_SCORING_COMPUTE_RESULTS.md`:
   (a) PR #27's own convergence gate did **not** certify `N_PERM=100` for
   every panel — it flagged 2/10 comparison pairs and forced 4 panels
   (`F_Colon-specific`, `F_Gut-specific`, `revCSC_extended28_minus_CLU`,
   `revCSC_primary27_minus_CLU_ASS1`) to `N_PERM=500`; (b) that
   convergence result was certified on the 665,473-cell primary atlas —
   since the null strata are built from expression-detectability in the
   population actually being scored, the same certification does not
   transparently transfer to HTAN's or CRLM's very different, much
   smaller populations. **Corrected: use `N_PERM=500` for all 13 panels
   on both extension datasets** — simpler than re-running a new
   100-vs-500 convergence gate per dataset, and the original motivation
   for `N_PERM=100` (compute cost at 665,473 cells) is far weaker at
   these datasets' scale (47,107 / 75,104 cells). Donor/patient-aware
   pooled + leave-one-donor-out correlation (HTAN: 29 patients, real
   robustness checks feasible; CRLM: 6 donors, reported as a descriptive
   sensitivity diagnostic only per §5, not a robustness classification).
5. **HTAN malignant-vs-normal contrast** (§3): joint null-calibrated
   scoring of the full epithelial dataset, then one malignant summary and
   one normal-epithelial summary per patient per panel, paired
   (within-patient) test across those per-patient summaries — patients,
   not cells, are the statistical unit (§3's corrected contract).
6. **Write-up**: two dataset-specific results sections (mirroring
   `STEP6_GUT_SCORING_COMPUTE_RESULTS.md`'s structure), each explicitly
   stating its own robustness/coverage caveats — not a single combined
   claim across datasets of very different size and population
   structure.

## What this design does not do

Does not re-derive or re-score any frozen D/F/P or revCSC gene set
(reused as-is). Does not re-run the secondary/tertiary composition
analyses on either dataset by default (§6, follow-up only if warranted).
Does not extend to `GSE178318` (no cell-type annotation exists — already
flagged as deferred in the original design doc, unchanged). Does not
touch Q5/Q6 of the original 6-question framework (explicitly out of
scope per the 2026-08-16 user confirmation, `docs/Q1_Q2_Q4_CROSS_REFERENCE.md`).

## Review history

- **Round 1 (REQUEST_CHANGES — 4 real technical/methodological issues,
  all independently verified against the committed source docs before
  fixing, no expansion of scientific scope requested)**: (1) the `.raw.X`
  loader fix was incomplete — it didn't account for
  `compute_detectability()`'s independent dependency on
  `layers["counts"]`, and didn't assert raw-feature-axis alignment
  explicitly rather than assuming positional correspondence; fixed to
  standardize extension datasets into the primary atlas's exact internal
  contract (§1). (2) The gene-ID version-suffix check needed to run
  against `raw.var_names` (the actual count-matrix feature axis per §1's
  fix), not necessarily `var_names`, plus a duplicate-collision assertion
  after any suffix-stripping (§2). (3) HTAN's malignant-vs-normal
  contrast description was ambiguous about the statistical unit — fixed
  to explicitly lock patient-level paired inference (one summary per
  patient per group) and joint null-calibration (not separately
  recalibrated malignant/normal scales), plus softened
  "matched-normal"/"own-normal" to "patient-matched" given HTAN's
  primary/metastasis/non-tumor specimen mix (§3). (4) A real
  factual error — the design claimed `N_PERM=100` was already certified
  for "this panel set," but PR #27's own convergence gate forced 4 of the
  13 panels to `N_PERM=500`, and that certification was population-
  specific to the primary atlas's detectability strata, not
  automatically transferable — confirmed directly against
  `docs/STEP6_GUT_SCORING_COMPUTE_RESULTS.md`'s convergence table; fixed
  to `N_PERM=500` for all panels on both extension datasets (compute
  plan item 4). Also fixed a real internal contradiction (§5 said CRLM
  gets "composition, same as primary/secondary," §6 said no default
  composition re-run — §6 is correct, §5 corrected to match) and
  clarified CRLM's exploratory rationale (6 donors is the real limiting
  factor, not "only 4,051 malignant cells") and LODO's descriptive-only
  status for CRLM specifically.

Submitting for round-2 review before compute, same discipline as every
prior step this session.
