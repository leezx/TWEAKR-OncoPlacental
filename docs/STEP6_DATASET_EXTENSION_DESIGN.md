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
style choice. **Fix**: a new, dataset-parameterized loader
(`load_extension_dataset(path, raw_counts_source="raw.X")`) that reads
counts from `adata.raw.X` for these two datasets, rebuilds `X` via the
identical `normalize_total(target_sum=1e4)` + `log1p` recipe, and is
verified against a handful of canonical epithelial/malignant markers
before being trusted (same "verify the mechanism" discipline as every
prior Step 6 loader change).

### 2. Gene-ID version-suffix check (verify, don't assume)

The inventory table records both datasets' gene IDs as "Ensembl" but
does not record whether `var_names` carry version suffixes (e.g.
`ENSG00000141510.4` vs. bare `ENSG00000141510`) — the frozen panels are
bare-ID. **Required first compute step, before any scoring**: confirm
`var_names` format directly (not assumed from the inventory doc's
shorthand); if versioned, strip suffixes with the same logic already
used elsewhere in this project, and report per-panel
`n_testable = |panel ∩ dataset var_names|` for all 13 panels on both
datasets before proceeding — if coverage is unexpectedly low for any
panel, that blocks compute until understood, not silently accepted.

### 3. HTAN's malignant-vs-matched-normal-epithelial contrast — a genuinely new analysis this dataset uniquely enables

Unlike `CRC_single_cell_atlas_2025` (malignant cells only,
`atlas_cell_type_middle`), `HTAN_CRC_progressive_plasticity` carries
**named, specific normal epithelial subtypes** (early colonocyte,
secretory, ISC, goblet, BEST4+, tuft, enteroendocrine) from the *same
patients* as the malignant cells, per the original inventory. This is
the direct within-patient malignant-vs-own-normal-epithelial contrast
flagged as an open item in the original design doc (item 4, never
executed). **New analysis, scoped here**: for each of the 5 D/F/P panels
plus the primary revCSC panel, compare each patient's malignant cells'
null-calibrated percentile distribution against that same patient's
normal epithelial cells (all subtypes pooled, and separately by named
subtype where n permits), using a paired (within-patient) test — not a
between-patient comparison, avoiding the exact confound flagged in the
original design's item 4. This is the most direct test in this project
of "oncofetal reactivation," and it only exists for `HTAN_CRC_
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
described as "replication" anywhere in this project's docs. Until it
runs, HTAN results are reported as "technical/external validation," per
the original design's already-adopted fallback label — this design does
not relax that requirement.

### 5. `CRLM_NMP_ATLAS` is explicitly scoped as exploratory only

Per the original inventory: only 4,051 of 75,104 cells are malignant
(the rest is TME — T cell, neutrophil, NK, macrophage), across 6 donors.
**Locked scope**: run the same revCSC↔D/F/P null-calibrated scoring on
its 4,051 malignant cells (correlation + composition, same construction
as the primary/secondary analyses), but report every result explicitly
labeled exploratory/small-n — no donor-cluster bootstrap CI, no
leave-one-donor-out robustness claim, no "replication" or "confirms"
language, matching this project's standing fallback-naming discipline
("if too few cells clear a sensible threshold, rename the deliverable
rather than claim resolved," carried over from the original design's
explicit fallback instruction). The `Before/After_NMP` timepoint
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
   uses "replication" language.
4. **Primary-analysis-equivalent scoring**: null-calibrated percentile +
   z-score, all 13 panels, both datasets, `N_PERM=100` (same convergence
   gate already certified for this panel set — no new convergence probe
   needed, same gene-set sizes as the primary atlas run). Donor/patient-
   aware pooled + leave-one-donor-out correlation (HTAN: 29 patients,
   real robustness checks feasible; CRLM: 6 donors, reported but
   explicitly flagged as underpowered per §5).
5. **HTAN malignant-vs-normal contrast** (§3): paired within-patient
   test, all patients with both malignant and normal-epithelial cells
   present.
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
