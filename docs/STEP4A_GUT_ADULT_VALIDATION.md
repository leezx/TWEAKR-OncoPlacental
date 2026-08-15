# Step 4a design: external adult-negative validation of the gut D/F/P gene sets

Per `docs/PROJECT_SUMMARY.md`'s "What's next": the frozen gut-specific gene
sets (`results/04a_dfp_gut/dfp_gut_gene_sets/`) — `F_Colon-developmental`,
`F_SI-developmental`, `D_Colon-shared`, `D_SI-shared`, `F_Gut-core`,
`P_Colon-specific`, `P_SI-specific` — were all built from **one** data
source: same-atlas fetal-vs-adult DE within the Gut Cell Atlas
(Elmentaite et al., *Nature* 2021). This design tests them against
**external** adult references that never touched construction, reusing
the two data sources and the exact methodology already reviewed and
approved for the pan-organ D/F/P's Step 5 Tier-2 validation
(`docs/STEP5_TIER2_VALIDATION.md`, PR #14/#15) rather than inventing a new
one. Design submitted before any qsub compute, same discipline as every
prior step.

## Why this is a genuinely independent check (stronger than Step 5's)

Step 5's original Tier-2 validation used GTEx as **Tier 1** (mandatory,
used to *define* the pan-organ signature via `adult_excluded`) and
Tabula Sapiens as **Tier 2** (held out specifically for independent
validation). For the gut gene sets, neither GTEx nor Tabula Sapiens nor
HPA was used anywhere in construction — `F_Colon-developmental`/
`F_SI-developmental` came entirely from the Gut Cell Atlas's own internal
fetal-vs-adult contrast (`~10X + source_family + Age_group`, PR #21).
**Both data sources below are fully external for the gut sets** — a
cleaner independence property than the original pan-organ check had.

## Data sources (both already on disk, verified in prior steps — not re-downloaded)

- **Tabula Sapiens `Large_Intestine`** (`datasets/TabulaSapiens/`, also on
  Argos): 13,680 cells, 2 donors (`TSP14`, `TSP2`), raw counts in
  `layers['raw_counts']` (already confirmed 100% integer-valued, Step 5
  inventory). `cell_ontology_class` includes real adult gut epithelial
  types matching the gut D/F/P's own epithelial-lineage-only construction
  population: `enterocyte of epithelium of large intestine`, `immature
  enterocyte`, `mature enterocyte`, `paneth cell of epithelium of large
  intestine`, `large intestine goblet cell` (plus non-epithelial types —
  T/B cells, fibroblast, monocyte — excluded from the epithelial-matched
  check below). **Real, stated gap**: Tabula Sapiens has no `Small_Intestine`
  organ file (confirmed, Step 5 inventory covers only Liver/Skin/Spleen/
  Thymus/Large_Intestine) — `F_SI-developmental` gets no single-cell
  organ-matched check from this source, only the GTEx bulk check below.
- **GTEx v11 median TPM** (`datasets/GTEx_v11_median_tpm/`): confirmed
  columns include `Colon_Sigmoid`, `Colon_Transverse` (+Mucosa/Muscularis/
  Mixed_Cell) and `Small_Intestine_Terminal_Ileum` (+Lymphoid_Aggregate/
  Mixed_Cell) — real organ-matched bulk tissue for **both** Colon and SI,
  filling the Tabula Sapiens SI gap at bulk (not single-cell) resolution.

## Method — reusing two already-approved contracts, not inventing new ones

### Check 1: Tabula Sapiens organ-matched, epithelial-lineage-restricted (new refinement)

Same donor-aware pseudobulk contract as `tier2_validation.py` (PR #14/#15,
2 review rounds, approved): primary unit is `(donor, cell_ontology_class)`
pseudobulk from raw counts, CPM computed per pseudobulk library's own
size, minimum 20-cell floor applied at the donor×cell-type level before
aggregation, cross-donor-consistent detection (CPM≥1 in **all** eligible
donors, ≥2 eligible donors) as the primary flag, single-donor hits
reported but not over-interpreted, 500-permutation expression-
detectability-matched null with global BH-FDR across all tests (exact
method from PR #15 round 1/2 fixes — reused, not redesigned).

**One refinement over the original**: restrict to `compartment ==
"epithelial"` cell types only (to be confirmed against the real column
values before running — `enterocyte`/`goblet cell`/`paneth cell` variants
above). The original whole-body check tested D-shared/P-specific across
*all* compartments in all 5 organs deliberately (those sets have no
lineage restriction in their own definition). `F_Colon-developmental`,
by contrast, was built exclusively from epithelial cells (PR #20's locked
population) — testing it against, say, Large_Intestine T cells would be
checking the wrong adult population. Only 2 donors here (`TSP14`,
`TSP2`), both must be eligible for a cell type to earn cross-donor-
consistent status — small but real, same "single-donor evidence" caveat
labeling as the original Thymus/Liver-etc. donor counts (2-3) in Step 5.

Target gene set: `F_Colon-developmental` (1,456 genes) against
Large_Intestine epithelial cell types. `F_SI-developmental` is not
tested here (no matching organ) — covered by Check 2 only, an explicit,
stated asymmetry rather than a silently-skipped gap.

### Check 2: GTEx bulk adult-exclusion, organ-matched (reuses `adult_excluded_percentile_audit.py`'s method exactly)

Same not-detected-floor + percentile-among-detected design already
built and approved for the pan-organ D/F/P (`docs/STEP4_STATISTICAL_DESIGN.md`
section 3): `not_detected(gene,tissue) = TPM<1`; for detected genes,
percentile rank computed only within the detected subset (avoids the
same zero-tie-block degeneracy already found and fixed there). Applied
here purely as **post-hoc validation reporting** (percentile distributions
at multiple candidate cutoffs, not a new frozen threshold gating anything)
— unlike the original use, nothing about the gut D/F/P gene sets changes
based on this result; it is reported honestly regardless of outcome.

- `F_Colon-developmental` vs. `Colon_Sigmoid` + `Colon_Transverse` (both
  columns, reported separately — no combination rule needed for reporting).
- `F_SI-developmental` vs. `Small_Intestine_Terminal_Ileum`.
- `F_Gut-core`, `D_Colon-shared`, `D_SI-shared` vs. all 68 GTEx tissues
  (whole-body-style, matching P-specific's original whole-body design —
  these are cross-region/summary sets, not organ-scoped by their own
  definition).

### Check 3: Tabula Sapiens whole-body-style, all 5 organs (unchanged from Step 5's approved design)

`D_Colon-shared` (5 genes), `D_SI-shared` (4 genes), `F_Gut-core` (712
genes), `P_Colon-specific`, `P_SI-specific` — tested exactly as Step 5
tested the pan-organ D-shared/P-specific: gene×cell-type×donor evidence
across all 5 Tabula Sapiens organs combined (Liver, Skin, Spleen, Thymus,
Large_Intestine), no lineage restriction (matches these sets' own
whole-body/cross-region definition), same permutation null + BH-FDR.

## What this does NOT do

- Does not re-open or re-run any part of the gut D/F/P's own construction
  (PR #19/#20/#21, frozen) or the mike_verzi external validation (PR #22,
  frozen) — this is a strictly additive, independent check layered on top.
- Does not propose a new cutoff or gate any gene's membership in the
  frozen sets — purely reports whether the already-frozen genes show
  adult single-cell/bulk signal externally, same "report honestly, don't
  optimize the story" posture as Step 5.
- `F_Colon-specific`/`F_SI-specific` (the F-minus-D/P-overlap partition
  sets, distinct from `F_Colon-developmental`/`F_SI-developmental`
  themselves) are not separately tested — `F_Colon-developmental`/
  `F_SI-developmental` are supersets of them and the more direct
  construction-level object; testing the partition sets too would be
  redundant with Check 1/2 without adding new information, unless review
  disagrees.

## What's next after this design is approved

Real compute: `scripts/04a_dfp_gut/gut_adult_validation_tabula_sapiens.py`
(Checks 1+3, adapting `tier2_validation.py`) and
`scripts/04a_dfp_gut/gut_adult_validation_gtex.py` (Check 2, adapting
`adult_excluded_percentile_audit.py`), run via qsub, results pulled back
byte-exact, written up honestly (including any null or contradictory
findings, same standing discipline as every step so far), submitted for
review, merged only after user confirmation.
