# Step 8: D/F/P/revCSC scoring for the Step 7 scRNA-seq cohorts — design

Follow-on to Step 7 (`docs/STEP7_CLIM_DATA_ACQUISITION_DESIGN.md`, PR #35;
`docs/STEP7_INVENTORY_RESULTS.md`, PR #36 — both merged), which
downloaded and structurally characterized all public CLiM/CLuM cohorts
but explicitly deferred any D/F/P/revCSC scoring to a separate PR. This
is that PR — design only, mirroring this project's own Step 6
dataset-extension precedent (`docs/STEP6_DATASET_EXTENSION_DESIGN.md`,
PR #33): locks the loader/coverage-check contract and pre-compute
investigation before any scoring compute runs.

## Scope (deliberately narrowed to keep this PR checkable in one pass)

**In scope**: scoring the 3 Step-7 **scRNA-seq** cohorts
(`GSE231559`, `GSE285990`, `GSE225857`) against the existing, frozen
13-panel D/F/P/revCSC gene-set inventory
(`scripts/06_crc_projection/crc_gut_scoring_core.py`), reusing that
module's scoring machinery unchanged — same null-calibration method,
same `N_PERM=500`, same empirical-percentile + z-score outputs, per this
project's Step 6 dataset-extension precedent.

**Explicitly out of scope, deferred to a separate future phase**:
scoring the 4 Step-7 **bulk** cohorts (`GSE131418`, `GSE17536`+
`GSE17537`, `GSE21510`, `TCGA-CRC`). Bulk scoring needs infrastructure
this project has not built: the 3 microarray cohorts are raw Affymetrix
CEL probe-level intensities requiring RMA normalization and probe→gene
summarization before any gene-level signature scoring is even possible,
and `scanpy.tl.score_genes`'s per-cell detectability-stratified null
calibration (this project's whole existing scoring method) has no bulk
equivalent — a bulk-appropriate single-sample scoring method (e.g.
ssGSEA-style rank scoring) would need its own design and review, not an
improvised adaptation folded into this PR. Not attempted here.

Also out of scope: any cross-cohort composition/contrast analysis
(mirrors patient groups, revCSC-high definitions, etc. from Step 6's
secondary/tertiary analyses) — this PR locks the scoring *compute*
contract only; composition analysis is a further follow-on once scoring
itself is reviewed and run, matching Step 6's own design→compute→
results 3-PR staging.

## Real pre-compute findings (this round — every claim checked directly
against the actual downloaded files, not assumed from the Step 7
inventory's summary-level characterization)

| Dataset | Loader contract | Real finding this round |
|---|---|---|
| **GSE231559** | 10x-style MTX triplets per GSM, bare/versioned Ensembl gene IDs, raw integer counts (both confirmed in Step 7 inventory) | **New finding**: the 26 samples split into two DIFFERENT gene-reference batches by sequencing batch — 18 `SC10`/`SC143`-prefixed samples use a 33,538-gene reference; 8 `SC173`/`SC216`-prefixed samples use a 36,601-gene reference. Confirmed directly from the committed `GSE231559_inventory.tsv`. This is a real cross-sample reference-annotation inconsistency within one GEO series, not assumed uniform. **Locked contract**: load each sample independently (native gene axis), then intersect to the common Ensembl-ID gene set across both references before any pooled/cross-sample step; per-panel `n_testable` genes must be computed against this intersected background, not either reference alone. Intersection computed this round directly from both references' real `features.tsv`: 33,538 ∩ 36,601 = **32,732 common genes** (bare Ensembl IDs). |
| **GSE285990** | Same as above | All 10 `P01_LM`-`P10_LM` samples confirmed uniform 37,487-gene reference in the Step 7 inventory (re-verified here) — no cross-sample reference split, unlike GSE231559. |
| **GSE225857** | TXT count matrices (genes × cells) per GSM, gene **symbols** (not Ensembl — confirmed `A1BG` etc., `first_row_looks_ensembl=False` in the Step 7 inventory), 2 pooled cross-patient matrices (immune/non-immune) | **New finding**: both matrices are already reduced to a filtered gene set — 17,066 genes (immune) / 17,515 genes (non-immune) — not the full ~20-36K-gene reference every other dataset in this project uses. Confirmed directly by counting the full gene list in both files (not assumed from a truncated sample). Symbol→Ensembl mapping reuses this project's own existing `DATA/1.Databases/HGNC_gene_id_mapping/processed/v0.1/hgnc_symbol_ensembl_map.tsv` (45,033 rows, built in Step 2) — not re-derived. **Locked contract**: map symbols to Ensembl via this table, drop any symbol with no unambiguous Ensembl mapping (logged, not silently dropped), then compute per-panel coverage against this already-reduced ~17K-gene background specifically — a materially smaller background than GSE231559/GSE285990's, so `n_testable` will be lower here and that is expected, not a bug, PROVIDED it is checked and reported explicitly rather than assumed comparable. |

## Coverage check (required pre-compute gate, per the Step 6 precedent's
"blocks compute until understood" contract)

Locked as: for each of the 13 panels × 3 datasets (`GSE225857` counted
as 2 populations, immune and non-immune, since they're pre-filtered to
different gene sets), count `n_testable` genes (panel genes present in
that dataset's own background gene set, after the axis/mapping steps
above) and flag any panel/dataset pair whose coverage falls materially
below that dataset's own median coverage (the same relative-deviation
rule Step 6's dataset extension used — `>15 points below that dataset's
own median` — not a flat floor). Any flagged pair must be genuinely
investigated (same discipline as Step 6's real CRLM `P_Gut-specific`
investigation) before compute proceeds — not waived.

**Real coverage numbers, computed this round** (not deferred to compute
time — the actual gene backgrounds for all 3 datasets, including
GSE231559's reference-intersected background from the finding above,
were fetched and checked directly against all 13 panels' real committed
Ensembl-ID files):

| Panel | n_genes | GSE231559 | GSE285990 | GSE225857 (immune) | GSE225857 (non-immune) |
|---|---|---|---|---|---|
| `D_Gut-shared` | 8 | 8 | 8 | 6 | 7 |
| `F_Gut-specific` | 2,192 | 2,182 | 2,185 | 1,838 | 1,914 |
| `F_Colon-specific` | 1,451 | 1,442 | 1,447 | 1,201 | 1,267 |
| `F_SI-specific` | 1,452 | 1,447 | 1,446 | 1,241 | 1,270 |
| `P_Gut-specific` | 76 | 76 | 76 | **39** | **44** |
| `revCSC_primary27_full` | 27 | 27 | 27 | 23 | 25 |
| `revCSC_primary27_minus_CLU` | 26 | 26 | 26 | 22 | 24 |
| `revCSC_primary27_minus_ASS1` | 26 | 26 | 26 | 22 | 24 |
| `revCSC_primary27_minus_CLU_ASS1` | 25 | 25 | 25 | 21 | 23 |
| `revCSC_extended28_full` | 28 | 27 | 27 | 23 | 25 |
| `revCSC_extended28_minus_CLU` | 27 | 26 | 26 | 22 | 24 |
| `revCSC_extended28_minus_ASS1` | 27 | 26 | 26 | 22 | 24 |
| `revCSC_extended28_minus_CLU_ASS1` | 26 | 25 | 25 | 21 | 23 |

`GSE231559` (33,538/36,601-gene references, intersected to a 32,732-gene
common background) and `GSE285990` (37,487 genes) show near-complete
coverage (≥96%) on every panel — no deviation, no gate triggered.
`GSE225857`'s two populations (16,028/16,616 genes successfully
symbol→Ensembl-mapped, out of 17,066/17,515 total) show broadly lower
coverage across the board (~82-89% for most panels, consistent with its
smaller pre-filtered background — expected, not itself a flag), **but
`P_Gut-specific` drops much further, to 51% (immune) / 58% (non-immune)
— clearly >15 points below this dataset's own ~82-89% median for every
other panel.** This is flagged as a required pre-compute investigation
item, **not waived**, and is structurally the same shape of finding as
Step 6's already-investigated CRLM `P_Gut-specific` coverage gap
(`docs/STEP6_DATASET_EXTENSION_RESULTS.md`) — worth checking during
compute whether the same underlying cause (a `P_Gut-specific`-specific
reference/annotation gap, not a data-quality problem) recurs here, or
whether this is a genuinely different, dataset-specific cause (e.g. the
pre-filtering step behind GSE225857's reduced gene list happening to
disproportionately exclude this panel's genes). Locked as a required
investigation during the compute PR, following exactly the same
discipline as the CRLM precedent — not assumed to be the same cause
without checking.

## Raw-counts / normalization contract

`GSE231559`/`GSE285990`: raw-vs-processed status (integer sampled
nonzeros) already confirmed in the Step 7 inventory — reused directly,
not re-verified. `GSE225857`: the Step 7 inventory only characterized
shape/orientation, not the raw-vs-processed status of the values — this
round directly sampled nonzero values from both count matrices (~1,600-
2,400 nonzero values each, first 200 genes × first 49 cells) and
confirmed **100% integer-valued** (immune: values up to 139, non-immune:
up to 301) — genuine raw UMI/read counts, not normalized. Both GSE225857
matrices are therefore in scope for this design's locked scoring method,
same as GSE231559/GSE285990.

## What this design does not do

Does not run any scoring compute (design-only, per this project's
standing discipline of separate design→compute PRs). Does not score the
4 bulk cohorts (explicitly deferred, see Scope). Does not run any
cross-cohort composition/contrast analysis. Does not modify the
existing, frozen `crc_gut_scoring_core.py` machinery or the 13-panel
gene-set inventory — reused exactly as-is, matching the Step 6
dataset-extension precedent's explicit non-goal of re-deriving already-
frozen scoring machinery.

Submitting for review before any compute, same discipline as every
prior step this session.
