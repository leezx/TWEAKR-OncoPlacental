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

## Round-1 review fixes (all independently re-verified against real
committed code/live sources before fixing, none disputed)

**Blocker 1 — the scoring population was never locked.** This scoring
method is population-dependent (detectability bins and null-calibration
controls are recomputed from the exact cells being scored — Step 6's
own extension design had to separate malignant-only vs. malignant+
normal joint scoring for exactly this reason), and the first draft only
named datasets and gene backgrounds, never which cells within them get
scored. **Fixed, locked per dataset**:
- `GSE231559`: primary scoring population is the paper's cited 9 CLiM +
  6 primary-tumor samples (already reconstructed exactly in Step 7,
  `results/07_clim_external_data/GSE231559_inventory.tsv`) — all cells
  within those 15 tumor samples, since no per-cell malignant/epithelial
  annotation exists in the deposited data (only sample-level tumor/
  normal labels from GEO's own `Sample_title` field). This sample-level
  vs. cell-level granularity mismatch is stated explicitly, not silently
  equated with the primary atlas's per-cell malignant-cell restriction.
  The 11 paired-normal samples are loaded and scored separately (not
  pooled into the primary population), available for a future normal-
  vs-tumor contrast but not part of this design's primary result.
- `GSE285990`: all 10 `P01_LM`-`P10_LM` samples (confirmed human liver-
  metastasis tumor tissue by GEO's own sample titles, no separate
  normal counterpart in this cohort) — single population, no split
  needed.
- `GSE225857`: **only the non-immune (CD45-) fraction** (`GSM7058755`,
  41,892 cells) is in scope for primary D/F/P/revCSC scoring — this
  panel set is built for epithelial/malignant developmental-program
  signatures, and CD45- is the closest available proxy to epithelial/
  tumor-compartment cells this dataset offers (still not malignant-only
  — CD45- includes stromal/endothelial cells too, no further cell-type
  annotation is available, stated as an explicit limitation, not
  smoothed over). **The 196,473-cell immune (CD45+) fraction
  (`GSM7058754`) is explicitly OUT OF SCOPE for D/F/P/revCSC scoring in
  this design** — an immune-cell population is not a sensible target for
  a developmental-epithelial-program gene panel, and scoring it would
  produce numbers with no established interpretation. If a future need
  arises to characterize the immune compartment, that requires its own
  separate design with panels appropriate to immune biology, not this
  panel set reused by default.

**Blocker 2 — the loader contract didn't reproduce `crc_gut_scoring_
core.py`'s exact required internal representation.** Confirmed directly
by re-reading the frozen module: `compute_detectability()` reads
`adata.layers["counts"]` directly (not `adata.X`), and `load_atlas()`
rebuilds `X` as `normalize_total(target_sum=1e4)` + `log1p` of that same
layer. The first draft only discussed generic "raw-vs-processed status,"
never locking this precise axis/layer contract — the same category of
bug the PR #33 review caught in the Step 6 extension design (missed
`compute_detectability()`'s independent dependency on `layers["counts"]`
). **Fixed, locked explicitly**: every dataset's loader must produce
`layers["counts"]` holding raw integer counts on the canonical gene
axis (see Blocker 3), then rebuild `X` via the identical
`normalize_total(1e4)`+`log1p` sequence — and this must hold true after
any sample concatenation (GSE231559's 15-sample primary population,
GSE231559's 11-sample normal population), not just per-sample before
merging.

**Blocker 3 — gene-axis canonicalization was underspecified relative to
the coverage numbers already reported.** The first draft's coverage
table used version-suffix-stripped ("bare") Ensembl IDs without ever
locking that stripping step as part of the design's actual contract, and
never addressed what happens if stripping (GSE231559/GSE285990) or
symbol→Ensembl mapping (GSE225857) produces duplicate/colliding gene
IDs, nor whether the coverage-check computation and the eventual
scoring loader would necessarily share the same axis-construction code.
**Fixed, locked explicitly**: (1) strip Ensembl version suffixes
(`.split(".")[0]`) as an explicit, shared canonicalization step; (2)
assert no duplicate collisions result from stripping/mapping — raise
loudly if any gene ID maps to >1 row after canonicalization, don't
silently sum/average/pick-first; (3) the coverage-check script and the
compute-time scoring loader must call the identical shared
canonicalization function (not two independent reimplementations that
could silently drift apart) — the coverage numbers already reported
below are only meaningful if the actual scoring object ends up on the
same axis.

**Blocker 4 — GSE225857's reduced gene universe was framed as a neutral
smaller-reference-background, when it's actually the product of the
original authors' own systematic, gene-class-biased preprocessing.**
Confirmed directly by fetching the real authors' analysis code
(`github.com/jalon9358/LianLab_CRCLM/data_import_and_filter.R`, the same
repository independently verified for GSE225857's patient-cohort
provenance in Step 7 PR #36): Seurat objects are created with
`min.cells = 50` (an expression-frequency floor, not a reference-
annotation limit), and a regex removes noncoding-gene-name patterns
(antisense transcripts, `LINC`-prefixed, `MIR`-prefixed genes, etc.)
plus ribosomal (`RP[LS]`) and mitochondrial (`MT.`) genes from a
separate integration-feature list. This means gene absence in the
deposited ~17K-gene matrices is **expression/preprocessing-dependent,
not a neutral smaller-reference phenomenon** — the first draft's
"expected, not itself a flag" framing for the general ~82-89% coverage
undersold this. **Investigated directly this round, not left as a
deferred TODO**: checked which of `D_Gut-shared`'s 2 missing genes
(`LGALS14`, `KISS1`) match the noncoding/RP/MT removal patterns —
neither does (both are protein-coding, correctly-named genes), so their
absence is more consistent with the `min.cells=50` expression floor
(both `LGALS14`/galectin-14 and `KISS1`/kisspeptin have narrow,
tissue-restricted expression profiles biologically plausible for
falling below a 50-cell floor in a large aggregate object). **Fixed,
locked**: the compute-time investigation required for any coverage
deviation is extended beyond `P_Gut-specific` to explicitly include
`D_Gut-shared` and the small revCSC panels too (proportionally large
losses, e.g. 2/8 = 25% for `D_Gut-shared`, don't trigger the aggregate
"−15 points vs. dataset median" rule but are still worth checking per-
gene, as done here) — check each missing panel gene's status against
the documented noncoding/RP/MT patterns AND the `min.cells=50`
expression-floor hypothesis directly (not assumed), for every flagged
or small-panel deviation, not only the largest one.

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
compute whether the same underlying cause recurs here, or whether this
is a genuinely different, dataset-specific cause. **Round-1 review fix**:
confirmed directly this round (see Blocker 4 above) that GSE225857's
gene reduction is NOT a neutral smaller-reference phenomenon but the
product of the original authors' own `min.cells=50` expression floor
plus systematic noncoding/RP/MT gene-class removal (both confirmed
directly from `github.com/jalon9358/LianLab_CRCLM`'s real preprocessing
code) — so the required investigation is per-gene (checked against both
mechanisms directly), not a generic "reference gap" assumption, and is
extended to `D_Gut-shared` and the small revCSC panels too (proportional
losses that don't trigger the aggregate deviation rule but are still
worth checking for a small panel), not only `P_Gut-specific`. Locked as
a required investigation during the compute PR.

## Cell/sample-scoring-population contract (round-1 review fix — see
Blocker 1 above for the full finding)

Locked per dataset: `GSE231559` primary population = the 9 CLiM + 6
primary-tumor samples' cells (15 samples, sample-level tumor label, no
per-cell malignant annotation available); 11 paired-normal samples
loaded/scored separately, not pooled into the primary result.
`GSE285990` = all 10 `P01_LM`-`P10_LM` samples, single population.
`GSE225857` = the non-immune (CD45-) fraction only, `GSM7058755`,
41,892 cells; the immune (CD45+) fraction, `GSM7058754`, 196,473 cells,
is explicitly OUT OF SCOPE for this panel set (not a sensible target for
an epithelial-developmental-program gene signature).

## Raw-counts / normalization contract (round-1 review fix — see
Blocker 2 above for the full finding)

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

**Locked internal-representation contract** (matches
`crc_gut_scoring_core.py`'s exact requirement, confirmed by directly
re-reading the frozen module — not the generic "raw-vs-processed" check
alone): every loaded dataset must place its raw integer counts, on the
canonicalized gene axis (see the gene-axis contract below), into
`adata.layers["counts"]` — the exact layer `compute_detectability()`
reads directly — then rebuild `adata.X` via
`sc.pp.normalize_total(adata, target_sum=1e4)` followed by
`sc.pp.log1p(adata)` applied to that same layer, identical to
`load_atlas()`'s own construction. This must hold after any sample
concatenation within a scoring population (e.g. GSE231559's 15-sample
primary population), not just per-sample before merging.

## Gene-axis canonicalization contract (round-1 review fix — see
Blocker 3 above for the full finding)

Locked, shared across the coverage check and the compute-time scoring
loader (one function, not two independent reimplementations that could
drift apart): (1) strip Ensembl version suffixes
(`gene_id.split(".")[0]`) as the canonical bare-ID axis for
`GSE231559`/`GSE285990`; (2) for `GSE225857`, map gene symbols to
Ensembl IDs via `hgnc_symbol_ensembl_map.tsv`, dropping any symbol with
no unambiguous mapping (logged with counts, not silently dropped); (3)
after either canonicalization path, assert no duplicate Ensembl IDs
result — raise loudly if two input genes collapse onto the same ID,
never silently sum/average/pick-first; (4) `GSE231559`'s two-reference
intersection (32,732 common genes, computed above) is itself part of
this canonical axis, applied before the version-stripping/mapping step
interacts with anything else.

## What this design does not do

Does not run any scoring compute (design-only, per this project's
standing discipline of separate design→compute PRs). Does not score the
4 bulk cohorts (explicitly deferred, see Scope). Does not score
GSE225857's immune (CD45+) fraction (explicitly out of scope, see the
cell/sample-scoring-population contract above). Does not run any
cross-cohort composition/contrast analysis. Does not modify the
existing, frozen `crc_gut_scoring_core.py` machinery or the 13-panel
gene-set inventory — reused exactly as-is, matching the Step 6
dataset-extension precedent's explicit non-goal of re-deriving already-
frozen scoring machinery.

## Review history

- **Round 1 (REQUEST_CHANGES — 4 real blockers, all independently
  re-verified against real committed code/live sources before fixing,
  none disputed)**: (1) the scoring population (which cells/samples get
  scored, per dataset) was never locked — fixed with an explicit
  per-dataset contract, including excluding GSE225857's immune fraction
  entirely as not a sensible target for this panel set. (2) the loader
  contract didn't reproduce `crc_gut_scoring_core.py`'s exact
  `layers["counts"]`/`normalize_total`+`log1p`-into-`X` requirement,
  the same category of gap PR #33's review caught in the Step 6
  extension design — fixed, locked explicitly. (3) gene-axis
  canonicalization (version-suffix stripping, duplicate-collision
  handling, shared code path between coverage-check and scoring loader)
  was underspecified relative to the coverage numbers already reported
  — fixed, locked as one shared canonicalization function. (4)
  GSE225857's reduced ~17K-gene universe was framed as a neutral smaller
  reference background — confirmed directly via the real authors' own
  preprocessing code (`github.com/jalon9358/LianLab_CRCLM`) that it's
  actually `min.cells=50` expression filtering plus systematic
  noncoding/RP/MT gene-class removal; investigated directly which of
  `D_Gut-shared`'s 2 missing genes matched either mechanism (neither
  matched the noncoding/RP/MT patterns, consistent with the
  `min.cells=50` expression floor instead) — fixed, required
  investigation extended beyond `P_Gut-specific` to `D_Gut-shared` and
  the small revCSC panels too.

Submitting for round-2 review before any compute, same discipline as
every prior step this session.
