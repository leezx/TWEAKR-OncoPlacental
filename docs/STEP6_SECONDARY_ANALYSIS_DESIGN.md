# Step 6 secondary analysis: revCSC-high developmental composition + M11 concordance — design

Implements the **secondary analysis** from the already-approved Step 6 design
(`docs/STEP6_CRC_PROJECTION_DESIGN.md`, "Revised analysis structure", item 2),
explicitly deferred out of scope in PR #27 (primary analysis, merged). Reuses
the primary compute's already-scored, already-verified per-cell percentiles
(`crc_gut_scoring_all_panels.parquet`, job 3621066) wherever possible —
this analysis does **not** re-score revCSC or gut D/F/P; it only scores the
one new gene set (M11) and analyzes already-computed data.

## What the approved design requires (context, not re-litigated here)

Two questions, both restricted to `CRC_single_cell_atlas_2025` (same dataset
bound as the primary analysis; per-dataset extension is separate, later work):

1. Are revCSC-high malignant cells preferentially P-specific, a particular
   F-lineage (especially GI/intestinal), D-shared, or split into multiple
   separable developmental substates — using a **calibrated, not
   post-hoc-picked** revCSC-high threshold.
2. **M11 concordance**, within the 297,307-cell M11 subset only: do
   revCSC-high cells enrich for high M11 scores (the atlas's own independent
   NMF finding)?

## 1. revCSC-high threshold (pre-registered, not outcome-dependent)

**Primary threshold: top decile (≥90th null-calibrated empirical percentile)
of `revCSC_primary27_full`**, the frozen 27-gene primary anchor, evaluated on
all 665,473 malignant cells — already computed in
`crc_gut_scoring_full/crc_gut_scoring_all_panels.parquet` (PR #27, job
3621066, byte-verified). Top-decile is a standard, generic convention for
"high-scoring" gene-set-score cells in the single-cell literature (not
selected after looking at any composition result), and is directly
interpretable against the null-calibration itself (a cell in the top decile
of its own signature's matched-null distribution is, by construction,
higher than ~90% of size/detectability-matched random draws).

**Two pre-registered sensitivity cutoffs, both reported alongside the
primary threshold, not substituted for it if they disagree**: top 5%
(percentile ≥95) and top quintile (percentile ≥80). If the composition
conclusion changes materially across these three cutoffs, that
threshold-sensitivity is itself reported as a finding, not resolved by
picking whichever cutoff gives the cleanest story.

**`revCSC_extended28_full`'s equivalent thresholds are also computed and
reported alongside the primary-variant result** (same convention as PR #27),
not treated as a separate primary threshold.

## 2. Developmental composition of revCSC-high cells

For each revCSC-high cell (at each of the 3 cutoffs), report:

- **Per-axis percentile distribution**: `D_Gut-shared`, `F_Colon-specific`
  (primary regional), `F_SI-specific` (secondary regional), `F_Gut-specific`
  (global/coarse), `P_Gut-specific` — all already computed, no re-scoring.
  Median + IQR per axis, reported for revCSC-high cells vs. the full cohort
  (not just revCSC-high in isolation, so a reader can see the actual shift).
- **Predominant-axis assignment**: per cell, the D/F(regional)/P axis with
  the highest percentile among {`D_Gut-shared`, `F_Colon-specific`,
  `F_SI-specific`, `P_Gut-specific`} (global `F_Gut-specific` excluded from
  this specific assignment step — it is the coarse/secondary summary, not a
  lineage-resolved axis, so including it would double-count with the
  regional F axes it's built from). Reports the distribution of predominant
  axis across revCSC-high cells (e.g., "62% predominantly F_Colon-specific,
  9% predominantly P_Gut-specific, ...") — this directly answers "split into
  multiple separable developmental substates" by showing whether one axis
  dominates or the assignment spreads across several.
- **Donor/study-aware aggregation (required, per the approved design's
  standing discipline)**: the composition breakdown above is computed once
  pooled across all revCSC-high cells, and again as an unweighted mean
  across donors (composite `donor_key`, PR #27's contract) and across
  studies — a composition conclusion driven by one large-donor/study is
  reported as such, not smoothed into a pooled number.

## 3. M11 scoring (new gene set, only new compute this analysis needs)

**297,307-cell M11 subset identification**: exact barcode match between
`NMF/viz_signature_MM_alt_clean_byscore_wardD2/addmodulescore.df.tsv`'s
index column (297,307 rows, confirmed) and the atlas's `obs_names` — same
provenance already verified in the approved design (100% overlap by exact
barcode string). No re-derivation of the subset; used as-is.

**M11 gene list**: `NMF/metamodule_fnmf/MM_alt_clean_byscore_wardD2/deliver.mm_top_genes.csv`'s
`M11` column (top50, 50 Ensembl IDs, already bare-Ensembl — no mapping
needed) is the primary version, matching the approved design's Jaccard
validation (top50 unique-gene version had M11's highest Jaccard match to
revCSC). `unique_top_versions/deliver.mm_top_genes.unique100.csv` (39 genes)
and `...unique200.csv` (62 genes) are scored as secondary/sensitivity
versions, not substituted for top50.

**M11 × gut-D/F/P overlap audit (redone against the gut-specific sets,
since the approved design's audit predates PR #25's gut re-anchor and used
pan-organ D/F/P)** — checked directly, real result:

| M11 version | D_Gut-shared | F_Gut-specific | F_Colon-specific | F_SI-specific | P_Gut-specific |
|---|---|---|---|---|---|
| top50 (n=50) | 0 | 4 | 2 | 3 | 0 |
| top100/unique (n=39) | 0 | 5 | 2 | 4 | 0 |
| top200/unique (n=62) | 0 | 8 | 4 | 7 | 0 |

Zero overlap with `D_Gut-shared`/`P_Gut-specific` throughout (consistent
with revCSC's own gut-overlap finding, PR #25) — no exclusion needed for
those two axes. Real (not negligible, unlike the pan-organ audit's
single-gene `MALAT1` finding) overlap with the F axes: top50 shares
`ENSG00000102265`/`ENSG00000135404`/`ENSG00000163191`/`ENSG00000213719`
with `F_Gut-specific`. **Overlap-exclusion contract**: M11 scored against
each F axis uses that axis's overlap genes dropped from the M11 gene list
(same pattern as revCSC's `_minus_CLU`/`_minus_ASS1` variants); M11 scored
against `D_Gut-shared`/`P_Gut-specific`/`revCSC` uses the full,
non-excluded M11 gene list (zero overlap, no exclusion needed).

**Scoring method**: identical null-calibrated empirical-percentile method as
the primary compute (`score_genes_fast`, `N_BINS=20`, same fixed seed
`20260815`, `N_PERM=100` — this subset is 297,307 cells, well within the
regime the convergence check already validated as adequate at 20,000-cells
scale; no separate convergence check re-run needed for a same-order-of-
magnitude cell count using the identical method already validated).
Detectability strata computed on the 297,307-cell subset itself (per the
core library's existing contract: strata always match the population being
scored), not reused from the full-665,473-cell run.

## 4. M11 concordance test

Within the 297,307-cell subset only:

1. **Continuous correlation**: M11 percentile (full + overlap-excluded
   variants) vs. `revCSC_primary27_full` percentile (already computed, no
   re-scoring) — Pearson + Spearman, same within-study-then-pooled and
   donor/study-aware leave-one-out validation as the primary analysis
   (reusing `crc_gut_scoring_primary_analysis.py`'s machinery directly, not
   reimplemented).
2. **Enrichment test (the design's actual stated question — "do
   revCSC-high cells enrich for high M11 scores")**: using the same 3
   revCSC-high thresholds (§1) restricted to this subset, test whether
   revCSC-high cells are enriched for M11-high status (M11 percentile ≥90,
   the matched threshold convention) via a 2×2 odds ratio, computed within
   each donor and pooled via Mantel-Haenszel (donor-stratified, avoiding
   the same single-large-donor pseudoreplication risk as the primary
   analysis) rather than a single pooled contingency table.

## Scope boundary (explicit)

- Still `CRC_single_cell_atlas_2025` only — no other dataset.
- Does not re-score or re-derive revCSC or any gut D/F/P set.
- Does not touch the tertiary analysis (full-atlas revCSC-independent D/F/P
  landscape) — separate, later work.
- M11's top50/100/200 gene lists and the 297,307-cell subset are used
  exactly as already independently derived (per the approved design's own
  "not re-derived by this project" discipline) — no re-clustering, no
  re-running the NMF pipeline.

Submitting for review before any qsub job runs, same discipline as every
prior step.
