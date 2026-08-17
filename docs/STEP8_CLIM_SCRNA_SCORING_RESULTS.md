# Step 8: CLiM/CLuM scRNA-seq D/F/P/revCSC scoring — real compute results

Executes `docs/STEP8_CLIM_SCRNA_SCORING_DESIGN.md` (PR #37, APPROVE after
3 review rounds) for real: loads the 4 locked scoring populations from
the Step 7 scRNA-seq cohorts (GSE231559, GSE285990, GSE225857), runs the
locked coverage check + required empirical missing-gene investigation,
then null-calibrates D/F/P/revCSC panel scores (`n_perm=500`, same
method as `crc_gut_scoring_core.py`'s existing full-atlas run) for every
cell. **Scoring only** — cross-cohort/cross-population statistical
comparison and any founding-question-relevant interpretation are
explicitly out of scope here, per the design's own scope boundary
("bulk cohorts + cross-cohort composition explicitly deferred"); this PR
reports the real, per-population scoring output and the real coverage/
investigation findings, nothing more.

All compute ran on Argos via `qsub` (job 3621331, `argos-codex` conda
env, `-pe pvm 8`, wallclock 15,739s ≈ 4.37h, exit status 0, 0 failures —
same qsub discipline as Step 6's full-atlas run). Scripts pushed and
results pulled back via `cat | ssh argos-qsub1 "cat > ..."` /
`ssh argos-qsub1 "cat ..." > local`, byte-exact md5 verified in both
directions (scp/rsync blocked in this sandbox).

## Scripts

- `scripts/08_clim_scrna_scoring/clim_scrna_scoring_core.py` — shared
  loader/coverage-check module implementing the design's locked
  contract: axis-verified `adata` construction (`layers["counts"]` +
  `normalize_total(1e4)`+`log1p` rebuild of `X`, identical to
  `crc_gut_scoring_core.load_atlas()`'s own construction), canonicalize-
  then-intersect gene-axis handling, symbol→Ensembl mapping for
  GSE225857 via this project's own Step-2 HGNC table.
- `scripts/08_clim_scrna_scoring/clim_scrna_scoring_driver.py` — loads
  all 4 populations, runs the coverage check + required empirical
  missing-gene investigation, then scores each population against all
  13 panels via `crc_gut_scoring_core.score_all_panels(fast=True)`.
- `scripts/08_clim_scrna_scoring/run_clim_scrna_scoring.sh` — qsub
  submission script (`-pe pvm 8`, same shape as
  `run_crc_gut_scoring_full.sh`).

All three loaders and the end-to-end scoring pipeline were smoke-tested
against real Argos data (not synthetic) before submitting the full qsub
run — see the loader-shape and 300-cell-subset scoring checks in the
prior commit on this branch.

## Scoring populations — real cell counts

| Population | Samples | Cells loaded | Genes (canonical axis) |
|---|---|---|---|
| GSE231559 primary (9 CLiM + 6 primary-tumor) | 15 | 88,792 | 32,732 (exact match to design doc's locked intersection count) |
| GSE231559 normal (separate calibration pass) | 11 | 45,174 | 32,732 |
| GSE285990 (P01_LM–P10_LM) | 10 | 123,677 | 37,487 (uniform reference, no cross-sample split) |
| GSE225857 non-immune (GSM7058755 only) | 1 | 41,892 | 16,616 (16,616/17,515 input symbols mapped to Ensembl; 899 dropped, no unambiguous HGNC mapping) |

GSE225857's immune fraction (GSM7058754, 196,473 cells) was never loaded
— confirmed by construction, not just by omission from the output file
list: `load_gse225857_nonimmune_population()` only ever opens
`GSM7058755_non_immune_counts.txt.gz`.

## Coverage check — all 13 panels × all 4 populations

Full table: `results/08_clim_scrna_scoring/coverage_check.tsv`.

- **GSE231559 primary/normal and GSE285990**: every panel ≥96.15% testable
  gene coverage (`revCSC_extended28_minus_CLU_ASS1`, 25/26 genes, is the
  minimum in all 3 populations), population-median coverage 99.7–99.8%,
  **0 panels flagged** (`>15pts below that population's own median`
  rule). Matches the design doc's "≥96% coverage everywhere" estimate
  closely (the revCSC 28-gene "extended" panels sit at 96.2–96.4%, the
  rest higher). **Round-1 review correction**: an earlier draft of this
  doc stated "≥92.5%" for these 3 populations — that number is real, but
  belongs to a different population entirely (GSE225857 non-immune's
  `revCSC_primary27_full`, 25/27 = 92.6%, see below), not these 3; fixed
  after independently re-deriving each population's true minimum
  directly from `coverage_check.tsv`.
- **GSE225857 non-immune**: general-panel coverage 87.3–92.6% (population
  median 88.9%), consistent with the design doc's "~82–89%" estimate.
  `P_Gut-specific` is the one flagged pair: **44/76 genes testable
  (57.9%)**, vs. the design doc's pre-compute estimate of 51%/58% —
  matches closely.

## Required empirical missing-gene investigation (GSE225857)

Per the design's round-2-corrected contract (the original round-1
mechanism attribution — `min.cells=50`/noncoding/RP-MT filtering — was
retracted in round 2 after being verified against the wrong R object),
this PR does **not** assign GSE225857's coverage gap to any specific
named preprocessing step. Instead, `investigate_missing_genes()` checks
each of the 53 gene×panel rows missing from GSE225857's non-immune gene
axis (D_Gut-shared: 1 gene; the 8 revCSC panels: 2–3 genes each; P_Gut-
specific: 32 genes) against the two reference populations that *do* have
near-complete coverage (GSE231559 primary, GSE285990) — full results in
`results/08_clim_scrna_scoring/gse225857_missing_gene_investigation.tsv`.

**Real finding**: 49/53 of these **gene×panel rows** — not 53 unique
genes; the same missing gene recurs across multiple revCSC
overlap-exclusion panel variants (e.g. `ENSG00000118523`,
`ENSG00000142871`, `ENSG00000291309` each appear in more than one row).
The 53 rows correspond to **36 unique missing genes** (35/36 present in
both reference populations). **Round-1 review correction**: an earlier
draft of this doc described this as "49/53 of these genes," implying 53
distinct genes — fixed to the correct row-vs-gene distinction after
independently re-deriving both counts directly from
`gse225857_missing_gene_investigation.tsv`.

Of the 49 (row-level) / 35 (unique-gene-level) present-in-both cases,
**median real-count detection fraction is 0.000 in both reference
populations** (i.e., at or near zero cells in GSE231559 primary and
GSE285990 actually detect these genes, even though the genes are present
on those axes). This is an honest, hedged conclusion, not a demonstrated
mechanism: it is consistent with most of GSE225857's missing genes being
inherently rare/very-low-expression genes for which the biological
detection signal is weak regardless of platform — but it does not by
itself prove GSE225857's absence isn't *also* driven by some GSE225857-
specific step, since a gene can be both rarely detected everywhere *and*
additionally excluded by a specific upstream filter. No further claim is
made; this remains a real, unresolved partial explanation, reported as
such.

## Null-calibrated scoring — real per-population, per-panel summary

Full per-cell scores: `results/08_clim_scrna_scoring/
{gse231559_primary,gse231559_normal,gse285990,gse225857_nonimmune}
_scores.parquet` (26 columns: `<panel>_percentile`/`<panel>_zscore` × 13
panels). Summary table:
`results/08_clim_scrna_scoring/summary_stats_by_population_panel.tsv`.
Percentile is on crc_gut_scoring_core.py's native 0–100 scale (empirical
percentile among each cell's own `n_perm=500` null draws); z-score is
`(observed − mean(null)) / std(null)`.

| Population | revCSC panels (8), median percentile range | D_Gut-shared median pct | F_Gut/Colon/SI-specific median pct | P_Gut-specific median pct |
|---|---|---|---|---|
| GSE231559 primary | 72.0–74.8 | 51.0 | 82.8 / 69.6 / 76.6 | 44.6 |
| GSE231559 normal | 77.8–82.0 | 55.9 | 77.8 / 62.6 / 67.4 | 49.7 |
| GSE285990 | 86.2–88.6 | 50.7 | 57.8 / 42.6 / 46.0 | 52.4 |
| GSE225857 non-immune | 74.4–76.4 | 62.8 | 84.8 / 50.6 / 79.4 | 61.4 |

Observed, reported without interpretation beyond the descriptive level
(cross-cohort/cross-population comparison is explicitly deferred, see
above):
- All 8 revCSC panels sit well above the 50th-percentile null
  expectation (median 72–89, `frac_cells_ge_p95` 0.15–0.34) in **every**
  population scored, including GSE231559's separately-calibrated
  paired-normal population — i.e., this descriptive elevation is not
  confined to the tumor/metastasis populations in this raw summary.
- `D_Gut-shared` and `P_Gut-specific` sit close to the 50th-percentile
  null expectation in all 4 populations (median percentile range
  44.6–62.8), unlike the revCSC panels. **Round-1 review correction**: an
  earlier draft of this doc stated "43.5–62.8" — 43.52 is
  `P_Gut-specific`'s *mean* percentile in `gse231559_primary`, not a
  median; the true minimum median across all 8 D/P population×panel
  combinations is 44.6 (`P_Gut-specific`, same population) — fixed after
  independently re-checking `summary_stats_by_population_panel.tsv`'s
  `median_percentile` column directly, not the `mean_percentile` column.
- The three F-panels (`F_Gut-specific`, `F_Colon-specific`,
  `F_SI-specific`) vary more by population than the revCSC or D/P
  panels — e.g. GSE285990 (liver-metastasis samples) sits lowest on
  `F_Colon-specific`/`F_SI-specific` (42.6/46.0) of the four
  populations, while GSE225857 non-immune sits highest on
  `F_Gut-specific`/`F_SI-specific` (84.8/79.4) but lowest on
  `F_Colon-specific` (50.6) of the four.

No claim beyond these raw descriptive numbers is made here — attributing
any of this to tumor biology, tissue-of-origin, or the project's
founding F-selectivity question requires the donor/study-aware
statistical treatment this design doc explicitly deferred, not a bare
median comparison across differently-sourced, differently-sized cell
populations.

## A real, expected NaN pattern in z-score columns (not a bug)

`D_Gut-shared_zscore` and `P_Gut-specific_zscore` — the two smallest
scored panels (8 and ≤76 testable genes) — contain NaN for a minority of
cells in every population (0.02–10.0% of cells, worst in
`gse231559_normal`'s `D_Gut-shared` at 10.0%). This is
`crc_gut_scoring_core.py`'s own existing, unchanged, documented behavior
(`zscore = np.where(null_std > 0, zscore, np.nan)`, line 336/379): when
all 500 null draws for a cell produce the identical score (typically 0,
for cells where none of a small panel's control-matched null gene sets
happen to be detected), `null_std` is exactly 0 and the z-score is
undefined by construction — not computable, not a missing-data artifact
of this PR's loaders. Occurs more on the smallest panels, as expected;
`percentile` (which does not divide by `null_std`) is never NaN.

## File manifest (all committed, md5-verified byte-exact Argos↔local)

- `results/08_clim_scrna_scoring/coverage_check.tsv` (52 rows: 13 panels
  × 4 populations)
- `results/08_clim_scrna_scoring/gse225857_missing_gene_investigation.tsv`
  (53 rows)
- `results/08_clim_scrna_scoring/{gse231559_primary,gse231559_normal,
  gse285990,gse225857_nonimmune}_scores.parquet` (per-cell scores)
- `results/08_clim_scrna_scoring/{gse231559_primary,gse231559_normal,
  gse285990,gse225857_nonimmune}_n_testable_genes_per_panel.tsv`
- `results/08_clim_scrna_scoring/summary_stats_by_population_panel.tsv`
  (this doc's summary table, full precision)

## Explicitly out of scope (per the locked design, unchanged here)

- Bulk RNA-seq/microarray cohorts (GSE17536, GSE17537, GSE21510,
  GSE131418, TCGA-CRC) — not scored in this PR.
- Cross-cohort/cross-population statistical comparison, donor/study-
  aware correlation analysis, or any founding-question interpretation —
  deferred, same as the design doc states.
- GSE231559 paired tumor-vs-normal contrast — the two populations were
  scored with **separate** null calibrations (per the design's locked
  contract) and are not directly comparable without an additional joint
  calibration pass, which this PR does not perform.

## Review history (PR #39)

**Round 1: REQUEST_CHANGES**, no compute blocker — the reviewer
confirmed the loader, canonicalize→intersect order, 4 populations,
`n_perm=500`, GSE225857 immune-fraction exclusion, and the
coverage-deviation→empirical-investigation-before-scoring sequence all
match the PR #37-locked contract. 4 real write-up/provenance errors
found, all independently re-verified against the committed artifacts
before fixing, no re-run of job 3621331 required:
1. The "≥92.5%" coverage claim for GSE231559/GSE285990 was real but
   misattributed — that number belongs to GSE225857 non-immune's
   `revCSC_primary27_full` (25/27=92.6%); the true minimum for
   GSE231559/GSE285990 is 96.15% — fixed above.
2. "49/53 of these genes" implied 53 unique genes; the file has 53
   gene×panel rows (36 unique genes, 35 present in both reference
   cohorts) — fixed above, both counts now stated explicitly.
3. "D/P median range 43.5–62.8" used 43.52, `P_Gut-specific`'s *mean*
   percentile in `gse231559_primary`, not a median (true minimum median
   is 44.6) — fixed above.
4. `docs/PROJECT_SUMMARY.md`'s Step 8 status was stale the moment this
   PR was opened (still said "has not been pushed... no open PR") —
   updated to reflect PR #39 actually being open.

**Round 2: APPROVE**, at head `fcd2058`. The reviewer independently
re-checked all 4 corrected numbers against the committed artifacts
before approving, and confirmed the compute itself continues to match
the PR #37 contract with no reason to re-run job 3621331.
