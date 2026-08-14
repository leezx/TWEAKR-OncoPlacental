# Archived scripts

These are the v1 implementations superseded during PR #19 round-1/round-2 review. Kept for audit history only — **do not run these**; the canonical, correct versions are `scripts/06a_normal_context/build_mike_verzi_human_final.py` and `scripts/06a_normal_context/mike_verzi_dfp_enrichment.py`.

- `build_mike_verzi_human_final.v1_buggy_case_mismatch.py` — mouse→human ortholog mapping had a case-sensitivity bug: the BioMart pull was filtered by the GMT's raw (inconsistently-cased) mouse symbols, so any capitalization mismatch against BioMart's canonical mouse gene name (e.g. `Col1A1` vs canonical `Col1a1`) silently dropped the row before any join happened, misclassifying real orthologs as `NOT_FOUND_IN_BIOMART`. 382/1,923 genes (19.9%) were wrongly excluded this way.
- `mike_verzi_dfp_enrichment.v1_buggy_universe.py` — hypergeometric background hardcoded to N=23,272 (all human protein-coding genes) without restricting `n`/`K`/`k` to a shared, meaningful gene universe — a human gene with no eligible mouse→human ortholog can never appear in a primary mike_verzi signature, yet was counted in the null denominator.

Both bugs are fixed in the canonical scripts. See commit history on branch `step06a-ortholog-overlap-compute-2026-08-14` and PR #19 review thread for the full diagnosis.
