# Step 5 Tier-2 Validation: Results and Interpretation

Design approved in PR #14 (`docs/STEP5_TIER2_VALIDATION.md`). This document
reports the compute results and the interpretation reached before submitting
for review.

**PR #15 review round 1 (REQUEST_CHANGES) fixed two things, both folded
into this version**: (1) the organ-matched background comparison originally
used a single unmatched random draw with no null distribution — replaced
with a 500-permutation, expression-detectability-matched null (see "Result
2" below); (2) a factual wording error ("81/84 genes show no adult
single-cell signal" — the correct split is 69/84 zero-hit, 15/84 with at
least one hit) — corrected in "Result 1" below. The whole-body D/P
pseudobulk methodology itself was confirmed consistent with PR #14's
approved contract and did not need to be rerun.

## What was run

- `scripts/05_tier2_validation/tier2_validation.py` (qsub job 3620542, Argos)
  — the two checks specified in the approved design:
  1. **Organ-matched check**: each of the 4 HDMA-organ-matching Tabula
     Sapiens files (Liver, Skin, Spleen, Thymus) against that organ's own
     frozen `F_developmental_<Organ>.txt` candidate list.
  2. **Whole-body-style check**: the 84 genes that actually made it into the
     frozen signature — D-shared (6) + P-specific (78) — checked as
     gene×cell-type×donor evidence across all 5 Tabula Sapiens organs
     (Liver, Skin, Spleen, Thymus, Large_Intestine) combined.
- `scripts/05_tier2_validation/background_detection_rate.py` (qsub job
  3620544) — for each of the 4 organ-matched organs, computes the identical
  cross-donor-consistent-hit rate for a random gene panel of the same size,
  to separate "these genes are elevated" from "this cell type is generically
  detected broadly regardless of gene identity."

All outputs pulled back from Argos and verified byte-exact (md5) against the
Argos copies before this write-up.

## Bugs found and fixed before these results were trustworthy

1. **`tier2_validation.py` had no `if __name__ == "__main__":` guard.** The
   first attempt at the background-rate comparison (`background_detection_rate.py`,
   job 3620540) imports `tier2_validation`'s helper functions — but importing
   a module with no main guard re-executes all of its top-level script code
   as an import side effect. The organ-matched check silently re-ran a second
   time (visible in the job log), then crashed when it reached the
   whole-body-style section. Fixed by wrapping both checks in
   `run_organ_matched_check()` / `run_whole_body_check()` functions called
   only under `if __name__ == "__main__":`.
2. **D-shared/P-specific/whole-body gene-list files were never pushed to
   Argos.** Only the 7 `F_developmental_<Organ>.txt` files had been
   re-pushed after the earlier OUT_DIR-clobbering incident; `D_shared_FINAL.txt`,
   `P_specific_FINAL.txt`, `P_developmental_primary84.txt`, etc. did not
   exist on Argos at all, so **the whole-body-style check — the part that
   actually matters for the frozen signature — had never successfully run
   before this round.** Fixed by pushing the full local `dfp_gene_sets/`
   directory to Argos and verifying md5-exact.

Both bugs are in my own harness/pipeline plumbing, not in the underlying
statistical design from PR #14.

## Result 1: whole-body D-shared + P-specific check (the signature that matters)

8,568 (organ, cell_type, gene) triples tested across 84 genes; **54 flagged
as cross-donor-consistent (0.63%)**.

Of the 84 genes tested, **69/84 show zero cross-donor-consistent hits
anywhere** across all 5 organs' cell types; **15/84 have at least one hit**
(counted directly from the 84-gene-by-gene breakdown, not estimated):

- **D-shared (6 genes): only 1 of 6 flagged** — `PCDH11X` in Thymus
  endothelial cell of artery (2/2 donors, CPM 3.55/4.97 — low).
- **P-specific (78 genes): 14 of 78 flagged** (53 flagged instances total,
  concentrated in these 14 genes): `RPA4`, `ZNF257`, `ZNF850`, `ZNF695`,
  `ZNF114`, `GJB7`, `ERVFRD-1`, `TMEM191C`, `GCM1`, `HTRA4`, `LAIR2`,
  `KISS1`, `TSKS`, `ZBTB8B`.

(An earlier draft of this document said "81/84 genes show no adult
single-cell signal" — that was a factual error caught in PR #15 review
round 1: the correct split is 69/84 zero-hit, 15/84 with at least one hit.
Corrected here; no change to the underlying compute or the flagged-gene
list itself.)

Two patterns worth flagging for manual review:

- **`ZNF257`, `ZNF850`, `RPA4` recur across nearly every organ and many
  unrelated cell types** (endothelial, macrophage, B/T/NK cells, monocytes,
  hepatic sinusoid endothelium, large-intestine transit-amplifying cells) at
  low, fairly uniform CPM (1.2–25). This pattern — same handful of genes,
  broad cross-tissue low-level hits — looks like generic low-level adult
  transcription these three genes specifically carry, not organ-specific
  biology. Candidates for exclusion from the frozen signature or at minimum
  flagging as lower-confidence P-specific markers.
- **`GCM1` is a genuine outlier**: Skin regulatory T cells, CPM 158–314
  (median/max), 2/2 donors consistent — far higher than any other flagged
  hit. GCM1 is the master syncytiotrophoblast transcription factor, so a
  hit this large in an adult non-placental cell type deserves a closer look
  before trusting it as placenta-specific, even though it's currently
  single-organ/single-cell-type.
- **`LAIR2`** recurs in NK cells specifically across two organs (Spleen CPM
  26–55, Thymus CPM 71–128) — moderately high and cross-organ-consistent
  within one cell lineage. LAIR2 is a known soluble decoy receptor secreted
  by activated/regulatory lymphocytes in normal adult immunology, which is a
  plausible innocent explanation independent of placental biology, but it's
  a real candidate false positive worth a literature check before the CRC
  projection step, since NK cells are relevant to a tumor microenvironment.

**Bottom line for the part of the signature that actually gets used
downstream: 69 of 84 genes (5 of 6 D-shared, 64 of 78 P-specific) show
zero cross-donor-consistent adult single-cell hits anywhere in this
Tier-2 check; 15 of 84 (1 D-shared, 14 P-specific) have at least one
hit**, and of those 15, most are low-CPM/single-cell-type occurrences —
this is still a clean result relative to the whole signature and does not
motivate re-opening the frozen Step 4 cutoffs, beyond flagging 4 genes
(`ZNF257`, `ZNF850`, `RPA4`, `ZNF695`, each appearing at low CPM in
multiple organs/cell types) as recurring low-confidence, `ZNF114` as a
single low-CPM (3.3-4.7) one-off occurrence, and 2 genes
(`GCM1`, `LAIR2`) as high-value manual-review items.

## Result 2: organ-matched F-specific check (the larger candidate lists)

79,628 (gene, cell_type) pairs tested; **12,187 flagged (15.3%)** — much
higher than the whole-body check, as previously observed. Two follow-ups
resolve most of this:

**(a) Highly organ-skewed**: Thymus alone accounts for 6,749 of 12,187
flags (55%), despite being only 1 of 4 organs and having the smallest
organ-specific F-list at 1,191 genes. Liver, by contrast, contributes only
296 flags despite 800 candidate genes.

**(b) Background comparison — corrected in PR #15 review round 1.** The
first version of this check (`background_detection_rate_comparison.tsv`,
job 3620544) drew a single random gene panel per organ (fixed seed, sampled
uniformly from the full transcriptome) and reported one fold-enrichment
number per cell type. Reviewer flagged two real problems: no null
distribution/uncertainty from a single random draw, and no
expression/detectability matching — F-specific genes already passed Step
4's fetal-expression selection, so an unmatched background pulls in large
numbers of low/never-expressed genes, mechanically deflating the
"background" rate and inflating the apparent enrichment.

**Fixed**: rewrote `background_detection_rate.py` to run a real
matched-permutation null (`background_permutation_null.tsv`, job 3620548):
genes are binned into detectability deciles by their own overall adult
detection frequency in that organ's Tabula Sapiens data (fraction of
eligible donor-celltype samples with CPM≥1, computed from the same data
being tested), and each of 500 permutations draws one random gene per
F-specific gene from its matching decile — giving a background panel
matched in size AND expression-level composition. Reports null median, 95%
interval, and an empirical one-sided p-value (fraction of the 500
permutation flag rates ≥ the observed rate) per (organ, cell_type).

**Result: the earlier ~2.2–3.2× Thymus enrichment was substantially an
artifact of the unmatched background.** Under the expression-matched null:

| Organ | Cell types tested | Significant at p<0.05 | Enrichment among significant hits |
|---|---|---|---|
| Liver | 12 | **0** | — (no organ-matched signal at all) |
| Skin | 18 | 4 (stromal cell, muscle cell, endothelial cell, T cell) | 1.11–1.27× |
| Spleen | 23 | 2 (endothelial cell, innate lymphoid cell) | 1.20–1.34× |
| Thymus | 29 | 7 (capillary/arterial/venous endothelium, fibroblast, macrophage, medullary thymic epithelial cell, vascular smooth muscle) | 1.06–1.29× |

Thymus still has the most cell types clearing significance (7/29, vs. 0/12
in Liver), and medullary thymic epithelial cell remains among the
significant hits (1.29×) — consistent with, though far more modest than
originally reported, the AIRE/FEZF2 "promiscuous gene expression"
candidate explanation. But **22 of Thymus's 29 cell types, and every one of
Liver's 12, show no significant difference from an expression-matched
random panel** — the naive uniform-background comparison had overstated
both the breadth (nearly every Thymus cell type) and the magnitude (2.2–3.2×
vs. the corrected 1.06–1.29×) of the effect. The real, corrected finding is
a narrower and much smaller signal in Thymus stromal/endothelial/mTEC
populations specifically, not an organ-wide 2–3× elevation.

**This does not motivate re-opening the frozen F-specific signature**,
for the same structural reasons as before (the frozen `F_specific_FINAL.txt`
already excludes anything overlapping `P_developmental_primary84` and
underwent its own adult-exclusion calibration in Step 4; this Tier-2 check
tests the pre-final per-organ *candidate* lists, a superset) — reinforced
now by the corrected effect sizes being modest (≤1.34×) rather than the
originally reported 2–3×. The organ-matched flagged list
(`organ_matched_F_specific_FLAGGED.tsv`) and the permutation-null table
(`background_permutation_null.tsv`) are retained as the audit trail;
Thymus's 7 significant cell types (and Spleen's/Skin's smaller sets) should
still get extra caution if the frozen signature is later projected onto
thymus- or spleen-adjacent contexts, but CRC is not thymus tissue, so this
remains a lower-priority finding for the immediate next step (CRC Oncofetal
projection).

## Files

- `results/05_tier2_validation/organ_matched_F_specific_validation.tsv` — full (organ,cell_type,gene) report, 79,628 rows.
- `results/05_tier2_validation/organ_matched_F_specific_FLAGGED.tsv` — cross-donor-consistent hits only, 12,187 rows.
- `results/05_tier2_validation/whole_body_DP_validation.tsv` — full (organ,cell_type,gene) report for D-shared+P-specific, 8,568 rows.
- `results/05_tier2_validation/whole_body_DP_FLAGGED.tsv` — cross-donor-consistent hits only, 54 rows.
- `results/05_tier2_validation/background_detection_rate_comparison.tsv` — **superseded**, single-draw unmatched-background comparison from before PR #15 review round 1; kept only as a record of what was corrected, not to be cited as a result.
- `results/05_tier2_validation/background_permutation_null.tsv` — the corrected 500-permutation, expression-matched null per (organ, cell_type): observed rate, null median/95% interval, empirical p-value, enrichment vs. null median.

## Recommendation for next step

1. Treat the whole-body D-shared/P-specific result as validation-passed,
   with `ZNF257`/`ZNF850`/`RPA4`/`ZNF695`/`ZNF114` downgraded to lower-confidence and
   `GCM1`/`LAIR2` flagged for a quick literature sanity check before the CRC
   projection (not a re-calibration of Step 4's frozen cutoffs).
2. Proceed to the next step in the reviewer's PR #13 guidance: project the
   frozen D/F/P signature onto real CRC Oncofetal single-cell/spatial data.
   Thymus's 7 significant cell types (mTEC, fibroblast, vascular smooth
   muscle, 3 endothelial subtypes) are noted as a modest (1.06-1.29x),
   expression-matched-null-confirmed audit finding but do not block this,
   since CRC is not thymus tissue.
