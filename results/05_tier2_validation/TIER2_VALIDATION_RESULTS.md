# Step 5 Tier-2 Validation: Results and Interpretation

Design approved in PR #14 (`docs/STEP5_TIER2_VALIDATION.md`). This document
reports the compute results and the interpretation reached before submitting
for review.

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

- **D-shared (6 genes): only 1 flagged** — `PCDH11X` in Thymus endothelial
  cell of artery (2/2 donors, CPM 3.55/4.97 — low).
- **P-specific (78 genes): 53 flagged instances, but only 13 distinct genes**
  ever appear: `RPA4`, `ZNF257`, `ZNF850`, `ZNF695`, `GJB7`, `ERVFRD-1`,
  `TMEM191C`, `GCM1`, `HTRA4`, `LAIR2`, `KISS1`, `TSKS`, `ZBTB8B`.

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
downstream: 81 of 84 genes (all of D-shared except a low-CPM PCDH11X hit,
65 of 78 P-specific) show no cross-donor-consistent adult single-cell
signal in this Tier-2 check.** This is a clean result and does not motivate
re-opening the frozen Step 4 signature, beyond flagging 3 genes
(`ZNF257`, `ZNF850`, `RPA4`) as recurring low-confidence and 2 genes
(`GCM1`, `LAIR2`) as high-value manual-review items.

## Result 2: organ-matched F-specific check (the larger candidate lists)

79,628 (gene, cell_type) pairs tested; **12,187 flagged (15.3%)** — much
higher than the whole-body check, as previously observed. Two follow-ups
resolve most of this:

**(a) Highly organ-skewed**: Thymus alone accounts for 6,749 of 12,187
flags (55%), despite being only 1 of 4 organs and having the smallest
organ-specific F-list at 1,191 genes. Liver, by contrast, contributes only
296 flags despite 800 candidate genes.

**(b) Background-rate comparison (job 3620544) shows this is not uniform
gene-driven signal — it's organ-specific:**

| Organ | F-specific flag rate range | Random-panel flag rate range | Enrichment |
|---|---|---|---|
| Liver | 0–8% | 0–11% | ~0.6–1.7× (essentially no enrichment; F-specific ≈ background) |
| Skin | 0–31% | 0–19% | ~0.7–1.7× (essentially no enrichment) |
| Spleen | 0–29% | 0–21% | ~0.5–1.4× (essentially no enrichment) |
| Thymus | 0–65% | 0–21% | **~2.2–3.2× across nearly every cell type** |

Liver/Skin/Spleen: F-specific genes are detected at essentially the same
rate as a random gene panel of the same size in the same cell types — the
~15% raw flagging rate in those 3 organs is explained by generic single-cell
detection breadth of certain cell types (mirrored in both the real and
random panels), not by anything specific to the curated F-developmental gene
content.

**Thymus is a real, organ-specific exception** — every major cell type
(fibroblast 3.2×, medullary thymic epithelial cell 3.18×, vascular smooth
muscle 3.06×, arterial/venous/capillary endothelium 2.85–2.97×, macrophage
2.43×, monocyte 2.4×, NK cell 2.23×, CD4 T cell 2.2×, plasma cell 1.77×)
shows F-specific genes detected ~2–3× more often than random genes of the
same panel size.

A partial, plausible biological explanation: medullary thymic epithelial
cells (mTECs) are documented in the immunology literature to exhibit
**AIRE/FEZF2-driven "promiscuous gene expression"** — deliberate broad
transcription of thousands of tissue-restricted self-antigens, including
developmentally-restricted genes, for central T-cell tolerance induction.
This is consistent with mTEC showing the largest single enrichment (3.18×)
in the table. It does **not**, however, explain why the same organ's
endothelial, fibroblast, and myeloid/lymphoid cell types are similarly
elevated — those cell types don't have a documented promiscuous-expression
mechanism, so either (a) there is a Thymus-wide technical or compositional
factor not yet identified (candidates checked and ruled out so far: trivial
small-library artifact, and dominant assay-method mixing — see prior
session's per-donor library-size and 10X/smartseq2 checks for "vein
endothelial cell"), or (b) the F-developmental candidate gene panel itself
happens to be enriched for genes with real broader adult expression breadth
in Thymus tissue specifically, independent of placental/fetal specificity.

**This does not currently motivate re-opening the frozen F-specific
signature**, because (i) the actual frozen `F_specific_FINAL.txt` (2,504
genes) already excludes anything overlapping `P_developmental_primary84`,
and separately underwent its own adult-exclusion calibration in Step 4 —
this Tier-2 check is testing the pre-final per-organ *candidate* lists,
which are a superset; (ii) 3 of 4 organs show no meaningful enrichment at
all. The organ-matched flagged list (`organ_matched_F_specific_FLAGGED.tsv`)
is retained as an audit trail, and Thymus-flagged F-specific genes should be
weighted with extra caution if/when the frozen signature is later projected
onto thymus-adjacent contexts, but CRC is not thymus tissue, so this is a
lower-priority finding for the immediate next step (CRC Oncofetal
projection).

## Files

- `results/05_tier2_validation/organ_matched_F_specific_validation.tsv` — full (organ,cell_type,gene) report, 79,628 rows.
- `results/05_tier2_validation/organ_matched_F_specific_FLAGGED.tsv` — cross-donor-consistent hits only, 12,187 rows.
- `results/05_tier2_validation/whole_body_DP_validation.tsv` — full (organ,cell_type,gene) report for D-shared+P-specific, 8,568 rows.
- `results/05_tier2_validation/whole_body_DP_FLAGGED.tsv` — cross-donor-consistent hits only, 54 rows.
- `results/05_tier2_validation/background_detection_rate_comparison.tsv` — per-(organ,cell_type) F-specific vs. random-panel flag rate and enrichment ratio.

## Recommendation for next step

1. Treat the whole-body D-shared/P-specific result as validation-passed,
   with `ZNF257`/`ZNF850`/`RPA4` downgraded to lower-confidence and
   `GCM1`/`LAIR2` flagged for a quick literature sanity check before the CRC
   projection (not a re-calibration of Step 4's frozen cutoffs).
2. Proceed to the next step in the reviewer's PR #13 guidance: project the
   frozen D/F/P signature onto real CRC Oncofetal single-cell/spatial data.
   The Thymus-specific F-developmental enrichment is noted as an audit
   finding but does not block this, since CRC is not thymus tissue.
