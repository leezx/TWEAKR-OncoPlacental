# Step 4 statistical design: the actual test/threshold for each evidence type

Builds on `docs/STEP4_DFP_DESIGN.md` (the F-developmental/P-developmental logic) and `results/04_dfp_signature/replicate_structure_audit.md` (which datasets/donors are eligible). This doc picks the actual statistical test and threshold logic for each of the three evidence types before any DE is run — per the reviewer's suggested next step after PR #7's APPROVE: "trophoblast pseudobulk DE 的 statistical-design/threshold audit."

## 1. P-developmental's positive evidence: trophoblast-vs-rest, paired by donor

The replicate-structure audit already established that all 4 eligible datasets have real donor-level pairing (each donor contributes both a trophoblast and a non-trophoblast population). That pairing should be used directly, not thrown away:

- **Pseudobulk construction**: per donor, per dataset — sum raw counts across all trophoblast-lineage cells → one trophoblast pseudobulk profile; sum raw counts across all non-trophoblast cells → one non-trophoblast pseudobulk profile. Normalize each pseudobulk profile to CPM, then log1p, **within that dataset only** (no cross-dataset raw-magnitude comparison — the 5 placental studies differ in protocol/depth just as much as the cross-platform boundary in `STEP3_METHOD_CONTRACT.md` addresses for HDMA/GTEx/HPA; the same discipline applies here even though all 5 are technically "the same platform type").
- **Test**: paired Wilcoxon signed-rank test per gene, per dataset, matched by donor (log-CPM trophoblast vs. log-CPM non-trophoblast) — leverages the donor pairing directly instead of an unpaired two-group test, which would throw away the fact that each donor's own baseline depth/composition is controlled for.
- **Effect size**: require both FDR-corrected significance *and* a minimum median paired log2FC across donors (exact cutoff for both — to be set from the real per-gene test-statistic distribution once run, not guessed in advance, same discipline as every other threshold in this project).

**Greenbaum's n=3 is a hard statistical ceiling, not just "weaker evidence"**: an exact paired Wilcoxon signed-rank test on 3 pairs has only 2³ = 8 possible sign arrangements, so the smallest two-sided p-value it can ever produce is 2/8 = 0.25 — it mathematically cannot clear a conventional significance threshold (e.g. FDR < 0.05) no matter how strong the true effect is. Treating it as a fourth equal vote in the quorum would be statistically wrong, not just imprecise. Instead:

- **`replicated_in_placenta(gene)`** = passes significance + effect-size in a quorum of the **3 adequately-powered datasets** (Arutyunyan n=17, Nature2026 n=23, VentoTormo n=12) — exact quorum (e.g. ≥2 of 3) to be set once the real per-gene result distributions are visible.
- **Greenbaum contributes only as an optional directional-concordance booster**: if all 3 of its donors show the same-sign log2FC (regardless of formal significance), that can be reported alongside a gene's call as supportive context, but it never substitutes for one of the 3 required votes and is never sufficient on its own.

## 2. F-developmental's positive evidence: consistency across HDMA's own individual samples

HDMA organs have no internal non-fetal-somatic population to contrast against (established in `STEP4_DFP_DESIGN.md` — the asymmetry is expected). But each organ *does* have real individual-level replicate structure, just discovered by checking `meta_value_counts['Sample']` in the Step 1 inventory JSONs rather than assumed:

| Organ | # individual samples |
|---|---|
| Adrenal | 4 |
| Thyroid | 3 |
| Spleen | 3 |
| Thymus | 3 |
| Skin | 3 |
| Liver | 7 |
| StomachEsophagus | 7 (mixes Stomach- and Esophagus-labeled samples — worth checking whether they should be split before pseudobulking, given `hdma_organ_to_gtex_tissue_map.tsv` already treats Stomach and Esophagus as separate GTEx/HPA tissues on the adult side) |

This lets F-developmental's positive evidence use the same spirit as P-developmental's, adapted to a single-population (no paired contrast) setting:

- **Pseudobulk construction**: per individual `Sample`, per organ — sum raw counts across all cells in that sample, CPM-normalize within that organ (same never-cross-dataset-raw-magnitude discipline).
- **"Elevated" criterion**: a gene's expression clears a threshold computed from **that organ's own distribution across its own samples** — e.g. median CPM across the organ's samples in the top-X-percentile of that organ's gene distribution, **and** detected above a floor in a majority of the organ's own samples (not driven by one outlier sample — directly uses the individual-level replicate structure, exact quorum again to be set empirically per organ given the 3–7 sample range).
- This keeps `elevated_in_fetal_somatic(gene, organ)` a real per-organ statement, consistent with `STEP4_DFP_DESIGN.md`'s framing that F-developmental is organ-specific, not a single fetal-somatic-wide number.

## 3. `adult_excluded`: within-dataset percentile, not an absolute number picked in advance

Per `STEP3_METHOD_CONTRACT.md`, `adult_excluded` is always evaluated as a rank/percentile/threshold **within** GTEx or HPA's own distribution, on `role == adult_negative_reference` rows only (excludes HPA's `placenta` row per the `role` column fix in PR #5).

- **Organ-matched** (feeds F-developmental): gene's expression in the matched GTEx/HPA tissue(s) falls below a percentile threshold computed from that tissue's own gene distribution.
- **Whole-body** (feeds P-developmental): gene must clear the "not elevated" bar **in every one of GTEx's 68 tissues** (or a high fraction, e.g. all but a small number of outlier tissues) — "whole-body adult-excluded" means genuinely absent everywhere in the adult body, not just low on average across tissues, which would let through a gene that's actually high in one specific adult tissue (e.g. testis, brain) as long as it's low elsewhere. This is the same "check every relevant reference, not just the average" spirit that motivated whole-body vs. organ-matched being two separate checks in `STEP4_DFP_DESIGN.md` originally.
- Exact percentile cutoff: not chosen here — the next sub-task before running the real DE is computing GTEx/HPA's real per-tissue expression-percentile distributions and reporting them, same discipline as every prior threshold decision in this project (collision count-mass, replicate-structure audit).

## What this design does NOT do yet

No DE actually run, no percentile cutoffs locked to numbers, no gene lists produced. This picks the test/statistic/pairing logic; the next sub-task computes real distributions against this logic before any cutoff is finalized.
