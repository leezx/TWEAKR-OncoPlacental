# Step 4 statistical design: the actual test/threshold for each evidence type

Builds on `docs/STEP4_DFP_DESIGN.md` (the F-developmental/P-developmental logic) and `results/04_dfp_signature/replicate_structure_audit.md` (which datasets/donors are eligible). This doc picks the actual statistical test and threshold logic for each of the three evidence types before any DE is run — per the reviewer's suggested next step after PR #7's APPROVE: "trophoblast pseudobulk DE 的 statistical-design/threshold audit."

## 1. P-developmental's positive evidence: trophoblast-vs-rest, paired by donor

The replicate-structure audit already established that all 4 eligible datasets have real donor-level pairing (each donor contributes both a trophoblast and a non-trophoblast population). That pairing should be used directly, not thrown away:

- **Pseudobulk construction**: per donor, per dataset — sum **raw** counts across all trophoblast-lineage cells → one trophoblast pseudobulk profile; sum raw counts across all non-trophoblast cells → one non-trophoblast pseudobulk profile. Kept as raw counts (not CPM/log-transformed) going into the primary model — see below.
- **Primary test (revised per PR #8 round-1 review)**: a **count-based paired RNA-seq DE model** — edgeR quasi-likelihood F-test or DESeq2, fit **independently per dataset**, with an explicit paired design (`~ donor + trophoblast_status`, donor as a blocking factor). This is the standard approach for paired pseudobulk RNA-seq DE and properly models mean-variance, library size, and dispersion from the actual count distribution — a paired Wilcoxon signed-rank test on log-CPM was the first draft's choice, but it throws away the count-model information the data actually has, and its behavior degrades under the heavy zero/tie structure typical of pseudobulk gene counts. Downgraded to a **secondary sensitivity/direction-consistency check**, not the primary DE engine.
- **Effect size**: the model's own logFC (edgeR/DESeq2), plus a pre-locked minimum effect-size threshold decided *before* looking at the real distribution (not tuned after the fact) — exact number still TBD, same discipline as every other threshold in this project. **Must be directional, not absolute**: `replicated_in_placenta`'s positive criterion is `logFC >= cutoff` (trophoblast-enriched), never `abs(logFC) >= cutoff` — since the edgeR contrast is defined as `troph vs nontroph`, a large *negative* logFC means strongly depleted-in-trophoblast (e.g. immune/endothelial markers), which is the opposite of the placenta-developmental signal P-developmental is trying to capture. Genes with strongly negative logFC can be kept as a separate `trophoblast_depleted` QC/informational set, but must never enter the P-developmental positive program through an absolute-value threshold (caught in PR #9 round 1 review: the first draft's `abs(logFC)>=1` proposal would have let PTPRC/PECAM1-like genes qualify).

**Greenbaum's n=3 exclusion reasoning revised to match the new primary model**: the original justification (paired Wilcoxon's minimum achievable two-sided p-value is 2/8 = 0.25 with 3 pairs) was mathematically correct but is now the wrong reason, since Wilcoxon is no longer the primary model. The real reason: **n=3 is too thin to reliably estimate dispersion/effect size in an independent per-dataset edgeR/DESeq2 fit** — dispersion estimation from 3 paired samples is unstable regardless of test choice. The conclusion is unchanged, just the justification now matches the actual primary model:

- **`replicated_in_placenta(gene)`** = passes significance + effect-size in a quorum of the datasets with real raw counts for the primary count-based model. **Revised per `results/04_dfp_signature/raw_counts_availability_audit.md`**: a direct check of the actual files (not assumed from the design stage) found VentoTormo's `decidua-v3.h5ad` has no raw counts anywhere (X is irreversibly normalized, no `.raw`, no layers, no separate raw-count file on disk) — it cannot supply a primary vote. Only **Arutyunyan (n=17) and Nature2026 (n=23)** are currently usable for the primary edgeR/DESeq2 model, so the quorum is **2-of-2** (interim; see audit doc's option 2 for a possible future 3rd-vote restoration via fresh ArrayExpress data, out of scope now).
  - **Effect-size cutoff proposal (per `results/04_dfp_signature/hpa_placenta_enrichment_audit.md`, pending sign-off)**: `logFC ≥ 0.75` per dataset. Calibrated against an independent HPA **placenta tissue-level** enrichment check (bulk RNA consensus across HPA's 40 tissues, not derived from the scRNA-seq DE itself and not a single-cell trophoblast reference) across candidates 0.5/0.75/1.0: enrichment purity rises monotonically with cutoff (2-of-2-overlap fold-enrichment 6.44×→11.15×→18.85×, all p≪1e-30), consistent whether measured against HPA's own unmodified fourfold-tissue-enrichment rule or that rule plus a project-defined nTPM≥1 noise floor (sensitivity check in the audit doc — the floor doesn't drive the conclusion). But 1.0 loses 2 of 9 canonical trophoblast markers in Nature2026 (ERVFRD-1, KRT7) while 0.75 retains all 9 in both datasets — 0.75 is the proposed balance point (real purity gain over 0.5, no marker-recall cost that 1.0 has).
- **VentoTormo and Greenbaum both contribute only as optional directional-concordance boosters** (paired Wilcoxon or simple sign-agreement on their own normalized/count data) — external supportive context, never required votes, never sufficient alone, and never part of the primary edgeR/DESeq2 model (VentoTormo: no usable raw counts; Greenbaum: n=3 too thin for dispersion/effect estimation in an independent fit).

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
