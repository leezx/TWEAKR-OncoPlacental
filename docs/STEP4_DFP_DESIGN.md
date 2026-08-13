# Step 4 design: D-shared / F-specific / P-specific gene set construction

Builds on `docs/STEP3_METHOD_CONTRACT.md` (locks the cross-platform comparison boundary) and directly addresses the open item from PR #5's APPROVE: "P-specific must reflect placenta/trophoblast independence relative to fetal-somatic, not just trophoblast-vs-other-placental-cells + adult-depletion." This doc defines the actual D/F/P split before any qsub job is written, so the statistical design gets reviewed once instead of after a wasted compute run.

## The asymmetry between the two developmental sides

The two "developmental" data sources are not symmetric, and the design has to reflect that rather than force them into identical treatment:

- **HDMA (fetal-somatic)**: each of the 7 organs is already a **pure fetal-somatic reference** — Step 1 confirmed zero trophoblast contamination (`annotv1`/`annotv2` have no placental labels). So "is this gene part of the fetal-somatic developmental program" is answered directly by **expression level within HDMA**, not by an internal DE against some other population.
- **Placental scRNA-seq (5 datasets)**: each dataset is a **mix** of trophoblast lineage cells (VCT/EVT/SCT + subtypes) and non-trophoblast cells (Hofbauer, fetal/decidual endothelial, decidual stromal/immune, maternal cells). So "is this gene part of the placental/trophoblast-specific program" requires an **internal DE**: trophoblast vs. every other cell type in the same dataset, to isolate the trophoblast signal from "just expressed somewhere in placental tissue" (which would leak in maternal/immune genes that have nothing to do with placental developmental biology).

This asymmetry is expected, not a design flaw — it falls directly out of what the two reference datasets actually are (Step 1 SUMMARY.md).

## Three evidence axes, symmetric where the data allows it

| Axis | Positive evidence source | Negative/exclusion evidence source |
|---|---|---|
| Fetal-somatic | HDMA pseudobulk expression level, per organ (no internal DE needed — already pure) | — |
| Trophoblast/placental | Trophoblast vs. other-placental-cells DE, computed **within each of the 5 placental datasets independently** (per `STEP3_METHOD_CONTRACT.md` — developmental evidence never touches GTEx/HPA) | — |
| Adult | — | GTEx (organ-matched + whole-body) / HPA (fills GTEx's Thymus gap; `placenta` row excluded — `role` column in the processed mapping tables), each dataset's own internal rank/percentile/threshold, never raw cross-platform magnitude |

The key structural fix from the reviewer's note: **the fetal-somatic and trophoblast axes are used as *exclusion* references for each other**, exactly the way GTEx/HPA are used as an exclusion reference for "adult":

- For **P-specific**, HDMA's own expression distribution is the fetal-somatic-exclusion check (is this gene *also* substantially expressed in fetal-somatic organs? if so, it's not placenta-specific, no matter how adult-depleted it is).
- For **F-specific**, the placental datasets' trophoblast-vs-rest DE result is the trophoblast-exclusion check (is this gene *also* elevated in trophoblast? if so, it's not fetal-somatic-specific).

## Definitions

- **elevated_in_fetal_somatic(gene, organ)**: gene's HDMA pseudobulk expression in that organ clears a threshold defined from HDMA's own distribution (e.g. top-X-percentile or a detection-rate + mean-expression floor — exact cutoff to be picked empirically against the real pseudobulk distribution, not assumed in advance).
- **elevated_in_trophoblast(gene, placental_dataset)**: gene is significantly higher in trophoblast-lineage pseudobulk vs. other-cell-type pseudobulk within that one placental dataset (Wilcoxon or equivalent on per-sample/per-donor pseudobulk values — requires checking each dataset actually has multiple donor/sample replicates for a valid test, not just multiple cells; to be verified per dataset before choosing the exact test).
- **adult_excluded(gene, ...)**: gene's expression in the relevant GTEx/HPA adult-negative-reference tissue(s) falls below that dataset's own internal threshold/percentile (`role == adult_negative_reference` rows only — `placenta` never counts here per the fixed `role` column).
- **replicated_in_placenta(gene)**: `elevated_in_trophoblast` holds in a **majority of the 5 placental datasets** it's tested in (not just one) — guards against one dataset's technical artifact driving a P-specific call. Exact quorum (e.g. ≥3 of 5, or ≥2 of however many actually have valid replicate structure) to be set once dataset-level DE is actually run and we see how many datasets clear the replicate-structure bar in the first place.

**D-shared** = `elevated_in_fetal_somatic` (any organ) AND `replicated_in_placenta` AND `adult_excluded` (both organ-matched *and* whole-body, since the gene is being claimed as shared across both domains).

**F-specific** = `elevated_in_fetal_somatic` (that organ) AND `adult_excluded` (organ-matched) AND NOT `replicated_in_placenta`.

**P-specific** = `replicated_in_placenta` AND `adult_excluded` (whole-body) AND NOT `elevated_in_fetal_somatic` (any organ) — **this NOT clause is the fix for the reviewer's flagged gap**; without it, a gene that's simply not adult-expressed anywhere (including being a normal fetal-somatic developmental gene) would incorrectly land in P-specific just because it also happens to be elevated in trophoblast.

## What's explicitly left open, to decide against real data rather than assume

- Exact percentile/threshold cutoffs for "elevated" and "adult_excluded" (placeholder language above, not numbers) — Step 4's first qsub pass should compute the real distributions per dataset and report them before locking a cutoff, same discipline as the collision-count-mass check in Step 3.
- Whether all 5 placental datasets have enough donor/sample-level replicate structure for within-dataset DE, or whether some need to be pseudobulk-per-cluster instead (weaker but usable) — needs a real inventory check, not an assumption.
- The exact quorum for `replicated_in_placenta`.
- Whether Nature2026's `snRNA_raw_counts` file (still missing usable annotation per Step 1 Finding #3) gets resolved before Step 4 runs or is excluded from this pass.

## What this design does NOT do yet

No qsub jobs, no gene lists, no thresholds chosen. This is the statistical/logical design only, submitted for review before any compute runs, matching the pattern that worked for `STEP3_METHOD_CONTRACT.md`.
