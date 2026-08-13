# Step 4 design: D-shared / F-specific / P-specific gene set construction

Builds on `docs/STEP3_METHOD_CONTRACT.md` (locks the cross-platform comparison boundary) and directly addresses the open item from PR #5's APPROVE: "P-specific must reflect placenta/trophoblast independence relative to fetal-somatic, not just trophoblast-vs-other-placental-cells + adult-depletion." This doc defines the actual D/F/P split before any qsub job is written, so the statistical design gets reviewed once instead of after a wasted compute run.

## The asymmetry between the two developmental sides

The two "developmental" data sources are not symmetric, and the design has to reflect that rather than force them into identical treatment:

- **HDMA (fetal-somatic)**: each of the 7 organs is already a **pure fetal-somatic reference** — Step 1 confirmed zero trophoblast contamination (`annotv1`/`annotv2` have no placental labels). So "is this gene part of the fetal-somatic developmental program" is answered directly by **expression level within HDMA**, not by an internal DE against some other population.
- **Placental scRNA-seq (5 datasets)**: each dataset is a **mix** of trophoblast lineage cells (VCT/EVT/SCT + subtypes) and non-trophoblast cells (Hofbauer, fetal/decidual endothelial, decidual stromal/immune, maternal cells). So "is this gene part of the placental/trophoblast-specific program" requires an **internal DE**: trophoblast vs. every other cell type in the same dataset, to isolate the trophoblast signal from "just expressed somewhere in placental tissue" (which would leak in maternal/immune genes that have nothing to do with placental developmental biology).

This asymmetry is expected, not a design flaw — it falls directly out of what the two reference datasets actually are (Step 1 SUMMARY.md).

## Building blocks

| Building block | Source | Role |
|---|---|---|
| Fetal-somatic expression evidence | HDMA pseudobulk expression level, per organ (no internal DE needed — already pure per Step 1) | Half of F-developmental (below) — on its own, NOT developmental evidence (see round-1 fix) |
| Trophoblast enrichment evidence | Trophoblast vs. other-placental-cells DE, computed **within each of the 5 placental datasets independently** (per `STEP3_METHOD_CONTRACT.md` — developmental evidence never touches GTEx/HPA), replicated across a quorum of datasets | Half of P-developmental (below) |
| Adult-exclusion evidence | GTEx (organ-matched + whole-body) / HPA (fills GTEx's Thymus gap; `placenta` row excluded — `role` column in the processed mapping tables), each dataset's own internal rank/percentile/threshold, never raw cross-platform magnitude | The other half of *both* F-developmental and P-developmental — folded into each program's definition, not a final filter applied after the fact (round-1 fix, see below) |

Once F-developmental and P-developmental are each independently adult-corrected this way, D/F/P falls out as simple set operations between the two (below) — no separate "exclude against the other developmental side" step is needed, because a gene that's genuinely elevated in *both* programs is D-shared by construction, not something requiring a bespoke cross-check.

## Definitions (revised per PR #6 round-1 review — see "Why adult-depletion has to be inside the developmental-evidence definition, not just a final filter" below)

- **elevated_in_fetal_somatic(gene, organ)**: gene's HDMA pseudobulk expression in that organ clears a threshold defined from HDMA's own distribution (e.g. top-X-percentile or a detection-rate + mean-expression floor — exact cutoff to be picked empirically against the real pseudobulk distribution, not assumed in advance). **On its own this is NOT developmental evidence** — see below.
- **elevated_in_trophoblast(gene, placental_dataset)**: gene is significantly higher in trophoblast-lineage pseudobulk vs. other-cell-type pseudobulk within that one placental dataset (Wilcoxon or equivalent on per-sample/per-donor pseudobulk values — requires checking each dataset actually has multiple donor/sample replicates for a valid test, not just multiple cells; to be verified per dataset before choosing the exact test).
- **replicated_in_placenta(gene)**: `elevated_in_trophoblast` holds in a **majority of the 5 placental datasets** it's tested in (not just one) — guards against one dataset's technical artifact driving a P-specific call. Exact quorum (e.g. ≥3 of 5, or ≥2 of however many actually have valid replicate structure) to be set once dataset-level DE is actually run and we see how many datasets clear the replicate-structure bar in the first place.
- **adult_excluded(gene, matched_organ | whole_body)**: gene's expression in the relevant GTEx/HPA adult-negative-reference tissue(s) falls below that dataset's own internal threshold/percentile (`role == adult_negative_reference` rows only — `placenta` never counts here per the fixed `role` column).

### Why adult-depletion has to be inside the developmental-evidence definition, not just a final filter

PR #6 round 1 caught a real conceptual error in the first draft: `elevated_in_fetal_somatic` alone was being treated as "fetal developmental evidence," but HDMA being trophoblast-free only proves it's a valid fetal-somatic *reference* — it says nothing about whether a gene highly expressed there is part of a genuine developmental program vs. a housekeeping gene, mature-organ-identity gene, or constitutive epithelial/metabolic gene (all of which would trivially satisfy `elevated_in_fetal_somatic`). Left uncorrected, this would have contaminated all three output sets: D-shared could collapse into "expressed in fetal tissue + trophoblast marker" rather than a real shared program; F-specific would absorb ordinary organ-identity genes; and P-specific's exclusion clause would over-prune — a genuine placenta gene with any normal expression in any one fetal organ would get wrongly excluded.

**Fix**: adult-depletion is folded into the definition of each developmental "program" itself, not applied only as a downstream filter:

- **F-developmental(gene, organ)** = `elevated_in_fetal_somatic(gene, organ)` AND `adult_excluded(gene, matched_organ)`
- **P-developmental(gene)** = `replicated_in_placenta(gene)` AND `adult_excluded(gene, whole_body)`

Both are now pre-corrected for adult expression before any set operation runs. This still respects `STEP3_METHOD_CONTRACT.md`'s cross-platform boundary — no HDMA-UMI-vs-GTEx-TPM fold-change is computed; `adult_excluded` is still evaluated via each bulk dataset's own internal rank/percentile.

### D/F/P as a clean three-way partition

With F-developmental and P-developmental defined as above (each already adult-corrected), D/F/P become simple set operations on the same two adult-corrected programs, not three ad hoc gene lists with different conditions:

- **D-shared** = F-developmental AND P-developmental
- **F-specific** = F-developmental AND NOT P-developmental
- **P-specific** = P-developmental AND NOT F-developmental

## What's explicitly left open, to decide against real data rather than assume

- Exact percentile/threshold cutoffs for "elevated" and "adult_excluded" (placeholder language above, not numbers) — Step 4's first qsub pass should compute the real distributions per dataset and report them before locking a cutoff, same discipline as the collision-count-mass check in Step 3.
- Whether all 5 placental datasets have enough donor/sample-level replicate structure for within-dataset DE, or whether some need to be pseudobulk-per-cluster instead (weaker but usable) — needs a real inventory check, not an assumption.
- The exact quorum for `replicated_in_placenta`.
- Whether Nature2026's `snRNA_raw_counts` file (still missing usable annotation per Step 1 Finding #3) gets resolved before Step 4 runs or is excluded from this pass.

## What this design does NOT do yet

No qsub jobs, no gene lists, no thresholds chosen. This is the statistical/logical design only, submitted for review before any compute runs, matching the pattern that worked for `STEP3_METHOD_CONTRACT.md`.
