# Step 4a design: re-anchoring F-developmental on real fetal gut/colon data

Builds on `docs/STEP4_DFP_DESIGN.md` / `docs/STEP4_STATISTICAL_DESIGN.md` (the frozen D/F/P logic and its statistical design) and directly addresses a scope mismatch the user flagged after reviewing the mike_verzi enrichment result: HDMA's 7-organ F-developmental reference (Adrenal/Liver/Skin/Spleen/StomachEsophagus/Thymus/Thyroid — `datasets/HumanDevelopmentMultiomicAtlas/dataset.md`) has **zero gut/colon representation**, which is scientifically mismatched for a CRC project. This doc defines `F_Colon-developmental` and `F_SI-developmental` — built from real human fetal gut single-cell data — as the new primary F-axis for CRC-facing work, before any qsub compute runs, matching the review-before-compute discipline that produced `STEP4_DFP_DESIGN.md`.

## What changes and what doesn't

- **P-developmental is unchanged.** It was already computed as a whole-body (not organ-scoped) program — `replicated_in_placenta(gene)` AND `adult_excluded(gene, whole_body)` — so it needs no gut-specific rebuild. It's reused as-is from the frozen Step 4 output (`results/04_dfp_signature/dfp_gene_sets/`).
- **F-developmental is rebuilt**, replacing HDMA's 7 pan-fetal organs with two real fetal gut regions as the primary reference:
  - `F_Colon-developmental` (large intestine) — **primary**, matches CRC's tissue of origin directly.
  - `F_SI-developmental` (small intestine) — **secondary, parallel**, run with the identical method for cross-region comparison.
  - A pooled `F_Gut-developmental` ("Gut-core" = Colon ∪ SI, or Colon ∩ SI as the stricter alternative) is reported only as a **tertiary summary**, and only if Colon and SI turn out concordant enough for pooling to be meaningful — not assumed up front.
- **D/F/P set operations are redefined** using the new gut-specific F, per the user's explicit formulas:
  - `D_Gut-shared = F_Gut-developmental ∩ P_developmental`
  - `F_Gut-specific = F_Gut-developmental \ P_developmental`
  - `P_Gut-specific = P_developmental \ F_Gut-developmental`
  - (Colon and SI versions computed the same way against the same unchanged `P_developmental`, before any Colon/SI pooling decision.)
- **The existing pan-organ HDMA-based D/F/P is not deleted.** It's demoted to a secondary "pan-fetal / non-gut developmental validation framework" — useful for asking "is this gene generically fetal, anywhere in the body" but no longer the primary CRC coordinate system.

## Reference data: Gut Cell Atlas fetal object (Elmentaite et al., Nature 2021)

Downloaded and verified this session (`scripts/04a_dfp_gut/download_gutcellatlas_fetal.sh`, `scripts/04a_dfp_gut/inventory_gutcellatlas_fetal.py`), both files byte-exact and md5-exact against source:

| File | Cells | Genes | Role |
|---|---|---|---|
| `fetal_RAWCOUNTS_cellxgene.h5ad` | 62,849 | 33,694 | Raw counts (verified integer, `X_integer_fraction=1.0`) — used for pseudobulk construction |
| `final_fetal_object_cellxgene.h5ad` | 62,849 | 26,757 (HVG subset) | Normalized, has `X_pca`/`X_umap` — metadata cross-check only, not used for counts |

Real structure confirmed by direct inventory, not assumed:

- **Gene IDs**: symbols in `var_names`, verified real Ensembl IDs in `var['gene_ids']` (33,694 unique `ENSG...` IDs, 1:1 with 33,694 unique symbols, no collisions) — same dual-ID discipline as HDMA/placental datasets, enables clean cross-dataset alignment.
- **Region** (`obs['Organ']`, 3 categories): `duojejunum` (duodenum+jejunum combined) 21,592 cells, `ileum` 20,110 cells, `colon` 21,147 cells. SI = duojejunum ∪ ileum; LI = colon.
- **Lineage** (`obs['cell_type_group']`, 4 categories): `epithelium` 16,937 / `mesenchymal` 40,671 / `vasculature` 3,274 / `immune` 1,967. Finer `cell_name_detailed` (28 categories) includes `Colonocytes`, `Small intestinal Epi`, `Enterocyte`, `Early/G2M/S enterocytes`, `LGR5 stem`, `BEST4+ enterocyte`, `Secretory Epi`.
- **Donor/sample structure** (real replicate structure, checked directly — same discipline as Step 4's HDMA per-organ-per-sample table):

  | Region | Epithelium-only cells | Donors contributing | Samples |
  |---|---|---|---|
  | duojejunum | 7,793 | 9/9 | 9 |
  | ileum | 6,302 | 9/9 | 10 |
  | colon | 2,842 | 9/9 | 9 |

  All 9 donors contribute epithelial cells to all 3 regions — no region is donor-confounded. Colon epithelium is the smallest (2,842 cells; thinnest donor contributes 64 cells) but still has full donor coverage, matching the Step 4 precedent of proceeding with organs that have 3–7 samples (Thyroid/Spleen/Thymus/Skin had only 3 each).
- **Age**: `PCW` (F6.1–F10.2) and `CRL` (crown-rump length, 9 categories) both present as parallel age proxies — same 9 donors, consistent with the previously-verified "6-11 PCW" description of this atlas.

## Primary population: epithelial lineage

Per the user's explicit direction, `cell_type_group == 'epithelium'` is the primary subset for the fetal-gut side of this analysis — it matches the CRC malignant-epithelial-cell-state focus the whole D/F/P signature exists to serve. Mesenchymal/vasculature/immune populations in this atlas are not used for `F_Colon-developmental`/`F_SI-developmental` construction (they may be useful later for stroma-context questions, out of scope here).

## Definitions (direct port of `STEP4_DFP_DESIGN.md` / `STEP4_STATISTICAL_DESIGN.md`'s F-developmental logic, gut regions substituted for HDMA organs)

Same asymmetry as before: this atlas, like HDMA, is a **single population per region** (no internal non-fetal-somatic contrast), so F-developmental's positive evidence uses the same "consistency across the region's own individual samples" design as Step 4's HDMA arm — not a fresh statistical framework.

- **`elevated_in_fetal_gut_epithelium(gene, region)`**: gene's expression in per-sample epithelial pseudobulk for that region clears a within-region percentile threshold, computed from the region's own sample distribution — same `elevated_pct` logic as HDMA (frozen at 75 in PR #12, to be re-validated against this dataset's own real distribution before reuse, not assumed to transfer unchanged).
  - **Pseudobulk construction**: per individual `Sample`, restricted to `cell_type_group == 'epithelium'`, sum raw counts, CPM-normalize within-region (never cross-region/cross-dataset raw magnitude, per `STEP3_METHOD_CONTRACT.md`'s standing boundary).
  - Region = `colon` for `F_Colon-developmental`; region = `duojejunum ∪ ileum` (pooled sample set) for `F_SI-developmental` — with duojejunum-vs-ileum internal concordance checked and reported (same "is a sub-split needed" question Step 4 resolved for StomachEsophagus, not assumed answered).
- **`adult_excluded(gene, matched_organ)`**: unchanged definition from `STEP4_STATISTICAL_DESIGN.md` §3 — gene's expression in the matched adult reference falls below that dataset's own internal percentile, `role == adult_negative_reference` rows only.
  - **Colon-matched adult reference**: GTEx v11 whole-body file (already in-project, `results/... GTEx_v11_median_tpm`) columns `Colon_Sigmoid`, `Colon_Transverse` (+`_Mucosa`/`_Muscularis`), plus HPA `colon`/`rectum` rows — same references already used for the CRC sanity check in Step 4.
  - **SI-matched adult reference**: GTEx v11's `Small_Intestine_Terminal_Ileum` (+`_Lymphoid_Aggregate`/`_Mixed_Cell`) column, plus HPA's `small intestine` row — **confirmed present in the already-downloaded GTEx v11 whole-body file** (verified this session via direct column-header check; this atlas contains all 68 GTEx tissues, not just the 4 columns previously pulled for the colon sanity check). **No new adult-reference download is needed for the SI arm** — this corrects an earlier assumption (surfaced when the gut-data gap was first discussed) that only colon-matched adult references were available in-project.
  - Tabula Sapiens `TS_Large_Intestine.h5ad` (already downloaded) remains available as a **Tier-2 independent validation** reference per the existing GTEx-v11 `link.md` design note — not used to define the signature, consistent with the project's standing adult-reference tiering.
- **`F_Colon-developmental(gene)` = `elevated_in_fetal_gut_epithelium(gene, colon)` AND `adult_excluded(gene, colon-matched)`**
- **`F_SI-developmental(gene)` = `elevated_in_fetal_gut_epithelium(gene, duojejunum ∪ ileum)` AND `adult_excluded(gene, SI-matched)`**

Both are pre-corrected for adult expression before any set operation, same as the HDMA arm — no bespoke "exclude against P" step, `D_Gut-shared` falls out of the set intersection by construction.

## Non-circularity boundary (explicit, per user directive)

The 5 `mike_verzi_fetal_signature.gmt` mouse gene sets (YAP_SIGNALING_GENES, REVIVAL_STEM_CELL_GENES, FETAL_SPHEROID_EPITHELIUM_GENES, REGENERATIVE_EPITHELIUM, FETAL_INTESTINE_GENES) **must not be used anywhere in `F_Colon-developmental`/`F_SI-developmental` construction** — no gene-selection, no threshold-tuning, no candidate filtering. Correct causal order: build `F_Colon-developmental`/`F_SI-developmental` independently from the Gut Cell Atlas data first; only afterward use the (fully independent, human-orthology-verified) mike_verzi signatures as external validation, per the three-layer statistical validation plan below. This mirrors the primary/extended-set-ordering lesson from PR #17 round 3 (Ly6a) — outcome-dependent evidence must never influence which genes make it into a primary definition.

## What this design does NOT do yet

No qsub jobs, no gene lists, no thresholds chosen — same discipline as `STEP4_DFP_DESIGN.md`. Explicitly left open, to decide against real data:

- Whether `elevated_pct=75` / `adult_excl_pct=25` (frozen for HDMA in PR #12) transfer unchanged to this dataset's own distribution, or need independent calibration — the marker-panel-check discipline from PR #12 (DLK1/IGF2/H19/LIN28B/PEG10 + gut-appropriate markers, e.g. LGR5/OLFM4 for crypt stem, MUC2 for goblet) should be re-run against this dataset before locking, not assumed to transfer.
- Whether duojejunum and ileum need to stay split (like HDMA's per-organ F-developmental) with a union/consensus reconciliation step, or can be pooled directly into one SI sample set — decide once the two regions' real elevated-gene overlap is visible (same "let real distribution decide" precedent as Step 4's F-specific union-vs-consensus resolution).
- Exact quorum for "detected in the region's own samples" (Step 4 used 0.5 vs 1.0, found it made almost no difference vs. the percentile cutoff) — to be re-checked, not assumed to transfer.
- Whether pooled `F_Gut-developmental` (Colon ∪ SI or ∩) is reported at all — contingent on Colon/SI concordance once both are built.
- The three-layer statistical validation (hypergeometric enrichment with CI/OR, preranked GSEA against a continuous fetal-vs-adult-colon differential statistic, size-and-expression-matched permutation nulls) against the 5 mike_verzi signatures — this is the next design doc after `F_Colon-developmental`/`F_SI-developmental` gene sets exist, not part of this one.
- Donor/dataset-provenance overlap audit between Gut Cell Atlas and the two other candidate fetal-gut datasets (GSE158702, GSE95630/GSE103239) — needed before any future "independent replication" claim across those atlases; not required to build `F_Colon-developmental` itself since only Gut Cell Atlas is used here.
