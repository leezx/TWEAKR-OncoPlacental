# HDMA RNA-assay gene-mapping gap: found and fixed

While building HDMA per-sample pseudobulk raw-count matrices (needed for
F-developmental's positive evidence, `docs/STEP4_STATISTICAL_DESIGN.md`
section 2), the first pseudobulk build (`build_hdma_pseudobulk.R`) failed
on all 7 organs.

## Finding

Step 2's already-merged gene-ID mapping deliverable
(`results/02_gene_id_mapping/canonical_feature_map/<organ>_canonical_feature_map.tsv`)
was built from `rownames(obj)` on each HDMA Seurat object — i.e. each
object's **default assay**. Direct inspection showed the default assay is
`decontX` (not `RNA`), which shares its (smaller) gene set with the `SCT`
assay — both are QC/normalization-filtered subsets of the raw gene space.

Pseudobulk must sum **true raw counts**, which live in the `RNA` assay's
`counts` layer — and that assay has *more* genes than the default assay:

| Organ | Default assay genes | RNA assay genes |
|---|---|---|
| Adrenal | 25,314 | 28,375 |
| Thyroid | — | 30,157 |
| Spleen | — | 34,351 |
| Thymus | — | 34,814 |
| Liver | — | 33,939 |
| Skin | — | 34,983 |
| StomachEsophagus | — | 34,423 |

(Default-assay counts for organs beyond Adrenal weren't re-derived here —
Adrenal's 25,314/28,375 relationship was checked directly and confirmed
the default set is a **strict subset** of the RNA-assay set, 0 genes lost;
same relationship assumed to hold structurally across organs since all 7
objects share the same processing pipeline, but only Adrenal was directly
verified before the fix.)

Every gene present in the RNA assay's raw counts but absent from the
existing canonical_feature_map caused the pseudobulk script's own
verify-before-trusting check to correctly abort rather than silently drop
genes (e.g. Adrenal: 3,061 genes in `RNA` counts had no map entry).

## Fix

1. **`scripts/02_gene_id_mapping/extract_rna_assay_gene_lists.sh`** — new
   script, mirrors the original `extract_all_organ_gene_lists.sh` but
   reads `rownames(GetAssayData(obj, assay="RNA", layer="counts"))`
   explicitly instead of the default-assay `rownames(obj)`. Run on Argos;
   found 186 ENSG IDs across the 7 organs not covered by the original
   gene-ID mapping's ENSG union.
2. Queried biotypes for those 186 new IDs via the existing
   `query_biotype.py` (Ensembl REST API), merged into
   `union_ensg_biotypes.tsv`.
3. **`scripts/02_gene_id_mapping/build_canonical_feature_map_rna_assay.py`**
   — new script, same mapping logic as the original
   `build_canonical_feature_map.py`, run against the RNA-assay gene lists.
   Output: `results/02_gene_id_mapping/canonical_feature_map_rna_assay/`
   — a **new, separate artifact**, keeping the original (already-reviewed
   and merged) `canonical_feature_map/` untouched rather than silently
   rewriting a prior deliverable.
4. **Consistency check, built into the new script**: for every
   `(organ, original_feature)` pair present in both the old and new maps,
   verified they resolve to the identical `canonical_symbol`. Result:
   **217,606 shared pairs checked, 0 mismatches** — confirms the new map
   is a clean superset extension of the old one, not a divergent
   remapping.
5. `build_hdma_pseudobulk.R` updated to read from the new
   `canonical_feature_map_rna_assay/` directory.

## Verification

- Raw-counts sanity check (already in `build_hdma_pseudobulk.R`, ran
  successfully before hitting the mapping gap): RNA assay's `counts` layer
  is confirmed integer-valued (fraction of sampled nonzero values
  integer = 1.0000) for every organ — true raw counts, not normalized.
- Gene-list extraction job (Argos job 3620434): all 7 organs OK, 0 failed.
- Mapping-build consistency check: 0/217,606 mismatches (see above).

This is the same "verify assumptions about which assay/layer a script
actually reads, don't assume alignment" discipline used earlier for the
Nature2026 raw-counts barcode-alignment check and the VentoTormo
raw-counts availability audit.
