# Canonical feature map (RNA assay) — build summary

Extends `../canonical_feature_map/` (built from each Seurat object's default assay, decontX/SCT) to cover the RNA assay's full raw-counts feature space, needed for pseudobulk. See module docstring for the gap this fixes.

Consistency check against the original map: 217606 shared (organ, original_feature) pairs checked, 0 mismatches (CLEAN — every gene the two maps share got the same canonical_symbol).

Per-organ mapping_status counts (summed across all 7 organs):

| mapping_status | count |
|---|---|
| native_symbol | 160440 |
| unmapped_kept_as_ensembl_id | 63814 |
| mapped_via_hgnc | 6743 |
| mapped_via_ensembl_display_name | 45 |

Total symbol collisions across all organs: 33 (see `collision_report.tsv`).

Collisions are reported, not resolved, by this script. Pseudobulk construction (build_hdma_pseudobulk.R) sums counts across colliding features by canonical_symbol.
