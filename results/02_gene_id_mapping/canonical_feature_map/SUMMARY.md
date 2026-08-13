# Canonical feature map — build summary

Per-organ mapping_status counts (summed across all 7 organs):

| mapping_status | count |
|---|---|
| native_symbol | 154046 |
| unmapped_kept_as_ensembl_id | 57283 |
| mapped_via_hgnc | 6236 |
| mapped_via_ensembl_display_name | 41 |

Total symbol collisions across all organs: 29 (see `collision_report.tsv` for the full list — one row per `(organ, canonical_symbol)` pair with >1 original feature).

Collisions are reported, not resolved, by this script. Step 3 pseudobulk construction must pick an explicit aggregation rule (e.g. sum counts across colliding features) before using canonical_symbol as the feature key.
