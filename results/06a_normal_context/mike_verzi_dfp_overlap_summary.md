# mike_verzi (normal-context fetal/revival) x D/F/P gene-overlap audit (v2, post case-fix)

v2: fixes PR #19 round-1 blocker -- v1's mouse->human mapping silently dropped genes whose GMT symbol case didn't exactly match BioMart's canonical mouse Gene name (e.g. Col1A1/Tnfrsf12A/Tmsb4X/Ly6A and bulk S100A*/Slc* families), misclassifying real orthologs as NOT_FOUND_IN_BIOMART. Fixed via a full unfiltered mouse-gene BioMart pull + explicit exact-then-case-insensitive local join (mike_verzi_symbol_resolution.tsv has the full raw->canonical->ensembl provenance for every gene).

## Signature -> human ortholog resolution

| Signature | n mouse genes | n primary (one2one) human | % resolved |
|---|---|---|---|
| `YAP_SIGNALING_GENES` | 444 | 340 | 76.6% |
| `REVIVAL_STEM_CELL_GENES` | 236 | 211 | 89.4% |
| `FETAL_SPHEROID_EPITHELIUM_GENES` | 270 | 231 | 85.6% |
| `REGENERATIVE_EPITHELIUM` | 203 | 153 | 75.4% |
| `FETAL_INTESTINE_GENES` | 1258 | 1147 | 91.2% |

## D/F/P overlap (primary human gene sets)

| Signature | D-shared (6) | F-specific global (2504) | P-specific (78) |
|---|---|---|---|
| `YAP_SIGNALING_GENES` | 0 | 30 | 0 |
| `REVIVAL_STEM_CELL_GENES` | 0 | 10 | 0 |
| `FETAL_SPHEROID_EPITHELIUM_GENES` | 0 | 24 | 0 |
| `REGENERATIVE_EPITHELIUM` | 0 | 9 | 0 |
| `FETAL_INTESTINE_GENES` | 0 | 135 | 1 |

## D/F/P overlap, F-lineage modules (organ ∩ F-specific)

| Signature | Adrenal (680) | Liver (797) | Skin (1123) | Spleen (1087) | Stomach (749) | Thymus (1189) | Thyroid (682) |
|---|---|---|---|---|---|---|---|
| `YAP_SIGNALING_GENES` | 6 | 3 | 7 | 7 | 4 | 21 | 6 |
| `REVIVAL_STEM_CELL_GENES` | 1 | 1 | 2 | 2 | 1 | 9 | 0 |
| `FETAL_SPHEROID_EPITHELIUM_GENES` | 8 | 4 | 6 | 10 | 4 | 18 | 5 |
| `REGENERATIVE_EPITHELIUM` | 1 | 0 | 2 | 3 | 1 | 5 | 0 |
| `FETAL_INTESTINE_GENES` | 23 | 19 | 31 | 44 | 22 | 96 | 11 |
