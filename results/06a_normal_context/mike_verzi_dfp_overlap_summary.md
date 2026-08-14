# mike_verzi (normal-context fetal/revival) x D/F/P gene-overlap audit

Per Step 6a design: before any scoring, check whether each of the 5 independently-published normal-tissue fetal/revival/regeneration gene sets (mouse->human orthology-mapped, primary = Compara-confirmed one-to-one only) shares genes with D-shared, F-specific (global + 7 lineage modules), or P-specific.

## Signature -> human ortholog resolution

| Signature | n mouse genes | n primary (one2one) human | % resolved |
|---|---|---|---|
| `YAP_SIGNALING_GENES` | 444 | 282 | 63.5% |
| `REVIVAL_STEM_CELL_GENES` | 236 | 169 | 71.6% |
| `FETAL_SPHEROID_EPITHELIUM_GENES` | 270 | 189 | 70.0% |
| `REGENERATIVE_EPITHELIUM` | 203 | 132 | 65.0% |
| `FETAL_INTESTINE_GENES` | 1258 | 927 | 73.7% |

## D/F/P overlap (primary human gene sets)

| Signature | D-shared (6) | F-specific global (2,504) | P-specific (78) |
|---|---|---|---|
| `YAP_SIGNALING_GENES` | 0 | 24 | 0 |
| `REVIVAL_STEM_CELL_GENES` | 0 | 6 | 0 |
| `FETAL_SPHEROID_EPITHELIUM_GENES` | 0 | 21 | 0 |
| `REGENERATIVE_EPITHELIUM` | 0 | 9 | 0 |
| `FETAL_INTESTINE_GENES` | 0 | 108 | 1 |

## D/F/P overlap, F-lineage modules (organ ∩ F-specific)

| Signature | Adrenal (680) | Liver (797) | Skin (1123) | Spleen (1087) | Stomach (749) | Thymus (1189) | Thyroid (682) |
|---|---|---|---|---|---|---|---|
| `YAP_SIGNALING_GENES` | 6 | 3 | 5 | 6 | 3 | 18 | 3 |
| `REVIVAL_STEM_CELL_GENES` | 1 | 0 | 1 | 1 | 0 | 6 | 0 |
| `FETAL_SPHEROID_EPITHELIUM_GENES` | 6 | 3 | 3 | 8 | 3 | 16 | 3 |
| `REGENERATIVE_EPITHELIUM` | 1 | 0 | 2 | 3 | 1 | 5 | 0 |
| `FETAL_INTESTINE_GENES` | 21 | 15 | 23 | 38 | 18 | 79 | 7 |
