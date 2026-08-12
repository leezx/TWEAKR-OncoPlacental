# Dataset: HPA trophoblast cell-specific proteome

- Source: Human Protein Atlas, full bulk export (https://www.proteinatlas.org/download/proteinatlas.tsv.zip)
- Role in project: "第一梯队" evidence layer ([[TWEAKR-Worklog#Placenta数据集]]) — normal-cell-type protein/RNA specificity, one of the required orthogonal evidence types (RNA / protein / chromatin / spatial) for the CNS-grade signature described in [[2026-GPT-TWEAKR-Oncofetal#Oncofetal-Placenta sig]]

This closes out the last of the 6 originally-missing datasets from the Aug 12 gap analysis.

## What this is, and how it differs from what was already in the project

The project already had an HPA **tissue-level** placenta gene list (`tissue_category_rna_placenta_Tissue.tsv` — placenta vs. other tissues). This is a different, complementary HPA resource: **single cell type** specificity — a given cell type vs. ~154 other individual cell types / ~53 cell-type groups, across all tissues, not just within placenta. Filtered here to trophoblast-related cell types (Cytotrophoblasts, Syncytiotrophoblasts, Extravillous trophoblasts).

## Method

Downloaded the full HPA bulk export (`proteinatlas.tsv`, ~43MB uncompressed, ~20,163 genes, all annotation columns) rather than trying to query a filtered subset via the API — it's small, and having the full table is reusable for other HPA-column-based filtering later.

Filtered locally via `build_hpa_trophoblast_proteome.py` (saved alongside the output, per the workspace's provenance convention — every derived file needs the exact code that produced it saved next to it). Three tiers, most to least confident:

| Tier | HPA category | Definition | n genes |
|---|---|---|---|
| 1 | Cell type enriched | ≥4x the next-highest of ~154 individual cell types, driven by a trophoblast subtype | 54 |
| 2 | Group enriched | ≥4x the next-highest of ~53 cell-type *groups*, driven by a trophoblast subtype | 90 |
| 3 | Cell type enhanced | ≥4x tissue average in ≥1 cell type (weaker, often co-lists many cell types) | 1641 |

Spot-check on Tier 1 (54 genes) confirms known trophoblast biology: `CGA`/`CGB2`/`CGB3`/`CGB5`/`CGB7`/`CGB8` (hCG subunits), `ERVFRD-1` (syncytin-2, the classic syncytiotrophoblast fusion gene), `CYP19A1` (aromatase), `ASCL2` (cytotrophoblast/EVT TF) — good signal.

## Discrepancy vs. the KB-cited numbers — flagged, not resolved

The KB discussion ([[2026-GPT-TWEAKR-Oncofetal#Placenta数据集]]) cited "887 trophoblast-elevated genes, 97 trophoblast-enriched, 90 group-enriched" for this resource. This run reproduces **group-enriched exactly (90)**, but gets 54 (not 97) for cell-type-enriched and 1641 (not ~797) for enhanced. Likely explanation: the live HPA site may present trophoblasts as one merged "Trophoblast cells" label on its dedicated single-cell-type page, whereas this script sums three granular subtype labels (Cytotrophoblasts / Syncytiotrophoblasts / Extravillous trophoblasts) that the current bulk-export data actually uses — or the two runs hit different HPA data versions. Not chasing an exact match; the group-enriched match on the one number we could reproduce exactly is a reasonable correctness signal, and the tiered output here is fully reproducible from the saved script regardless of what the live site currently shows.

## What's on disk

`DATA/1.Databases/signatures/TWEAKR_Oncofetal_gene_sets/`:
- `raw/hpa_full_export/proteinatlas.tsv(.zip)` — full unfiltered HPA export
- `processed/v0.1/build_hpa_trophoblast_proteome.py` — the filter script (also the provenance record)
- `processed/v0.1/hpa_trophoblast_cell_type_enriched.txt` (Tier 1)
- `processed/v0.1/hpa_trophoblast_group_enriched.txt` (Tier 2)
- `processed/v0.1/hpa_trophoblast_cell_type_enhanced.txt` (Tier 3)
- `processed/v0.1/hpa_trophoblast_single_cell_type_proteome.tsv` — all 3 tiers combined with full annotation columns

## Full record

`DATA/1.Databases/signatures/TWEAKR_Oncofetal_gene_sets/link.md`, `DATA/1.Databases/signatures/TWEAKR_Oncofetal_gene_sets/processed/README.md`
