# Dataset: HPA RNA tissue consensus

- Source: Human Protein Atlas (proteinatlas.org), reused from an existing local copy in `DATA/CRC-Atlas` (same public file, not re-downloaded)
- Role in project: second half of Tier-1 adult normal reference for Step 3 D/F/P signatures, alongside `[[GTEx_v11_median_tpm]]`. Fills GTEx's one gap — GTEx has no Thymus tissue at all (donor population skews older, thymus involutes with age); HPA's tissue panel does include thymus.

## What's on disk

`DATA/1.Databases/HPA_RNA_tissue_consensus/`:
- `raw/rna_tissue_hpa.tsv.zip` (6,653,963 bytes; md5-verified byte-exact copy of `DATA/CRC-Atlas/phase2/03_data/raw/HPA_normal_tissue/rna_tissue_hpa.tsv.zip`)
- `processed/v0.1/hdma_organ_to_hpa_tissue_map.tsv` — HDMA-organ ↔ HPA-tissue-name mapping (all 7 organs + colon/rectum + a placenta cross-check row confirmed present)
- Also pushed to Argos (`~/DATA/1.Databases/HPA_RNA_tissue_consensus/`)

## Verification

- md5 of the copy matches the CRC-Atlas project's original file exactly.
- `unzip -t` passes; long-format table is 20,162 genes × 40 tissues = 806,480 rows, matches the file's actual row count exactly.
- All 7 HDMA organ names (+ colon, rectum, placenta) confirmed present among the 40 tissue values by direct lookup.

## Full record

`DATA/1.Databases/HPA_RNA_tissue_consensus/link.md`
