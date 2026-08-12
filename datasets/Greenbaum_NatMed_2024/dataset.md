# Dataset: Greenbaum et al. 2024, Nature Medicine — Spatial multiomic landscape of the human placenta at molecular resolution

- Journal: Nature Medicine, 2024
- Role in project: "第一梯队" evidence layer ([[TWEAKR-Worklog#Placenta数据集]]) — RNA + chromatin accessibility + spatial evidence for the same genes, feeds the "Placenta evidence score" concept from `2026-GPT-TWEAKR-Oncofetal.md`

**How this differs from the other 3 entries in this repo:** downloaded manually by the user (not investigated/pulled by this session's automated workflow like `Arutyunyan2023_MFI`). This file documents what's there after the fact.

## Source

- Broad Single Cell Portal: https://singlecell.broadinstitute.org/single_cell/study/SCP2601
- 8 donors, 6–11 weeks post-conceptional age
- Paired snRNA + snATAC (10x multiome) + 3 complementary spatial approaches: Slide-tags, STARmap-ISS, STARmap-ISH

## What's on disk

`DATA/scRNAseq/Greenbaum_NatMed_2024/raw/SCP2601/` (~7.5GB), original Broad SCP download structure preserved:

| Path | Content |
|---|---|
| `10x_genes_file/`, `10x_barcodes_file/` | RNA gene/barcode lists |
| `expression/665be9bd512b6f7d793fa94b/` (2.6GB) | 10x RNA matrix: `matrix.mtx` + `genes_rna2.tsv` + `barcodes_rna2.tsv` |
| `expression/665bea28aa3336975ba927f2/` (4.9GB) | 10x ATAC matrix: `matrix_atac.mtx` + `peakAnnotations_atac.tsv` |
| `metadata/metadata.csv` | cell-level metadata |
| `other/humanplacenta_spatial.csv` | spatial coordinates |
| `other/humanplacenta_cluster.csv` | cluster assignments |
| `other/humanplacenta_expression.csv.gz` / `_raw.csv.gz` | normalized / raw expression |
| `other/humanplacenta_genescore.csv.gz` | gene activity scores (from ATAC) |
| `other/humanplacenta_motifzscore.csv.gz` | TF motif z-scores (chromVAR-style) — useful for the TEAD/YAP regulon check mentioned in the KB discussion's "第三层" evidence requirement |
| `file_supplemental_info.tsv` | SCP's own file-type manifest |

## Design note

Kept as one `scRNAseq` dataset_id rather than split across `scATACseq`/`SpatialTranscriptomics` like `Arutyunyan2023_MFI` — this came as a single unified Broad SCP delivery (one accession, one directory tree from the portal), not separate per-modality accessions from the source, so splitting it would fight the source structure rather than reflect it.

## Full record

`DATA/scRNAseq/Greenbaum_NatMed_2024/link.md`
