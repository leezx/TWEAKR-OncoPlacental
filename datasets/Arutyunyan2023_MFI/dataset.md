# Dataset: Arutyunyan et al. 2023, Nature — Spatial multiomics map of trophoblast development in early pregnancy

- Paper: *Spatial multiomics map of trophoblast development in early pregnancy*
- Journal: Nature, 2023
- DOI: `10.1038/s41586-023-05869-0`
- PubMed: 36991123
- bioRxiv preprint: `10.1101/2022.11.06.515326`
- Lab: Vento-Tormo lab (Wellcome Sanger Institute)
- Role in project: independent early-placenta reference for `P1/P2/P3/D` module construction ([[TWEAKR-Worklog#Placenta数据集]] evidence layer, "第一梯队" tier — Nature 2023 Arutyunyan)

This single study deposited 4 ArrayExpress accessions covering 4 different modalities. This file is the master record for all 4.

## Accessions overview

| Accession | Modality | Samples | ENA project | Status here |
|---|---|---|---|---|
| E-MTAB-12421 | scRNA-seq / snRNA-seq, primary tissue (implantation sites, decidua, placenta) | 3 implantation-site snRNA + 5 decidual + 3 placental (scRNA/snRNA), 6–13 PCW | ERP144668 | ✅ downloaded (processed h5ad) |
| E-MTAB-12650 | scRNA-seq / snRNA-seq of primary trophoblast organoids (PTO) and trophoblast stem cells (TSC), incl. EVT differentiation | PTO + TSC + unstimulated organoid subsets | ERP144781 (+ EGAD00001010017, controlled-access raw) | ✅ downloaded (processed h5ad) |
| E-MTAB-12698 | 10x Visium spatial transcriptomics of maternal-fetal interface incl. myometrium | 8 Visium sections | (ArrayExpress Files, no separate ENA raw needed for processed output) | ✅ downloaded (Space Ranger processed output) |
| E-MTAB-12595 | 10x multiome snRNA-seq + snATAC-seq, primary tissue (implantation sites, decidua, placenta), 6–9 PCW | 3 samples × (RNA + ATAC) = 6 runs | ERP144790 | ⚠️ NOT downloaded — see "E-MTAB-12595 decision" below |

## Where the files are

### E-MTAB-12421 — primary tissue scRNA/snRNA (processed)

No processed matrix is attached to the ArrayExpress record itself (only `.idf.txt`/`.sdrf.txt` metadata there; raw FASTQ lives at ENA `ERP144668`). The paper's own data portal — **Reproductive Cell Atlas** (`reproductivecellatlas.org/mfi.html`) — hosts author-processed AnnData objects (`raw_counts` in `.raw`, normalized log counts in `.X`):

| File | Scope | Size | URL |
|---|---|---|---|
| `adata_all_donors_all_cell_states_..._UPD_20230307.h5ad` | **all donors, all cell states** (superset — includes trophoblast + decidual + immune + stromal + endothelial) | 12.38 GB | `https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/mfi/adata_all_donors_all_cell_states_raw_counts_in_raw_normlog_counts_in_X_for_download_UPD_20230307.h5ad` |
| `adata_all_donors_trophoblast_..._UPD_20230307.h5ad` | all donors, trophoblast-only subset | 3.69 GB | *(not downloaded — derivable from the file above; see decision note)* |
| `adata_P13_trophoblast_....h5ad` | single donor (P13) trophoblast only | 1.69 GB | *(not downloaded — subset of the above)* |

**Decision:** only the "all donors, all cell states" superset was pulled into `raw/`. The trophoblast-only and single-donor files are strict subsets of it and would just duplicate storage; if a trophoblast-only working copy is needed later, derive it locally into `processed/` (same pattern as `GSE131696`'s `placenta_trophoblast_filtered_list.csv`).

### E-MTAB-12650 — PTO / TSC organoid scRNA/snRNA (processed)

Also served from the Reproductive Cell Atlas portal (raw FASTQ is at ENA `ERP144781`, and some runs are controlled-access at EGA `EGAD00001010017` — not needed since processed matrices are open):

| File | Scope | Size | URL |
|---|---|---|---|
| `Organoid_PTO_cellxgene.h5ad` | primary trophoblast organoids (PTO) | 1.56 GB | `https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/mfi/Organoid_PTO_cellxgene.h5ad` |
| `Organoid_TSC_cellxgene.h5ad` | trophoblast stem cells (TSC) | 0.79 GB | `https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/mfi/Organoid_TSC_cellxgene.h5ad` |
| `adata_Fig3_trophoblast_organoids_unstimulated_....h5ad` | unstimulated organoids, Fig. 3 subset | 2.01 GB | `https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/mfi/adata_Fig3_trophoblast_organoids_unstimulated_raw_counts_in_rawX_normlog_counts_in_X_for_download.h5ad` |

All three downloaded — they are not simple subsets of one another (PTO vs TSC are different in vitro systems; Fig3 is a specific stimulation-state slice used in the paper).

### E-MTAB-12698 — Visium spatial transcriptomics (processed)

Unlike 12421/12650, ArrayExpress itself hosts the processed Space Ranger output directly under `Files/` (no portal needed):

Base URL: `https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/698/E-MTAB-12698/Files/`

| Sample | `*_spaceranger_output.tar.gz` | Raw slide image (`.ndpi`/`.tif`) |
|---|---|---|
| Pla_Camb9518737 | 278 MB | 216 MB `.ndpi` + 20 MB `.tiff` |
| Pla_HDBR9518710 | 189 MB | 115 MB `.ndpi` + 19 MB `.tiff` |
| WS_PLA_S9101764 | 324 MB | 40 MB `.ndpi` + 20 MB `.tif` |
| WS_PLA_S9101765 | 301 MB | 44 MB `.ndpi` + 21 MB `.tif` |
| WS_PLA_S9101766 | 212 MB | 38 MB `.ndpi` + 20 MB `.tif` |
| WS_PLA_S9101767 | 325 MB | 45 MB `.ndpi` + 21 MB `.tif` |
| WS_PLA_S9101769 | 271 MB | 41 MB `.ndpi` + 21 MB `.tif` |
| WS_PLA_S9101770 | 272 MB | 39 MB `.ndpi` + 20 MB `.tif` |

**Decision:** downloaded all 8 `*_spaceranger_output.tar.gz` (filtered/raw feature-barcode matrices + spatial coordinates + hi-res image already embedded in Space Ranger output; ≈2.17 GB total). Skipped the full-resolution `.ndpi`/`.tif` whole-slide scans (≈700 MB total) — not needed for expression/signature analysis, can be pulled later only if a publication-quality histology figure is needed.

The Reproductive Cell Atlas portal (`invivo-spatial.html`) also serves smaller pre-packaged `Visium_spatial_ID_*_raw.h5ad` per sample (30–64 MB each, ≈420 MB total) if a quick AnnData load without re-parsing Space Ranger output is preferred later.

### E-MTAB-12595 — multiome snRNA-seq + snATAC-seq — decision needed

No processed matrix for this accession is on ArrayExpress `Files/` (only `.idf.txt`/`.sdrf.txt`) and none is on the Reproductive Cell Atlas portal either (checked `mfi.html`, `invivo-spatial.html`, and the `ventolab/MFI` GitHub repo — no multiome download link found anywhere public).

The only available copy is **raw FASTQ at ENA `ERP144790`**, 6 runs (3 snATAC + 3 snRNA, one triplet per donor: `Pla_HDBR10084192/3/4` for ATAC, `Pla_HDBR10142863/4/5` for RNA):

| Run | Assay | Sample | Approx. size |
|---|---|---|---|
| ERR10855815/16/17 | ATAC-seq | Pla_HDBR10084192/3/4 | ~31.5 / 31.5 / 33.1 GB each |
| ERR10855818/19/20 | RNA-seq | Pla_HDBR10142863/4/5 | ~66.8 / 61.1 / 75.2 GB each |

**Total ≈ 299 GB raw FASTQ**, and it would still need local `cellranger-arc count` processing (needs a reference genome + GTF + a machine with substantial RAM/CPU) before it's usable — this is a much bigger commitment than every other file in this dataset combined (which totals ~19 GB).

**Not downloaded yet.** Flagging for your decision:
1. Skip E-MTAB-12595 for now — the other 3 accessions already give RNA-level primary tissue, organoid, and spatial coverage; multiome ATAC mainly adds chromatin-accessibility evidence (useful for the "regulatory support" evidence layer, not required for the RNA-level P1/P2/P3/D signature construction in Aim 1).
2. Download the ~300 GB raw FASTQ and run `cellranger-arc` locally (needs GRCh38 arc-reference + compute planning — a separate task, not a quick follow-on).
3. Keep checking periodically — authors sometimes add processed multiome matrices to the portal after initial publication (the RNA-only files were tagged `UPD_20230307`, i.e. updated after the original release).

Default assumption until told otherwise: **option 1 (skip)**. The dataset directory has been scaffolded as a skeleton so it's easy to resume.

## Cell type / annotation summary (from `all_donors_all_cell_states` obs metadata — to be filled in after first load)

*(TODO: once the h5ad is loaded once, paste the `obs` column names and `celltype` value counts here, same pattern as `DATA/scRNAseq/2026_human_maternal_fetal_Nature/raw/README.md`.)*

## Storage locations (local)

- `DATA/scRNAseq/Arutyunyan2023_MFI/raw/` — E-MTAB-12421 (`all_donors_all_cell_states`) + E-MTAB-12650 (PTO/TSC/Fig3) h5ad files
- `DATA/SpatialTranscriptomics/Arutyunyan2023_MFI_Visium/raw/` — E-MTAB-12698 Space Ranger outputs (8 samples)
- `DATA/scATACseq/Arutyunyan2023_MFI_multiome/` — skeleton only, E-MTAB-12595 not downloaded (see decision above)

## Source links

- Reproductive Cell Atlas portal: https://www.reproductivecellatlas.org/mfi.html
- Spatial data page: https://www.reproductivecellatlas.org/invivo-spatial.html
- ArrayExpress E-MTAB-12421: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12421
- ArrayExpress E-MTAB-12595: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12595
- ArrayExpress E-MTAB-12698: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12698
- ArrayExpress E-MTAB-12650: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12650
- Analysis code (paper's own): https://github.com/ventolab/MFI
