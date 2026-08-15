#!/bin/bash
#$ -N tweakr_gutatlas_fetal_dl
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -uo pipefail

OUT_DIR=/home/zz950/DATA/scRNAseq/GutCellAtlas_Elmentaite2021/raw
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "Downloading fetal_RAWCOUNTS_cellxgene.h5ad..."
curl -L -o fetal_RAWCOUNTS_cellxgene.h5ad "https://cellgeni.cog.sanger.ac.uk/gutcellatlas/fetal_RAWCOUNTS_cellxgene.h5ad"
echo "Downloading final_fetal_object_cellxgene.h5ad (normalized, for metadata/annotation cross-check)..."
curl -L -o final_fetal_object_cellxgene.h5ad "https://cellgeni.cog.sanger.ac.uk/gutcellatlas/final_fetal_object_cellxgene.h5ad"

echo "--- md5sums ---"
md5sum fetal_RAWCOUNTS_cellxgene.h5ad final_fetal_object_cellxgene.h5ad

echo "--- file sizes ---"
ls -la fetal_RAWCOUNTS_cellxgene.h5ad final_fetal_object_cellxgene.h5ad
