#!/bin/bash
#$ -N tweakr_gutatlas_epi_lifespan_dl
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -uo pipefail

# Per PR #20 round-1 REQUEST_CHANGES: F_Colon-developmental/F_SI-developmental
# should be built as a same-atlas fetal-vs-adult epithelial DE, not a
# mechanical port of HDMA's single-population + external-adult-percentile
# design. epi_raw_counts02_v2.h5ad is the Space-Time Gut Cell Atlas's
# full-lifespan (fetal + pediatric + adult) epithelium-only raw-counts
# object -- the file this primary construction needs. Real size/md5
# verified via source HEAD request this session (not assumed):
#   content-length: 666604440
#   md5 (x-amz-meta-s3cmd-attrs): 2a149b8cf04567569707e9d1fab27209

OUT_DIR=/home/zz950/DATA/scRNAseq/GutCellAtlas_Elmentaite2021/raw
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "Downloading epi_raw_counts02_v2.h5ad (full-lifespan epithelium, fetal+pediatric+adult)..."
curl -sS -o epi_raw_counts02_v2.h5ad "https://cellgeni.cog.sanger.ac.uk/gutcellatlas/epi_raw_counts02_v2.h5ad"

echo "--- md5sum ---"
md5sum epi_raw_counts02_v2.h5ad

echo "--- file size ---"
ls -la epi_raw_counts02_v2.h5ad

echo "--- expected ---"
echo "size=666604440 md5=2a149b8cf04567569707e9d1fab27209"
