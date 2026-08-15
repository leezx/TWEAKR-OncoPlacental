#!/bin/bash
#$ -N tweakr_gut_scoring_gene_sets
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== Build 13-panel gut scoring gene-set inventory (bare-Ensembl-ID) ==="
python3 scripts/06_crc_projection/build_gut_scoring_gene_sets.py results/06_crc_projection/gut_scoring
