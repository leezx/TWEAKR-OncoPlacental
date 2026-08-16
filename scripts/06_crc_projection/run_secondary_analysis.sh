#!/bin/bash
#$ -N tweakr_secondary_analysis
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 4
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== Step 1: M11 scoring, 297,307-cell subset, N_PERM=500, 6 panels ==="
python3 scripts/06_crc_projection/m11_scoring_full.py \
    results/06_crc_projection/m11_scoring_full

echo "=== Step 2: revCSC-high cohort definition + developmental composition ==="
python3 scripts/06_crc_projection/secondary_analysis_composition.py \
    results/06_crc_projection/gut_scoring_full \
    results/06_crc_projection/secondary_analysis_composition

echo "=== Step 3: M11 x revCSC concordance (continuous correlation + enrichment) ==="
python3 scripts/06_crc_projection/secondary_analysis_concordance.py \
    results/06_crc_projection/m11_scoring_full \
    results/06_crc_projection/gut_scoring_full \
    results/06_crc_projection/secondary_analysis_composition \
    results/06_crc_projection/secondary_analysis_concordance
