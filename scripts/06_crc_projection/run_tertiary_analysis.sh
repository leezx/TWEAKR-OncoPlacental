#!/bin/bash
#$ -N tweakr_tertiary_analysis
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== Step 6 tertiary analysis: full-atlas revCSC-independent D/F/P composition ==="
python3 scripts/06_crc_projection/tertiary_analysis_composition.py \
    results/06_crc_projection/gut_scoring_full \
    results/06_crc_projection/tertiary_analysis_composition
