#!/bin/bash
#$ -N tweakr_gut_primary_analysis_only
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 4
set -euo pipefail
source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex
cd /home/zz950/TWEAKR-OncoPlacental
echo "=== Re-run primary analysis only (robustness definition fix, PR #27 round 2) ==="
python3 scripts/06_crc_projection/crc_gut_scoring_primary_analysis.py \
    results/06_crc_projection/gut_scoring_full \
    results/06_crc_projection/gut_scoring_primary_analysis
