#!/bin/bash
#$ -N tweakr_gut_scoring_full
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 8
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== Full 665,473-cell null-calibrated scoring, 13 panels ==="
echo "Reads results/06_crc_projection/gut_scoring_convergence_check/nperm500_required_panels.txt as the N_PERM gate"
python3 scripts/06_crc_projection/crc_gut_scoring_full.py \
    results/06_crc_projection/gut_scoring_full \
    results/06_crc_projection/gut_scoring_convergence_check

echo ""
echo "=== Primary analysis: donor/study-aware revCSC<->D/F/P correlation ==="
python3 scripts/06_crc_projection/crc_gut_scoring_primary_analysis.py \
    results/06_crc_projection/gut_scoring_full \
    results/06_crc_projection/gut_scoring_primary_analysis
