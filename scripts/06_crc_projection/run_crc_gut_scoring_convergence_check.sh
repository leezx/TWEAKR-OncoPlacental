#!/bin/bash
#$ -N tweakr_gut_scoring_convergence
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 4
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== N_PERM=100 vs 500 convergence check (20,000-cell stratified subset, 13 panels) ==="
echo "MUST complete and be reviewed before run_crc_gut_scoring_full.sh"
python3 scripts/06_crc_projection/crc_gut_scoring_convergence_check.py \
    results/06_crc_projection/gut_scoring_convergence_check
