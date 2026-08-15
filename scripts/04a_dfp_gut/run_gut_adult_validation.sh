#!/bin/bash
#$ -N tweakr_gut_adult_validation
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 4
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== Check 2: GTEx bulk adult-expression/consistency audit ==="
python3 scripts/04a_dfp_gut/gut_adult_validation_gtex.py results/04a_dfp_gut/adult_validation

echo ""
echo "=== Checks 1+3: Tabula Sapiens adult-expression/consistency audit ==="
python3 scripts/04a_dfp_gut/gut_adult_validation_tabula_sapiens.py results/04a_dfp_gut/adult_validation
