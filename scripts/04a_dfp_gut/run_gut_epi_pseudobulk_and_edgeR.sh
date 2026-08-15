#!/bin/bash
#$ -N tweakr_gutatlas_epi_edgeR
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 4
set -uo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex
export R_LIBS_USER=/home/zz950/softwares/miniforge3/envs/argos-codex/lib/R/user-library

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== Building pseudobulk ==="
python3 scripts/04a_dfp_gut/build_gut_epi_pseudobulk.py

echo ""
echo "=== edgeR: LargeInt ==="
Rscript scripts/04a_dfp_gut/run_gut_epi_edgeR.R LargeInt

echo ""
echo "=== edgeR: SmallInt ==="
Rscript scripts/04a_dfp_gut/run_gut_epi_edgeR.R SmallInt
