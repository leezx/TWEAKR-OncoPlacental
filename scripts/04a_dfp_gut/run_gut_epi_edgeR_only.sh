#!/bin/bash
#$ -N tweakr_gutatlas_epi_edgeR_retry
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex
export R_LIBS_USER=/home/zz950/softwares/miniforge3/envs/argos-codex/lib/R/user-library

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== edgeR: LargeInt ==="
Rscript scripts/04a_dfp_gut/run_gut_epi_edgeR.R LargeInt

echo ""
echo "=== edgeR: SmallInt ==="
Rscript scripts/04a_dfp_gut/run_gut_epi_edgeR.R SmallInt
