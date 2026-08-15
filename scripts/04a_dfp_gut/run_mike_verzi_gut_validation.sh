#!/bin/bash
#$ -N tweakr_mike_verzi_gut_validation
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 4
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex
export R_LIBS_USER=/home/zz950/softwares/miniforge3/envs/argos-codex/lib/R/user-library

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== Layer 1 + 3: hypergeometric enrichment + permutation nulls ==="
python3 scripts/04a_dfp_gut/mike_verzi_gut_enrichment_permutation.py results/04a_dfp_gut/mike_verzi_validation

echo ""
echo "=== Layer 2: preranked GSEA (primary evidence) ==="
Rscript scripts/04a_dfp_gut/mike_verzi_gut_gsea.R results/04a_dfp_gut/mike_verzi_validation
