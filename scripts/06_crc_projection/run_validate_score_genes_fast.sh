#!/bin/bash
#$ -N tweakr_validate_score_genes_fast
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 4
set -euo pipefail
source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex
cd /home/zz950/TWEAKR-OncoPlacental
python3 scripts/06_crc_projection/validate_score_genes_fast.py
