#!/bin/bash
#$ -N tweakr_gutatlas_fetal_inventory
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 4
set -uo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental
python3 scripts/04a_dfp_gut/inventory_gutcellatlas_fetal.py results/04a_dfp_gut/inventory
