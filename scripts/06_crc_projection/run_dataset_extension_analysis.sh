#!/bin/bash
#$ -N tweakr_dataset_extension_analysis
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 1
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== Step 6 dataset extension: analysis (correlation + patient-matched contrast) ==="
python3 scripts/06_crc_projection/dataset_extension_analysis.py \
    results/06_crc_projection/dataset_extension \
    results/06_crc_projection/dataset_extension_analysis
