#!/bin/bash
#$ -N tweakr_dataset_extension_scoring
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== Step 6 dataset extension: HTAN_CRC_progressive_plasticity + CRLM_NMP_ATLAS scoring ==="
python3 scripts/06_crc_projection/dataset_extension_scoring.py \
    results/06_crc_projection/dataset_extension
