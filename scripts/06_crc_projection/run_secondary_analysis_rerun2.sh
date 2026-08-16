#!/bin/bash
#$ -N tweakr_secondary_analysis_rerun2
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

# PR #29 round 2 review fixes -- composition (Step0xStepA cross-tab) was
# confirmed correct in round 2, not rerun. Only concordance's zero-cell
# donor-exclusion rule and the signed/max-abs column naming changed.

echo "=== Step 3 (rerun 2): M11 x revCSC concordance, literal zero-cell exclusion rule ==="
python3 scripts/06_crc_projection/secondary_analysis_concordance.py \
    results/06_crc_projection/m11_scoring_full \
    results/06_crc_projection/gut_scoring_full \
    results/06_crc_projection/secondary_analysis_composition \
    results/06_crc_projection/secondary_analysis_concordance
