#!/bin/bash
#$ -N tweakr_step8_clim_scrna_scoring
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 8
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== Step 8: D/F/P/revCSC scoring, GSE231559 (primary+normal) / GSE285990 / GSE225857 non-immune ==="
echo "Per docs/STEP8_CLIM_SCRNA_SCORING_DESIGN.md (PR #37, APPROVE at d06edf1)"
python3 scripts/08_clim_scrna_scoring/clim_scrna_scoring_driver.py \
    results/08_clim_scrna_scoring
