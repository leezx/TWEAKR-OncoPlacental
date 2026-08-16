#!/bin/bash
#$ -N tweakr_secondary_analysis_rerun
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

# PR #29 round 1 review fixes -- M11 scoring (job 3621118, N_PERM=500) is
# NOT rerun here (reviewer's own note: "the expensive M11 scoring does
# not need to be rerun" -- only the composition Step0xStepA cross-tab and
# the concordance MH zero-cell-donor exclusion logic changed).

echo "=== Step 2 (rerun): revCSC-high cohort definition + developmental composition ==="
python3 scripts/06_crc_projection/secondary_analysis_composition.py \
    results/06_crc_projection/gut_scoring_full \
    results/06_crc_projection/secondary_analysis_composition

echo "=== Step 3 (rerun): M11 x revCSC concordance (continuous correlation + enrichment) ==="
python3 scripts/06_crc_projection/secondary_analysis_concordance.py \
    results/06_crc_projection/m11_scoring_full \
    results/06_crc_projection/gut_scoring_full \
    results/06_crc_projection/secondary_analysis_composition \
    results/06_crc_projection/secondary_analysis_concordance
