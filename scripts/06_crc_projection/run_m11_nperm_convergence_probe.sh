#!/bin/bash
#$ -N tweakr_m11_nperm_probe
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 4
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

echo "=== Step 1: build M11 Ensembl gene sets + M11 x revCSC overlap audit ==="
python3 scripts/06_crc_projection/build_m11_gene_sets.py \
    results/06_crc_projection/gut_scoring

echo "=== Step 2: M11 N_PERM=500 vs 1000 convergence probe (20,000-cell M11-subset sample) ==="
python3 scripts/06_crc_projection/m11_nperm_convergence_probe.py \
    results/06_crc_projection/m11_nperm_convergence_probe
