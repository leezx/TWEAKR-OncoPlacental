#!/bin/bash
#$ -N tweakr_04_hdma_pseudobulk
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 4
set -uo pipefail
# Per-organ loop, not `-e`: one organ's rds failing to load shouldn't stop
# the other 6 -- same discipline as run_inventory.sh (PR #3 review), exit
# code tracked explicitly per organ so a partial failure can't look like a
# clean run.

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

cd /home/zz950/TWEAKR-OncoPlacental

declare -a FAILED=()
declare -a OK=()

for organ in Adrenal Thyroid Spleen Thymus Liver Skin StomachEsophagus; do
  echo "=== $organ — $(date -Iseconds) ==="
  if Rscript scripts/04_dfp_signature/build_hdma_pseudobulk.R "$organ"; then
    OK+=("$organ")
  else
    echo "*** FAILED: $organ ***"
    FAILED+=("$organ")
  fi
done

echo
echo "=== Summary — $(date -Iseconds) ==="
echo "OK (${#OK[@]}): ${OK[*]:-none}"
echo "FAILED (${#FAILED[@]}): ${FAILED[*]:-none}"

if [ "${#FAILED[@]}" -gt 0 ]; then
  exit 1
fi
