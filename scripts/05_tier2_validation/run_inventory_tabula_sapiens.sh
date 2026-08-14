#!/bin/bash
#$ -N tweakr_05_ts_inventory
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 4
set -uo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

ROOT=/home/zz950/TWEAKR-OncoPlacental
DATA=/home/zz950/DATA
OUT="$ROOT/results/05_tier2_validation/inventory"
SCRIPTS="$ROOT/scripts/05_tier2_validation"
mkdir -p "$OUT"

declare -a FAILED=()
declare -a OK=()

for organ in Liver Skin Spleen Thymus Large_Intestine; do
  echo "=== $organ — $(date -Iseconds) ==="
  if python3 "$SCRIPTS/inventory_tabula_sapiens.py" \
      "$DATA/1.Databases/TabulaSapiens/raw/TS_${organ}.h5ad.zip" \
      "TabulaSapiens_${organ}" "$OUT"; then
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
