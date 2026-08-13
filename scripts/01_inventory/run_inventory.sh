#!/bin/bash
#$ -N tweakr_01_inventory
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -uo pipefail
# Deliberately not `-e`: this is a batch of 15 independent per-dataset
# inventory checks, and one failing shouldn't stop the other 14 from
# running. Instead, every invocation's exit code is tracked explicitly and
# the job exits non-zero at the end if anything failed, with a clear
# FAILED/OK summary — a caught-and-logged per-file error must not look like
# a clean "=== Done ===" job success (see PR #3 review).

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

ROOT=/home/zz950/TWEAKR-OncoPlacental
DATA=/home/zz950/DATA
OUT="$ROOT/results/01_inventory"
SCRIPTS="$ROOT/scripts/01_inventory"
mkdir -p "$OUT"

echo "=== Step 1: Inventory pass — $(date -Iseconds) ==="

declare -a FAILED=()
declare -a OK=()

run_step() {
  local label="$1"
  shift
  if "$@"; then
    OK+=("$label")
  else
    echo "*** FAILED: $label (exit $?) ***"
    FAILED+=("$label")
  fi
}

# --- Placental / trophoblast reference (h5ad) ---
run_step Arutyunyan_primary_tissue python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/Arutyunyan2023_MFI/raw/primary_tissue/adata_all_donors_all_cell_states_UPD_20230307.h5ad" \
  Arutyunyan_primary_tissue "$OUT"

run_step Arutyunyan_organoid_PTO python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/Arutyunyan2023_MFI/raw/organoid/Organoid_PTO_cellxgene.h5ad" \
  Arutyunyan_organoid_PTO "$OUT"

run_step Arutyunyan_organoid_TSC python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/Arutyunyan2023_MFI/raw/organoid/Organoid_TSC_cellxgene.h5ad" \
  Arutyunyan_organoid_TSC "$OUT"

run_step Arutyunyan_organoid_Fig3 python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/Arutyunyan2023_MFI/raw/organoid/adata_Fig3_trophoblast_organoids_unstimulated.h5ad" \
  Arutyunyan_organoid_Fig3 "$OUT"

run_step Nature2026_scPlacenta_host python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/2026_human_maternal_fetal_Nature/raw/scPlacenta_host.h5ad" \
  Nature2026_scPlacenta_host "$OUT"

run_step Nature2026_snRNA_raw_counts python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/2026_human_maternal_fetal_Nature/raw/snRNA_raw_counts.h5ad" \
  Nature2026_snRNA_raw_counts "$OUT"

run_step VentoTormo_decidua_v3 python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/VentoTormo_Nature_2018/raw/decidua-v3.h5ad" \
  VentoTormo_decidua_v3 "$OUT"

run_step Greenbaum_NatMed_2024 python3 "$SCRIPTS/inventory_greenbaum_mtx.py" \
  "$DATA/scRNAseq/Greenbaum_NatMed_2024/raw/SCP2601" "$OUT"

# --- Fetal-somatic reference (HDMA RDS, 7 organs) ---
for organ in Adrenal Thyroid Spleen Thymus Liver Skin StomachEsophagus; do
  run_step "HDMA_${organ}" Rscript "$SCRIPTS/inventory_seurat_rds.R" \
    "$DATA/scRNAseq/HumanDevelopmentMultiomicAtlas/raw/per_organ_RNA_seurat/${organ}_RNA_obj_clustered_final.rds" \
    "HDMA_${organ}" "$OUT"
done

echo
echo "=== Summary — $(date -Iseconds) ==="
echo "OK (${#OK[@]}): ${OK[*]:-none}"
echo "FAILED (${#FAILED[@]}): ${FAILED[*]:-none}"
echo "Results in: $OUT"
ls -la "$OUT"

if [ "${#FAILED[@]}" -gt 0 ]; then
  echo "=== Done WITH FAILURES — $(date -Iseconds) ==="
  exit 1
fi
echo "=== Done — all $(( ${#OK[@]} )) datasets OK — $(date -Iseconds) ==="
