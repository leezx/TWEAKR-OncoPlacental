#!/bin/bash
#$ -N tweakr_01_inventory
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -uo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

ROOT=/home/zz950/TWEAKR-OncoPlacental
DATA=/home/zz950/DATA
OUT="$ROOT/results/01_inventory"
SCRIPTS="$ROOT/scripts/01_inventory"
mkdir -p "$OUT"

echo "=== Step 1: Inventory pass — $(date -Iseconds) ==="

# --- Placental / trophoblast reference (h5ad) ---
python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/Arutyunyan2023_MFI/raw/primary_tissue/adata_all_donors_all_cell_states_UPD_20230307.h5ad" \
  Arutyunyan_primary_tissue "$OUT"

python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/Arutyunyan2023_MFI/raw/organoid/Organoid_PTO_cellxgene.h5ad" \
  Arutyunyan_organoid_PTO "$OUT"

python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/Arutyunyan2023_MFI/raw/organoid/Organoid_TSC_cellxgene.h5ad" \
  Arutyunyan_organoid_TSC "$OUT"

python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/Arutyunyan2023_MFI/raw/organoid/adata_Fig3_trophoblast_organoids_unstimulated.h5ad" \
  Arutyunyan_organoid_Fig3 "$OUT"

python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/2026_human_maternal_fetal_Nature/raw/scPlacenta_host.h5ad" \
  Nature2026_scPlacenta_host "$OUT"

python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/2026_human_maternal_fetal_Nature/raw/snRNA_raw_counts.h5ad" \
  Nature2026_snRNA_raw_counts "$OUT"

python3 "$SCRIPTS/inventory_h5ad.py" \
  "$DATA/scRNAseq/VentoTormo_Nature_2018/raw/decidua-v3.h5ad" \
  VentoTormo_decidua_v3 "$OUT"

python3 "$SCRIPTS/inventory_greenbaum_mtx.py" \
  "$DATA/scRNAseq/Greenbaum_NatMed_2024/raw/SCP2601" "$OUT"

# --- Fetal-somatic reference (HDMA RDS, 7 organs) ---
for organ in Adrenal Thyroid Spleen Thymus Liver Skin StomachEsophagus; do
  Rscript "$SCRIPTS/inventory_seurat_rds.R" \
    "$DATA/scRNAseq/HumanDevelopmentMultiomicAtlas/raw/per_organ_RNA_seurat/${organ}_RNA_obj_clustered_final.rds" \
    "HDMA_${organ}" "$OUT"
done

echo "=== Done — $(date -Iseconds) ==="
echo "Results in: $OUT"
ls -la "$OUT"
