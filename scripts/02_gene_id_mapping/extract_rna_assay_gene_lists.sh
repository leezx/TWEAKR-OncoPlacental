#!/bin/bash
#$ -N tweakr_02_rna_assay_gene_lists
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -uo pipefail
# Extends Step 2's gene-ID mapping to cover the RNA assay's raw-counts
# feature space, per a real gap found while building HDMA per-sample
# pseudobulk (Step 4): the existing canonical_feature_map/<organ>*.tsv
# files were built from rownames(obj) (the object's DEFAULT assay --
# decontX/SCT, a QC-filtered reduced gene set), NOT the RNA assay's raw
# "counts" layer, which has MORE genes (e.g. Adrenal: 25,314 default vs
# 28,375 RNA -- confirmed the default set is a strict subset, 0 genes in
# default but missing from RNA, so this only ADDS coverage, doesn't
# change any existing mapping). Pseudobulk must sum true raw counts
# (RNA assay), so it needs a mapping covering RNA's full feature space.

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

ROOT=/home/zz950/TWEAKR-OncoPlacental
DATA=/home/zz950/DATA
OUT="$ROOT/results/02_gene_id_mapping/gene_lists"
mkdir -p "$OUT"

declare -a FAILED=()
declare -a OK=()

for organ in Adrenal Thyroid Spleen Thymus Liver Skin StomachEsophagus; do
  echo "=== $organ — $(date -Iseconds) ==="
  if Rscript -e '
    suppressMessages(library(SeuratObject))
    args <- commandArgs(trailingOnly = TRUE)
    organ <- args[1]; out_dir <- args[2]
    obj <- readRDS(sprintf("'"$DATA"'/scRNAseq/HumanDevelopmentMultiomicAtlas/raw/per_organ_RNA_seurat/%s_RNA_obj_clustered_final.rds", organ))
    g <- rownames(GetAssayData(obj, assay = "RNA", layer = "counts"))
    is_ensg <- grepl("^ENSG", g)
    writeLines(g, file.path(out_dir, paste0(organ, "_rna_assay_all_genes.txt")))
    writeLines(g[is_ensg], file.path(out_dir, paste0(organ, "_rna_assay_ensg_only.txt")))
    cat(organ, "RNA-assay total:", length(g), "ensg:", sum(is_ensg), "symbol:", sum(!is_ensg), "\n")
  ' "$organ" "$OUT"; then
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

# Union of RNA-assay ENSG IDs across all 7 organs, and the delta vs the
# already-covered union from the original (default-assay) mapping -- this
# is the actual new set that needs biotype/HGNC lookups, if any beyond
# what union_all_organs_ensg.txt already covers.
cat "$OUT"/*_rna_assay_ensg_only.txt | sort -u > "$OUT/union_all_organs_rna_assay_ensg.txt"
echo "Union unique RNA-assay ENSG IDs across all organs: $(wc -l < "$OUT/union_all_organs_rna_assay_ensg.txt")"
if [ -f "$OUT/union_all_organs_ensg.txt" ]; then
  comm -23 "$OUT/union_all_organs_rna_assay_ensg.txt" <(sort -u "$OUT/union_all_organs_ensg.txt") > "$OUT/new_ensg_not_in_original_union.txt"
  echo "New ENSG IDs not covered by the original (default-assay) union: $(wc -l < "$OUT/new_ensg_not_in_original_union.txt")"
fi

if [ "${#FAILED[@]}" -gt 0 ]; then
  exit 1
fi
