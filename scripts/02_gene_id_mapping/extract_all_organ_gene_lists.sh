#!/bin/bash
#$ -N tweakr_02_gene_lists
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -uo pipefail
# Full gene-name list per HDMA organ (all 7, not just Adrenal — per PR #4
# review). Splits into already-symbol vs ENSG-fallback per organ so the
# ENSG set feeding the biotype/collision QC downstream is verifiably
# complete, not a single-organ sample.

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
    g <- rownames(obj)
    is_ensg <- grepl("^ENSG", g)
    writeLines(g, file.path(out_dir, paste0(organ, "_all_genes.txt")))
    writeLines(g[is_ensg], file.path(out_dir, paste0(organ, "_ensg_only.txt")))
    cat(organ, "total:", length(g), "ensg:", sum(is_ensg), "symbol:", sum(!is_ensg), "\n")
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

# Union of ENSG IDs across all 7 organs (dedupe) — this is the actual set
# the biotype QC needs to cover, not any single organ's list.
cat "$OUT"/*_ensg_only.txt | sort -u > "$OUT/union_all_organs_ensg.txt"
echo "Union unique ENSG IDs across all organs: $(wc -l < "$OUT/union_all_organs_ensg.txt")"

if [ "${#FAILED[@]}" -gt 0 ]; then
  exit 1
fi
