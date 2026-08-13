#!/bin/bash
#$ -N tweakr_03_collision_mass
#$ -cwd
#$ -o /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -e /home/zz950/TWEAKR-OncoPlacental/_qsub_logs/
#$ -pe pvm 2
set -uo pipefail
# Step 3 prep: decide the aggregation rule for the 29 symbol collisions found
# in PR #4's canonical_feature_map (results/02_gene_id_mapping/canonical_feature_map/
# collision_report.tsv). "Sum counts across colliding features" is the standard
# rule (same one Cell Ranger/STARsolo use for duplicate gene symbols), but we
# quantify rather than assume: what fraction of an organ's total RNA counts
# actually sit on the colliding original_features? If negligible, the choice
# of aggregation rule (sum vs pick-one) can't materially change downstream
# pseudobulk results, so "sum" is safe to adopt without further debate.

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

ROOT=/home/zz950/TWEAKR-OncoPlacental
DATA=/home/zz950/DATA
COLLISIONS="$ROOT/results/02_gene_id_mapping/canonical_feature_map/collision_report.tsv"
OUT="$ROOT/results/03_pseudobulk_prep"
mkdir -p "$OUT"

declare -a FAILED=()
declare -a OK=()

for organ in Adrenal Thyroid Spleen Thymus Liver Skin StomachEsophagus; do
  echo "=== $organ — $(date -Iseconds) ==="
  if Rscript -e '
    suppressMessages(library(SeuratObject))
    suppressMessages(library(Matrix))
    args <- commandArgs(trailingOnly = TRUE)
    organ <- args[1]; collisions_path <- args[2]; out_dir <- args[3]

    obj <- readRDS(sprintf("'"$DATA"'/scRNAseq/HumanDevelopmentMultiomicAtlas/raw/per_organ_RNA_seurat/%s_RNA_obj_clustered_final.rds", organ))
    counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
    total_counts <- sum(counts)

    coll <- read.delim(collisions_path, stringsAsFactors = FALSE)
    coll_organ <- coll[coll$organ == organ, ]
    feats <- unique(unlist(strsplit(coll_organ$original_features, ",")))
    feats_present <- intersect(feats, rownames(counts))
    collision_counts <- if (length(feats_present) > 0) sum(counts[feats_present, , drop = FALSE]) else 0

    frac <- collision_counts / total_counts
    cat(organ, "total_counts:", total_counts,
        "collision_feature_counts:", collision_counts,
        "n_collision_features:", length(feats_present),
        "fraction:", format(frac, scientific = TRUE), "\n")

    out <- data.frame(organ = organ, total_counts = total_counts,
                       collision_counts = collision_counts,
                       n_collision_features = length(feats_present),
                       fraction = frac)
    write.table(out, file.path(out_dir, paste0(organ, "_collision_mass.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
  ' "$organ" "$COLLISIONS" "$OUT"; then
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

# Combine per-organ TSVs into one table
head -1 "$OUT/Adrenal_collision_mass.tsv" > "$OUT/all_organs_collision_mass.tsv" 2>/dev/null || true
for organ in Adrenal Thyroid Spleen Thymus Liver Skin StomachEsophagus; do
  f="$OUT/${organ}_collision_mass.tsv"
  [ -f "$f" ] && tail -n +2 "$f" >> "$OUT/all_organs_collision_mass.tsv"
done

if [ "${#FAILED[@]}" -gt 0 ]; then
  exit 1
fi
