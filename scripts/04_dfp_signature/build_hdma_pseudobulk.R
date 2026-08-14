#!/usr/bin/env Rscript
# Per-organ, per-individual-sample pseudobulk raw-count matrices for
# F-developmental's positive evidence, per docs/STEP4_STATISTICAL_DESIGN.md
# section 2: "per individual Sample, per organ -- sum raw counts across all
# cells in that sample". No cell-type filtering -- F-developmental's
# "elevated_in_fetal_somatic" is a whole-organ statement, not restricted to
# a specific fetal cell type (HDMA has no internal non-fetal-somatic
# contrast population to subset against; see STEP4_DFP_DESIGN.md).
#
# Genes collapsed to Step 2's canonical_symbol via
# results/02_gene_id_mapping/canonical_feature_map_rna_assay/<organ>_canonical_feature_map.tsv
# -- the RNA-assay-scoped map (build_canonical_feature_map_rna_assay.py),
# not the original canonical_feature_map/ (built from each object's
# DEFAULT assay -- decontX/SCT -- which covers fewer genes than the RNA
# assay's raw "counts" layer this script actually reads: e.g. Adrenal
# 25,314 default vs 28,375 RNA. Confirmed default is a strict subset of
# RNA in every organ, and the RNA-assay map gives identical
# canonical_symbol for every gene the two maps share (0/217,606
# mismatches) -- a superset extension, not a divergent remapping.
# Collisions (multiple original features mapping to one canonical_symbol)
# are summed, per that map's own SUMMARY.md note.
#
# Usage: Rscript build_hdma_pseudobulk.R <organ>
# Writes results/04_dfp_signature/hdma_pseudobulk/<organ>_pseudobulk_{counts,meta}.tsv

suppressMessages(library(SeuratObject))
suppressMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
organ <- args[1]

ROOT <- "/home/zz950/TWEAKR-OncoPlacental"
DATA <- "/home/zz950/DATA"
OUT_DIR <- file.path(ROOT, "results/04_dfp_signature/hdma_pseudobulk")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

rds_path <- sprintf("%s/scRNAseq/HumanDevelopmentMultiomicAtlas/raw/per_organ_RNA_seurat/%s_RNA_obj_clustered_final.rds", DATA, organ)
map_path <- file.path(ROOT, sprintf("results/02_gene_id_mapping/canonical_feature_map_rna_assay/%s_canonical_feature_map.tsv", organ))

cat(sprintf("%s: reading %s\n", organ, rds_path))
obj <- readRDS(rds_path)

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
cat(sprintf("%s: raw RNA counts %d genes x %d cells\n", organ, nrow(counts), ncol(counts)))

# Verify raw-count-ness directly rather than trusting the layer name (same
# discipline as the h5ad raw-counts audit in build_trophoblast_pseudobulk.py).
sample_vals <- counts@x[seq_len(min(10000, length(counts@x)))]
frac_integer <- mean(abs(sample_vals - round(sample_vals)) < 1e-6)
cat(sprintf("%s: fraction of sampled nonzero values that are integer-valued: %.4f\n", organ, frac_integer))
if (frac_integer < 0.99) {
  stop(sprintf("%s: RNA counts layer does not look like raw integer counts (only %.2f%% integer-valued) -- aborting rather than silently building pseudobulk from normalized data", organ, 100 * frac_integer))
}

# ---- Collapse to canonical_symbol (name-matched, not positional) ----
map <- read.delim(map_path, stringsAsFactors = FALSE)
rownames(map) <- map$original_feature
missing <- setdiff(rownames(counts), map$original_feature)
if (length(missing) > 0) {
  stop(sprintf("%s: %d genes in counts have no canonical_feature_map entry (e.g. %s) -- map is stale relative to this rds", organ, length(missing), paste(head(missing, 5), collapse = ", ")))
}
canonical <- map[rownames(counts), "canonical_symbol"]

uniq_canon <- sort(unique(canonical))
gene_indicator <- sparseMatrix(
  i = seq_along(canonical), j = match(canonical, uniq_canon), x = 1,
  dims = c(length(canonical), length(uniq_canon))
)
collapsed <- Matrix::t(gene_indicator) %*% counts  # n_canonical x n_cells
rownames(collapsed) <- uniq_canon
n_collisions <- length(canonical) - length(uniq_canon)
cat(sprintf("%s: %d original features -> %d canonical genes (%d collapsed by collision-summing)\n",
            organ, length(canonical), length(uniq_canon), n_collisions))

# ---- Per-Sample pseudobulk (name-matched to meta.data, not positional) ----
meta <- obj@meta.data
stopifnot(all(colnames(collapsed) == rownames(meta)))  # Seurat invariant, verified not assumed
samples <- as.character(meta$Sample)
uniq_samples <- sort(unique(samples))
cell_indicator <- sparseMatrix(
  i = seq_along(samples), j = match(samples, uniq_samples), x = 1,
  dims = c(length(samples), length(uniq_samples))
)
pseudobulk <- collapsed %*% cell_indicator  # n_canonical x n_samples
colnames(pseudobulk) <- uniq_samples

pb_df <- as.data.frame(as.matrix(pseudobulk))
pb_df <- cbind(gene = rownames(pb_df), pb_df)
counts_out <- file.path(OUT_DIR, paste0(organ, "_pseudobulk_counts.tsv"))
write.table(pb_df, counts_out, sep = "\t", row.names = FALSE, quote = FALSE)

n_cells_per_sample <- as.integer(table(samples)[uniq_samples])
meta_out_df <- data.frame(sample = uniq_samples, organ = organ, n_cells = n_cells_per_sample)
meta_out <- file.path(OUT_DIR, paste0(organ, "_pseudobulk_meta.tsv"))
write.table(meta_out_df, meta_out, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("%s: wrote %d canonical genes x %d samples to %s\n", organ, nrow(pb_df), ncol(pseudobulk), counts_out))
cat(sprintf("%s: sample sizes (cells): %s\n", organ, paste(sprintf("%s=%d", uniq_samples, n_cells_per_sample), collapse = ", ")))
