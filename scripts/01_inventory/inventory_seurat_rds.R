#!/usr/bin/env Rscript
# Inventory pass for HDMA per-organ Seurat RDS objects (Aim 1, fetal-somatic
# reference side). For each object: class, assay names + dims, gene-naming
# convention, meta.data columns + value counts for low-cardinality columns,
# reductions, graphs.
#
# Usage: Rscript inventory_seurat_rds.R <path/to/organ.rds> <label> <output_dir>

suppressMessages(library(SeuratObject))
suppressMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: inventory_seurat_rds.R <path> <label> <output_dir>")
}
path <- args[1]
label <- args[2]
out_dir <- args[3]

info <- list(label = label, path = path)

obj <- tryCatch(readRDS(path), error = function(e) e)
if (inherits(obj, "error")) {
  info$error <- conditionMessage(obj)
} else {
  info$class <- class(obj)[1]
  info$assays <- names(obj@assays)
  info$dims <- list(genes = nrow(obj), cells = ncol(obj))
  info$gene_sample <- head(rownames(obj), 10)
  info$looks_like_ensembl <- grepl("^ENSG", rownames(obj)[1])
  info$reductions <- names(obj@reductions)
  info$graphs <- names(obj@graphs)
  info$meta_columns <- colnames(obj@meta.data)

  value_counts <- list()
  for (col in colnames(obj@meta.data)) {
    v <- obj@meta.data[[col]]
    if (is.character(v) || is.factor(v)) {
      nun <- length(unique(v))
      if (nun > 1 && nun <= 60) {
        tab <- table(v)
        value_counts[[col]] <- as.list(setNames(as.integer(tab), names(tab)))
      }
    }
  }
  info$meta_value_counts <- value_counts

  for (a in names(obj@assays)) {
    layer_names <- tryCatch(SeuratObject::Layers(obj[[a]]), error = function(e) character(0))
    info[[paste0("assay_", a, "_layers")]] <- layer_names
  }
}

out_path <- file.path(out_dir, paste0(label, ".json"))
write(toJSON(info, auto_unbox = TRUE, pretty = TRUE, force = TRUE), out_path)

cat("===", label, "===\n")
cat(toJSON(info, auto_unbox = TRUE, pretty = TRUE, force = TRUE), "\n")
cat("Wrote:", out_path, "\n")
