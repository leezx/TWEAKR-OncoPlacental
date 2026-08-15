#!/usr/bin/env Rscript
# Three-layer statistical validation, layer 2 (PRIMARY evidence per the
# user's explicit directive): preranked GSEA using the real edgeR logFC
# (T_g, the continuous fetal-vs-adult-colon/SI differential statistic) as
# the ranking statistic -- not a binary overlap cut. Tests each of the 5
# independent mike_verzi signatures against the ranked gene list from
# each region's primary edgeR fit (LargeInt_edgeR_primary.tsv /
# SmallInt_edgeR_primary.tsv), using fgsea (Bioconductor's canonical
# preranked-GSEA implementation).
#
# PR #22 round-1 fix (reviewer blocker 1): the ranked list is now
# restricted to the SAME universe U = one2one-ortholog-eligible ^
# region-filterByExpr-tested locked for layers 1/3 (mike_verzi_gut_
# enrichment_permutation.py), not the full ~15k filterByExpr-tested
# gene list. Genes without a one-to-one mouse ortholog can structurally
# never be a signature member (the mike_verzi human_primary.txt files
# are themselves built only from one2one orthologs), so including them
# in the ranked list as guaranteed misses shifts the ES/NES/p reference
# space away from layers 1/3's universe -- not a cosmetic difference.
#
# Usage: Rscript mike_verzi_gut_gsea.R <out_dir>
# Reads results/06a_normal_context/<SIG>_human_primary.txt (5 files)
#       results/06a_normal_context/mouse_biomart_full.tsv (one2one filter)
#       results/04a_dfp_gut/edgeR/{LargeInt,SmallInt}_edgeR_primary.tsv
# Writes <out_dir>/mike_verzi_gut_gsea_results.tsv

suppressMessages(library(fgsea))

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1) args[1] else "results/04a_dfp_gut/mike_verzi_validation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ROOT <- "/home/zz950/TWEAKR-OncoPlacental"
MIKE_VERZI_DIR <- file.path(ROOT, "results/06a_normal_context")
GUT_DIR <- file.path(ROOT, "results/04a_dfp_gut")

SIGNATURES <- c("YAP_SIGNALING_GENES", "REVIVAL_STEM_CELL_GENES", "FETAL_SPHEROID_EPITHELIUM_GENES",
                "REGENERATIVE_EPITHELIUM", "FETAL_INTESTINE_GENES")
REGIONS <- list(Colon = "LargeInt", SI = "SmallInt")

sig_lists <- lapply(SIGNATURES, function(s) {
  readLines(file.path(MIKE_VERZI_DIR, paste0(s, "_human_primary.txt")))
})
names(sig_lists) <- SIGNATURES

# one2one ortholog-eligible universe, same source/filter as layers 1/3
biomart <- read.delim(file.path(MIKE_VERZI_DIR, "mouse_biomart_full.tsv"))
one2one <- unique(biomart$Human.gene.name[biomart$Human.homology.type == "ortholog_one2one" &
                                           nzchar(biomart$Human.gene.name)])
cat(sprintf("Human genes with >=1 Compara one2one mouse ortholog: %d\n", length(one2one)))

all_results <- list()

for (label in names(REGIONS)) {
  region <- REGIONS[[label]]
  edger <- read.delim(file.path(GUT_DIR, "edgeR", paste0(region, "_edgeR_primary.tsv")))

  # ranking statistic: real edgeR logFC (T_g), per the design doc's stated
  # preference -- the exact continuous fetal-vs-adult differential
  # statistic, not a p-value-derived or binarized substitute.
  # Restricted to U = one2one ^ region-filterByExpr-tested (fix above).
  edger_U <- edger[edger$gene %in% one2one, ]
  ranks <- edger_U$logFC
  names(ranks) <- edger_U$gene
  ranks <- sort(ranks, decreasing = TRUE)
  cat(sprintf("=== %s (%s): %d genes ranked by real edgeR logFC (U = one2one ^ filterByExpr-tested) ===\n",
              label, region, length(ranks)))

  set.seed(20260815)  # fixed seed for fgsea's permutation-based p-value estimation, stated explicitly
  fgsea_res <- fgsea(pathways = sig_lists, stats = ranks, minSize = 5, maxSize = 5000, eps = 0)
  fgsea_res$region <- label
  fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, function(x) paste(x, collapse = ";"))
  all_results[[label]] <- fgsea_res

  print(fgsea_res[order(fgsea_res$pval), c("pathway", "NES", "pval", "padj", "size")])
}

combined <- do.call(rbind, all_results)
combined <- combined[, c("region", "pathway", "size", "ES", "NES", "pval", "padj", "leadingEdge")]
colnames(combined)[colnames(combined) == "pathway"] <- "signature"

out_path <- file.path(OUT_DIR, "mike_verzi_gut_gsea_results.tsv")
write.table(combined, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nWrote %d rows to %s\n", nrow(combined), out_path))
