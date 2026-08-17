#!/usr/bin/env Rscript
# Step 9 (developmental ternary map) -- Track C: real one-vs-rest edgeR DE
# within the single HCL atlas (GSE134355), the one candidate that has
# Fetal intestine + Placenta + Adult intestine on one uniformly-processed
# platform -- avoids the cross-dataset/cross-platform batch-as-biology
# risk Tracks A/B cannot fully escape (Placenta from Arutyunyan/
# Nature2026 vs Fetal/Adult gut from Gut Cell Atlas/HDMA are genuinely
# different studies).
#
# Design: 3 groups (Fetal n=5, Adult n=8, Placenta n=1 -- single donor,
# no replication, reported honestly as a real limitation, not smoothed
# over). ~0+group design, three one-vs-rest contrasts (group - mean of
# other two). Placenta's n=1 means its dispersion is borrowed entirely
# from the common/trended dispersion estimate (edgeR's standard behavior
# for an unreplicated group) -- p-values for Placenta contrasts are not
# meaningful and are not used; only the logFC effect size is used
# downstream (matches the design's own principle: use effect size, not
# p-value, per the KB note this analysis is based on).
#
# Usage: Rscript hcl_edgeR_onevsrest.R
suppressMessages(library(edgeR))

OUT_DIR <- "/Volumes/Stelligen_SSD/Stelligen/PR/TWEAKR-OncoPlacental/results/09_developmental_ternary"

counts <- read.delim(file.path(OUT_DIR, "hcl_pseudobulk_counts.tsv"), row.names = 1, check.names = FALSE)
meta <- read.delim(file.path(OUT_DIR, "hcl_pseudobulk_meta.tsv"))
stopifnot(all(colnames(counts) == meta$sample))

group <- factor(meta$group, levels = c("Fetal", "Placenta", "Adult"))
cat("Group sizes:\n"); print(table(group))

y <- DGEList(counts = counts, group = group)
keep <- filterByExpr(y, group = group)
cat(sprintf("filterByExpr: %d / %d genes kept\n", sum(keep), length(keep)))
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
y <- estimateDisp(y, design)
cat(sprintf("common dispersion: %.4f\n", y$common.dispersion))
fit <- glmQLFit(y, design)

# One-vs-rest contrasts: group - mean(other two)
contrasts <- makeContrasts(
  Fetal_vs_rest    = Fetal - (Placenta + Adult) / 2,
  Placenta_vs_rest = Placenta - (Fetal + Adult) / 2,
  Adult_vs_rest    = Adult - (Fetal + Placenta) / 2,
  levels = design
)

for (cn in colnames(contrasts)) {
  qlf <- glmQLFTest(fit, contrast = contrasts[, cn])
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  tt$gene <- rownames(tt)
  out_path <- file.path(OUT_DIR, sprintf("hcl_edgeR_%s.tsv", cn))
  write.table(tt[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")],
              out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("Wrote %s (%d genes)\n", out_path, nrow(tt)))
}

cat("\nDone.\n")
