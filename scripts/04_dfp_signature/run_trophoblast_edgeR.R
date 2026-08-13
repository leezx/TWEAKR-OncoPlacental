#!/usr/bin/env Rscript
# Primary trophoblast-vs-rest DE per docs/STEP4_STATISTICAL_DESIGN.md (PR #8,
# revised after finding VentoTormo has no raw counts): edgeR quasi-likelihood
# F-test, paired design (~ donor + status), fit independently per dataset.
#
# Usage: Rscript run_trophoblast_edgeR.R <dataset_label>
# Reads results/04_dfp_signature/pseudobulk/<label>_pseudobulk_{counts,meta}.tsv
# Writes results/04_dfp_signature/edgeR/<label>_edgeR_results.tsv

suppressMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
label <- args[1]

ROOT <- "/home/zz950/TWEAKR-OncoPlacental"
IN_DIR <- file.path(ROOT, "results/04_dfp_signature/pseudobulk")
OUT_DIR <- file.path(ROOT, "results/04_dfp_signature/edgeR")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

counts <- read.delim(file.path(IN_DIR, paste0(label, "_pseudobulk_counts.tsv")),
                      row.names = 1, check.names = FALSE)
meta <- read.delim(file.path(IN_DIR, paste0(label, "_pseudobulk_meta.tsv")))

stopifnot(all(colnames(counts) == meta$sample))

# Paired design: keep only donors that actually have both groups (should be
# all of them per the replicate-structure audit, but re-verify here rather
# than assume the upstream script's filtering was exhaustive).
donor_counts <- table(meta$donor)
paired_donors <- names(donor_counts[donor_counts == 2])
keep <- meta$donor %in% paired_donors
cat(sprintf("%s: %d/%d donors have both groups (used for paired DE)\n",
            label, length(paired_donors), length(donor_counts)))

counts <- counts[, keep]
meta <- meta[keep, ]
meta$donor <- factor(meta$donor)
meta$status <- factor(meta$status, levels = c("nontroph", "troph"))
# order by donor so paired design is unambiguous
ord <- order(meta$donor, meta$status)
counts <- counts[, ord]
meta <- meta[ord, ]

y <- DGEList(counts = counts, group = meta$status)
keep_genes <- filterByExpr(y, group = meta$status)
cat(sprintf("%s: %d/%d genes pass filterByExpr\n", label, sum(keep_genes), length(keep_genes)))
y <- y[keep_genes, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~ donor + status, data = meta)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "statustroph")

res <- topTags(qlf, n = Inf)$table
res$gene <- rownames(res)
res <- res[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")]

out_path <- file.path(OUT_DIR, paste0(label, "_edgeR_results.tsv"))
write.table(res, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("%s: wrote %d gene results to %s\n", label, nrow(res), out_path))

cat(sprintf("\n%s dispersion summary: common.dispersion=%.4f\n", label, y$common.dispersion))
cat(sprintf("%s: %d genes FDR<0.05, %d genes FDR<0.05 & |logFC|>=1\n",
            label, sum(res$FDR < 0.05, na.rm = TRUE),
            sum(res$FDR < 0.05 & abs(res$logFC) >= 1, na.rm = TRUE)))
