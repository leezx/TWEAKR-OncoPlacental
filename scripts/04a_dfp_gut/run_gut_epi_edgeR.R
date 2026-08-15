#!/usr/bin/env Rscript
# Step 4a primary DE (docs/STEP4A_GUT_FDEV_DESIGN.md, PR #20 APPROVE round 4):
# Second-trimester-fetal vs. Adult epithelium, edgeR quasi-likelihood F-test,
# unpaired design ~ tenX + source_family + Age_group (donors are disjoint
# individuals between fetal/adult, confirmed in the design audit -- no
# donor-blocking term).
#
# Two robustness checks, per PR #21 REQUEST_CHANGES round 1 (both are
# SUBSETS of the same primary donor set -- NOT independent data -- worded
# that way throughout, not as "independent replication"):
#   1. "5'-only subset sensitivity" (~ source_family + Age_group, tenX
#      dropped since constant within the subset): all 5'-chemistry donors,
#      a subset of the primary donor set (LargeInt 8/12, SmallInt 7/10).
#   2. "Human_colon_16S + 5'-only exact-matched subset" (~ Age_group only,
#      both tenX and source_family constant within this subset): the
#      cleanest possible matched-stratum check (LargeInt 2 fetal vs 5
#      adult, SmallInt 2 vs 4) -- per PR #20's APPROVE round-4 non-blocking
#      implementation note, not run in the first PR #21 submission,
#      added here. Purely descriptive: not expected/required to reach
#      independent significance at this sample size, only checked for
#      effect direction/rough magnitude concordance with the primary
#      model.
#
# Usage: Rscript run_gut_epi_edgeR.R <region_label>   # LargeInt | SmallInt
# Reads results/04a_dfp_gut/pseudobulk/<label>_pseudobulk_{counts,meta}.tsv
# Writes results/04a_dfp_gut/edgeR/<label>_edgeR_primary.tsv
#        results/04a_dfp_gut/edgeR/<label>_edgeR_5prime_subset.tsv
#        results/04a_dfp_gut/edgeR/<label>_edgeR_Human_colon_16S_5prime_exact_matched.tsv
#        results/04a_dfp_gut/edgeR/<label>_concordance_summary.txt

suppressMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
label <- args[1]

ROOT <- "/home/zz950/TWEAKR-OncoPlacental"
IN_DIR <- file.path(ROOT, "results/04a_dfp_gut/pseudobulk")
OUT_DIR <- file.path(ROOT, "results/04a_dfp_gut/edgeR")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

counts <- read.delim(file.path(IN_DIR, paste0(label, "_pseudobulk_counts.tsv")),
                      row.names = 1, check.names = FALSE)
meta <- read.delim(file.path(IN_DIR, paste0(label, "_pseudobulk_meta.tsv")))
stopifnot(all(colnames(counts) == meta$sample))

meta$Age_group <- factor(meta$Age_group, levels = c("Adult", "Second trim"))
meta$tenX <- factor(meta$tenX)
meta$source_family <- factor(meta$source_family)

cat(sprintf("=== %s: donor-level design ===\n", label))
print(table(meta$Age_group, meta$tenX, meta$source_family))

run_edgeR <- function(counts_sub, meta_sub, design_formula, coef_name, tag) {
  # Drop unused factor levels after subsetting -- same bug class caught in the
  # Python design audit (pandas category dtype retaining unused levels
  # silently inflates the design matrix). model.matrix() has the identical
  # failure mode with an un-dropped R factor: a level with zero rows in this
  # subset still gets a dummy column, making a genuinely full-rank design
  # look rank-deficient (or, for a formula term whose only surviving level
  # is confounded with the intercept, actually cause a mismatched column
  # count against qr()$rank). Always re-drop within this function, not just
  # once at load time, since callers pass different row subsets.
  meta_sub$Age_group <- droplevels(meta_sub$Age_group)
  meta_sub$tenX <- droplevels(meta_sub$tenX)
  meta_sub$source_family <- droplevels(meta_sub$source_family)

  y <- DGEList(counts = counts_sub, group = meta_sub$Age_group)
  keep_genes <- filterByExpr(y, group = meta_sub$Age_group)
  cat(sprintf("%s [%s]: %d/%d genes pass filterByExpr\n", label, tag, sum(keep_genes), length(keep_genes)))
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  design <- model.matrix(design_formula, data = meta_sub)
  stopifnot(qr(design)$rank == ncol(design))  # re-verify full rank on the actual data going into the fit, not just trust the earlier audit
  cat(sprintf("%s [%s]: design matrix %d x %d, rank %d (full rank confirmed)\n",
              label, tag, nrow(design), ncol(design), qr(design)$rank))

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = coef_name)

  res <- topTags(qlf, n = Inf)$table
  res$gene <- rownames(res)
  res <- res[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")]
  list(res = res, y = y, design = design)
}

# --- primary model: all donors, ~ tenX + source_family + Age_group ---
primary <- run_edgeR(counts, meta, ~ tenX + source_family + Age_group,
                      "Age_groupSecond trim", "primary")
out_primary <- file.path(OUT_DIR, paste0(label, "_edgeR_primary.tsv"))
write.table(primary$res, out_primary, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("%s: wrote %d gene results to %s\n", label, nrow(primary$res), out_primary))
cat(sprintf("%s [primary]: %d genes FDR<0.05, %d genes FDR<0.05 & |logFC|>=1\n",
            label, sum(primary$res$FDR < 0.05, na.rm = TRUE),
            sum(primary$res$FDR < 0.05 & abs(primary$res$logFC) >= 1, na.rm = TRUE)))

# --- mandatory check 1: 5'-only SUBSET (not independent data) ---
fp <- meta$tenX == "5'"
cat(sprintf("\n%s: 5'-only subset = %d/%d primary donors (%s)\n",
            label, sum(fp), nrow(meta), paste(table(meta$Age_group[fp]), collapse = " vs ")))
sens <- run_edgeR(counts[, fp], meta[fp, ], ~ source_family + Age_group,
                   "Age_groupSecond trim", "5prime_subset")
out_sens <- file.path(OUT_DIR, paste0(label, "_edgeR_5prime_subset.tsv"))
write.table(sens$res, out_sens, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("%s: wrote %d gene results to %s\n", label, nrow(sens$res), out_sens))

# --- mandatory check 2: Human_colon_16S + 5'-only EXACT-MATCHED subset ---
# (PR #20 APPROVE round-4 non-blocking note, not run in PR #21's first
# submission -- added here.) Both tenX and source_family are constant
# within this subset by construction, so the model is ~ Age_group alone.
em <- meta$tenX == "5'" & meta$source_family == "Human_colon_16S"
cat(sprintf("\n%s: Human_colon_16S+5'-only exact-matched subset = %d/%d primary donors (%s)\n",
            label, sum(em), nrow(meta), paste(table(meta$Age_group[em]), collapse = " vs ")))
exact <- run_edgeR(counts[, em], meta[em, ], ~ Age_group,
                    "Age_groupSecond trim", "exact_matched")
out_exact <- file.path(OUT_DIR, paste0(label, "_edgeR_Human_colon_16S_5prime_exact_matched.tsv"))
write.table(exact$res, out_exact, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("%s: wrote %d gene results to %s\n", label, nrow(exact$res), out_exact))

# --- concordance checks (mandatory per design doc, not contingent on diagnostics) ---
# Both checks below compare against SUBSETS of the primary donor set, not
# independent data -- reported as subset/matched-stratum robustness, never
# as "independent replication".
concordance <- function(primary_res, other_res, tag) {
  shared <- merge(primary_res, other_res, by = "gene", suffixes = c(".primary", ".other"))
  cor_pearson <- cor(shared$logFC.primary, shared$logFC.other, method = "pearson")
  cor_spearman <- cor(shared$logFC.primary, shared$logFC.other, method = "spearman")
  sign_concordant <- mean(sign(shared$logFC.primary) == sign(shared$logFC.other), na.rm = TRUE)
  sig_primary <- shared[shared$FDR.primary < 0.05, ]
  sign_concordant_sig <- if (nrow(sig_primary) > 0) {
    mean(sign(sig_primary$logFC.primary) == sign(sig_primary$logFC.other), na.rm = TRUE)
  } else NA
  c(
    sprintf("%s primary-vs-%s subset concordance (%d shared genes after filterByExpr in both -- %s is a SUBSET of the primary donors, not independent data):",
            label, tag, nrow(shared), tag),
    sprintf("  Pearson r (logFC): %.4f", cor_pearson),
    sprintf("  Spearman rho (logFC): %.4f", cor_spearman),
    sprintf("  Sign-concordance (all shared genes): %.4f", sign_concordant),
    sprintf("  Sign-concordance (primary FDR<0.05 genes, n=%d): %s",
            nrow(sig_primary), ifelse(is.na(sign_concordant_sig), "NA (no primary FDR<0.05 genes)", sprintf("%.4f", sign_concordant_sig)))
  )
}

summary_lines <- c(
  concordance(primary$res, sens$res, "5prime_subset"),
  "",
  concordance(primary$res, exact$res, "Human_colon_16S_5prime_exact_matched")
)
writeLines(summary_lines, file.path(OUT_DIR, paste0(label, "_concordance_summary.txt")))
cat(paste(summary_lines, collapse = "\n"), "\n")
