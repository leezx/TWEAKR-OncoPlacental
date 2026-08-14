#!/usr/bin/env python3
"""
Background detection-rate check, run BEFORE reporting the Tier-2
validation results: the F-specific organ-matched check flagged ~15% of
(gene, cell_type) pairs as "cross-donor-consistent" -- before treating
that as a finding about the frozen signature, need to know whether that
rate is actually elevated relative to a random gene panel, or whether
some cell types (e.g. endothelial cells, known in the literature for
broad low-level transcriptional activity) are just generically
"promiscuous" at the single-cell level regardless of which genes are
tested.

For each organ, samples a random gene panel (same size as that organ's
F-specific gene list, for a fair comparison) from the full transcriptome
and computes the same cross-donor-consistent-detection rate per cell
type, using the identical pseudobulk/CPM/floor logic as
tier2_validation.py.

Usage: python3 background_detection_rate.py
(run on Argos, argos-codex env)
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/zz950/TWEAKR-OncoPlacental/scripts/05_tier2_validation")
from tier2_validation import load_organ, build_donor_celltype_pseudobulk, gene_celltype_report, load_gene_list, F_MATCHED_ORGANS, OUT_DIR

rng = np.random.default_rng(0)
rows = []

for organ in F_MATCHED_ORGANS:
    print(f"=== {organ} ===")
    adata = load_organ(organ)
    pb_counts, meta = build_donor_celltype_pseudobulk(adata)
    lib_size = pb_counts.sum(axis=0)
    pb_cpm = pb_counts.div(lib_size, axis=1) * 1e6

    f_genes = load_gene_list(f"F_developmental_{organ}")
    n_f = len(f_genes)
    all_genes = list(pb_cpm.index)
    random_genes = list(rng.choice(all_genes, size=min(n_f, len(all_genes)), replace=False))

    f_report = gene_celltype_report(pb_cpm, meta, f_genes)
    rand_report = gene_celltype_report(pb_cpm, meta, random_genes)

    f_rate = f_report.cross_donor_consistent_hit.mean() if len(f_report) else float("nan")
    rand_rate = rand_report.cross_donor_consistent_hit.mean() if len(rand_report) else float("nan")

    # Per-cell-type comparison
    f_by_ct = f_report.groupby("cell_type").cross_donor_consistent_hit.mean()
    rand_by_ct = rand_report.groupby("cell_type").cross_donor_consistent_hit.mean()
    for ct in f_by_ct.index:
        rows.append({
            "organ": organ, "cell_type": ct,
            "F_specific_flag_rate": round(float(f_by_ct.get(ct, np.nan)), 4),
            "random_gene_flag_rate": round(float(rand_by_ct.get(ct, np.nan)), 4),
            "n_F_specific_genes": n_f,
        })

    print(f"{organ}: F-specific overall flag rate = {f_rate:.4f}, random-gene-panel flag rate = {rand_rate:.4f} "
          f"(panel size {len(random_genes)})")

result = pd.DataFrame(rows)
result["enrichment_over_background"] = (result["F_specific_flag_rate"] / result["random_gene_flag_rate"].replace(0, np.nan)).round(2)
result.to_csv(f"{OUT_DIR}/background_detection_rate_comparison.tsv", sep="\t", index=False)
print("\n=== Per-cell-type comparison ===")
print(result.to_string(index=False))
print(f"\nWrote {OUT_DIR}/background_detection_rate_comparison.tsv")
