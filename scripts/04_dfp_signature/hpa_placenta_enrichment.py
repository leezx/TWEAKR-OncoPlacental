#!/usr/bin/env python3
"""
HPA trophoblast/placenta known-gene enrichment calibration for the
P-developmental effect-size cutoff, per the PR #9 reviewer's explicit next
step: "先把 0.75 作为 leading candidate、再用独立的 HPA trophoblast/placenta
enrichment 决定是否冻结".

Independent of the edgeR pseudobulk DE itself: uses HPA's own bulk RNA
consensus data (rna_tissue_hpa.tsv, per STEP3_METHOD_CONTRACT.md's Tier-1
adult reference; placenta row used here only as a POSITIVE cross-check,
never as adult-negative reference, per the role-column fix in PR #5).

Defines "HPA placenta-enriched" using HPA's own published definition of
"Tissue enriched": >=4-fold higher nTPM in one tissue vs every other
tissue, with a floor (nTPM>=1) to avoid noise from near-zero genes.

For each dataset and each candidate effect-size cutoff (0.5/0.75/1.0),
computes: universe = filterByExpr-tested genes; pass_set = logFC>=cutoff &
FDR<0.05 (directional, per the round-1 fix); background/pass-set HPA-
placenta-enriched rates; fold enrichment; hypergeometric p-value. Also
reports the 2-of-2 overlap set at each cutoff.

Usage: python3 hpa_placenta_enrichment.py
(run on Argos, argos-codex env)
"""
import zipfile
import pandas as pd
import numpy as np
from scipy.stats import hypergeom

ROOT = "/home/zz950/TWEAKR-OncoPlacental"
HPA_ZIP = "/home/zz950/DATA/1.Databases/HPA_RNA_tissue_consensus/raw/rna_tissue_hpa.tsv.zip"
EDGE_DIR = f"{ROOT}/results/04_dfp_signature/edgeR"
OUT = f"{ROOT}/results/04_dfp_signature/hpa_placenta_enrichment.tsv"

# ---- Load HPA long-format table, pivot to gene x tissue nTPM ----
with zipfile.ZipFile(HPA_ZIP) as z:
    inner = z.namelist()[0]
    with z.open(inner) as f:
        hpa = pd.read_csv(f, sep="\t")

wide = hpa.pivot_table(index="Gene name", columns="Tissue", values="nTPM", aggfunc="max")
print(f"HPA wide table: {wide.shape[0]} genes x {wide.shape[1]} tissues")
assert "placenta" in wide.columns

placenta = wide["placenta"]
other = wide.drop(columns=["placenta"])
other_max = other.max(axis=1)

# HPA's own "Tissue enriched" definition: >=4x every other tissue, plus a
# floor so near-zero-everywhere genes don't count as "enriched" on noise.
placenta_enriched = (placenta >= 1.0) & (placenta >= 4 * other_max)
hpa_placenta_set = set(placenta.index[placenta_enriched])
print(f"HPA placenta-enriched gene set (>=4x every other tissue, nTPM>=1): {len(hpa_placenta_set)} genes")

# ---- Load edgeR results ----
datasets = {}
for label in ["Arutyunyan", "Nature2026"]:
    df = pd.read_csv(f"{EDGE_DIR}/{label}_edgeR_results.tsv", sep="\t")
    datasets[label] = df

cutoffs = [0.5, 0.75, 1.0]
rows = []


def enrichment_row(name, universe_genes, pass_genes):
    universe_genes = set(universe_genes)
    pass_genes = set(pass_genes) & universe_genes
    N = len(universe_genes)
    K = len(universe_genes & hpa_placenta_set)  # HPA-enriched genes in universe
    n = len(pass_genes)
    k = len(pass_genes & hpa_placenta_set)
    bg_rate = K / N if N else float("nan")
    pass_rate = k / n if n else float("nan")
    fold = pass_rate / bg_rate if bg_rate else float("nan")
    # hypergeometric: P(X >= k) drawing n from population N with K successes
    pval = hypergeom.sf(k - 1, N, K, n) if n else float("nan")
    rows.append({
        "set": name, "universe_N": N, "HPA_enriched_in_universe_K": K,
        "pass_n": n, "pass_HPA_enriched_k": k,
        "background_rate": round(bg_rate, 4) if N else np.nan,
        "pass_rate": round(pass_rate, 4) if n else np.nan,
        "fold_enrichment": round(fold, 2) if bg_rate else np.nan,
        "hypergeom_pval": pval,
    })


for label, df in datasets.items():
    universe = df["gene"].values
    for c in cutoffs:
        pass_genes = df.loc[(df.logFC >= c) & (df.FDR < 0.05), "gene"].values
        enrichment_row(f"{label}_cutoff{c}", universe, pass_genes)

# 2-of-2 overlap sets (universe = intersection of both datasets' tested genes)
a_df, n_df = datasets["Arutyunyan"], datasets["Nature2026"]
shared_universe = set(a_df["gene"]) & set(n_df["gene"])
for c in cutoffs:
    a_pass = set(a_df.loc[(a_df.logFC >= c) & (a_df.FDR < 0.05), "gene"])
    n_pass = set(n_df.loc[(n_df.logFC >= c) & (n_df.FDR < 0.05), "gene"])
    overlap = a_pass & n_pass
    enrichment_row(f"2of2_overlap_cutoff{c}", shared_universe, overlap)

result = pd.DataFrame(rows)
result.to_csv(OUT, sep="\t", index=False)
print(result.to_string(index=False))
print(f"\nWrote {OUT}")
