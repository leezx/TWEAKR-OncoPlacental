#!/usr/bin/env python3
"""
HPA placenta tissue-level enrichment calibration for the P-developmental
effect-size cutoff, per the PR #9 reviewer's explicit next step: "先把 0.75
作为 leading candidate、再用独立的 HPA trophoblast/placenta enrichment 决定
是否冻结".

Independent of the edgeR pseudobulk DE itself: uses HPA's own bulk RNA
consensus data (rna_tissue_hpa.tsv, per STEP3_METHOD_CONTRACT.md's Tier-1
adult reference; placenta row used here only as a POSITIVE cross-check,
never as adult-negative reference, per the role-column fix in PR #5).

Three placenta reference sets, checked against three cutoff candidates:

1. `official_65`: HPA's own published "Tissue enriched (placenta)"
   classification, fetched directly from proteinatlas.org's search API
   (query `tissue_category_rna:placenta;Tissue enriched`) -- exactly 65
   genes, matching HPA's own reported count. This is the real HPA
   classification, not a project recomputation. Saved locally as
   `results/04_dfp_signature/hpa_placenta_official_tissue_enriched_65genes.tsv`.
2. `no_floor`: a **project-recomputed** fourfold-ratio set from the raw HPA
   nTPM matrix (placenta >= 4x every other tissue's max), WITHOUT an
   expression floor. PR #10 round 2 review caught a real bug in this
   recomputation: when placenta and every other tissue are both exactly 0,
   `0 >= 4*0` is trivially true, so near-zero/zero-everywhere genes flood
   into this set on a division-by-zero-shaped artifact -- this is why it
   doesn't reproduce HPA's actual 65-gene classification (604 vs 65). It is
   NOT "HPA's rule exactly as published" and must never be called that --
   it's the project's own recomputation attempt from the matrix, kept only
   as a robustness/sensitivity check alongside the real official set.
3. `floored`: the same project recomputation, but with an explicit
   project-defined nTPM>=1 floor to exclude the zero/zero artifact above.
   Still a project recomputation, not HPA's official classification.

For each dataset and each candidate effect-size cutoff (0.5/0.75/1.0),
computes: universe = filterByExpr-tested genes; pass_set = logFC>=cutoff &
FDR<0.05 (directional, per the PR #9 round-1 fix); background/pass-set
enriched-gene rates against each of the 3 reference sets; fold enrichment;
hypergeometric p-value. Also reports the 2-of-2 overlap set at each cutoff.

This is an HPA **placenta tissue-level** external calibration (bulk RNA
consensus data), not an HPA single-cell trophoblast cell-type calibration
-- flagged per reviewer note not to conflate the two in downstream naming
unless an HPA single-cell trophoblast reference is separately incorporated.

Usage: python3 hpa_placenta_enrichment.py
(run on Argos, argos-codex env)
"""
import zipfile
import pandas as pd
import numpy as np
from scipy.stats import hypergeom

ROOT = "/home/zz950/TWEAKR-OncoPlacental"
HPA_ZIP = "/home/zz950/DATA/1.Databases/HPA_RNA_tissue_consensus/raw/rna_tissue_hpa.tsv.zip"
OFFICIAL_65 = f"{ROOT}/results/04_dfp_signature/hpa_placenta_official_tissue_enriched_65genes.tsv"
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

# Project recomputation of a fourfold ratio from the matrix, unmodified:
# placenta >= 4x every other tissue. NOT HPA's official classification --
# when placenta and all other tissues are 0, 0>=4*0 is trivially true, so
# zero-everywhere genes leak in (this is why it doesn't match the real
# 65-gene official set -- see module docstring, PR #10 round 2 review).
fourfold_only = placenta >= 4 * other_max
hpa_set_no_floor = set(placenta.index[fourfold_only])
print(f"Project-recomputed fourfold-ratio set, no floor (NOT HPA official): {len(hpa_set_no_floor)} genes")

# Project-defined addition: also require nTPM>=1, to exclude the zero/zero
# artifact above. Still a project recomputation, not HPA's own definition.
placenta_enriched_floored = fourfold_only & (placenta >= 1.0)
hpa_set_floored = set(placenta.index[placenta_enriched_floored])
print(f"Project-recomputed fourfold-ratio set, with nTPM>=1 floor (NOT HPA official): {len(hpa_set_floored)} genes")

# The real HPA official classification -- fetched from proteinatlas.org's
# search API, not recomputed. Ground truth for this calibration.
official_genes = pd.read_csv(OFFICIAL_65, sep="\t")["Gene"]
hpa_set_official = set(official_genes)
print(f"HPA official 'Tissue enriched (placenta)' set (proteinatlas.org, ground truth): {len(hpa_set_official)} genes")

# ---- Load edgeR results ----
datasets = {}
for label in ["Arutyunyan", "Nature2026"]:
    df = pd.read_csv(f"{EDGE_DIR}/{label}_edgeR_results.tsv", sep="\t")
    datasets[label] = df

cutoffs = [0.5, 0.75, 1.0]
rows = []


def enrichment_row(name, universe_genes, pass_genes, ref_set, ref_label):
    universe_genes = set(universe_genes)
    pass_genes = set(pass_genes) & universe_genes
    N = len(universe_genes)
    K = len(universe_genes & ref_set)
    n = len(pass_genes)
    k = len(pass_genes & ref_set)
    bg_rate = K / N if N else float("nan")
    pass_rate = k / n if n else float("nan")
    fold = pass_rate / bg_rate if bg_rate else float("nan")
    # hypergeometric: P(X >= k) drawing n from population N with K successes
    pval = hypergeom.sf(k - 1, N, K, n) if n else float("nan")
    rows.append({
        "set": name, "reference": ref_label, "universe_N": N, "ref_in_universe_K": K,
        "pass_n": n, "pass_ref_k": k,
        "background_rate": round(bg_rate, 4) if N else np.nan,
        "pass_rate": round(pass_rate, 4) if n else np.nan,
        "fold_enrichment": round(fold, 2) if bg_rate else np.nan,
        "hypergeom_pval": pval,
    })


references = [
    ("official_65", hpa_set_official),
    ("floored", hpa_set_floored),
    ("no_floor", hpa_set_no_floor),
]

for label, df in datasets.items():
    universe = df["gene"].values
    for c in cutoffs:
        pass_genes = df.loc[(df.logFC >= c) & (df.FDR < 0.05), "gene"].values
        for ref_label, ref_set in references:
            enrichment_row(f"{label}_cutoff{c}", universe, pass_genes, ref_set, ref_label)

# 2-of-2 overlap sets (universe = intersection of both datasets' tested genes)
a_df, n_df = datasets["Arutyunyan"], datasets["Nature2026"]
shared_universe = set(a_df["gene"]) & set(n_df["gene"])
for c in cutoffs:
    a_pass = set(a_df.loc[(a_df.logFC >= c) & (a_df.FDR < 0.05), "gene"])
    n_pass = set(n_df.loc[(n_df.logFC >= c) & (n_df.FDR < 0.05), "gene"])
    overlap = a_pass & n_pass
    for ref_label, ref_set in references:
        enrichment_row(f"2of2_overlap_cutoff{c}", shared_universe, overlap, ref_set, ref_label)

result = pd.DataFrame(rows)
result.to_csv(OUT, sep="\t", index=False)
print(result.to_string(index=False))
print(f"\nWrote {OUT}")
