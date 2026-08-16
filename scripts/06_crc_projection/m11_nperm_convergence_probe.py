#!/usr/bin/env python3
"""
Step 6 secondary analysis (PR #28), round-2 review requirement: PR #27's
N_PERM=100-vs-500 convergence check only validated *continuous* percentile
correlation, never binary hard-threshold cohort membership, and it did so
for the *old* percentile-threshold construction, not the cross-cell
z-score-rank construction this design actually uses. Two things follow:

  1. The existing 20,000-cell convergence-check parquet files
     (results/06_crc_projection/gut_scoring_convergence_check/) already
     contain both N_PERM=100 and N_PERM=500 scores for
     revCSC_primary27_minus_CLU_ASS1 -- recomputed locally (no new compute
     needed, see docs/STEP6_SECONDARY_ANALYSIS_DESIGN.md's updated §1) at
     the real z-rank construction: Jaccard(100,500) = 0.7748/0.8199/0.8788
     at top5%/10%/20%. This *replaces* the design doc's old (wrong-
     construction) "~90% Jaccard" claim -- lower than previously stated,
     which if anything strengthens (not weakens) the case for N_PERM=500.
  2. That 100-vs-500 number tells us 100 is insufficient; it does NOT tell
     us 500 has converged. This script runs the actual new probe the
     reviewer asked for: score the 45-gene M11_top50_minus_revCSC_overlap
     panel at N_PERM=500 AND N_PERM=1000 on a fixed 20,000-cell
     stratified-by-study_id subset of the 297,307-cell M11 population
     (not the full-atlas convergence-check subset -- M11 is only defined
     on its own barcode-matched subset), then reports z-rank cohort
     membership Jaccard(500,1000) at top5/10/20%. If 500 already agrees
     well with 1000, N_PERM=500 is certified for the real 297,307-cell run;
     if not, the real run uses whichever draw count this probe shows is
     adequate.

Usage: python3 m11_nperm_convergence_probe.py <out_dir>
(run on Argos via qsub, BEFORE the 297,307-cell M11 run)
"""
import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crc_gut_scoring_core import (
    ATLAS_H5AD, compute_detectability, bin_detectability, score_panel,
    N_BINS, SEED,
)
import scanpy as sc
import anndata as ad

REPO = "/home/zz950/TWEAKR-OncoPlacental"
ADDMODULESCORE_TSV = ("/home/zz950/DATA/scRNAseq/meta_study/CRC_single_cell_atlas_2025/"
                       "NMF/viz_signature_MM_alt_clean_byscore_wardD2/addmodulescore.df.tsv")
GENE_SET_DIR = f"{REPO}/results/06_crc_projection/gut_scoring"
M11_PANEL = "M11_top50_minus_revCSC_overlap"
N_CELLS_SUBSET = 20000


def load_m11_subset_atlas(n_cells_subset, seed):
    """Full atlas -> restricted to the 297,307-cell M11 barcode subset ->
    stratified-by-study_id n_cells_subset sample of THAT population (not
    of the full 665,473-cell atlas)."""
    print(f"Loading {ATLAS_H5AD} ...", flush=True)
    adata = ad.read_h5ad(ATLAS_H5AD)
    print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

    m11_barcodes = pd.read_csv(ADDMODULESCORE_TSV, sep="\t", index_col=0, usecols=[0]).index
    m11_barcodes = pd.Index(m11_barcodes.astype(str))
    overlap = adata.obs_names.intersection(m11_barcodes)
    print(f"M11 subset barcode match: {len(overlap)}/{len(m11_barcodes)} "
          f"addmodulescore.df.tsv rows found in atlas obs_names "
          f"(expected 297,307 exact match per approved design)", flush=True)
    assert len(overlap) == 297307, f"M11 barcode match count changed: got {len(overlap)}, expected 297307"

    adata = adata[adata.obs_names.isin(overlap)].copy()
    print(f"Restricted to M11 subset: {adata.n_obs} cells", flush=True)

    rng = np.random.default_rng(seed)
    study = adata.obs["study_id"].astype(str)
    frac = n_cells_subset / adata.n_obs
    keep_idx = []
    for s, idx in study.groupby(study).groups.items():
        idx = np.asarray(idx)
        n_take = max(1, int(round(len(idx) * frac)))
        n_take = min(n_take, len(idx))
        take = rng.choice(idx, size=n_take, replace=False)
        keep_idx.extend(take.tolist())
    keep_idx = np.asarray(sorted(set(keep_idx)))
    pos = adata.obs.index.get_indexer(keep_idx)
    adata = adata[pos].copy()
    print(f"Stratified subset of M11 population (by study_id, seed={seed}): "
          f"{adata.n_obs} cells (target {n_cells_subset})", flush=True)

    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def zrank_jaccard(scores_a, scores_b, index, cutoffs=(0.95, 0.90, 0.80)):
    rank_a = pd.Series(scores_a, index=index).rank(pct=True, method="average")
    rank_b = pd.Series(scores_b, index=index).rank(pct=True, method="average")
    rows = []
    for cutoff in cutoffs:
        set_a = set(index[rank_a >= cutoff])
        set_b = set(index[rank_b >= cutoff])
        inter = len(set_a & set_b)
        union = len(set_a | set_b)
        rows.append({
            "top_pct": round((1 - cutoff) * 100, 1),
            "n_a": len(set_a), "n_b": len(set_b),
            "intersection": inter, "union": union,
            "jaccard": round(inter / union, 4) if union else float("nan"),
        })
    return pd.DataFrame(rows)


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else \
        f"{REPO}/results/06_crc_projection/m11_nperm_convergence_probe"
    os.makedirs(out_dir, exist_ok=True)

    adata = load_m11_subset_atlas(N_CELLS_SUBSET, SEED)

    panel_path = f"{GENE_SET_DIR}/{M11_PANEL}.ensembl.txt"
    with open(panel_path) as f:
        panel_ens = [l.strip() for l in f if l.strip()]
    print(f"M11 panel: {M11_PANEL}, {len(panel_ens)} genes", flush=True)

    detectability = compute_detectability(adata)
    bin_labels = bin_detectability(detectability, N_BINS)

    print("\n=== Scoring M11 panel at N_PERM=500 ===", flush=True)
    pct500, z500, n_testable_500 = score_panel(adata, M11_PANEL, panel_ens, bin_labels, n_perm=500, seed=SEED)
    print(f"n_testable={n_testable_500}", flush=True)

    print("\n=== Scoring M11 panel at N_PERM=1000 ===", flush=True)
    pct1000, z1000, n_testable_1000 = score_panel(adata, M11_PANEL, panel_ens, bin_labels, n_perm=1000, seed=SEED)
    print(f"n_testable={n_testable_1000}", flush=True)

    out_df = pd.DataFrame({
        f"{M11_PANEL}_percentile_n500": pct500, f"{M11_PANEL}_zscore_n500": z500,
        f"{M11_PANEL}_percentile_n1000": pct1000, f"{M11_PANEL}_zscore_n1000": z1000,
    }, index=adata.obs_names)
    scores_path = f"{out_dir}/m11_scores_nperm500_vs_1000.parquet"
    out_df.to_parquet(scores_path)
    print(f"\nWrote {scores_path}", flush=True)

    jacc_df = zrank_jaccard(z500, z1000, adata.obs_names)
    jacc_path = f"{out_dir}/m11_nperm500_vs_1000_zrank_jaccard.tsv"
    jacc_df.to_csv(jacc_path, sep="\t", index=False)
    print(f"\nWrote {jacc_path}")
    print(jacc_df.to_string(index=False), flush=True)

    from scipy.stats import pearsonr
    r, _ = pearsonr(pct500, pct1000)
    print(f"\nContinuous percentile Pearson r (500 vs 1000): {r:.4f}", flush=True)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
