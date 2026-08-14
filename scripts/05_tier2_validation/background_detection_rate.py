#!/usr/bin/env python3
"""
Background detection-rate check, PR #15 review round 1 fix.

Round-1 version drew ONE random gene panel (fixed seed) per organ from the
FULL transcriptome and reported a single fold-enrichment number. Reviewer
(REQUEST_CHANGES) identified two real problems with that:

1. No null distribution / uncertainty: a single random-panel realization
   can't tell you whether an observed 2-3x "enrichment" would still hold
   with a different random draw, or whether it's noise from that one draw.
2. No expression/detectability matching: F-specific genes already passed
   Step 4's fetal-expression selection. Sampling uniformly from the full
   transcriptome pulls in large numbers of genes that are lowly or never
   expressed anywhere in adult tissue, which mechanically deflates the
   "background" detection rate and inflates the apparent enrichment ratio
   -- an artifact of the comparison, not necessarily real biology.

Fix: a real matched-permutation null.
  - For each organ, compute each gene's overall adult detectability score
    (fraction of eligible donor-celltype pseudobulk samples, across ALL
    cell types, where CPM >= floor) -- an expression-level proxy computed
    from the SAME Tabula Sapiens data being tested, not an external prior.
  - Bin genes into detectability deciles.
  - Run N=500 permutations: for each permutation, draw one random gene
    per F-specific gene from the SAME decile bin (without replacement
    within a permutation), giving a background panel matched in size AND
    detectability composition to the real F-specific list.
  - For each (organ, cell_type), report the observed F-specific flag rate,
    the null median, null 95% interval (2.5/97.5 percentiles) from the 500
    permutations, and an empirical one-sided p-value (fraction of null
    permutations with a flag rate >= the observed rate).

Cross-donor-consistent-hit status is precomputed ONCE per organ for every
gene in the transcriptome (a genes x cell_types boolean matrix), so that
each of the 500 permutations is just an index + mean() -- avoids
recomputing the expensive per-gene donor-consistency logic 500x.

Usage: python3 background_detection_rate.py
(run on Argos, argos-codex env)
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/zz950/TWEAKR-OncoPlacental/scripts/05_tier2_validation")
from tier2_validation import load_organ, build_donor_celltype_pseudobulk, load_gene_list, F_MATCHED_ORGANS, OUT_DIR, MIN_CELLS, NOT_DETECTED_FLOOR

N_PERM = 500
N_BINS = 10
SEED = 0


def compute_hit_matrix_and_detectability(pb_counts, meta, min_cells=MIN_CELLS, floor=NOT_DETECTED_FLOOR):
    """Vectorized, once-per-organ: returns
    hit_matrix (genes x eligible cell_types, bool: cross-donor-consistent hit)
    detectability (genes,): fraction of ALL eligible donor-celltype samples
    (across every cell type) where CPM >= floor -- used for stratified
    background sampling, not itself part of the hit criterion.
    """
    lib_size = pb_counts.sum(axis=0)
    pb_cpm = pb_counts.div(lib_size, axis=1) * 1e6
    eligible_meta = meta[meta.n_cells >= min_cells]
    eligible_cols = eligible_meta["sample"].values
    detected = pb_cpm[eligible_cols] >= floor  # genes x eligible samples, bool

    detectability = detected.mean(axis=1)  # genes,

    hit_cols = {}
    for ct, sub in eligible_meta.groupby("cell_type"):
        cols = sub["sample"].values
        n_donors = len(cols)
        if n_donors == 0:
            continue
        sub_detected = detected[cols]
        n_detected = sub_detected.sum(axis=1)
        hit = (n_detected == n_donors) & (n_donors >= 2) & (n_detected > 0)
        hit_cols[ct] = hit
    hit_matrix = pd.DataFrame(hit_cols, index=pb_cpm.index)
    return hit_matrix, detectability


def matched_permutation_null(hit_matrix, detectability, f_genes, n_perm=N_PERM, n_bins=N_BINS, seed=SEED):
    """Returns a DataFrame: cell_type x [observed, null_median, null_p2_5,
    null_p97_5, empirical_p, n_perm_used]."""
    rng = np.random.default_rng(seed)
    f_genes_present = [g for g in f_genes if g in hit_matrix.index]

    # Bin ALL genes by detectability decile (matching pool), using the same
    # bin edges the F-specific genes themselves fall into.
    try:
        bin_labels, bin_edges = pd.qcut(detectability, q=n_bins, labels=False, retbins=True, duplicates="drop")
    except ValueError:
        # detectability has too few unique values for n_bins -- fall back to fewer bins
        bin_labels, bin_edges = pd.qcut(detectability, q=min(n_bins, detectability.nunique()), labels=False, retbins=True, duplicates="drop")
    bin_labels = pd.Series(bin_labels, index=detectability.index)

    f_bins = bin_labels.loc[f_genes_present]
    bin_pools = {b: bin_labels.index[bin_labels == b].tolist() for b in bin_labels.dropna().unique()}

    observed = hit_matrix.loc[f_genes_present].mean(axis=0)  # per cell_type

    null_rates = np.zeros((n_perm, hit_matrix.shape[1]))
    for i in range(n_perm):
        panel = []
        for b, genes_in_bin in f_bins.groupby(f_bins).groups.items():
            n_needed = len(genes_in_bin)
            pool = bin_pools.get(b, [])
            if len(pool) >= n_needed:
                draw = rng.choice(pool, size=n_needed, replace=False)
            else:
                # bin too small (shouldn't normally happen) -- fall back to sampling with replacement
                draw = rng.choice(pool if len(pool) > 0 else hit_matrix.index.values, size=n_needed, replace=True)
            panel.extend(draw)
        null_rates[i, :] = hit_matrix.loc[panel].mean(axis=0).values

    null_median = np.median(null_rates, axis=0)
    null_p2_5 = np.percentile(null_rates, 2.5, axis=0)
    null_p97_5 = np.percentile(null_rates, 97.5, axis=0)
    # one-sided empirical p: fraction of permutations with rate >= observed
    empirical_p = (null_rates >= observed.values[None, :]).mean(axis=0)

    result = pd.DataFrame({
        "cell_type": hit_matrix.columns,
        "observed_F_flag_rate": observed.values,
        "null_median": null_median,
        "null_p2_5": null_p2_5,
        "null_p97_5": null_p97_5,
        "empirical_p_value": empirical_p,
        "n_permutations": n_perm,
        "n_F_specific_genes_matched": len(f_genes_present),
    })
    result["enrichment_vs_null_median"] = (result["observed_F_flag_rate"] / result["null_median"].replace(0, np.nan)).round(2)
    result["outside_null_95pct_interval"] = (result["observed_F_flag_rate"] > result["null_p97_5"]) | (result["observed_F_flag_rate"] < result["null_p2_5"])
    return result.round(4)


rows = []
for organ in F_MATCHED_ORGANS:
    print(f"=== {organ} ===")
    adata = load_organ(organ)
    pb_counts, meta = build_donor_celltype_pseudobulk(adata)
    del adata

    hit_matrix, detectability = compute_hit_matrix_and_detectability(pb_counts, meta)
    f_genes = load_gene_list(f"F_developmental_{organ}")
    print(f"{organ}: {len(f_genes)} F-specific candidate genes, {hit_matrix.shape[1]} eligible cell types, "
          f"{(detectability > 0).sum()}/{len(detectability)} genes with any adult detection")

    result = matched_permutation_null(hit_matrix, detectability, f_genes)
    result.insert(0, "organ", organ)
    rows.append(result)

    n_sig = int((result.empirical_p_value < 0.05).sum())
    print(f"{organ}: {n_sig}/{len(result)} cell types have observed F-specific rate significantly "
          f"above the {N_PERM}-permutation expression-matched null (p<0.05, one-sided)")

final = pd.concat(rows, ignore_index=True)
final.to_csv(f"{OUT_DIR}/background_permutation_null.tsv", sep="\t", index=False)
print(f"\n=== Full per-(organ,cell_type) permutation-null comparison ===")
print(final.to_string(index=False))
print(f"\nWrote {OUT_DIR}/background_permutation_null.tsv")

sig = final[final.empirical_p_value < 0.05].sort_values(["organ", "empirical_p_value"])
print(f"\n=== Cell types with F-specific rate significantly above expression-matched null (p<0.05) ===")
print(sig.to_string(index=False) if len(sig) else "(none)")
