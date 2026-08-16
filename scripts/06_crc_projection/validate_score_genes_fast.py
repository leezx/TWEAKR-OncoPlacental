#!/usr/bin/env python3
"""
Step 6 gut re-anchor: numerical validation that score_genes_fast (the
precomputed-bin reimplementation used for the full 665,473-cell run, see
crc_gut_scoring_core.py) reproduces real scanpy.tl.score_genes exactly.

Runs both implementations on a small (5,000-cell) subset for a spread of
panel sizes (smallest: D_Gut-shared, 8 genes; largest: F_Gut-specific,
~2,192 genes; one mid-size), for both the real signature and 3 null draws
each. Checks TWO things per case, in order of importance:
  1. The selected control-gene SET is byte-identical between
     scanpy.tl.score_genes' internal selection and score_genes_fast's --
     this is the real correctness proof (same algorithm, same inputs,
     same random_state=0 reseed -> must pick the same genes).
  2. The resulting score itself is numerically close (relative tolerance
     1e-4) -- expected to NOT be bit-identical even with an identical
     control-gene set, since scanpy's sparse mean
     (`_sparse_nanmean`: explicit NaN-masking then `sum/n_elements`) and
     score_genes_fast's plain `.mean(axis=1)` accumulate floating-point
     sums in a different order; residual float64 rounding noise scaling
     with gene-list size (observed empirically: ~1e-8 for an 8-gene
     panel, ~1e-6 for a 2,192-gene panel) is expected and not a defect.
Both checks must pass for score_genes_fast to be trusted for the full run.

Usage: python3 validate_score_genes_fast.py
(run on Argos, argos-codex env)
"""
import sys
import os
import time
import numpy as np
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crc_gut_scoring_core import (
    load_atlas, compute_detectability, bin_detectability, load_panel_ensembl_ids,
    testable_genes, draw_null_gene_set, precompute_score_genes_bins, score_genes_fast,
    N_BINS, SEED,
)
from scanpy.tools._score_genes import _score_genes_bins, _check_score_genes_args

CHECK_PANELS = ["D_Gut-shared", "F_Colon-specific", "F_Gut-specific"]
N_CELLS = 5000
N_NULL_CHECKS = 3


def main():
    adata = load_atlas(n_cells_subset=N_CELLS, seed=SEED)
    detectability = compute_detectability(adata)
    bin_labels = bin_detectability(detectability, N_BINS)

    t0 = time.time()
    _, obs_cut = precompute_score_genes_bins(adata)
    print(f"Precomputed bins in {time.time()-t0:.1f}s", flush=True)

    rng = np.random.default_rng(SEED)
    all_ok = True

    for panel in CHECK_PANELS:
        panel_ens = load_panel_ensembl_ids(panel)
        real_genes = testable_genes(panel_ens, set(adata.var_names))
        print(f"\n=== {panel} (n_testable={len(real_genes)}) ===", flush=True)

        gene_sets_to_check = [("real_signature", real_genes)]
        bin_pools = {b: idx.tolist() for b, idx in bin_labels.groupby(bin_labels).groups.items()}
        for i in range(N_NULL_CHECKS):
            null_genes = draw_null_gene_set(real_genes, bin_labels, bin_pools, rng)
            gene_sets_to_check.append((f"null_draw_{i}", null_genes))

        for label, genes in gene_sets_to_check:
            # --- check 1: control-gene SET identity against scanpy's own
            # internal _score_genes_bins (the real correctness proof) ---
            np.random.seed(0)
            gene_list_idx, gene_pool_idx, get_subset = _check_score_genes_args(
                adata, genes, None, use_raw=False, layer=None)
            ref_control_genes = set()
            for r_genes in _score_genes_bins(
                gene_list_idx, gene_pool_idx, ctrl_as_ref=True, ctrl_size=50,
                n_bins=25, get_subset=get_subset,
            ):
                ref_control_genes |= set(r_genes)

            _, fast_control_genes = score_genes_fast(adata, genes, obs_cut, random_state=0,
                                                       return_control_genes=True)
            fast_control_genes = set(fast_control_genes)
            sets_match = ref_control_genes == fast_control_genes

            # --- check 2: numerical closeness of the resulting score ---
            sc.tl.score_genes(adata, gene_list=genes, score_name="_tmp_ref", use_raw=False, random_state=0)
            ref_score = adata.obs["_tmp_ref"].values.copy()
            del adata.obs["_tmp_ref"]
            fast_score = score_genes_fast(adata, genes, obs_cut, random_state=0)
            close = bool(np.allclose(ref_score, fast_score, rtol=1e-4, atol=1e-6))
            max_abs_diff = float(np.max(np.abs(ref_score - fast_score)))

            ok = sets_match and close
            all_ok = all_ok and ok
            status = "OK" if ok else "MISMATCH"
            print(f"  {label}: control_gene_sets_match={sets_match} "
                  f"(n_ref={len(ref_control_genes)}, n_fast={len(fast_control_genes)}), "
                  f"max_abs_diff={max_abs_diff:.3e}, allclose(rtol=1e-4)={close} [{status}]", flush=True)

    print(f"\n=== {'ALL CHECKS PASSED' if all_ok else 'AT LEAST ONE MISMATCH -- DO NOT TRUST score_genes_fast'} ===", flush=True)
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
