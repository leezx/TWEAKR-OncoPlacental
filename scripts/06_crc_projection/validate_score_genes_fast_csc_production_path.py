#!/usr/bin/env python3
"""
Step 6 gut re-anchor: CSC production-path parity validation (PR #27
round 2 — reviewer correctly noted that `validate_score_genes_fast.py`
only validated `score_genes_fast` on the default CSR-loaded object,
never on the CSC-converted object the real full-scale run actually
scores against, and that `precompute_score_genes_bins` was computed
*before* CSC conversion in that first validation, whereas the real
production path (`crc_gut_scoring_full.py`) converts to CSC first via
`load_atlas(to_csc=True)` and only then precomputes bins on the
already-CSC object. A speed probe showing CSC is faster is not proof the
CSC-path scores are numerically identical to the reference implementation
-- this script is that proof, using the exact `load_atlas(to_csc=True)`
call and call order production uses, not a reconstruction of it.

Also directly checks the float32-vs-float64 concern the reviewer raised:
`X` is confirmed float32 after `normalize_total`+`log1p` (checked
directly, not assumed), and `scanpy.tl.score_genes`'s own internal
`_sparse_nanmean` explicitly upcasts to float64 for the sum before
dividing -- `_nan_means_dense_or_sparse` (crc_gut_scoring_core.py) was
fixed to do the same rather than relying on sparse `.mean()`'s
unspecified accumulation dtype. This script re-validates that fix too.

Usage: python3 validate_score_genes_fast_csc_production_path.py
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
    print("Loading atlas via the EXACT production call: "
          "load_atlas(n_cells_subset=N_CELLS, to_csc=True)", flush=True)
    adata = load_atlas(n_cells_subset=N_CELLS, seed=SEED, to_csc=True)
    print(f"X dtype after production-path load: {adata.X.dtype}, format: {adata.X.format}", flush=True)
    assert adata.X.format == "csc", "Expected CSC -- production-path load_atlas(to_csc=True) did not convert"

    detectability = compute_detectability(adata)
    bin_labels = bin_detectability(detectability, N_BINS)

    t0 = time.time()
    _, obs_cut = precompute_score_genes_bins(adata)
    print(f"Precomputed bins on the CSC object (production order) in {time.time()-t0:.1f}s", flush=True)

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
            # --- reference: real scanpy.tl.score_genes on the SAME
            # CSC-converted adata (not a fresh CSR reload) ---
            np.random.seed(0)
            gene_list_idx, gene_pool_idx, get_subset = _check_score_genes_args(
                adata, genes, None, use_raw=False, layer=None)
            ref_control_genes = set()
            for r_genes in _score_genes_bins(
                gene_list_idx, gene_pool_idx, ctrl_as_ref=True, ctrl_size=50,
                n_bins=25, get_subset=get_subset,
            ):
                ref_control_genes |= set(r_genes)

            sc.tl.score_genes(adata, gene_list=genes, score_name="_tmp_ref", use_raw=False, random_state=0)
            ref_score = adata.obs["_tmp_ref"].values.copy()
            del adata.obs["_tmp_ref"]

            # --- fast path, same CSC object, production code path ---
            fast_score, fast_control_genes = score_genes_fast(
                adata, genes, obs_cut, random_state=0, return_control_genes=True)
            fast_control_genes = set(fast_control_genes)

            sets_match = ref_control_genes == fast_control_genes
            close = bool(np.allclose(ref_score, fast_score, rtol=1e-4, atol=1e-6))
            max_abs_diff = float(np.max(np.abs(ref_score - fast_score)))

            ok = sets_match and close
            all_ok = all_ok and ok
            status = "OK" if ok else "MISMATCH"
            print(f"  {label}: control_gene_sets_match={sets_match} "
                  f"(n_ref={len(ref_control_genes)}, n_fast={len(fast_control_genes)}), "
                  f"max_abs_diff={max_abs_diff:.3e}, allclose(rtol=1e-4)={close} [{status}]", flush=True)

    # Also: does the percentile/z-score computation (not just raw score)
    # agree between CSC-production-path fast scoring and a CSR-reference
    # full null-calibration run, for one small panel? This is the
    # end-to-end check the reviewer asked for ("ideally also compare
    # resulting empirical percentiles").
    print("\n=== End-to-end percentile check: D_Gut-shared, N_PERM=20 (small, for speed) ===", flush=True)
    from crc_gut_scoring_core import score_panel, score_panel_fast
    panel = "D_Gut-shared"
    panel_ens = load_panel_ensembl_ids(panel)
    N_PERM_CHECK = 20

    adata_csr = load_atlas(n_cells_subset=N_CELLS, seed=SEED, to_csc=False)
    detectability_csr = compute_detectability(adata_csr)
    bin_labels_csr = bin_detectability(detectability_csr, N_BINS)
    pct_ref, z_ref, n_ref = score_panel(adata_csr, panel, panel_ens, bin_labels_csr, N_PERM_CHECK, SEED)

    pct_fast, z_fast, n_fast = score_panel_fast(adata, panel, panel_ens, bin_labels, obs_cut, N_PERM_CHECK, SEED)

    from scipy.stats import pearsonr
    pct_r, _ = pearsonr(pct_ref, pct_fast)
    z_corr = np.corrcoef(z_ref[~np.isnan(z_ref) & ~np.isnan(z_fast)],
                          z_fast[~np.isnan(z_ref) & ~np.isnan(z_fast)])[0, 1]
    print(f"  n_testable: ref={n_ref}, fast={n_fast}", flush=True)
    print(f"  percentile Pearson r (CSR-reference vs CSC-production-fast): {pct_r:.4f}", flush=True)
    print(f"  z-score Pearson r (CSR-reference vs CSC-production-fast): {z_corr:.4f}", flush=True)
    e2e_ok = bool(pct_r > 0.99 and z_corr > 0.99)
    all_ok = all_ok and e2e_ok
    print(f"  End-to-end check: {'OK' if e2e_ok else 'MISMATCH'} (threshold r>0.99, N_PERM=20 so exact "
          f"agreement isn't expected -- this checks the same qualitative signal, not bit-identity)", flush=True)

    print(f"\n=== {'ALL CHECKS PASSED' if all_ok else 'AT LEAST ONE MISMATCH -- DO NOT TRUST score_genes_fast on CSC'} ===", flush=True)
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
