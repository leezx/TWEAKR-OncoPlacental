#!/usr/bin/env python3
"""Throwaway timing probe for score_genes_fast at scale (not part of the
reviewed pipeline). Companion to _timing_probe.py."""
import sys
import os
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crc_gut_scoring_core import (
    load_atlas, load_panel_ensembl_ids, testable_genes,
    precompute_score_genes_bins, score_genes_fast,
)

for n_cells in [300000, 665473]:
    adata = load_atlas(n_cells_subset=(None if n_cells >= 665473 else n_cells), seed=20260815)
    genes = load_panel_ensembl_ids("F_Gut-specific")
    real_genes = testable_genes(genes, set(adata.var_names))
    t0 = time.time()
    _, obs_cut = precompute_score_genes_bins(adata)
    t_precompute = time.time() - t0
    t0 = time.time()
    for i in range(5):
        score_genes_fast(adata, real_genes, obs_cut)
    dt = (time.time() - t0) / 5
    print(f"n_cells={adata.n_obs}: precompute={t_precompute:.2f}s (one-time), "
          f"{dt:.3f}s per score_genes_fast call (F_Gut-specific, {len(real_genes)} genes)", flush=True)
    del adata
