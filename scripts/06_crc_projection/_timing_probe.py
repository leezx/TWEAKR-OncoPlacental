#!/usr/bin/env python3
"""Throwaway timing probe (not part of the reviewed pipeline): measures
score_genes per-call cost at several cell-count scales to check whether
per-draw cost scales with n_cells before committing to a full-job resource
request. Deleted / not referenced by any qsub script once the real timing
question is answered."""
import sys
import os
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crc_gut_scoring_core import load_atlas, compute_detectability, bin_detectability, load_panel_ensembl_ids, testable_genes
import scanpy as sc
import numpy as np

for n_cells in [20000, 100000, 300000]:
    adata = load_atlas(n_cells_subset=n_cells, seed=20260815)
    genes = load_panel_ensembl_ids("F_Gut-specific")
    real_genes = testable_genes(genes, set(adata.var_names))
    t0 = time.time()
    for i in range(3):
        sc.tl.score_genes(adata, gene_list=real_genes, score_name="_tmp", use_raw=False)
    dt = (time.time() - t0) / 3
    print(f"n_cells={n_cells}: {dt:.2f}s per score_genes call (F_Gut-specific, {len(real_genes)} genes)", flush=True)
    del adata
