#!/usr/bin/env python3
"""Throwaway timing probe: is CSR->CSC conversion of adata.X (before
repeated column-slicing calls) a real further speedup for score_genes_fast
at full scale? layers['counts'] confirmed CSR (row-oriented) -- column
slicing on CSR is known to be inefficient regardless of column count."""
import sys
import os
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crc_gut_scoring_core import load_atlas, load_panel_ensembl_ids, testable_genes, precompute_score_genes_bins, score_genes_fast
from scipy import sparse

adata = load_atlas(n_cells_subset=None)
genes27 = testable_genes(load_panel_ensembl_ids("revCSC_primary27_full"), set(adata.var_names))

_, obs_cut = precompute_score_genes_bins(adata)

t0 = time.time()
for i in range(3):
    score_genes_fast(adata, genes27, obs_cut)
print(f"CSR: {(time.time()-t0)/3:.2f}s per call (27-gene panel)", flush=True)

t0 = time.time()
adata.X = adata.X.tocsc()
print(f"CSR->CSC conversion: {time.time()-t0:.2f}s (one-time)", flush=True)

t0 = time.time()
for i in range(3):
    score_genes_fast(adata, genes27, obs_cut)
print(f"CSC: {(time.time()-t0)/3:.2f}s per call (27-gene panel)", flush=True)

genes_big = testable_genes(load_panel_ensembl_ids("F_Gut-specific"), set(adata.var_names))
t0 = time.time()
for i in range(3):
    score_genes_fast(adata, genes_big, obs_cut)
print(f"CSC: {(time.time()-t0)/3:.2f}s per call (F_Gut-specific, {len(genes_big)} genes)", flush=True)
