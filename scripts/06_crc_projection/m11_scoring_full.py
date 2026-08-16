#!/usr/bin/env python3
"""
Step 6 secondary analysis compute (per docs/STEP6_SECONDARY_ANALYSIS_DESIGN.md,
PR #28, APPROVE after 4 review rounds): scores all 6 M11 gene-set panels
(top50/100/200 x full/minus_revCSC_overlap) on the 297,307-cell M11
subset at N_PERM=500 -- certified adequate by the real 500-vs-1000
z-rank convergence probe run during PR #28's review (job 3621108,
Jaccard 0.9437/0.9666/0.9817 at top5/10/20%, Pearson r=0.9997).

Barcode-match to the M11 subset re-verified exact (297,307/297,307) here
on the full population (m11_nperm_convergence_probe.py did the same
check on a 20k-cell sample of it).

Uses score_genes_fast (fast production path, CSC-converted X) -- same
numerically-validated machinery as the primary compute's full run
(crc_gut_scoring_full.py), not the real-scanpy path (that path is
reserved for convergence/validation checks, per that script's docstring).

Usage: python3 m11_scoring_full.py <out_dir>
(run on Argos via qsub)
"""
import sys
import os
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crc_gut_scoring_core import ATLAS_H5AD, score_all_panels, SEED

REPO = "/home/zz950/TWEAKR-OncoPlacental"
ADDMODULESCORE_TSV = ("/home/zz950/DATA/scRNAseq/meta_study/CRC_single_cell_atlas_2025/"
                       "NMF/viz_signature_MM_alt_clean_byscore_wardD2/addmodulescore.df.tsv")
M11_PANELS = [
    "M11_top50_full", "M11_top50_minus_revCSC_overlap",
    "M11_top100_full", "M11_top100_minus_revCSC_overlap",
    "M11_top200_full", "M11_top200_minus_revCSC_overlap",
]
N_PERM = 500  # certified via job 3621108, see STEP6_SECONDARY_ANALYSIS_DESIGN.md sec 3


def load_m11_subset_atlas():
    print(f"Loading {ATLAS_H5AD} ...", flush=True)
    adata = ad.read_h5ad(ATLAS_H5AD)
    print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

    m11_barcodes = pd.read_csv(ADDMODULESCORE_TSV, sep="\t", index_col=0, usecols=[0]).index
    m11_barcodes = pd.Index(m11_barcodes.astype(str))
    overlap = adata.obs_names.intersection(m11_barcodes)
    print(f"M11 subset barcode match: {len(overlap)}/{len(m11_barcodes)} "
          f"(expected 297,307 exact match per approved design)", flush=True)
    assert len(overlap) == 297307, f"M11 barcode match count changed: got {len(overlap)}, expected 297307"

    adata = adata[adata.obs_names.isin(overlap)].copy()
    print(f"Restricted to M11 subset: {adata.n_obs} cells", flush=True)

    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.X = adata.X.tocsc()
    return adata


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else f"{REPO}/results/06_crc_projection/m11_scoring_full"
    os.makedirs(out_dir, exist_ok=True)
    checkpoint_dir = f"{out_dir}/_checkpoints"

    adata = load_m11_subset_atlas()

    scores, n_testable = score_all_panels(
        adata, n_perm=N_PERM, panels=M11_PANELS, seed=SEED,
        checkpoint_dir=checkpoint_dir, fast=True,
    )
    scores_path = f"{out_dir}/m11_scores.parquet"
    scores.to_parquet(scores_path)
    print(f"Wrote {scores_path} ({scores.shape[0]} cells x {scores.shape[1]} cols)", flush=True)

    meta_cols = ["study_id", "patient_id", "platform", "atlas_cell_type_middle",
                 "donor_id", "sample_id"]
    meta_cols = [c for c in meta_cols if c in adata.obs.columns]
    meta = adata.obs[meta_cols].copy()
    meta["donor_key"] = meta["study_id"].astype(str) + "||" + meta["patient_id"].astype(str)
    meta_path = f"{out_dir}/m11_cell_metadata.parquet"
    meta.to_parquet(meta_path)
    print(f"Wrote {meta_path}", flush=True)

    n_testable_path = f"{out_dir}/n_testable_genes_per_m11_panel.tsv"
    pd.DataFrame([
        {"panel": p, "n_testable": n_testable[p], "n_perm_used": N_PERM} for p in M11_PANELS
    ]).to_csv(n_testable_path, sep="\t", index=False)
    print(f"Wrote {n_testable_path}", flush=True)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
