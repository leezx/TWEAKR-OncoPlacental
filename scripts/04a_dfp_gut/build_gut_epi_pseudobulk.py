#!/usr/bin/env python3
"""
Build donor x region pseudobulk raw-count matrices for the locked Step 4a
primary DE contrast (docs/STEP4A_GUT_FDEV_DESIGN.md, PR #20 APPROVE,
round 3/4): Second-trimester-fetal vs. Adult epithelium, LargeInt primary /
SmallInt secondary, within epi_raw_counts02_v2.h5ad.

Includes ALL donors in {Second trim, Adult} x {LargeInt, SmallInt} x
Epithelial, regardless of 10X chemistry -- the R step (run_gut_epi_edgeR.R)
subsets to the mandatory 5'-only sensitivity model itself, so this script
only needs to build the pseudobulk once per region and carry the full
donor-level covariate table (tenX, source_family, Age_group) alongside it,
same separation-of-concerns pattern as Step 4's
build_trophoblast_pseudobulk.py / run_trophoblast_edgeR.R.

Usage: python3 build_gut_epi_pseudobulk.py
(run on Argos, argos-codex env; writes to results/04a_dfp_gut/pseudobulk/)
"""
import os
import re
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

PATH = "/home/zz950/DATA/scRNAseq/GutCellAtlas_Elmentaite2021/raw/epi_raw_counts02_v2.h5ad"
OUT_DIR = "/home/zz950/TWEAKR-OncoPlacental/results/04a_dfp_gut/pseudobulk"
os.makedirs(OUT_DIR, exist_ok=True)

DONOR_COL = "Sample name"


def source_family(batch):
    return re.sub(r"\d+$", "", str(batch))


def build_pseudobulk(adata_sub, region_label):
    obs = adata_sub.obs
    donor = obs[DONOR_COL].astype(str).str.strip().values
    counts = adata_sub.X
    gene_names = adata_sub.var_names.values

    uniq_donors = sorted(set(donor))
    donor_idx = {d: i for i, d in enumerate(uniq_donors)}
    row_idx = np.arange(len(donor))
    col_idx = np.array([donor_idx[d] for d in donor])
    indicator = sparse.csr_matrix(
        (np.ones(len(donor)), (row_idx, col_idx)),
        shape=(len(donor), len(uniq_donors)),
    )
    pseudobulk = indicator.T @ counts
    pseudobulk = np.asarray(pseudobulk.todense()) if sparse.issparse(pseudobulk) else np.asarray(pseudobulk)

    pb_df = pd.DataFrame(pseudobulk.T, index=gene_names, columns=uniq_donors)

    meta = obs.groupby(DONOR_COL, observed=True).agg(
        Age_group=("Age_group", lambda s: s.astype(str).mode().iloc[0]),
        tenX=("10X", lambda s: s.astype(str).mode().iloc[0]),
        batch=("batch", lambda s: s.astype(str).mode().iloc[0]),
        n_cells=("Region", "size"),
    ).reset_index().rename(columns={DONOR_COL: "sample"})
    meta["donor"] = meta["sample"]
    meta["source_family"] = meta["batch"].apply(source_family)
    meta = meta.set_index("sample").loc[uniq_donors].reset_index()

    pb_df.to_csv(f"{OUT_DIR}/{region_label}_pseudobulk_counts.tsv", sep="\t")
    meta.to_csv(f"{OUT_DIR}/{region_label}_pseudobulk_meta.tsv", sep="\t", index=False)
    print(f"{region_label}: {pb_df.shape[0]} genes x {pb_df.shape[1]} donors written", flush=True)
    print(meta.groupby(["Age_group", "tenX", "source_family"], observed=True).size(), flush=True)
    return pb_df, meta


def main():
    print(f"Reading {PATH}...", flush=True)
    adata = ad.read_h5ad(PATH)  # need real X in memory for the indicator-matrix sum, not backed
    adata = adata[adata.obs["category"] == "Epithelial"]
    adata = adata[adata.obs["Age_group"].isin(["Second trim", "Adult"])]

    for region in ["LargeInt", "SmallInt"]:
        sub = adata[adata.obs["Region"] == region].copy()
        build_pseudobulk(sub, region)


if __name__ == "__main__":
    main()
