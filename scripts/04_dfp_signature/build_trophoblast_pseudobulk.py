#!/usr/bin/env python3
"""
Build donor x trophoblast-status pseudobulk raw-count matrices for the
placental datasets confirmed to have real raw counts, per the Step 4
statistical design (docs/STEP4_STATISTICAL_DESIGN.md).

Raw-counts availability audit (done directly against the actual files,
not assumed from Step 1's inventory which didn't capture this):
  - Arutyunyan primary_tissue: X is normalized (non-integer); .raw.X
    holds true integer raw counts. USABLE.
  - Nature2026 scPlacenta_host: X is normalized, no .raw slot. But the
    sibling snRNA_raw_counts.h5ad holds integer raw counts, and its
    obs_names match scPlacenta_host's obs_names 100% (191,735/191,735)
    -- no complex join needed, direct index alignment. USABLE. This
    also resolves Step 1 Finding #3 ("snRNA_raw_counts has no usable
    annotation") -- it does now, borrowed from scPlacenta_host.obs.
  - VentoTormo decidua-v3: X is normalize_total+log1p'd (confirmed by
    checking expm1(X) is NOT close to integers, only ~1.4% match); no
    .raw, no layers, no separate raw-counts file on disk. NOT USABLE
    for the primary count-based model -- would need fresh raw data
    from ArrayExpress E-MTAB-6701, out of scope here. Demoted to the
    same secondary/sensitivity role as Greenbaum.

Usage: python build_trophoblast_pseudobulk.py
(run on Argos, argos-codex env; writes to results/04_dfp_signature/pseudobulk/)
"""
import os
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

OUT_DIR = "/home/zz950/TWEAKR-OncoPlacental/results/04_dfp_signature/pseudobulk"
os.makedirs(OUT_DIR, exist_ok=True)


def build_pseudobulk(counts, gene_names, donor, status, dataset_label):
    """counts: cells x genes sparse matrix (raw). donor/status: arrays aligned to counts' rows."""
    donor = pd.Series(donor).astype(str).str.strip().values
    status = np.asarray(status)
    groups = pd.DataFrame({"donor": donor, "status": status})
    groups["group"] = groups["donor"] + "__" + groups["status"]
    uniq_groups = sorted(groups["group"].unique())

    # indicator matrix (n_cells x n_groups), then pseudobulk = indicator.T @ counts
    group_idx = {g: i for i, g in enumerate(uniq_groups)}
    row_idx = np.arange(len(groups))
    col_idx = groups["group"].map(group_idx).values
    indicator = sparse.csr_matrix(
        (np.ones(len(groups)), (row_idx, col_idx)),
        shape=(len(groups), len(uniq_groups)),
    )
    pseudobulk = indicator.T @ counts  # n_groups x n_genes
    pseudobulk = np.asarray(pseudobulk.todense()) if sparse.issparse(pseudobulk) else np.asarray(pseudobulk)

    pb_df = pd.DataFrame(pseudobulk.T, index=gene_names, columns=uniq_groups)
    meta = pd.DataFrame({"sample": uniq_groups})
    meta["donor"] = [g.split("__")[0] for g in uniq_groups]
    meta["status"] = [g.split("__")[1] for g in uniq_groups]
    meta["n_cells"] = groups.groupby("group").size().reindex(uniq_groups).values

    pb_df.to_csv(f"{OUT_DIR}/{dataset_label}_pseudobulk_counts.tsv", sep="\t")
    meta.to_csv(f"{OUT_DIR}/{dataset_label}_pseudobulk_meta.tsv", sep="\t", index=False)
    print(f"{dataset_label}: {pb_df.shape[0]} genes x {pb_df.shape[1]} samples "
          f"({meta['donor'].nunique()} donors) written")
    return pb_df, meta


def build_arutyunyan():
    a = ad.read_h5ad(
        "/home/zz950/DATA/scRNAseq/Arutyunyan2023_MFI/raw/primary_tissue/"
        "adata_all_donors_all_cell_states_UPD_20230307.h5ad"
    )
    counts = a.raw.X  # true raw counts
    gene_names = a.raw.var_names.values
    donor = a.obs["donor"].values
    status = np.where(a.obs["coarse_annot"].values == "Trophoblast", "troph", "nontroph")
    return build_pseudobulk(counts, gene_names, donor, status, "Arutyunyan")


def build_nature2026():
    host = ad.read_h5ad(
        "/home/zz950/DATA/scRNAseq/2026_human_maternal_fetal_Nature/raw/scPlacenta_host.h5ad"
    )
    raw = ad.read_h5ad(
        "/home/zz950/DATA/scRNAseq/2026_human_maternal_fetal_Nature/raw/snRNA_raw_counts.h5ad"
    )
    # Direct index alignment -- verified 100% obs_names overlap already.
    raw_aligned = raw[host.obs_names, :]
    assert (raw_aligned.obs_names == host.obs_names).all(), "obs_names alignment failed"

    counts = raw_aligned.X
    gene_names = raw_aligned.var_names.values
    donor = host.obs["sample_id"].values
    status = np.where(host.obs["major_class"].isin(["SCT", "VCT", "EVT"]).values, "troph", "nontroph")
    return build_pseudobulk(counts, gene_names, donor, status, "Nature2026")


def main():
    build_arutyunyan()
    build_nature2026()
    print("\nVentoTormo skipped: no raw counts available in decidua-v3.h5ad "
          "(X is normalize_total+log1p'd, not reversible; no .raw; no separate "
          "raw-counts file on disk). Not usable for the primary count-based model.")


if __name__ == "__main__":
    main()
