#!/usr/bin/env python3
"""
Phase I inventory pass -- before designing how to project the frozen D/F/P
signature onto real CRC data, check the actual structure of every CRC
single-cell/spatial dataset candidate directly. Same discipline as Step
1's inventory and Step 5's Tabula Sapiens inventory: never assume raw
counts, gene-ID convention, or an "Oncofetal" annotation exist -- open the
real file and check.

Candidates found on Argos (via /home/zz950/DATA), none previously vetted
inside TWEAKR-OncoPlacental:
  1. GSE178318 -- 10x mtx triplet, primary CRC + liver metastasis, 9 patients
  2. HTAN_CRC_progressive_plasticity -- h5ad, epithelial-only export, 29 patients
  3. CRLM_NMP_ATLAS -- h5ad, CRC liver metastasis atlas, 6 donors
  4. CRC_single_cell_atlas_2025 (adata_nmf.h5ad) -- 77-dataset meta-atlas with
     pre-computed NMF meta-programs (candidate "Oncofetal"-like M11/revCSC
     program already computed per the KB strategy doc, NOT re-derived here)

For each: report shape, gene-ID convention (symbol vs Ensembl), obs columns,
low-cardinality obs columns (candidate cell-type/patient columns), raw-counts
integer-check on X/.raw.X/each layer, and a case-insensitive search of every
low-cardinality obs column's values for "oncofetal"/"fetal"/"placent"/"trophoblast"
substrings (to directly check, not assume, whether any pre-existing label
already exists).

Usage: python3 inventory_crc_datasets.py <out_dir>
(run on Argos, argos-codex env)
"""
import sys
import os
import gzip
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import io as sio
from scipy import sparse

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "results/06_crc_projection/inventory"
os.makedirs(OUT_DIR, exist_ok=True)

CRC_ATLAS_RAW = "/home/zz950/DATA/CRC-Atlas/phase2/03_data/raw"
META_STUDY = "/home/zz950/DATA/scRNAseq/meta_study/CRC_single_cell_atlas_2025"

KEYWORDS = ["oncofetal", "fetal", "placent", "trophoblast", "revcsc", "m11"]


def integer_check(mat, n_sample=2000):
    if mat is None:
        return None
    sub = mat[:min(n_sample, mat.shape[0])]
    arr = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
    nz = arr[arr != 0]
    if len(nz) == 0:
        return None
    return float(np.mean(np.abs(nz - np.round(nz)) < 1e-6))


def scan_obs_for_keywords(obs):
    hits = {}
    for col in obs.columns:
        try:
            if obs[col].dtype == object or str(obs[col].dtype) == "category":
                vals = obs[col].astype(str).unique()
                for v in vals:
                    vl = v.lower()
                    for kw in KEYWORDS:
                        if kw in vl:
                            hits.setdefault(col, []).append(v)
        except Exception:
            pass
    return hits


def inventory_adata(adata, label, patient_col_guess=None):
    info = {"label": label}
    info["shape"] = {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)}
    info["var_names_sample"] = list(adata.var_names[:10])
    info["looks_like_ensembl"] = bool(str(adata.var_names[0]).startswith("ENSG"))
    info["obs_columns"] = list(adata.obs.columns)
    info["layers"] = list(adata.layers.keys())
    info["has_raw"] = adata.raw is not None
    if adata.raw is not None:
        info["raw_shape"] = {"n_obs": int(adata.raw.n_obs), "n_vars": int(adata.raw.n_vars)}

    value_counts = {}
    for col in adata.obs.columns:
        try:
            nun = adata.obs[col].nunique()
            if 1 < nun <= 200:
                vc = adata.obs[col].value_counts()
                value_counts[col] = {"n_unique": int(nun), "top10": vc.head(10).to_dict()}
        except Exception:
            pass
    info["low_cardinality_obs_columns"] = value_counts

    info["X_integer_fraction"] = integer_check(adata.X)
    if adata.raw is not None:
        info["raw_X_integer_fraction"] = integer_check(adata.raw.X)
    for layer_name in adata.layers.keys():
        info[f"layer_{layer_name}_integer_fraction"] = integer_check(adata.layers[layer_name])

    info["keyword_hits_in_obs"] = scan_obs_for_keywords(adata.obs)
    return info


results = {}

# 1. GSE178318 -- mtx triplet
print("=== GSE178318 ===")
try:
    d = f"{CRC_ATLAS_RAW}/GSE178318"
    mat = sio.mmread(f"{d}/GSE178318_matrix.mtx.gz").tocsr()
    genes = pd.read_csv(f"{d}/GSE178318_genes.tsv.gz", header=None, sep="\t")
    barcodes = pd.read_csv(f"{d}/GSE178318_barcodes.tsv.gz", header=None, sep="\t")
    info = {"label": "GSE178318"}
    info["matrix_shape_as_stored"] = list(mat.shape)
    info["n_genes_file"] = len(genes)
    info["n_barcodes_file"] = len(barcodes)
    # mtx convention is usually genes x cells -- confirm by matching dims
    info["orientation_matches_genes_rows"] = mat.shape[0] == len(genes)
    info["gene_id_columns_sample"] = genes.head(5).astype(str).values.tolist()
    info["barcode_sample"] = barcodes[0].head(5).astype(str).tolist()
    info["X_integer_fraction"] = integer_check(mat.T if mat.shape[0] == len(genes) else mat)
    # sample-id parsing: barcodes often prefixed/suffixed with a sample tag
    bc_series = barcodes[0].astype(str)
    info["barcode_suffix_sample"] = bc_series.head(20).tolist()
    results["GSE178318"] = info
    print(json.dumps(info, indent=2, default=str)[:3000])
except Exception as e:
    results["GSE178318"] = {"error": str(e)}
    print(f"ERROR: {e}")

# 2. HTAN_CRC_progressive_plasticity -- h5ad
print("\n=== HTAN_CRC_progressive_plasticity ===")
try:
    adata = ad.read_h5ad(f"{CRC_ATLAS_RAW}/HTAN_CRC_progressive_plasticity/epithelial.h5ad")
    info = inventory_adata(adata, "HTAN_CRC_progressive_plasticity")
    results["HTAN_CRC_progressive_plasticity"] = info
    print(json.dumps(info, indent=2, default=str)[:3000])
    del adata
except Exception as e:
    results["HTAN_CRC_progressive_plasticity"] = {"error": str(e)}
    print(f"ERROR: {e}")

# 3. CRLM_NMP_ATLAS -- h5ad
print("\n=== CRLM_NMP_ATLAS ===")
try:
    adata = ad.read_h5ad(f"{CRC_ATLAS_RAW}/CRLM_NMP_ATLAS/crlm_nmp_atlas.h5ad")
    info = inventory_adata(adata, "CRLM_NMP_ATLAS")
    results["CRLM_NMP_ATLAS"] = info
    print(json.dumps(info, indent=2, default=str)[:3000])
    del adata
except Exception as e:
    results["CRLM_NMP_ATLAS"] = {"error": str(e)}
    print(f"ERROR: {e}")

# 4. CRC_single_cell_atlas_2025 meta-atlas -- h5ad with precomputed NMF meta-programs
print("\n=== CRC_single_cell_atlas_2025 (adata_nmf.h5ad) ===")
try:
    adata = ad.read_h5ad(f"{META_STUDY}/adata_nmf.h5ad")
    info = inventory_adata(adata, "CRC_single_cell_atlas_2025")
    # Also specifically check for an "M11"/meta-program score column or obsm entry
    info["obsm_keys"] = list(adata.obsm.keys())
    info["uns_keys"] = list(adata.uns.keys())
    results["CRC_single_cell_atlas_2025"] = info
    print(json.dumps(info, indent=2, default=str)[:3000])
    del adata
except Exception as e:
    results["CRC_single_cell_atlas_2025"] = {"error": str(e)}
    print(f"ERROR: {e}")

with open(f"{OUT_DIR}/crc_dataset_inventory.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nWrote {OUT_DIR}/crc_dataset_inventory.json")
