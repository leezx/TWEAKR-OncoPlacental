#!/usr/bin/env python3
"""
Real inventory of epi_raw_counts02_v2.h5ad (Space-Time Gut Cell Atlas,
epithelium-only, full lifespan) -- per PR #20 round-1 REQUEST_CHANGES:
this file becomes the primary construction data for F_Colon-developmental/
F_SI-developmental (same-atlas fetal-vs-adult epithelial DE), replacing the
mechanical HDMA-percentile port. Before any DE compute: check real
age/region/donor structure, and specifically check the reviewer's factual
claim that second-trimester (12-17 PCW) fetal + adult (29-69 yr) samples
share 10x 5' v2 chemistry while first-trimester (6-11 PCW) is a separately-
processed earlier batch -- verify, don't assume.

Usage: python3 inventory_gutcellatlas_epi_lifespan.py <out_dir>
"""
import sys
import os
import json
import numpy as np
import anndata as ad

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "results/04a_dfp_gut/inventory"
os.makedirs(OUT_DIR, exist_ok=True)

PATH = "/home/zz950/DATA/scRNAseq/GutCellAtlas_Elmentaite2021/raw/epi_raw_counts02_v2.h5ad"


def integer_check(mat, n_sample=2000):
    if mat is None:
        return None
    try:
        sub = mat[:min(n_sample, mat.shape[0])]
        arr = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
        nz = arr[arr != 0]
        if len(nz) == 0:
            return None
        return float(np.mean(np.abs(nz - np.round(nz)) < 1e-6))
    except Exception as exc:
        return f"ERROR: {exc}"


print(f"Reading {PATH} (backed='r')...", flush=True)
adata = ad.read_h5ad(PATH, backed="r")

info = {"path": PATH}
info["shape"] = {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)}
info["var_names_sample"] = list(adata.var_names[:10])
info["looks_like_ensembl"] = bool(str(adata.var_names[0]).startswith("ENSG"))
info["var_columns"] = list(adata.var.columns)
info["obs_columns"] = list(adata.obs.columns)
info["layers"] = list(adata.layers.keys())
info["has_raw"] = adata.raw is not None
info["obsm_keys"] = list(adata.obsm.keys())
info["uns_keys"] = list(adata.uns.keys())
info["X_integer_fraction"] = integer_check(adata.X)

print(json.dumps(info, indent=2, default=str), flush=True)

# low-cardinality obs columns -- looking specifically for age/region/donor/
# chemistry/assay/batch columns
value_counts = {}
for col in adata.obs.columns:
    try:
        nun = adata.obs[col].nunique()
        if 1 < nun <= 250:
            vc = adata.obs[col].value_counts()
            value_counts[col] = {"n_unique": int(nun), "top30": vc.head(30).to_dict()}
    except Exception:
        pass
info["low_cardinality_obs_columns"] = value_counts

print("\n--- low-cardinality obs columns ---", flush=True)
for col, v in value_counts.items():
    print(f"\n{col} (n_unique={v['n_unique']}):")
    for k, n in list(v["top30"].items())[:15]:
        print(f"  {k}: {n}")

# keyword scan for age/chemistry/assay/region/donor columns specifically
AGE_KW = ["age", "pcw", "week", "stage", "trimester", "development"]
CHEM_KW = ["chemistry", "assay", "kit", "platform", "10x", "protocol", "chem"]
REGION_KW = ["region", "organ", "tissue", "site", "gut", "intestin", "colon", "duoden", "jejun", "ileum"]
DONOR_KW = ["donor", "sample", "patient", "individual", "subject"]

def scan_cols(keywords):
    return [c for c in adata.obs.columns if any(kw.lower() in c.lower() for kw in keywords)]

info["age_related_columns"] = scan_cols(AGE_KW)
info["chemistry_related_columns"] = scan_cols(CHEM_KW)
info["region_related_columns"] = scan_cols(REGION_KW)
info["donor_related_columns"] = scan_cols(DONOR_KW)

print("\n--- keyword-matched column names ---")
print("age-related:", info["age_related_columns"])
print("chemistry-related:", info["chemistry_related_columns"])
print("region-related:", info["region_related_columns"])
print("donor-related:", info["donor_related_columns"])

out_path = f"{OUT_DIR}/gutcellatlas_epi_lifespan_inventory.json"
with open(out_path, "w") as f:
    json.dump(info, f, indent=2, default=str)
print(f"\nWrote {out_path}")
