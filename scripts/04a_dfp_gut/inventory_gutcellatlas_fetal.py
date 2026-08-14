#!/usr/bin/env python3
"""
Real, direct inventory pass on the newly-downloaded Gut Cell Atlas fetal
object (Elmentaite et al., Nature 2021) -- before any F_Colon-developmental /
F_SI-developmental construction work begins. Same discipline as Step 1's
inventory_h5ad.py and Step 6's inventory_crc_datasets.py: never assume raw
counts, gene-ID convention, region annotation, or epithelial-lineage labels
exist -- open the real files and check.

Two files inventoried:
  1. fetal_RAWCOUNTS_cellxgene.h5ad  -- raw counts, candidate for
     F_Colon-developmental / F_SI-developmental construction
  2. final_fetal_object_cellxgene.h5ad -- normalized, for metadata/annotation
     cross-check only (e.g. if raw-counts file has thinner obs metadata)

For each: shape, gene-ID convention (symbol vs Ensembl), obs columns,
low-cardinality obs columns (top values), raw-counts integer-check on
X/.raw.X/each layer, and a case-insensitive keyword scan of every
low-cardinality obs column's values for region (duodenum/jejunum/ileum/
colon/large intestine/small intestine), lineage/cell-type (epithel/stem/
enterocyte/goblet/paneth/tuft/EEC), and donor/sample-structure columns.

Usage: python3 inventory_gutcellatlas_fetal.py <out_dir>
(run on Argos, backed='r' throughout -- inventory only, not the real analysis)
"""
import sys
import os
import json
import numpy as np
import anndata as ad

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "results/04a_dfp_gut/inventory"
os.makedirs(OUT_DIR, exist_ok=True)

RAW_DIR = "/home/zz950/DATA/scRNAseq/GutCellAtlas_Elmentaite2021/raw"

REGION_KEYWORDS = ["duoden", "jejun", "ileum", "ileal", "colon", "cecum", "caecum",
                    "large intestine", "small intestine", "rectum", "sigmoid",
                    "transverse", "ascending", "descending", "large_intestine",
                    "small_intestine", "LargeInt", "SmallInt"]
LINEAGE_KEYWORDS = ["epithel", "stem", "enterocyte", "goblet", "paneth", "tuft",
                     "eec", "enteroendocrine", "progenitor", "crypt", "tuft"]
DONOR_KEYWORDS = ["donor", "sample", "patient", "individual", "subject"]
AGE_KEYWORDS = ["age", "pcw", "week", "stage", "gestation"]


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


def scan_obs_for_keywords(obs, keyword_groups):
    hits = {}
    for group_name, keywords in keyword_groups.items():
        group_hits = {}
        for col in obs.columns:
            try:
                if obs[col].dtype == object or str(obs[col].dtype) == "category":
                    vals = obs[col].astype(str).unique()
                    matched = []
                    for v in vals:
                        vl = v.lower()
                        if any(kw.lower() in vl for kw in keywords):
                            matched.append(v)
                    if matched:
                        group_hits[col] = matched[:30]
            except Exception:
                pass
        if group_hits:
            hits[group_name] = group_hits
    return hits


def inventory_one(path, label):
    info = {"label": label, "path": path}
    try:
        adata = ad.read_h5ad(path, backed="r")
    except Exception as exc:
        info["error"] = f"{exc.__class__.__name__}: {exc}"
        return info

    info["shape"] = {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)}
    info["var_names_sample"] = list(adata.var_names[:10])
    info["looks_like_ensembl"] = bool(str(adata.var_names[0]).startswith("ENSG"))
    info["var_columns"] = list(adata.var.columns)
    if "feature_name" in adata.var.columns:
        info["feature_name_sample"] = list(adata.var["feature_name"].astype(str)[:10])
    info["obs_columns"] = list(adata.obs.columns)
    info["layers"] = list(adata.layers.keys())
    info["has_raw"] = adata.raw is not None
    if adata.raw is not None:
        info["raw_shape"] = {"n_obs": int(adata.raw.n_obs), "n_vars": int(adata.raw.n_vars)}
    info["obsm_keys"] = list(adata.obsm.keys())
    info["uns_keys"] = list(adata.uns.keys())

    value_counts = {}
    for col in adata.obs.columns:
        try:
            nun = adata.obs[col].nunique()
            if 1 < nun <= 200:
                vc = adata.obs[col].value_counts()
                value_counts[col] = {"n_unique": int(nun), "top20": vc.head(20).to_dict()}
        except Exception:
            pass
    info["low_cardinality_obs_columns"] = value_counts

    info["X_integer_fraction"] = integer_check(adata.X)
    if adata.raw is not None:
        info["raw_X_integer_fraction"] = integer_check(adata.raw.X)
    for layer_name in adata.layers.keys():
        info[f"layer_{layer_name}_integer_fraction"] = integer_check(adata.layers[layer_name])

    info["keyword_hits_in_obs"] = scan_obs_for_keywords(adata.obs, {
        "region": REGION_KEYWORDS,
        "epithelial_lineage": LINEAGE_KEYWORDS,
        "donor_sample": DONOR_KEYWORDS,
        "age_stage": AGE_KEYWORDS,
    })

    return info


results = {}

print("=== fetal_RAWCOUNTS_cellxgene.h5ad ===")
info = inventory_one(f"{RAW_DIR}/fetal_RAWCOUNTS_cellxgene.h5ad", "fetal_RAWCOUNTS_cellxgene")
results["fetal_RAWCOUNTS_cellxgene"] = info
print(json.dumps(info, indent=2, default=str)[:6000])

print("\n=== final_fetal_object_cellxgene.h5ad ===")
info = inventory_one(f"{RAW_DIR}/final_fetal_object_cellxgene.h5ad", "final_fetal_object_cellxgene")
results["final_fetal_object_cellxgene"] = info
print(json.dumps(info, indent=2, default=str)[:6000])

with open(f"{OUT_DIR}/gutcellatlas_fetal_inventory.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nWrote {OUT_DIR}/gutcellatlas_fetal_inventory.json")
