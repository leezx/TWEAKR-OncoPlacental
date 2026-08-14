#!/usr/bin/env python3
"""
Tabula Sapiens inventory pass -- checks actual file structure (assay/
layer names, gene ID convention, cell-type annotation columns, raw-counts
availability) directly, before designing the Tier-2 validation compute.
Same discipline as Step 1's inventory pass: never assume structure from
prior notes, open the real file.

Usage: python3 inventory_tabula_sapiens.py <path/to/TS_Organ.h5ad> <label> <output_dir>
"""
import sys
import json
import zipfile
import tempfile
import os
import anndata as ad
import numpy as np

path, label, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out_dir, exist_ok=True)

info = {"label": label, "path": path}

# Files are zipped -- extract to a temp dir first
if path.endswith(".zip"):
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        h5ad_names = [n for n in names if n.endswith(".h5ad")]
        info["zip_contents"] = names
        if not h5ad_names:
            info["error"] = f"no .h5ad file found in zip, contents: {names}"
            with open(f"{out_dir}/{label}.json", "w") as f:
                json.dump(info, f, indent=2)
            print(json.dumps(info, indent=2))
            sys.exit(1)
        with tempfile.TemporaryDirectory() as tmpdir:
            z.extract(h5ad_names[0], tmpdir)
            real_path = os.path.join(tmpdir, h5ad_names[0])
            adata = ad.read_h5ad(real_path)
else:
    adata = ad.read_h5ad(path)

info["shape"] = {"n_obs": adata.n_obs, "n_vars": adata.n_vars}
info["var_names_sample"] = list(adata.var_names[:10])
info["looks_like_ensembl"] = str(adata.var_names[0]).startswith("ENSG")
info["obs_columns"] = list(adata.obs.columns)
info["layers"] = list(adata.layers.keys())
info["has_raw"] = adata.raw is not None
if adata.raw is not None:
    info["raw_shape"] = {"n_obs": adata.raw.n_obs, "n_vars": adata.raw.n_vars}

# Cell-type-like columns (low cardinality, likely annotation)
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

# Raw-counts check: sample X and .raw.X (if present) for integer-valued-ness
def integer_check(mat, n_sample=2000):
    if mat is None:
        return None
    sub = mat[:min(n_sample, mat.shape[0])]
    arr = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
    nz = arr[arr != 0]
    if len(nz) == 0:
        return None
    return float(np.mean(np.abs(nz - np.round(nz)) < 1e-6))

info["X_integer_fraction"] = integer_check(adata.X)
if adata.raw is not None:
    info["raw_X_integer_fraction"] = integer_check(adata.raw.X)
for layer_name in adata.layers.keys():
    info[f"layer_{layer_name}_integer_fraction"] = integer_check(adata.layers[layer_name])

with open(f"{out_dir}/{label}.json", "w") as f:
    json.dump(info, f, indent=2, default=str)

print(f"=== {label} ===")
print(json.dumps(info, indent=2, default=str)[:4000])
print(f"\nWrote {out_dir}/{label}.json")
