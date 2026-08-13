#!/usr/bin/env python3
"""
Inventory pass for h5ad datasets in the normal-development reference (Aim 1).

For each file: shape, gene-ID convention (symbol vs Ensembl), obs columns,
value_counts for low-cardinality obs columns (candidate cell-type/annotation
fields), presence of a duplicate `raw` layer, and layers/obsm/obsp keys.

Output: one JSON file per dataset under results/, plus a combined summary
printed to stdout (captured by the qsub log).

Usage: python inventory_h5ad.py <path/to/file.h5ad> <label> <output_dir>
"""
import sys
import json
import anndata


def inventory_one(path: str, label: str) -> dict:
    info: dict = {"label": label, "path": path}
    try:
        # backed='r' — this is a metadata inventory pass, not the actual
        # analysis; avoid materializing the full sparse X unnecessarily.
        adata = anndata.read_h5ad(path, backed="r")
    except Exception as exc:
        info["error"] = f"{exc.__class__.__name__}: {exc}"
        return info

    info["shape"] = list(adata.shape)
    info["var_index_sample"] = list(adata.var_names[:10])
    info["looks_like_ensembl"] = bool(adata.var_names[0].startswith("ENSG"))
    info["obs_columns"] = list(adata.obs.columns)
    info["has_raw"] = adata.raw is not None
    if adata.raw is not None:
        info["raw_nnz"] = int(adata.raw.X.nnz) if hasattr(adata.raw.X, "nnz") else None
    info["X_nnz"] = int(adata.X.nnz) if hasattr(adata.X, "nnz") else None
    info["X_dtype"] = str(adata.X.dtype)
    info["layers"] = list(adata.layers.keys())
    info["obsm"] = list(adata.obsm.keys())
    info["obsp"] = list(adata.obsp.keys())

    value_counts = {}
    for col in adata.obs.columns:
        try:
            nunique = adata.obs[col].nunique()
        except TypeError:
            continue
        if 1 < nunique <= 60:
            vc = adata.obs[col].value_counts()
            value_counts[col] = {str(k): int(v) for k, v in vc.items()}
    info["obs_value_counts"] = value_counts

    return info


def main():
    if len(sys.argv) != 4:
        print("Usage: inventory_h5ad.py <path> <label> <output_dir>", file=sys.stderr)
        sys.exit(2)
    path, label, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]

    result = inventory_one(path, label)

    out_path = f"{out_dir}/{label}.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=str)

    print(f"=== {label} ===")
    print(json.dumps(result, indent=2, default=str))
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
