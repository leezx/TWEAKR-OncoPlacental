#!/usr/bin/env python3
"""
Inventory pass for the Greenbaum_NatMed_2024 SCP2601 multiome bundle
(mtx-format RNA + ATAC matrices, plus metadata/other CSVs) — the other
half of the placental/trophoblast reference side of Aim 1.

Usage: python inventory_greenbaum_mtx.py <SCP2601_dir> <output_dir>
"""
import sys
import json
import gzip
import os

import scipy.io as sio
import pandas as pd


def inspect_mtx(path: str) -> dict:
    m = sio.mmread(path)
    return {
        "shape": list(m.shape),
        "nnz": int(m.nnz),
        "dtype": str(m.dtype),
    }


def main():
    if len(sys.argv) != 3:
        print("Usage: inventory_greenbaum_mtx.py <SCP2601_dir> <output_dir>", file=sys.stderr)
        sys.exit(2)
    base, out_dir = sys.argv[1], sys.argv[2]

    info: dict = {}

    rna_dir = f"{base}/expression/665be9bd512b6f7d793fa94b"
    atac_dir = f"{base}/expression/665bea28aa3336975ba927f2"
    info["rna_matrix"] = inspect_mtx(f"{rna_dir}/matrix.mtx")
    info["atac_matrix"] = inspect_mtx(f"{atac_dir}/matrix_atac.mtx")

    with open(f"{rna_dir}/genes_rna2.tsv") as f:
        genes = [line.split("\t")[0].strip() for line in f.readlines()[:10]]
    info["rna_gene_sample"] = genes
    info["rna_looks_like_ensembl"] = genes[0].startswith("ENSG") if genes else None

    meta = pd.read_csv(f"{base}/metadata/metadata.csv", low_memory=False)
    info["metadata_columns"] = list(meta.columns)
    info["metadata_nrows"] = len(meta)
    value_counts = {}
    for col in meta.columns:
        try:
            nunique = meta[col].nunique()
        except TypeError:
            continue
        if 1 < nunique <= 60:
            vc = meta[col].value_counts()
            value_counts[col] = {str(k): int(v) for k, v in vc.items()}
    info["metadata_value_counts"] = value_counts

    cluster_path = f"{base}/other/humanplacenta_cluster.csv"
    if os.path.exists(cluster_path):
        cl = pd.read_csv(cluster_path, low_memory=False)
        info["cluster_file_columns"] = list(cl.columns)
        info["cluster_file_nrows"] = len(cl)
        for col in cl.columns:
            try:
                nunique = cl[col].nunique()
            except TypeError:
                continue
            if 1 < nunique <= 60:
                vc = cl[col].value_counts()
                info.setdefault("cluster_value_counts", {})[col] = {str(k): int(v) for k, v in vc.items()}

    spatial_path = f"{base}/other/humanplacenta_spatial.csv"
    if os.path.exists(spatial_path):
        sp = pd.read_csv(spatial_path, nrows=5)
        info["spatial_file_columns"] = list(sp.columns)

    out_path = f"{out_dir}/Greenbaum_NatMed_2024.json"
    with open(out_path, "w") as f:
        json.dump(info, f, indent=2, default=str)

    print("=== Greenbaum_NatMed_2024 ===")
    print(json.dumps(info, indent=2, default=str))
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
