#!/usr/bin/env python3
"""
Gene-ID / duplicate-symbol audit, per PR #21 REQUEST_CHANGES blocker #2:
"gene_id -> source var_name -> canonical symbol" must be preserved and
audited before D/F/P set arithmetic or marker calibration trusts raw
var_names symbol strings.

Real findings this audit exists to document (verified via BioMart, not
assumed):
  - 506 var_names (104 base symbols) in epi_raw_counts02_v2.h5ad are
    anndata-uniquified duplicates (e.g. "IGF2" and "IGF2-1" both present).
  - 0 duplicate gene_ids (Ensembl IDs) -- every var_name row is a genuinely
    distinct Ensembl gene, NOT a case of one biological gene split into
    multiple pseudobulk features. No DE refit is required.
  - IGF2 specifically: var_name "IGF2" -> ENSG00000284779 ("novel protein",
    BioMart external_gene_name is blank -- NOT actually HGNC-curated
    "IGF2" despite the symbol string in this dataset's var_names).
    var_name "IGF2-1" -> ENSG00000167244, HGNC:5466, "insulin like growth
    factor 2" -- the real, canonical IGF2. The marker-calibration script
    was looking up "IGF2" and silently finding the wrong (non-canonical)
    record, missing IGF2-1's massive real signal (logFC=12.53 LargeInt,
    6.59 SmallInt).
  - Checked against results/04_dfp_signature/dfp_gene_sets/
    P_developmental_primary84.txt (frozen D/F/P input): 0 of the 84 genes
    collide with any of the 104 affected base symbols -- the existing
    D_Colon-shared/D_SI-shared/P_Colon-specific/P_SI-specific gene-set
    arithmetic is NOT affected by this bug.

Usage: python3 gene_id_audit.py <out_dir>
(run on Argos, argos-codex env)
"""
import sys
import os
import re
import json
import pandas as pd
import anndata as ad

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "results/04a_dfp_gut/gene_id_audit"
os.makedirs(OUT_DIR, exist_ok=True)

PATH = "/home/zz950/DATA/scRNAseq/GutCellAtlas_Elmentaite2021/raw/epi_raw_counts02_v2.h5ad"

print(f"Reading {PATH} (backed='r')...", flush=True)
a = ad.read_h5ad(PATH, backed="r")
vn = a.var_names
gene_ids = a.var["gene_ids"]

base = pd.Series([re.sub(r"-\d+$", "", v) for v in vn], index=vn)
var_id_map = pd.DataFrame({
    "var_name": vn,
    "gene_id": gene_ids.values,
    "base_symbol": base.values,
})
var_id_map["is_suffixed"] = var_id_map["var_name"] != var_id_map["base_symbol"]
var_id_map["base_symbol_is_duplicated"] = var_id_map["base_symbol"].isin(
    base[base.duplicated(keep=False)].unique()
)

out_path = f"{OUT_DIR}/var_id_map.tsv"
var_id_map.to_csv(out_path, sep="\t", index=False)
print(f"Wrote {out_path} ({len(var_id_map)} rows)", flush=True)

n_dup_gene_id = gene_ids.duplicated(keep=False).sum()
n_dup_base = var_id_map["base_symbol_is_duplicated"].sum()
n_dup_base_symbols = var_id_map.loc[var_id_map["base_symbol_is_duplicated"], "base_symbol"].nunique()

print(f"\nDuplicate gene_ids (would indicate a real feature-collapse bug): {n_dup_gene_id}", flush=True)
print(f"var_names affected by duplicate-symbol suffixing: {n_dup_base} ({n_dup_base_symbols} base symbols)", flush=True)

summary = {
    "n_vars_total": int(len(vn)),
    "n_duplicate_gene_ids": int(n_dup_gene_id),
    "n_suffix_affected_var_names": int(n_dup_base),
    "n_suffix_affected_base_symbols": int(n_dup_base_symbols),
}
with open(f"{OUT_DIR}/gene_id_audit_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(f"Wrote {OUT_DIR}/gene_id_audit_summary.json")
