#!/usr/bin/env python3
"""
Gene-ID / duplicate-symbol audit, v2 -- fixes PR #21 round-3 REQUEST_CHANGES:
the v1 heuristic (strip a trailing "-\\d+" from var_name, treat matches as
anndata-uniquification artifacts) is unreliable, because many real HGNC
symbols legitimately end in a dash-number as part of their own canonical
name -- immunoglobulin/TCR variable-region genes like "IGHVIII-5" are
exactly this case, and are exactly the gene family v1's heuristic flagged
as "mostly Ig/TCR families" -- which the reviewer correctly identified as
the heuristic's own likely failure mode, not evidence of anything real.

Correct method: compare each var_name against BioMart's AUTHORITATIVE
external_gene_name for that var_name's gene_id (Ensembl ID) --
human_biomart_full.tsv, a single bulk query, same pattern already used
successfully elsewhere in this project (mouse_biomart_full.tsv). A
var_name is "renamed" (potentially anndata-uniquified or otherwise not
matching current canonical nomenclature) if and only if it differs from
BioMart's real external_gene_name for that exact gene_id -- not based on
string-pattern guessing.

Also fixes a real gap the reviewer caught: v1's docstring claimed a
P_developmental_primary84.txt collision check that the code never actually
performed. This version reads the file and computes the real check.

Usage: python3 gene_id_audit.py <out_dir> <biomart_tsv>
(run on Argos, argos-codex env; biomart_tsv = human_biomart_full.tsv,
pulled in from the local bulk BioMart query, see query_human_biomart_full.py)
"""
import sys
import os
import json
import pandas as pd
import anndata as ad

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "results/04a_dfp_gut/gene_id_audit"
BIOMART_TSV = sys.argv[2] if len(sys.argv) > 2 else "human_biomart_full.tsv"
os.makedirs(OUT_DIR, exist_ok=True)

PATH = "/home/zz950/DATA/scRNAseq/GutCellAtlas_Elmentaite2021/raw/epi_raw_counts02_v2.h5ad"

print(f"Reading {PATH} (backed='r')...", flush=True)
a = ad.read_h5ad(PATH, backed="r")
vn = a.var_names
gene_ids = a.var["gene_ids"]

print(f"Reading {BIOMART_TSV} (authoritative gene_id -> external_gene_name)...", flush=True)
biomart = pd.read_csv(BIOMART_TSV, sep="\t")
biomart = biomart.rename(columns={"Gene stable ID": "gene_id", "Gene name": "authoritative_symbol"})
biomart_map = dict(zip(biomart["gene_id"], biomart["authoritative_symbol"]))

var_id_map = pd.DataFrame({
    "var_name": vn,
    "gene_id": gene_ids.values,
})
var_id_map["authoritative_symbol"] = var_id_map["gene_id"].map(biomart_map)
var_id_map["found_in_biomart"] = var_id_map["authoritative_symbol"].notna()
# treat "not found in biomart" separately from "found but renamed" -- don't
# conflate the two (a gene absent from this BioMart archive snapshot isn't
# necessarily mis-annotated, just unresolved)
var_id_map["renamed_or_mismatched"] = (
    var_id_map["found_in_biomart"]
    & (var_id_map["authoritative_symbol"].fillna("") != "")
    & (var_id_map["var_name"] != var_id_map["authoritative_symbol"])
)

out_path = f"{OUT_DIR}/var_id_map.tsv"
var_id_map.to_csv(out_path, sep="\t", index=False)
print(f"Wrote {out_path} ({len(var_id_map)} rows)", flush=True)

n_dup_gene_id = gene_ids.duplicated(keep=False).sum()
n_not_found = (~var_id_map["found_in_biomart"]).sum()
n_renamed = var_id_map["renamed_or_mismatched"].sum()

print(f"\nDuplicate gene_ids (would indicate a real feature-collapse bug): {n_dup_gene_id}", flush=True)
print(f"var_names not found in this BioMart archive snapshot: {n_not_found}", flush=True)
print(f"var_names that differ from BioMart's authoritative external_gene_name "
      f"(real renamed/mismatched genes, NOT a string-pattern guess): {n_renamed}", flush=True)

renamed_examples = var_id_map[var_id_map["renamed_or_mismatched"]].head(20)
print("\nExamples of renamed/mismatched genes (var_name -> authoritative_symbol):")
for _, row in renamed_examples.iterrows():
    print(f"  {row['var_name']} -> {row['authoritative_symbol']} ({row['gene_id']})")

# --- P_developmental collision check (real code, matching every gene's
# actual symbol string against the set of var_names known to be
# renamed/mismatched by the authoritative BioMart comparison) ---
P_DEV_PATH = "/home/zz950/TWEAKR-OncoPlacental/results/04_dfp_signature/dfp_gene_sets/P_developmental_primary84.txt"
renamed_var_names = set(var_id_map.loc[var_id_map["renamed_or_mismatched"], "var_name"])
renamed_authoritative_symbols = set(var_id_map.loc[var_id_map["renamed_or_mismatched"], "authoritative_symbol"].dropna())

with open(P_DEV_PATH) as f:
    p_dev_genes = {l.strip() for l in f if l.strip()}

# a P_developmental gene is at risk if its symbol matches either side of a
# renamed pair: the h5ad's (wrong/old) var_name, or the authoritative
# symbol that a renamed var_name should have been -- either way, a naive
# symbol-based overlap check between P_developmental and this dataset's
# F-developmental sets could silently miss or misassign it.
p_dev_collisions_with_old_varname = sorted(p_dev_genes & renamed_var_names)
p_dev_collisions_with_authoritative = sorted(p_dev_genes & renamed_authoritative_symbols)
p_dev_collisions = sorted(set(p_dev_collisions_with_old_varname) | set(p_dev_collisions_with_authoritative))

print(f"\nP_developmental_primary84.txt genes: {len(p_dev_genes)}", flush=True)
print(f"P_developmental genes matching a renamed var_name (old symbol side): {p_dev_collisions_with_old_varname}", flush=True)
print(f"P_developmental genes matching a renamed gene's authoritative symbol: {p_dev_collisions_with_authoritative}", flush=True)
print(f"P_developmental genes at risk overall: {len(p_dev_collisions)} {p_dev_collisions}", flush=True)

with open(f"{OUT_DIR}/p_developmental_collision_check.tsv", "w") as f:
    f.write("p_developmental_gene\tmatches_renamed_old_varname\tmatches_renamed_authoritative_symbol\n")
    for g in sorted(p_dev_genes):
        f.write(f"{g}\t{g in renamed_var_names}\t{g in renamed_authoritative_symbols}\n")
print(f"Wrote {OUT_DIR}/p_developmental_collision_check.tsv", flush=True)

summary = {
    "n_vars_total": int(len(vn)),
    "n_duplicate_gene_ids": int(n_dup_gene_id),
    "n_not_found_in_biomart": int(n_not_found),
    "n_renamed_or_mismatched_by_authoritative_symbol": int(n_renamed),
    "p_developmental_primary84_path": P_DEV_PATH,
    "n_p_developmental_genes": int(len(p_dev_genes)),
    "n_p_developmental_collisions": int(len(p_dev_collisions)),
    "p_developmental_collisions": p_dev_collisions,
    "method": "authoritative BioMart external_gene_name comparison per gene_id, "
              "NOT a string-pattern (-\\d+ suffix) heuristic",
}
with open(f"{OUT_DIR}/gene_id_audit_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(f"Wrote {OUT_DIR}/gene_id_audit_summary.json")
