#!/usr/bin/env python3
"""
Step 6 secondary analysis (PR #28): build M11's Ensembl-ID gene-set
inventory directly against CRC_single_cell_atlas_2025's var_names, and
commit the M11 x revCSC gene-overlap audit as a script/result artifact
(per round-2 review: "I would also commit the overlap audit as a
script/result artifact rather than leaving exact overlap counts only in
design prose").

M11's gene list is *already Ensembl-ID-native* in the source NMF deliver
files (confirmed directly on Argos, not assumed -- `deliver.mm_top_genes.csv`'s
M11 column and `unique_top_versions/deliver.mm_top_genes.unique{50,100,200}.csv`
both contain bare ENSG IDs), so no symbol->Ensembl mapping step is needed
for M11 itself (unlike the D/F/P panels, which started as var_name lists).

For each of the 3 M11 cutoffs (top50 = deliver.mm_top_genes.csv's first 50
non-null M11 values; top100/200 = the corresponding unique_top_versions
file's M11 column), writes both the full panel and the
revCSC_primary27-overlap-excluded panel (the primary variant for every
M11<->revCSC comparison per round-1 review).

Usage: python3 build_m11_gene_sets.py <out_dir>
(run on Argos; no scRNA data touched, only gene-ID text files -- but reads
files that only exist on Argos, so this does not run locally.)
"""
import sys
import os
import csv
import argparse

REPO = "/home/zz950/TWEAKR-OncoPlacental"
NMF_DIR = "/home/zz950/DATA/scRNAseq/meta_study/CRC_single_cell_atlas_2025/NMF"
TOP50_SRC = f"{NMF_DIR}/metamodule_fnmf/MM_alt_clean_byscore_wardD2/deliver.mm_top_genes.csv"
TOP100_SRC = f"{NMF_DIR}/unique_top_versions/deliver.mm_top_genes.unique100.csv"
TOP200_SRC = f"{NMF_DIR}/unique_top_versions/deliver.mm_top_genes.unique200.csv"
REVCSC_DIR = f"{REPO}/results/06_crc_projection/revcsc_overlap_audit"

# revCSC_primary27_full genes M11 shares -- confirmed directly (round 1
# review + independently re-verified against revCSC_human_FINAL.tsv here):
# all 5 are include_primary rows, i.e. members of revCSC_primary27_full
# (and every _minus_CLU/_minus_ASS1/_minus_CLU_ASS1 variant, since none of
# the 5 is CLU or ASS1). PMEPA1 (6th, top200-only) is also include_primary.
M11_REVCSC_OVERLAP_SYMBOLS = {"ANXA1", "KRT18", "SFN", "TMSB4X", "TNFRSF12A"}
M11_REVCSC_OVERLAP_SYMBOLS_TOP200 = M11_REVCSC_OVERLAP_SYMBOLS | {"PMEPA1"}


def load_m11_column(path, n_top=None):
    """Reads the M11 column of a deliver.mm_top_genes*.csv (Ensembl IDs,
    already ranked top-to-bottom by construction), optionally truncated to
    the first n_top non-empty values."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        vals = []
        for row in reader:
            v = row.get("M11", "").strip()
            if v:
                vals.append(v)
    if n_top is not None:
        vals = vals[:n_top]
    return vals


def load_revcsc_symbol_ensembl_map():
    """human_ensembl_id -> human_ortholog_symbol, primary27 rows only."""
    m = {}
    with open(f"{REVCSC_DIR}/revCSC_human_FINAL.tsv", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["inclusion_decision"] == "include_primary":
                m[row["human_ensembl_id"].strip()] = row["human_ortholog_symbol"].strip()
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", nargs="?", default=f"{REPO}/results/06_crc_projection/gut_scoring")
    ap.add_argument("--audit-out", default=f"{REPO}/results/06_crc_projection/m11_revcsc_overlap_audit")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.audit_out, exist_ok=True)

    revcsc_map = load_revcsc_symbol_ensembl_map()  # ens -> symbol, primary27
    revcsc_sym_to_ens = {s: e for e, s in revcsc_map.items()}
    overlap_ens_top50_100 = {revcsc_sym_to_ens[s] for s in M11_REVCSC_OVERLAP_SYMBOLS}
    overlap_ens_top200 = {revcsc_sym_to_ens[s] for s in M11_REVCSC_OVERLAP_SYMBOLS_TOP200}
    print(f"revCSC primary27 overlap genes (top50/100 M11): "
          f"{sorted(M11_REVCSC_OVERLAP_SYMBOLS)} -> {sorted(overlap_ens_top50_100)}")
    print(f"revCSC primary27 overlap genes (top200 M11, adds PMEPA1): "
          f"{sorted(M11_REVCSC_OVERLAP_SYMBOLS_TOP200)} -> {sorted(overlap_ens_top200)}")

    cutoffs = {
        "M11_top50": (load_m11_column(TOP50_SRC, n_top=50), overlap_ens_top50_100),
        "M11_top100": (load_m11_column(TOP100_SRC), overlap_ens_top50_100),
        "M11_top200": (load_m11_column(TOP200_SRC), overlap_ens_top200),
    }

    audit_rows = []
    for name, (genes, overlap_ens) in cutoffs.items():
        genes_set = set(genes)
        assert len(genes_set) == len(genes), f"{name}: duplicate Ensembl IDs in source column"
        found_overlap = genes_set & (overlap_ens_top50_100 | overlap_ens_top200)
        excl_genes = genes_set - overlap_ens
        full_path = f"{args.out_dir}/{name}_full.ensembl.txt"
        excl_path = f"{args.out_dir}/{name}_minus_revCSC_overlap.ensembl.txt"
        with open(full_path, "w") as fh:
            fh.write("\n".join(sorted(genes_set)) + "\n")
        with open(excl_path, "w") as fh:
            fh.write("\n".join(sorted(excl_genes)) + "\n")
        print(f"{name}: n_full={len(genes_set)}, revCSC_primary27-overlap="
              f"{len(found_overlap)} {sorted(revcsc_map.get(e, e) for e in found_overlap)}, "
              f"n_minus_revCSC_overlap={len(excl_genes)}")
        print(f"  Wrote {full_path}, {excl_path}")
        audit_rows.append({
            "m11_cutoff": name, "n_full": len(genes_set),
            "n_overlap_with_revCSC_primary27": len(found_overlap),
            "overlap_genes": ",".join(sorted(revcsc_map.get(e, e) for e in found_overlap)),
            "n_minus_revCSC_overlap": len(excl_genes),
        })

    audit_path = f"{args.audit_out}/m11_revcsc_overlap_audit.tsv"
    with open(audit_path, "w", newline="") as f:
        cols = ["m11_cutoff", "n_full", "n_overlap_with_revCSC_primary27",
                "overlap_genes", "n_minus_revCSC_overlap"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in audit_rows:
            w.writerow(r)
    print(f"\nWrote {audit_path}")
    print("=== DONE ===")


if __name__ == "__main__":
    main()
