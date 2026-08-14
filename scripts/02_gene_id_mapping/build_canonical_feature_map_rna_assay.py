#!/usr/bin/env python3
"""
RNA-assay-scoped canonical_feature_map, built to fix a real gap found while
building HDMA per-sample pseudobulk (Step 4): the original
canonical_feature_map/<organ>_canonical_feature_map.tsv (build_canonical_
feature_map.py) was built from rownames(obj) -- the Seurat object's
DEFAULT assay (decontX/SCT, a QC-filtered reduced gene set) -- NOT the RNA
assay's raw "counts" layer, which has more genes (e.g. Adrenal: 25,314
default vs 28,375 RNA). Pseudobulk needs true raw counts (RNA assay), so
it needs a mapping covering RNA's full feature space.

Confirmed the default-assay set is a strict subset of the RNA-assay set in
every organ (0 genes in default missing from RNA) -- this script's output
is a superset of the original canonical_feature_map's mappings, not a
divergent one; same original_feature always maps to the same
canonical_symbol where both cover it (verified below, not assumed).

Reuses the same HGNC/biotype/ensembl-display-name lookups as Step 2's
original mapping, extended with results/02_gene_id_mapping/
new_ensg_biotypes.tsv for the 186 ENSG IDs present only in the RNA assay
(queried via query_biotype.py, merged into union_ensg_biotypes.tsv).

Inputs:
  results/02_gene_id_mapping/gene_lists/<organ>_rna_assay_all_genes.txt
  results/02_gene_id_mapping/union_ensg_biotypes.tsv (now includes the 186 new IDs)
  results/02_gene_id_mapping/unresolved_protein_coding_display_names.tsv
  DATA/1.Databases/HGNC_gene_id_mapping/processed/v0.1/hgnc_symbol_ensembl_map.tsv

Output:
  results/02_gene_id_mapping/canonical_feature_map_rna_assay/<organ>_canonical_feature_map.tsv
  results/02_gene_id_mapping/canonical_feature_map_rna_assay/collision_report.tsv
  results/02_gene_id_mapping/canonical_feature_map_rna_assay/SUMMARY.md
"""
import csv
import os
from collections import defaultdict, Counter

REPO = "/Volumes/Stelligen_SSD/Stelligen/PR/TWEAKR-OncoPlacental"
DATA = "/Volumes/Stelligen_SSD/Stelligen/DATA"
GENE_LISTS = f"{REPO}/results/02_gene_id_mapping/gene_lists"
OUT_DIR = f"{REPO}/results/02_gene_id_mapping/canonical_feature_map_rna_assay"
OLD_MAP_DIR = f"{REPO}/results/02_gene_id_mapping/canonical_feature_map"
os.makedirs(OUT_DIR, exist_ok=True)

ORGANS = ["Adrenal", "Thyroid", "Spleen", "Thymus", "Liver", "Skin", "StomachEsophagus"]


def load_hgnc_ensg_to_symbol():
    m = {}
    path = f"{DATA}/1.Databases/HGNC_gene_id_mapping/processed/v0.1/hgnc_symbol_ensembl_map.tsv"
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            ens = row["ensembl_id"]
            if ens:
                m[ens] = row["symbol"]
    return m


def load_biotype():
    m = {}
    path = f"{REPO}/results/02_gene_id_mapping/union_ensg_biotypes.tsv"
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            m[row["ensembl_id"]] = row["biotype"]
    return m


def load_ensembl_display_name():
    m = {}
    path = f"{REPO}/results/02_gene_id_mapping/unresolved_protein_coding_display_names.tsv"
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row["ensembl_display_name"]:
                m[row["ensembl_id"]] = row["ensembl_display_name"]
    return m


def load_old_map(organ):
    """The original default-assay map, for a direct consistency check."""
    path = f"{OLD_MAP_DIR}/{organ}_canonical_feature_map.tsv"
    m = {}
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            m[row["original_feature"]] = row["canonical_symbol"]
    return m


def main():
    hgnc_map = load_hgnc_ensg_to_symbol()
    biotype_map = load_biotype()
    ensembl_name_map = load_ensembl_display_name()

    status_counts_overall = Counter()
    collision_rows = []
    n_collisions_total = 0
    n_consistency_checked = 0
    n_consistency_mismatches = 0

    for organ in ORGANS:
        genes_path = f"{GENE_LISTS}/{organ}_rna_assay_all_genes.txt"
        with open(genes_path) as f:
            genes = [l.strip() for l in f if l.strip()]

        old_map = load_old_map(organ)

        rows = []
        by_symbol = defaultdict(list)
        for g in genes:
            is_ensg = g.startswith("ENSG")
            if not is_ensg:
                canonical = g
                status = "native_symbol"
            elif g in hgnc_map:
                canonical = hgnc_map[g]
                status = "mapped_via_hgnc"
            elif g in ensembl_name_map:
                canonical = ensembl_name_map[g]
                status = "mapped_via_ensembl_display_name"
            else:
                canonical = g  # kept as its own ENSG ID, not dropped
                status = "unmapped_kept_as_ensembl_id"

            bt = biotype_map.get(g, "") if is_ensg else ""
            rows.append((g, is_ensg, canonical, status, bt))
            by_symbol[canonical].append(g)
            status_counts_overall[status] += 1

            if g in old_map:
                n_consistency_checked += 1
                if old_map[g] != canonical:
                    n_consistency_mismatches += 1
                    print(f"  CONSISTENCY MISMATCH {organ} {g}: old={old_map[g]} new={canonical}")

        out_path = f"{OUT_DIR}/{organ}_canonical_feature_map.tsv"
        with open(out_path, "w", newline="") as fo:
            w = csv.writer(fo, delimiter="\t")
            w.writerow(["original_feature", "is_ensg", "canonical_symbol", "mapping_status", "biotype"])
            for row in rows:
                w.writerow(row)

        organ_collisions = {sym: feats for sym, feats in by_symbol.items() if len(feats) > 1}
        n_collisions_total += len(organ_collisions)
        for sym, feats in sorted(organ_collisions.items(), key=lambda kv: -len(kv[1])):
            collision_rows.append((organ, sym, len(feats), ",".join(feats)))

        print(f"{organ}: {len(genes)} RNA-assay features -> {len(by_symbol)} unique canonical_symbols "
              f"({len(organ_collisions)} collisions)")

    with open(f"{OUT_DIR}/collision_report.tsv", "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["organ", "canonical_symbol", "n_original_features", "original_features"])
        for row in collision_rows:
            w.writerow(row)

    with open(f"{OUT_DIR}/SUMMARY.md", "w") as fo:
        fo.write("# Canonical feature map (RNA assay) — build summary\n\n")
        fo.write("Extends `../canonical_feature_map/` (built from each Seurat object's "
                  "default assay, decontX/SCT) to cover the RNA assay's full raw-counts "
                  "feature space, needed for pseudobulk. See module docstring for the "
                  "gap this fixes.\n\n")
        fo.write(f"Consistency check against the original map: {n_consistency_checked} "
                  f"shared (organ, original_feature) pairs checked, "
                  f"{n_consistency_mismatches} mismatches "
                  f"({'CLEAN — every gene the two maps share got the same canonical_symbol' if n_consistency_mismatches == 0 else 'SEE MISMATCHES ABOVE, NOT CLEAN'}).\n\n")
        fo.write("Per-organ mapping_status counts (summed across all 7 organs):\n\n")
        fo.write("| mapping_status | count |\n|---|---|\n")
        for status, n in status_counts_overall.most_common():
            fo.write(f"| {status} | {n} |\n")
        fo.write(f"\nTotal symbol collisions across all organs: {n_collisions_total} "
                  f"(see `collision_report.tsv`).\n")
        fo.write("\nCollisions are reported, not resolved, by this script. Pseudobulk "
                  "construction (build_hdma_pseudobulk.R) sums counts across colliding "
                  "features by canonical_symbol.\n")

    print(f"\nConsistency check: {n_consistency_checked} shared pairs, {n_consistency_mismatches} mismatches")
    print(f"Total collisions across all organs: {n_collisions_total}")
    print(f"Wrote: {OUT_DIR}/*_canonical_feature_map.tsv, collision_report.tsv, SUMMARY.md")


if __name__ == "__main__":
    main()
