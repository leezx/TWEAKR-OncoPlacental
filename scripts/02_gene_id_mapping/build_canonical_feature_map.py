#!/usr/bin/env python3
"""
Build the per-organ canonical_feature_map for HDMA (PR #4 review item e),
and detect symbol collisions introduced by the ENSG->symbol mapping step
(review item d).

For each of the 7 HDMA organs, every original rowname (feature) is one of:
  - already a gene symbol (not ENSG-format)              -> kept as-is
  - an ENSG ID that resolves to an HGNC-approved symbol   -> mapped to that symbol
  - an ENSG ID with an Ensembl-native display_name        -> mapped to that name
  - an ENSG ID with no symbol anywhere we checked         -> kept as its own
                                                              ENSG ID (not dropped)

A "collision" is >1 distinct original_feature in the same organ landing on
the same canonical_symbol after this mapping (e.g. an ENSG's HGNC symbol
already exists as a separate native-symbol feature in that organ's matrix,
or two different ENSG IDs mapping to the same symbol). These need an
explicit aggregation/drop rule before Step 3 pseudobulk can sum counts by
symbol; this script only detects and reports them, does not resolve them.

Inputs (all already pulled from Argos / downloaded locally):
  results/02_gene_id_mapping/gene_lists/<organ>_all_genes.txt
  results/02_gene_id_mapping/union_ensg_biotypes.tsv
  results/02_gene_id_mapping/unresolved_protein_coding_display_names.tsv
  DATA/1.Databases/HGNC_gene_id_mapping/processed/v0.1/hgnc_symbol_ensembl_map.tsv

Output:
  results/02_gene_id_mapping/canonical_feature_map/<organ>_canonical_feature_map.tsv
    columns: original_feature, is_ensg, canonical_symbol, mapping_status, biotype
  results/02_gene_id_mapping/canonical_feature_map/collision_report.tsv
    columns: organ, canonical_symbol, n_original_features, original_features (comma-joined)
  results/02_gene_id_mapping/canonical_feature_map/SUMMARY.md
"""
import csv
import os
from collections import defaultdict, Counter

REPO = "/Volumes/Stelligen_SSD/Stelligen/PR/TWEAKR-OncoPlacental"
DATA = "/Volumes/Stelligen_SSD/Stelligen/DATA"
GENE_LISTS = f"{REPO}/results/02_gene_id_mapping/gene_lists"
OUT_DIR = f"{REPO}/results/02_gene_id_mapping/canonical_feature_map"
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


def main():
    hgnc_map = load_hgnc_ensg_to_symbol()
    biotype_map = load_biotype()
    ensembl_name_map = load_ensembl_display_name()

    status_counts_overall = Counter()
    collision_rows = []
    n_collisions_total = 0

    for organ in ORGANS:
        genes_path = f"{GENE_LISTS}/{organ}_all_genes.txt"
        with open(genes_path) as f:
            genes = [l.strip() for l in f if l.strip()]

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

        print(f"{organ}: {len(genes)} features -> {len(by_symbol)} unique canonical_symbols "
              f"({len(organ_collisions)} collisions)")

    with open(f"{OUT_DIR}/collision_report.tsv", "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["organ", "canonical_symbol", "n_original_features", "original_features"])
        for row in collision_rows:
            w.writerow(row)

    with open(f"{OUT_DIR}/SUMMARY.md", "w") as fo:
        fo.write("# Canonical feature map — build summary\n\n")
        fo.write("Per-organ mapping_status counts (summed across all 7 organs):\n\n")
        fo.write("| mapping_status | count |\n|---|---|\n")
        for status, n in status_counts_overall.most_common():
            fo.write(f"| {status} | {n} |\n")
        fo.write(f"\nTotal symbol collisions across all organs: {n_collisions_total} "
                  f"(see `collision_report.tsv` for the full list — one row per "
                  f"`(organ, canonical_symbol)` pair with >1 original feature).\n")
        fo.write("\nCollisions are reported, not resolved, by this script. Step 3 "
                  "pseudobulk construction must pick an explicit aggregation rule "
                  "(e.g. sum counts across colliding features) before using "
                  "canonical_symbol as the feature key.\n")

    print(f"\nTotal collisions across all organs: {n_collisions_total}")
    print(f"Wrote: {OUT_DIR}/*_canonical_feature_map.tsv, collision_report.tsv, SUMMARY.md")


if __name__ == "__main__":
    main()
