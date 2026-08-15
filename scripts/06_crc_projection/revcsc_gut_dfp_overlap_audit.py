#!/usr/bin/env python3
"""
Step 6 re-anchor: revCSC x GUT-D/F/P gene-overlap audit.

Per Worklog.md's "what's next": re-run the M11/revCSC CRC decomposition
(Step 6, PR #16/#17, design fully approved but compute never run) against
the new Gut-specific D/F/P instead of the pan-organ one. This is the first
concrete step -- same discipline as the original revcsc_dfp_overlap_audit.py
(PR #16/#17).

**PR #25 round-1 correction**: the F-arm input must be `F_{region}-specific`
(F ^ NOT P, the same P-deduplicated role `F_specific_FINAL` played in the
original design -- the first version of this script used the raw
`F_{region}-developmental` instead, which still contains the D-overlap
genes and would have broken the D/F/P mutual-exclusivity the original
Step 6 design explicitly required).

**PR #25 round-2 correction**: the original design scored BOTH a global
`F_specific_FINAL` (coarse/secondary) AND 7 per-organ lineage modules
(primary, region-resolved) -- round-1's fix mapped only the lineage-module
role. Fixed: `F_Gut-specific` computed here as `(F_Colon-developmental U
F_SI-developmental) \\ P_developmental` -- the gut-scoped analogue of
`F_specific_FINAL`'s own definition.

**PR #25 round-3 correction (two issues)**:
1. The approved Step 6 design has exactly ONE global D axis
   (`D_shared_FINAL`) and ONE global P axis (`P_specific_FINAL`) -- only F
   was ever organ/lineage-resolved. Mapping the global D/P axes directly to
   two REGIONAL axes each (`D_Colon-shared`+`D_SI-shared`,
   `P_Colon-specific`+`P_SI-specific`) silently turned one global axis into
   two regional ones, which the original design never had for D/P. Fixed:
   `D_Gut-shared` = `D_Colon-shared UNION D_SI-shared` (8 unique genes,
   computed here) and `P_Gut-specific` = `P_Colon-specific INTERSECT
   P_SI-specific` (76 genes, computed here -- equivalent to `P_developmental
   \\ (F_Colon-developmental U F_SI-developmental)`, i.e. genes that are
   P-developmental but fetal-developmental in NEITHER region) are now the
   correct global D/P axes; the 4 regional D/P sets remain available as
   secondary descriptive views, explicitly not independent replication
   panels (per PR #24's own finding that P_Colon-specific/P_SI-specific
   already share 76 of 79/80 genes).
2. The overlap test was gated on BioMart `authoritative_symbol` resolution,
   leaving 44-65 genes per set NOT_TESTABLE even though the underlying GCA
   `gene_id` (Ensembl) is known and unique (zero duplicate gene_ids across
   all 33,538 GCA variables, per PR #21's audit) and revCSC's own frozen
   provenance table (`revCSC_human_FINAL.tsv`) already carries a verified
   human Ensembl ID for every primary/extended member. Fixed: overlap is
   now tested primarily via GCA `gene_id` == revCSC `human_ensembl_id`
   (Ensembl-ID identity, not symbol-string matching), with
   `authoritative_symbol` retained only for human-readable reporting, not
   as the testability gate -- same lesson already applied to GTEx matching
   in PR #23/#24.

Uses the FINAL frozen, ortholog-audited revCSC provenance table
(`revCSC_human_FINAL.tsv`) and its primary27/extended28 symbol lists (for
reporting) -- not re-derived from the raw CSC table, per PR #17 round 2's
frozen labeling.

Usage: python3 revcsc_gut_dfp_overlap_audit.py <out_dir>
(run on Argos, argos-codex env)
"""
import sys
import os
import csv

REPO = "/home/zz950/TWEAKR-OncoPlacental"
REVCSC_DIR = f"{REPO}/results/06_crc_projection/revcsc_overlap_audit"
GUT_DFP_DIR = f"{REPO}/results/04a_dfp_gut/dfp_gut_gene_sets"
VAR_ID_MAP = f"{REPO}/results/04a_dfp_gut/gene_id_audit/var_id_map.tsv"
DFP_DIR = f"{REPO}/results/04_dfp_signature/dfp_gene_sets"

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else f"{REPO}/results/06_crc_projection/revcsc_gut_overlap_audit"
os.makedirs(OUT_DIR, exist_ok=True)

# Pre-existing gut D/F/P files loaded directly.
GUT_FILE_SETS = ["F_Colon-specific", "F_SI-specific", "D_Colon-shared", "D_SI-shared",
                  "P_Colon-specific", "P_SI-specific", "F_Gut-core"]
# Computed sets (global axes + F_Gut-core's companions), built in main().
GUT_COMPUTED_SETS = ["F_Gut-specific", "D_Gut-shared", "P_Gut-specific"]


def load_gene_set(path):
    with open(path) as f:
        return {l.strip() for l in f if l.strip()}


def load_var_id_map():
    """var_name -> (gene_id, authoritative_symbol_or_None)."""
    m = {}
    with open(VAR_ID_MAP, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sym = row["authoritative_symbol"] if row["found_in_biomart"] == "True" else None
            m[row["var_name"]] = (row["gene_id"], sym)
    return m


def load_revcsc_ensembl_map():
    """human_ensembl_id -> human_ortholog_symbol, restricted to include_primary
    rows (the frozen primary27 set) plus the one extended-only row (Ly6a)."""
    primary_map, extended_only_map = {}, {}
    with open(f"{REVCSC_DIR}/revCSC_human_FINAL.tsv", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            ens = row["human_ensembl_id"].strip()
            sym = row["human_ortholog_symbol"].strip()
            if not ens:
                continue
            if row["inclusion_decision"] == "include_primary":
                primary_map[ens] = sym
            elif row["inclusion_decision"] == "include_extended_sensitivity_only":
                extended_only_map[ens] = sym
    return primary_map, extended_only_map


def resolve_gut_set(var_names, var_map):
    """Returns (ensembl_ids: set, n_ensembl_resolved, n_no_gene_id,
    symbol_report: dict ensembl_id->best-known symbol for readability).
    Ensembl gene_id is the PRIMARY identity test (round-3 fix); a var_name
    is only truly untestable if var_id_map has no gene_id for it at all
    (essentially never, per PR #21's audit -- reported explicitly if found)."""
    ensembl_ids = set()
    symbol_report = {}
    n_no_gene_id = 0
    for vn in var_names:
        gene_id, sym = var_map.get(vn, (None, None))
        if not gene_id:
            n_no_gene_id += 1
            continue
        ensembl_ids.add(gene_id)
        if sym:
            symbol_report[gene_id] = sym
    return ensembl_ids, len(ensembl_ids), n_no_gene_id, symbol_report


def main():
    revcsc_primary_map, revcsc_extended_only_map = load_revcsc_ensembl_map()
    revcsc_primary_ens = set(revcsc_primary_map)
    revcsc_extended_only_ens = set(revcsc_extended_only_map)
    print(f"revCSC primary (Ensembl-ID identity, {len(revcsc_primary_ens)} genes): "
          f"{sorted(revcsc_primary_map.values())}")
    print(f"revCSC extended-only (adds Ly6a/LY6A, {len(revcsc_extended_only_ens)} genes): "
          f"{sorted(revcsc_extended_only_map.values())}")

    var_map = load_var_id_map()
    print(f"var_id_map: {len(var_map)} GCA var_names with a gene_id")

    gut_sets = {name: load_gene_set(f"{GUT_DFP_DIR}/{name}.txt") for name in GUT_FILE_SETS}

    # ---- Computed global axes (round-2/round-3 fixes) ----
    f_colon_dev = load_gene_set(f"{GUT_DFP_DIR}/F_Colon-developmental.txt")
    f_si_dev = load_gene_set(f"{GUT_DFP_DIR}/F_SI-developmental.txt")
    p_developmental = load_gene_set(f"{DFP_DIR}/P_developmental_primary84.txt")

    f_gut_specific = (f_colon_dev | f_si_dev) - p_developmental
    d_gut_shared = gut_sets["D_Colon-shared"] | gut_sets["D_SI-shared"]
    p_gut_specific = gut_sets["P_Colon-specific"] & gut_sets["P_SI-specific"]

    print(f"F_Gut-specific = (F_Colon-developmental[{len(f_colon_dev)}] u "
          f"F_SI-developmental[{len(f_si_dev)}]) \\ P_developmental[{len(p_developmental)}] "
          f"= {len(f_gut_specific)} genes")
    print(f"D_Gut-shared = D_Colon-shared[{len(gut_sets['D_Colon-shared'])}] u "
          f"D_SI-shared[{len(gut_sets['D_SI-shared'])}] = {len(d_gut_shared)} genes")
    print(f"P_Gut-specific = P_Colon-specific[{len(gut_sets['P_Colon-specific'])}] n "
          f"P_SI-specific[{len(gut_sets['P_SI-specific'])}] = {len(p_gut_specific)} genes")

    for name, geneset in [("F_Gut-specific", f_gut_specific),
                           ("D_Gut-shared", d_gut_shared),
                           ("P_Gut-specific", p_gut_specific)]:
        path = f"{OUT_DIR}/{name}.txt"
        with open(path, "w") as fh:
            fh.write("\n".join(sorted(geneset)) + "\n")
        print(f"Wrote {path}")
        gut_sets[name] = geneset

    # ---- Overlap test, Ensembl-ID-primary (round-3 fix #2) ----
    all_set_names = GUT_FILE_SETS + GUT_COMPUTED_SETS
    rows = []
    for set_name in all_set_names:
        var_names = gut_sets[set_name]
        ensembl_ids, n_resolved, n_no_gene_id, symbol_report = resolve_gut_set(var_names, var_map)

        overlap_primary_ens = sorted(ensembl_ids & revcsc_primary_ens)
        overlap_extended_only_ens = sorted(ensembl_ids & revcsc_extended_only_ens)
        overlap_primary_syms = [revcsc_primary_map[e] for e in overlap_primary_ens]
        overlap_extended_only_syms = [revcsc_extended_only_map[e] for e in overlap_extended_only_ens]

        rows.append({
            "gut_set": set_name, "n_var_names": len(var_names),
            "n_ensembl_resolved": n_resolved, "n_no_gene_id": n_no_gene_id,
            "overlap_primary27_n": len(overlap_primary_ens),
            "overlap_primary27_genes": ";".join(overlap_primary_syms),
            "overlap_extended_only_n": len(overlap_extended_only_ens),
            "overlap_extended_only_genes": ";".join(overlap_extended_only_syms),
        })
        print(f"{set_name}: n_var_names={len(var_names)}, n_ensembl_resolved={n_resolved}, "
              f"n_no_gene_id={n_no_gene_id}, overlap_primary27={overlap_primary_syms}, "
              f"overlap_extended_only={overlap_extended_only_syms}")

    out_path = f"{OUT_DIR}/revcsc_gut_dfp_overlap_audit.tsv"
    with open(out_path, "w", newline="") as f:
        cols = ["gut_set", "n_var_names", "n_ensembl_resolved", "n_no_gene_id",
                "overlap_primary27_n", "overlap_primary27_genes",
                "overlap_extended_only_n", "overlap_extended_only_genes"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
