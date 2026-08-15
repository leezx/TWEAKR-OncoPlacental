#!/usr/bin/env python3
"""
Step 6 re-anchor: revCSC x GUT-D/F/P gene-overlap audit.

Per Worklog.md's "what's next": re-run the M11/revCSC CRC decomposition
(Step 6, PR #16/#17, design fully approved but compute never run) against
the new Gut-specific D/F/P instead of the pan-organ one. This is the first
concrete step -- same discipline as the original revcsc_dfp_overlap_audit.py
(PR #16/#17), but against F_Colon-developmental/F_SI-developmental/
D_Colon-shared/D_SI-shared/P_Colon-specific/P_SI-specific instead of
D_shared_FINAL/F_specific_FINAL/P_specific_FINAL/F_developmental_<Organ>.

Uses the FINAL frozen, ortholog-audited revCSC sets (not re-derived from
the raw CSC table): revCSC_symbols.primary27.txt (primary, Compara
one2one-confirmed) and revCSC_symbols.extended28.txt (adds Ly6a, sensitivity
only) -- per PR #17 round 2's frozen labeling.

Cross-dataset gene-ID contract: same as the Step 4a adult-validation work
(PR #23/#24) -- the gut F/D/P gene lists are raw Gut Cell Atlas var_names
(may include anndata duplicate-symbol artifacts like IGF2-1), so overlap
against revCSC's plain human gene symbols is computed via each gut gene's
authoritative_symbol (PR #21's var_id_map.tsv), not the raw var_name
string -- avoiding exactly the IGF2-1-vs-IGF2 mismatch class.

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

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else f"{REPO}/results/06_crc_projection/revcsc_gut_overlap_audit"
os.makedirs(OUT_DIR, exist_ok=True)

GUT_SETS = [
    "F_Colon-developmental", "F_SI-developmental",
    "D_Colon-shared", "D_SI-shared",
    "P_Colon-specific", "P_SI-specific",
    "F_Gut-core",
]


def load_gene_set(path):
    with open(path) as f:
        return {l.strip() for l in f if l.strip()}


def load_var_id_map():
    """var_name -> authoritative_symbol (only where BioMart-resolved)."""
    m = {}
    with open(VAR_ID_MAP, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["found_in_biomart"] == "True":
                m[row["var_name"]] = row["authoritative_symbol"]
    return m


def to_symbols(var_names, var_map):
    """Map a gut D/F/P var_name set to its authoritative_symbol set.
    var_names with no BioMart resolution are dropped (reported as
    unresolved, not silently treated as non-overlapping -- consistent
    with the NOT_TESTABLE discipline from PR #23/#24)."""
    resolved = set()
    unresolved = []
    for vn in var_names:
        if vn in var_map:
            resolved.add(var_map[vn])
        else:
            unresolved.append(vn)
    return resolved, unresolved


def main():
    primary = load_gene_set(f"{REVCSC_DIR}/revCSC_symbols.primary27.txt")
    extended = load_gene_set(f"{REVCSC_DIR}/revCSC_symbols.extended28.txt")
    print(f"revCSC primary (Compara one2one, 27 genes): {sorted(primary)}")
    print(f"revCSC extended (28 genes, adds Ly6a/LY6A): {sorted(extended)}")

    var_map = load_var_id_map()
    print(f"var_id_map: {len(var_map)} GCA var_names with a BioMart-resolved authoritative_symbol")

    rows = []
    for set_name in GUT_SETS:
        var_names = load_gene_set(f"{GUT_DFP_DIR}/{set_name}.txt")
        symbols, unresolved = to_symbols(var_names, var_map)
        overlap_primary = sorted(primary & symbols)
        overlap_extended_only = sorted((extended - primary) & symbols)  # Ly6a/LY6A only, if present
        rows.append({
            "gut_set": set_name, "n_var_names": len(var_names),
            "n_resolved_symbols": len(symbols), "n_unresolved": len(unresolved),
            "overlap_primary27_n": len(overlap_primary),
            "overlap_primary27_genes": ";".join(overlap_primary),
            "overlap_extended_only_n": len(overlap_extended_only),
            "overlap_extended_only_genes": ";".join(overlap_extended_only),
        })
        print(f"{set_name}: n_var_names={len(var_names)}, n_resolved={len(symbols)}, "
              f"n_unresolved={len(unresolved)}, overlap_primary27={overlap_primary}, "
              f"overlap_extended_only={overlap_extended_only}")

    out_path = f"{OUT_DIR}/revcsc_gut_dfp_overlap_audit.tsv"
    with open(out_path, "w", newline="") as f:
        cols = ["gut_set", "n_var_names", "n_resolved_symbols", "n_unresolved",
                "overlap_primary27_n", "overlap_primary27_genes",
                "overlap_extended_only_n", "overlap_extended_only_genes"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
