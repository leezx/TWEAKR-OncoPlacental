#!/usr/bin/env python3
"""
Step 6 gut re-anchor: build the 13-panel scoring gene-set inventory as
bare-Ensembl-ID lists directly usable against CRC_single_cell_atlas_2025's
var_names. Per docs/STEP6_GUT_SCORING_COMPUTE_DESIGN.md (PR #26, APPROVE
after 2 review rounds).

8 revCSC panels (27-gene primary + 28-gene extended/sensitivity, each in
up to 4 overlap-exclusion forms per the locked contract) + 5 D/F/P panels
(D_Gut-shared, F_Gut-specific, F_Colon-specific, F_SI-specific,
P_Gut-specific). Gene-ID mapping reuses PR #25's contract exactly:
var_name -> gene_id via var_id_map.tsv for D/F/P sets; revCSC used
directly via its human_ensembl_id (already Ensembl). Every panel's final
testable list is {signature Ensembl IDs} ∩ {atlas var_names}, with
n_input/n_mapped_to_ensembl/n_present_in_atlas reported per panel --
same auditable-denominator discipline as every prior gene-ID contract
this project has used.

Usage: python3 build_gut_scoring_gene_sets.py <out_dir>
(run on Argos, argos-codex env; also runs fine locally -- no scRNA data
touched here, only gene-ID text files. CRC_single_cell_atlas_2025's
var_names are only used if --atlas-h5ad is given; otherwise the
"n_present_in_atlas" column is left blank and must be filled in by the
scoring script itself against the real atlas.)
"""
import sys
import os
import csv
import argparse

REPO = "/home/zz950/TWEAKR-OncoPlacental"
REVCSC_DIR = f"{REPO}/results/06_crc_projection/revcsc_overlap_audit"
GUT_OVERLAP_DIR = f"{REPO}/results/06_crc_projection/revcsc_gut_overlap_audit"
GUT_DFP_DIR = f"{REPO}/results/04a_dfp_gut/dfp_gut_gene_sets"
VAR_ID_MAP = f"{REPO}/results/04a_dfp_gut/gene_id_audit/var_id_map.tsv"

# CLU/ASS1 are the only two revCSC-primary overlap genes found anywhere in
# the gut D/F/P sets (PR #25's audit); LY6A (extended-only) overlaps none.
CLU_ENS = None  # resolved below from revCSC_human_FINAL.tsv, not hardcoded
ASS1_ENS = None


def load_gene_set(path):
    with open(path) as f:
        return {l.strip() for l in f if l.strip()}


def load_var_id_map():
    """var_name -> gene_id (bare Ensembl, no version suffix)."""
    m = {}
    with open(VAR_ID_MAP, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["gene_id"]:
                m[row["var_name"]] = row["gene_id"]
    return m


def load_revcsc_ensembl_map():
    """Returns (primary_map, extended_only_map): human_ensembl_id -> symbol."""
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


def map_varnames_to_ensembl(var_names, var_map):
    """var_name set -> (ensembl_id set, n_mapped, n_unmapped)."""
    ens_ids = set()
    n_unmapped = 0
    for vn in var_names:
        gid = var_map.get(vn)
        if gid:
            ens_ids.add(gid)
        else:
            n_unmapped += 1
    return ens_ids, len(ens_ids), n_unmapped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", nargs="?", default=f"{REPO}/results/06_crc_projection/gut_scoring")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    var_map = load_var_id_map()
    revcsc_primary_map, revcsc_extended_only_map = load_revcsc_ensembl_map()
    revcsc_primary_ens = set(revcsc_primary_map)
    revcsc_extended_ens = revcsc_primary_ens | set(revcsc_extended_only_map)
    print(f"revCSC primary: {len(revcsc_primary_ens)} Ensembl IDs")
    print(f"revCSC extended: {len(revcsc_extended_ens)} Ensembl IDs "
          f"(adds {sorted(revcsc_extended_only_map.values())})")

    # symbol -> ensembl for CLU/ASS1 exclusion, resolved from the primary map
    sym_to_ens = {s: e for e, s in revcsc_primary_map.items()}
    clu_ens = sym_to_ens.get("CLU")
    ass1_ens = sym_to_ens.get("ASS1")
    assert clu_ens and ass1_ens, "CLU/ASS1 must resolve from revCSC_human_FINAL.tsv primary rows"
    print(f"CLU -> {clu_ens}, ASS1 -> {ass1_ens} (overlap-exclusion targets)")

    # ---- 8 revCSC panels ----
    revcsc_panels = {
        "revCSC_primary27_full": revcsc_primary_ens,
        "revCSC_primary27_minus_CLU": revcsc_primary_ens - {clu_ens},
        "revCSC_primary27_minus_ASS1": revcsc_primary_ens - {ass1_ens},
        "revCSC_primary27_minus_CLU_ASS1": revcsc_primary_ens - {clu_ens, ass1_ens},
        "revCSC_extended28_full": revcsc_extended_ens,
        "revCSC_extended28_minus_CLU": revcsc_extended_ens - {clu_ens},
        "revCSC_extended28_minus_ASS1": revcsc_extended_ens - {ass1_ens},
        "revCSC_extended28_minus_CLU_ASS1": revcsc_extended_ens - {clu_ens, ass1_ens},
    }
    expected_sizes = {
        "revCSC_primary27_full": 27, "revCSC_primary27_minus_CLU": 26,
        "revCSC_primary27_minus_ASS1": 26, "revCSC_primary27_minus_CLU_ASS1": 25,
        "revCSC_extended28_full": 28, "revCSC_extended28_minus_CLU": 27,
        "revCSC_extended28_minus_ASS1": 27, "revCSC_extended28_minus_CLU_ASS1": 26,
    }
    for name, s in revcsc_panels.items():
        assert len(s) == expected_sizes[name], f"{name}: expected {expected_sizes[name]}, got {len(s)}"

    # ---- 5 D/F/P panels: var_name -> Ensembl via var_id_map.tsv ----
    # D_Gut-shared / F_Gut-specific / P_Gut-specific were computed+written by
    # PR #25's revcsc_gut_dfp_overlap_audit.py; F_Colon-specific/F_SI-specific
    # are direct Step 4a frozen files. All are var_name-space until mapped here.
    dfp_sources = {
        "D_Gut-shared": f"{GUT_OVERLAP_DIR}/D_Gut-shared.txt",
        "F_Gut-specific": f"{GUT_OVERLAP_DIR}/F_Gut-specific.txt",
        "F_Colon-specific": f"{GUT_DFP_DIR}/F_Colon-specific.txt",
        "F_SI-specific": f"{GUT_DFP_DIR}/F_SI-specific.txt",
        "P_Gut-specific": f"{GUT_OVERLAP_DIR}/P_Gut-specific.txt",
    }

    all_panels_ens = dict(revcsc_panels)
    mapping_rows = []
    for name, path in dfp_sources.items():
        var_names = load_gene_set(path)
        ens_ids, n_mapped, n_unmapped = map_varnames_to_ensembl(var_names, var_map)
        all_panels_ens[name] = ens_ids
        mapping_rows.append({
            "panel": name, "n_input_var_names": len(var_names),
            "n_mapped_to_ensembl": n_mapped, "n_unmapped": n_unmapped,
        })
        print(f"{name}: n_input={len(var_names)}, n_mapped_to_ensembl={n_mapped}, "
              f"n_unmapped={n_unmapped}")

    for name, s in revcsc_panels.items():
        mapping_rows.append({
            "panel": name, "n_input_var_names": len(s),
            "n_mapped_to_ensembl": len(s), "n_unmapped": 0,
        })

    # ---- write out every panel as a bare-Ensembl-ID .txt ----
    for name, ens_ids in all_panels_ens.items():
        path = f"{args.out_dir}/{name}.ensembl.txt"
        with open(path, "w") as fh:
            fh.write("\n".join(sorted(ens_ids)) + "\n")
        print(f"Wrote {path} ({len(ens_ids)} genes)")

    mapping_path = f"{args.out_dir}/gut_scoring_gene_id_mapping_summary.tsv"
    with open(mapping_path, "w", newline="") as f:
        cols = ["panel", "n_input_var_names", "n_mapped_to_ensembl", "n_unmapped"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in mapping_rows:
            w.writerow(r)
    print(f"\nWrote {mapping_path}")
    print(f"\nTotal panels: {len(all_panels_ens)} (8 revCSC + 5 D/F/P)")
    print("=== DONE ===")


if __name__ == "__main__":
    main()
