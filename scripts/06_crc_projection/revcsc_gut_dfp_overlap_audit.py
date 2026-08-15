#!/usr/bin/env python3
"""
Step 6 re-anchor: revCSC x GUT-D/F/P gene-overlap audit.

Per Worklog.md's "what's next": re-run the M11/revCSC CRC decomposition
(Step 6, PR #16/#17, design fully approved but compute never run) against
the new Gut-specific D/F/P instead of the pan-organ one. This is the first
concrete step -- same discipline as the original revcsc_dfp_overlap_audit.py
(PR #16/#17), against D_Colon-shared/D_SI-shared/P_Colon-specific/
P_SI-specific/F_Colon-specific/F_SI-specific instead of
D_shared_FINAL/F_specific_FINAL/P_specific_FINAL/F_developmental_<Organ>.

**PR #25 round-1 correction**: the F-arm input must be `F_{region}-specific`
(F ^ NOT P, the same P-deduplicated role `F_specific_FINAL` played in the
original design -- reviewer caught that the first version of this script
used the raw `F_{region}-developmental` instead, which still contains the
D-overlap genes and would have broken the D/F/P mutual-exclusivity the
original Step 6 design explicitly required). `F_{region}-developmental`
is intentionally NOT used here for that reason.

**PR #25 round-2 correction**: the original Step 6 design scored BOTH a
global `F_specific_FINAL` set (used as a coarse, secondary summary) AND
the 7 per-organ lineage modules (the primary, region-resolved
interpretation) -- round-1's fix mapped only the lineage-module role
(to `F_Colon-specific`/`F_SI-specific`) and dropped the global-summary
role entirely. Fixed: `F_Gut-specific` is now computed here (not
hardcoded) as `(F_Colon-developmental ∪ F_SI-developmental) \\
P_developmental` -- the exact gut-scoped analogue of
`F_specific_FINAL`'s own definition (`F-specific = F-developmental AND
NOT P-developmental`, just unioned across both gut regions first). This
plays the coarse/global role; `F_Colon-specific`/`F_SI-specific` remain
the primary, region-resolved interpretation; `F_Gut-core`
(Colon∩SI, NOT P-deduplicated) remains a separate tertiary
concordance/core summary, not a substitute for either.

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

DFP_DIR = f"{REPO}/results/04_dfp_signature/dfp_gene_sets"

GUT_SETS = [
    "F_Colon-specific", "F_SI-specific",  # P-deduplicated Layer 2 F input (see docstring correction)
    "D_Colon-shared", "D_SI-shared",
    "P_Colon-specific", "P_SI-specific",
    "F_Gut-core",  # NOT P-deduplicated (F_Colon-developmental ^ F_SI-developmental, no P_developmental
                   # filter at all) -- reported for completeness only, not a Layer 2 scoring input
                   # (per docs/STEP4A_GUT_ADULT_VALIDATION.md's own estimand-mismatch finding for
                   # this same set)
    # F_Gut-specific is computed in main() (not a pre-existing file) and appended below --
    # the coarse/global P-deduplicated analogue of the original F_specific_FINAL.
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

    # F_Gut-specific = (F_Colon-developmental UNION F_SI-developmental) MINUS P_developmental --
    # computed here, not hardcoded, the gut-scoped analogue of F_specific_FINAL
    # (round-2 fix: the original delta dropped this coarse/global role entirely).
    f_colon_dev = load_gene_set(f"{GUT_DFP_DIR}/F_Colon-developmental.txt")
    f_si_dev = load_gene_set(f"{GUT_DFP_DIR}/F_SI-developmental.txt")
    p_developmental = load_gene_set(f"{DFP_DIR}/P_developmental_primary84.txt")
    f_gut_specific = (f_colon_dev | f_si_dev) - p_developmental
    print(f"F_Gut-specific = (F_Colon-developmental[{len(f_colon_dev)}] u "
          f"F_SI-developmental[{len(f_si_dev)}]) \\ P_developmental[{len(p_developmental)}] "
          f"= {len(f_gut_specific)} genes")
    f_gut_specific_path = f"{OUT_DIR}/F_Gut-specific.txt"
    with open(f_gut_specific_path, "w") as fh:
        fh.write("\n".join(sorted(f_gut_specific)) + "\n")
    print(f"Wrote {f_gut_specific_path}")

    gut_set_gene_lists = {name: load_gene_set(f"{GUT_DFP_DIR}/{name}.txt") for name in GUT_SETS}
    gut_set_gene_lists["F_Gut-specific"] = f_gut_specific

    rows = []
    for set_name in GUT_SETS + ["F_Gut-specific"]:
        var_names = gut_set_gene_lists[set_name]
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
