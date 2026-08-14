#!/usr/bin/env python3
"""
Freeze mike_verzi_fetal_signature.gmt's 5 mouse gene sets to human-mapped,
orthology-verified signatures. v2 -- fixes PR #19 round-1 REQUEST_CHANGES:

Blocker 1 (this file): the v1 join matched the GMT's raw mouse symbol
against BioMart's canonical `Gene name` via exact string equality. BioMart's
bulk pull in v1 was itself filtered by that same raw symbol list, so any
case mismatch (GMT's "Col1A1" vs canonical "Col1a1"; "Tnfrsf12A" vs
"Tnfrsf12a"; "Tmsb4X" vs "Tmsb4x"; "Ly6A" vs "Ly6a"; bulk S100A*/Slc*
families) caused the row to be silently dropped before any join happened --
these were misclassified NOT_FOUND_IN_BIOMART when they're really just
capitalization-normalization failures, not missing orthologs. Fixed by
querying the FULL unfiltered mouse gene table (query_mouse_biomart_full.py,
83,845 rows, no symbol-list filter possible) and joining locally with an
explicit two-stage exact-then-case-insensitive lookup, per the reviewer's
prescribed mapping:

    GMT raw symbol -> canonical mouse symbol -> Ensembl Compara human ortholog

  1. exact match against canonical `Gene name` (preferred)
  2. exact match fails -> case-insensitive unique match against canonical
     `Gene name` (recovers Col1A1 -> Col1a1 style GMT-formatting drift)
  3. case-insensitive match hits >1 distinct canonical symbol -> AMBIGUOUS_SYMBOL,
     excluded from primary (matches Ly6a-style outcome-independence discipline
     from PR #17 round 3 -- ambiguity is resolved by data provenance, never by
     which candidate "looks more correct" downstream)
  4. no match at all (exact or case-insensitive) -> NOT_FOUND_IN_BIOMART (now a
     real absence, not a formatting artifact)

Every gene's raw_mouse_symbol -> canonical_mouse_symbol -> mouse_ensembl_id
provenance is preserved in the output (mike_verzi_symbol_resolution.tsv), so
the recovery step is auditable, not a silent rewrite.

Inputs:
  results/06a_normal_context/mike_verzi_sets_raw.json
  results/06a_normal_context/mouse_biomart_full.tsv (v2, full unfiltered pull)
  results/04_dfp_signature/dfp_gene_sets/{D_shared,F_specific,P_specific}_FINAL.txt
  results/04_dfp_signature/dfp_gene_sets/F_developmental_<Organ>.txt (7 organs)

Outputs:
  results/06a_normal_context/mike_verzi_symbol_resolution.tsv -- new: raw->canonical->ensembl provenance
  results/06a_normal_context/mike_verzi_human_FINAL.tsv
  results/06a_normal_context/mike_verzi_signature_gene_counts.tsv
  results/06a_normal_context/mike_verzi_dfp_overlap.tsv
  results/06a_normal_context/mike_verzi_dfp_overlap_summary.md
"""
import csv
import json
import os
from collections import defaultdict, Counter

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT_DIR = os.path.join(REPO, "results/06a_normal_context")
DFP_DIR = os.path.join(REPO, "results/04_dfp_signature/dfp_gene_sets")
ORGANS = ["Adrenal", "Liver", "Skin", "Spleen", "Stomach", "Thymus", "Thyroid"]


def load_gene_set(path):
    with open(path) as f:
        return {l.strip() for l in f if l.strip()}


def main():
    sets = json.load(open(f"{OUT_DIR}/mike_verzi_sets_raw.json"))
    all_mouse_genes = sorted({g for genes in sets.values() for g in genes})

    # full, unfiltered mouse gene table -- rows keyed by canonical Gene name
    rows = list(csv.DictReader(open(f"{OUT_DIR}/mouse_biomart_full.tsv"), delimiter="\t"))
    by_canonical_name = defaultdict(list)   # exact canonical symbol -> rows (one gene can have multiple ortholog rows)
    by_upper_name = defaultdict(set)        # uppercased symbol -> set of distinct canonical symbols sharing that uppercase form
    for r in rows:
        name = r["Gene name"]
        if not name:
            continue
        by_canonical_name[name].append(r)
        by_upper_name[name.upper()].add(name)

    # --- Stage 1: resolve each raw GMT symbol to a canonical mouse symbol ---
    resolution = {}  # raw_symbol -> (canonical_symbol or None, method, mouse_ensembl_ids)
    for raw in all_mouse_genes:
        if raw in by_canonical_name:
            canonical_candidates = {raw}
            method = "exact"
        else:
            canonical_candidates = by_upper_name.get(raw.upper(), set())
            method = "case_insensitive" if len(canonical_candidates) == 1 else (
                "ambiguous_case_insensitive" if len(canonical_candidates) > 1 else "not_found")

        if method in ("exact", "case_insensitive"):
            canonical = next(iter(canonical_candidates))
            ensembl_ids = sorted({r["Gene stable ID"] for r in by_canonical_name[canonical]})
            resolution[raw] = (canonical, method, ensembl_ids)
        elif method == "ambiguous_case_insensitive":
            resolution[raw] = (None, method, sorted(canonical_candidates))
        else:
            resolution[raw] = (None, "not_found", [])

    # --- Stage 2: classify orthology using the resolved canonical symbol's rows ---
    classification = {}  # raw_symbol -> (class, [(human_symbol, human_ensembl)])
    for raw in all_mouse_genes:
        canonical, method, extra = resolution[raw]
        if method == "ambiguous_case_insensitive":
            classification[raw] = ("AMBIGUOUS_SYMBOL", [])
            continue
        if method == "not_found":
            classification[raw] = ("NOT_FOUND_IN_BIOMART", [])
            continue
        entries = by_canonical_name[canonical]
        with_ortholog = [e for e in entries if e["Human gene name"]]
        if not with_ortholog:
            classification[raw] = ("NO_ORTHOLOG", [])
        else:
            types = {e["Human homology type"] for e in with_ortholog}
            if len(with_ortholog) == 1 and types == {"ortholog_one2one"}:
                classification[raw] = ("ortholog_one2one",
                                        [(with_ortholog[0]["Human gene name"], with_ortholog[0]["Human gene stable ID"])])
            else:
                classification[raw] = ("ambiguous",
                                        [(e["Human gene name"], e["Human gene stable ID"]) for e in with_ortholog])

    # symbol-resolution provenance artifact (new -- makes the case-recovery step auditable)
    with open(f"{OUT_DIR}/mike_verzi_symbol_resolution.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["raw_mouse_symbol", "canonical_mouse_symbol", "resolution_method", "mouse_ensembl_ids"])
        for raw in all_mouse_genes:
            canonical, method, extra = resolution[raw]
            ensembl_str = ";".join(extra) if method in ("exact", "case_insensitive") else (
                "AMBIGUOUS:" + ";".join(extra) if method == "ambiguous_case_insensitive" else "")
            w.writerow([raw, canonical or "", method, ensembl_str])

    n_recovered = sum(1 for raw in all_mouse_genes if resolution[raw][1] == "case_insensitive")
    n_ambiguous_symbol = sum(1 for raw in all_mouse_genes if resolution[raw][1] == "ambiguous_case_insensitive")
    n_not_found = sum(1 for raw in all_mouse_genes if resolution[raw][1] == "not_found")
    print(f"Symbol resolution: {len(all_mouse_genes)} raw GMT symbols -> "
          f"{len(all_mouse_genes) - n_recovered - n_ambiguous_symbol - n_not_found} exact, "
          f"{n_recovered} case-insensitive-recovered, "
          f"{n_ambiguous_symbol} ambiguous_case_insensitive (excluded), "
          f"{n_not_found} genuinely not_found", flush=True)

    # frozen artifact: one row per mouse gene, with which signature(s) it's in
    gene_to_signatures = defaultdict(list)
    for sig, genes in sets.items():
        for g in genes:
            gene_to_signatures[g].append(sig)

    with open(f"{OUT_DIR}/mike_verzi_human_FINAL.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["mouse_symbol", "canonical_mouse_symbol", "resolution_method", "signatures", "class",
                     "n_human_targets", "human_targets", "inclusion_decision"])
        for g in all_mouse_genes:
            cls, targets = classification[g]
            canonical, method, extra = resolution[g]
            sigs = ";".join(sorted(set(gene_to_signatures[g])))
            target_str = ";".join(f"{s}:{e}" for s, e in targets)
            decision = "include_primary" if cls == "ortholog_one2one" else (
                "extended_ambiguous_only" if cls == "ambiguous" else "exclude")
            w.writerow([g, canonical or "", method, sigs, cls, len(targets), target_str, decision])

    # per-signature primary (one2one) human symbol sets
    primary_human_sets = {}
    extended_human_sets = {}
    for sig, mouse_genes in sets.items():
        primary = set()
        extended = set()
        for g in mouse_genes:
            cls, targets = classification[g]
            if cls == "ortholog_one2one":
                primary.add(targets[0][0])
                extended.add(targets[0][0])
            elif cls == "ambiguous":
                for s, e in targets:
                    extended.add(s)
        primary_human_sets[sig] = primary
        extended_human_sets[sig] = extended
        with open(f"{OUT_DIR}/{sig}_human_primary.txt", "w") as f:
            f.write("\n".join(sorted(primary)) + "\n")
        with open(f"{OUT_DIR}/{sig}_human_extended.txt", "w") as f:
            f.write("\n".join(sorted(extended)) + "\n")

    # gene-count summary
    with open(f"{OUT_DIR}/mike_verzi_signature_gene_counts.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["signature", "n_mouse_genes", "n_primary_one2one_human",
                     "n_extended_human", "pct_resolved_primary"])
        for sig, mouse_genes in sets.items():
            n_mouse = len(mouse_genes)
            n_primary = len(primary_human_sets[sig])
            n_ext = len(extended_human_sets[sig])
            w.writerow([sig, n_mouse, n_primary, n_ext, f"{100*n_primary/n_mouse:.1f}%"])

    # D/F/P overlap audit (primary human sets only, per pre-registered rule)
    d_shared = load_gene_set(f"{DFP_DIR}/D_shared_FINAL.txt")
    f_specific = load_gene_set(f"{DFP_DIR}/F_specific_FINAL.txt")
    p_specific = load_gene_set(f"{DFP_DIR}/P_specific_FINAL.txt")
    f_lineage = {}
    for organ in ORGANS:
        organ_set = load_gene_set(f"{DFP_DIR}/F_developmental_{organ}.txt")
        f_lineage[organ] = organ_set & f_specific

    overlap_rows = []
    md_lines = ["# mike_verzi (normal-context fetal/revival) x D/F/P gene-overlap audit (v2, post case-fix)\n"]
    md_lines.append(
        "v2: fixes PR #19 round-1 blocker -- v1's mouse->human mapping silently "
        "dropped genes whose GMT symbol case didn't exactly match BioMart's "
        "canonical mouse Gene name (e.g. Col1A1/Tnfrsf12A/Tmsb4X/Ly6A and bulk "
        "S100A*/Slc* families), misclassifying real orthologs as "
        "NOT_FOUND_IN_BIOMART. Fixed via a full unfiltered mouse-gene BioMart "
        "pull + explicit exact-then-case-insensitive local join "
        "(mike_verzi_symbol_resolution.tsv has the full raw->canonical->ensembl "
        "provenance for every gene).\n"
    )
    md_lines.append("## Signature -> human ortholog resolution\n")
    md_lines.append("| Signature | n mouse genes | n primary (one2one) human | % resolved |")
    md_lines.append("|---|---|---|---|")
    for sig, mouse_genes in sets.items():
        n_mouse = len(mouse_genes)
        n_primary = len(primary_human_sets[sig])
        md_lines.append(f"| `{sig}` | {n_mouse} | {n_primary} | {100*n_primary/n_mouse:.1f}% |")

    md_lines.append(f"\n## D/F/P overlap (primary human gene sets)\n")
    md_lines.append(f"| Signature | D-shared ({len(d_shared)}) | F-specific global ({len(f_specific)}) | P-specific ({len(p_specific)}) |")
    md_lines.append("|---|---|---|---|")
    for sig in sets:
        hg = primary_human_sets[sig]
        ov_d = sorted(hg & d_shared)
        ov_f = sorted(hg & f_specific)
        ov_p = sorted(hg & p_specific)
        overlap_rows.append((sig, "D-shared", len(ov_d), ";".join(ov_d)))
        overlap_rows.append((sig, "F-specific-global", len(ov_f), ";".join(ov_f)))
        overlap_rows.append((sig, "P-specific", len(ov_p), ";".join(ov_p)))
        md_lines.append(f"| `{sig}` | {len(ov_d)} | {len(ov_f)} | {len(ov_p)} |")

    md_lines.append("\n## D/F/P overlap, F-lineage modules (organ ∩ F-specific)\n")
    header = "| Signature | " + " | ".join(f"{o} ({len(f_lineage[o])})" for o in ORGANS) + " |"
    md_lines.append(header)
    md_lines.append("|---" * (len(ORGANS) + 1) + "|")
    for sig in sets:
        hg = primary_human_sets[sig]
        cells = []
        for organ in ORGANS:
            ov = sorted(hg & f_lineage[organ])
            overlap_rows.append((sig, f"F-lineage-{organ}", len(ov), ";".join(ov)))
            cells.append(str(len(ov)))
        md_lines.append(f"| `{sig}` | " + " | ".join(cells) + " |")

    with open(f"{OUT_DIR}/mike_verzi_dfp_overlap.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["signature", "target", "n_overlap", "overlap_genes"])
        for row in overlap_rows:
            w.writerow(row)

    with open(f"{OUT_DIR}/mike_verzi_dfp_overlap_summary.md", "w") as f:
        f.write("\n".join(md_lines) + "\n")

    # console summary
    print("\nOrtholog classification (all unique mouse genes, v2):")
    print(Counter(v[0] for v in classification.values()))
    print("\nPer-signature resolution (v2):")
    for sig, mouse_genes in sets.items():
        print(f"  {sig}: {len(mouse_genes)} mouse -> {len(primary_human_sets[sig])} primary human "
              f"({100*len(primary_human_sets[sig])/len(mouse_genes):.1f}%)")
    print(f"\nWrote {OUT_DIR}/mike_verzi_symbol_resolution.tsv")
    print(f"Wrote {OUT_DIR}/mike_verzi_human_FINAL.tsv")
    print(f"Wrote {OUT_DIR}/mike_verzi_dfp_overlap.tsv")
    print(f"Wrote {OUT_DIR}/mike_verzi_dfp_overlap_summary.md")


if __name__ == "__main__":
    main()
