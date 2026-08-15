#!/usr/bin/env python3
"""
Freeze mike_verzi_fetal_signature.gmt's 5 mouse gene sets to human-mapped,
orthology-verified signatures, per the Step 6a design
(docs/STEP6A_NORMAL_CONTEXT_FETAL_DECOMPOSITION.md) and the same discipline
established for revCSC in PR #17: primary = Compara-confirmed one-to-one
orthologs only; ambiguous (one-to-many) calls reported separately, never
silently used to pick a primary human target. Unlike revCSC's 32 genes
(one-by-one Ensembl REST queries), this covers 1,923 unique mouse gene
symbols across 5 signatures -- queried via a single bulk Ensembl BioMart
call (jun2026.archive.ensembl.org/biomart/martservice), not one-by-one.

Inputs:
  results/06a_normal_context/mike_verzi_sets_raw.json
    -- per-signature mouse gene membership (5 sets, from the .gmt file)
  results/06a_normal_context/mike_verzi_biomart_raw.tsv
    -- raw BioMart bulk-query output (mouse gene, human ortholog, type, confidence)
  results/04_dfp_signature/dfp_gene_sets/{D_shared,F_specific,P_specific}_FINAL.txt
  results/04_dfp_signature/dfp_gene_sets/F_developmental_<Organ>.txt (7 organs)

Outputs:
  results/06a_normal_context/mike_verzi_human_FINAL.tsv
    -- one row per queried mouse gene: class, human target(s), inclusion decision
  results/06a_normal_context/mike_verzi_signature_gene_counts.tsv
    -- per signature: n mouse genes, n resolved to primary human symbol
  results/06a_normal_context/mike_verzi_dfp_overlap.tsv
    -- per signature x per D/F/P target: overlap gene count and gene list
  results/06a_normal_context/mike_verzi_dfp_overlap_summary.md
    -- human-readable summary table
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
    biomart_rows = list(csv.DictReader(open(f"{OUT_DIR}/mike_verzi_biomart_raw.tsv"), delimiter="\t"))

    by_mouse_gene = defaultdict(list)
    for r in biomart_rows:
        by_mouse_gene[r["Gene name"]].append(r)

    all_mouse_genes = sorted({g for genes in sets.values() for g in genes})

    # classify each mouse gene
    classification = {}  # mouse_symbol -> (class, [ (human_symbol, human_ensembl) ])
    for g in all_mouse_genes:
        entries = by_mouse_gene.get(g, [])
        with_ortholog = [e for e in entries if e["Human gene name"]]
        if not entries:
            classification[g] = ("NOT_FOUND_IN_BIOMART", [])
        elif not with_ortholog:
            classification[g] = ("NO_ORTHOLOG", [])
        else:
            types = {e["Human homology type"] for e in with_ortholog}
            if len(with_ortholog) == 1 and types == {"ortholog_one2one"}:
                classification[g] = ("ortholog_one2one",
                                      [(with_ortholog[0]["Human gene name"], with_ortholog[0]["Human gene stable ID"])])
            else:
                classification[g] = ("ambiguous",
                                      [(e["Human gene name"], e["Human gene stable ID"]) for e in with_ortholog])

    # frozen artifact: one row per mouse gene, with which signature(s) it's in
    gene_to_signatures = defaultdict(list)
    for sig, genes in sets.items():
        for g in genes:
            gene_to_signatures[g].append(sig)

    with open(f"{OUT_DIR}/mike_verzi_human_FINAL.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["mouse_symbol", "signatures", "class", "n_human_targets",
                     "human_targets", "inclusion_decision"])
        for g in all_mouse_genes:
            cls, targets = classification[g]
            sigs = ";".join(sorted(set(gene_to_signatures[g])))
            target_str = ";".join(f"{s}:{e}" for s, e in targets)
            decision = "include_primary" if cls == "ortholog_one2one" else (
                "extended_ambiguous_only" if cls == "ambiguous" else "exclude")
            w.writerow([g, sigs, cls, len(targets), target_str, decision])

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
    md_lines = ["# mike_verzi (normal-context fetal/revival) x D/F/P gene-overlap audit\n"]
    md_lines.append(
        "Per Step 6a design: before any scoring, check whether each of the 5 "
        "independently-published normal-tissue fetal/revival/regeneration gene "
        "sets (mouse->human orthology-mapped, primary = Compara-confirmed "
        "one-to-one only) shares genes with D-shared, F-specific (global + 7 "
        "lineage modules), or P-specific.\n"
    )
    md_lines.append("## Signature -> human ortholog resolution\n")
    md_lines.append("| Signature | n mouse genes | n primary (one2one) human | % resolved |")
    md_lines.append("|---|---|---|---|")
    for sig, mouse_genes in sets.items():
        n_mouse = len(mouse_genes)
        n_primary = len(primary_human_sets[sig])
        md_lines.append(f"| `{sig}` | {n_mouse} | {n_primary} | {100*n_primary/n_mouse:.1f}% |")

    md_lines.append("\n## D/F/P overlap (primary human gene sets)\n")
    md_lines.append("| Signature | D-shared (6) | F-specific global (2,504) | P-specific (78) |")
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
    print("Ortholog classification (all 1,923 unique mouse genes):")
    print(Counter(v[0] for v in classification.values()))
    print("\nPer-signature resolution:")
    for sig, mouse_genes in sets.items():
        print(f"  {sig}: {len(mouse_genes)} mouse -> {len(primary_human_sets[sig])} primary human "
              f"({100*len(primary_human_sets[sig])/len(mouse_genes):.1f}%)")
    print(f"\nWrote {OUT_DIR}/mike_verzi_human_FINAL.tsv")
    print(f"Wrote {OUT_DIR}/mike_verzi_dfp_overlap.tsv")
    print(f"Wrote {OUT_DIR}/mike_verzi_dfp_overlap_summary.md")


if __name__ == "__main__":
    main()
