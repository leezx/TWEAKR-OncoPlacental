#!/usr/bin/env python3
"""
Phase I design pivot: M11 is an NMF meta-program merely *annotated* as
containing Oncofetal-like features (via a Jaccard match to revCSC) -- it is
a proxy, not a definition, and should not have been treated as the primary
Oncofetal anchor. The actual published, independently-defined
Oncofetal(-like) signature is revCSC itself (Colorectal cancer revival
cancer-stem-cell signature, from CSC_subtype_signatures.ensembl_mapping.tsv,
mouse-derived gene symbols mapped to human Ensembl IDs). This script:

1. Extracts revCSC's unique, successfully-mapped human gene set from that
   table (dropping duplicate rows and rows with no ensembl_id, i.e. genes
   that failed mouse->human mapping).
2. Computes revCSC's gene overlap against the frozen D-shared, F-specific
   (global + each of the 7 per-organ lineage modules, each intersected with
   F_specific_FINAL), and P-specific signatures -- same discipline as the
   earlier M11 x D/F/P audit (results/06_crc_projection/m11_overlap_audit/),
   needed before any null-calibrated revCSC score can be correlated against
   D/F/P without risking mechanical (gene-sharing) correlation.

Usage: python3 revcsc_dfp_overlap_audit.py <out_dir>
(run on Argos, argos-codex env -- reads the CSC signature table from
/home/zz950/DATA/scRNAseq/meta_study/CRC_single_cell_atlas_2025/NMF/
csc_signature_jaccard/CSC_subtype_signatures.ensembl_mapping.tsv, and the
frozen D/F/P gene lists from this repo's own results/04_dfp_signature/)
"""
import sys
import os
import csv

REPO = "/home/zz950/TWEAKR-OncoPlacental"
CSC_TABLE = (
    "/home/zz950/DATA/scRNAseq/meta_study/CRC_single_cell_atlas_2025/"
    "NMF/csc_signature_jaccard/CSC_subtype_signatures.ensembl_mapping.tsv"
)
DFP_DIR = f"{REPO}/results/04_dfp_signature/dfp_gene_sets"
ORGANS = ["Adrenal", "Liver", "Skin", "Spleen", "Stomach", "Thymus", "Thyroid"]

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "results/06_crc_projection/revcsc_overlap_audit"
os.makedirs(OUT_DIR, exist_ok=True)


def load_gene_set(path):
    with open(path) as f:
        return {l.strip() for l in f if l.strip()}


def load_revcsc():
    rows = []
    unmapped = []
    with open(CSC_TABLE, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row["cluster"] != "revCSC":
                continue
            rows.append(row)
    n_raw = len(rows)
    symbols_mapped = set()
    for row in rows:
        ens = row["ensembl_id"].strip()
        sym = row["gene_symbol_norm"].strip()
        if ens:
            symbols_mapped.add(sym)
        else:
            unmapped.append(sym)
    return n_raw, sorted(symbols_mapped), sorted(set(unmapped))


def main():
    n_raw, revcsc_symbols, unmapped = load_revcsc()
    n_raw_unique_symbols = len({row for row in revcsc_symbols} | set(unmapped))

    d_shared = load_gene_set(f"{DFP_DIR}/D_shared_FINAL.txt")
    f_specific = load_gene_set(f"{DFP_DIR}/F_specific_FINAL.txt")
    p_specific = load_gene_set(f"{DFP_DIR}/P_specific_FINAL.txt")
    f_lineage = {}
    for organ in ORGANS:
        organ_set = load_gene_set(f"{DFP_DIR}/F_developmental_{organ}.txt")
        f_lineage[organ] = organ_set & f_specific  # intersected with frozen F_specific_FINAL

    revcsc_set = set(revcsc_symbols)

    overlap_d = sorted(revcsc_set & d_shared)
    overlap_f = sorted(revcsc_set & f_specific)
    overlap_p = sorted(revcsc_set & p_specific)
    overlap_lineage = {organ: sorted(revcsc_set & genes) for organ, genes in f_lineage.items()}

    lines = []
    lines.append("# revCSC x D/F/P gene-overlap audit\n")
    lines.append(
        "Per Phase I design pivot: revCSC (not M11) is the primary Oncofetal "
        "anchor. Checked directly whether revCSC's own mapped human gene set "
        "shares genes with D-shared, F-specific (global and each of the 7 "
        "per-organ lineage modules, intersected with the frozen "
        "F_specific_FINAL.txt), or P-specific, before any null-calibrated "
        "revCSC score is correlated against D/F/P scores.\n"
    )
    lines.append(f"## revCSC signature extraction\n")
    lines.append(
        f"- Raw rows in `CSC_subtype_signatures.ensembl_mapping.tsv` for "
        f"cluster=revCSC: {n_raw} (includes duplicate rows -- the source "
        f"table repeats some genes).\n"
        f"- Distinct gene symbols (raw, before mapping filter): "
        f"{n_raw_unique_symbols}.\n"
        f"- Successfully mapped to a human Ensembl ID: {len(revcsc_symbols)}.\n"
        f"- Failed mouse->human mapping (dropped from scoring set): "
        f"{len(unmapped)} ({', '.join(unmapped) if unmapped else 'none'}).\n"
    )
    lines.append("## Overlap result\n")
    lines.append("| Target signature (n genes) | revCSC overlap (n) | Overlapping genes |")
    lines.append("|---|---|---|")
    lines.append(f"| D-shared ({len(d_shared)}) | {len(overlap_d)} | {', '.join(overlap_d) or '-'} |")
    lines.append(f"| F-specific global ({len(f_specific)}) | {len(overlap_f)} | {', '.join(overlap_f) or '-'} |")
    lines.append(f"| P-specific ({len(p_specific)}) | {len(overlap_p)} | {', '.join(overlap_p) or '-'} |")
    for organ in ORGANS:
        genes = overlap_lineage[organ]
        lines.append(
            f"| F-lineage {organ} ({len(f_lineage[organ])}) | {len(genes)} | "
            f"{', '.join(genes) or '-'} |"
        )
    lines.append("")
    lines.append(f"revCSC scoring-set size used going forward: {len(revcsc_set)} genes "
                  f"({', '.join(sorted(revcsc_set))}).")

    md_path = f"{OUT_DIR}/revcsc_dfp_gene_overlap_audit.md"
    with open(md_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    with open(f"{OUT_DIR}/revCSC_symbols.txt", "w") as f:
        for s in sorted(revcsc_set):
            f.write(s + "\n")

    print(f"revCSC: {n_raw} raw rows -> {len(revcsc_set)} unique mapped genes "
          f"({len(unmapped)} unmapped: {unmapped})")
    print(f"Overlap: D-shared={len(overlap_d)}, F-specific={len(overlap_f)}, "
          f"P-specific={len(overlap_p)}")
    for organ in ORGANS:
        print(f"  F-lineage {organ}: {len(overlap_lineage[organ])}")
    print(f"Wrote {md_path}")
    print(f"Wrote {OUT_DIR}/revCSC_symbols.txt")


if __name__ == "__main__":
    main()
