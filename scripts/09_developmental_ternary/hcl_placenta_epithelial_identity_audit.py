#!/usr/bin/env python3
"""
Step 9 (developmental ternary map) -- Track C round-2 review requirement:
transcriptomic identity audit for the 1,095 HCL `Placenta1` barcodes
labeled "Epithelial cell" (the population Track C's pseudobulk uses as
its trophoblast proxy, since HCL's own annotation has no "Trophoblast"
label at all). The reviewer required empirical validation that this
proxy really is trophoblast-like, not just biological-reasoning
("placental epithelium = trophoblast") without a check.

Method: per-cell CP10K-normalized log1p expression (standard scRNA-seq
per-cell library-size normalization, independent of the pseudobulk
scores used elsewhere in this analysis), comparing the 1,095
"Epithelial cell" barcodes against the other 8,500 Placenta1 barcodes
(dominated by Fibroblast 72.7%, Macrophage 13.4% -- see
docs/STEP9_DEVELOPMENTAL_TERNARY_MAP.md). Two marker panels:
- trophoblast identity markers (classic CTB/STB/EVT markers spanning
  general trophoblast, syncytiotrophoblast, and extravillous
  trophoblast biology): KRT7, GATA3, TFAP2C, GCM1, ERVW-1, ERVFRD-1,
  CGA, CGB1/2/3/5/7/8, PSG1/3/4/5/6/7/8/9/11, HLA-G
- contamination markers (should be LOW if this is a clean epithelial
  compartment, not a coarsely-mislabeled mix): fibroblast (COL1A1,
  COL1A2, DCN, LUM, PDGFRB), immune/macrophage (PTPRC, CD68, CD14,
  CD163, LYZ)

Usage: python3 hcl_placenta_epithelial_identity_audit.py
(runs locally -- reuses the already-downloaded, byte-verified GSE134355
Placenta1 matrix + HCL annotation, lightweight, not heavy compute)
"""
import gzip
import numpy as np
import pandas as pd

RAW_DIR = "/Volumes/Stelligen_SSD/Stelligen/DATA/scRNAseq/GSE134355/raw/extracted"
ANNOT_PATH = "/Volumes/Stelligen_SSD/Stelligen/DATA/scRNAseq/GSE134355/raw/annotation/HCL_Fig1_cell_Info.xlsx"
OUT_DIR = "/Volumes/Stelligen_SSD/Stelligen/PR/TWEAKR-OncoPlacental/results/09_developmental_ternary"

TROPHOBLAST_MARKERS = [
    "KRT7", "GATA3", "TFAP2C", "GCM1", "ERVW-1", "ERVFRD-1",
    "CGA", "CGB1", "CGB2", "CGB3", "CGB5", "CGB7", "CGB8",
    "PSG1", "PSG3", "PSG4", "PSG5", "PSG6", "PSG7", "PSG8", "PSG9", "PSG11",
    "HLA-G",
]
CONTAMINATION_MARKERS = [
    "COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB",  # fibroblast
    "PTPRC", "CD68", "CD14", "CD163", "LYZ",       # immune/macrophage
]


def main():
    annot = pd.read_excel(ANNOT_PATH)
    annot["celltype"] = annot["celltype"].astype(str).str.strip()
    sub = annot[annot["sample"] == "Placenta"]

    epi_barcodes = set(x.split(".", 1)[1] for x in sub[sub["celltype"] == "Epithelial cell"]["cellnames"])
    rest_barcodes = set(x.split(".", 1)[1] for x in sub[sub["celltype"] != "Epithelial cell"]["cellnames"])

    with gzip.open(f"{RAW_DIR}/GSM4008722_Placenta1_dge.txt.gz", "rt") as f:
        mat = pd.read_csv(f, sep="\t", index_col=0)

    epi = [b for b in epi_barcodes if b in mat.columns]
    rest = [b for b in rest_barcodes if b in mat.columns]
    print(f"Epithelial cell (trophoblast-proxy): {len(epi)} barcodes")
    print(f"Rest of Placenta1 (Fibroblast/Macrophage/other): {len(rest)} barcodes")

    lib_size = mat.sum(axis=0)
    cpm10k = mat.div(lib_size, axis=1) * 1e4
    log_cpm = np.log1p(cpm10k)

    rows = []
    for gene in TROPHOBLAST_MARKERS + CONTAMINATION_MARKERS:
        category = "trophoblast" if gene in TROPHOBLAST_MARKERS else "contamination"
        if gene not in log_cpm.index:
            rows.append({"gene": gene, "category": category, "present": False})
            continue
        epi_vals, rest_vals = log_cpm.loc[gene, epi], log_cpm.loc[gene, rest]
        rows.append({
            "gene": gene, "category": category, "present": True,
            "epi_mean_logCP10K": round(epi_vals.mean(), 4),
            "rest_mean_logCP10K": round(rest_vals.mean(), 4),
            "epi_detect_frac": round((mat.loc[gene, epi] > 0).mean(), 4),
            "rest_detect_frac": round((mat.loc[gene, rest] > 0).mean(), 4),
            "log2fc_epi_vs_rest": round((epi_vals.mean() - rest_vals.mean()) / np.log(2), 4),
        })

    df = pd.DataFrame(rows)
    out_path = f"{OUT_DIR}/hcl_placenta_epithelial_identity_audit.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\nWrote {out_path}")

    present = df[df["present"]]
    troph = present[present["category"] == "trophoblast"]
    contam = present[present["category"] == "contamination"]
    n_troph_up = (troph["log2fc_epi_vs_rest"] > 0).sum()
    n_contam_down = (contam["log2fc_epi_vs_rest"] < 0).sum()
    print(f"\nTrophoblast markers positively enriched in Epithelial cell group: "
          f"{n_troph_up}/{len(troph)}")
    print(f"Contamination markers depleted in Epithelial cell group: "
          f"{n_contam_down}/{len(contam)}")
    print(f"\nStrongest trophoblast signals (log2FC):")
    print(troph.sort_values("log2fc_epi_vs_rest", ascending=False)[
        ["gene", "epi_detect_frac", "rest_detect_frac", "log2fc_epi_vs_rest"]
    ].to_string(index=False))


if __name__ == "__main__":
    main()
