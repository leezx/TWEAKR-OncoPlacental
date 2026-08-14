#!/usr/bin/env python3
"""
Tier-2 validation compute, per docs/STEP5_TIER2_VALIDATION.md (approved
after 1 review round: donor-aware, not organ-pooled).

Primary unit: (organ, donor, cell_ontology_class) pseudobulk from raw
counts (layers['raw_counts']), CPM per donor-celltype library's own size.
Minimum cell-count floor (>=20 cells) applied at the donor x cell-type
level. Aggregated to cell-type level for reporting, keeping donor
structure visible (eligible donor count, donor-level detection fraction,
median/range CPM); cell types backed by only 1 donor labeled explicitly.

- Organ-matched check: each of the 4 HDMA-organ-matching Tabula Sapiens
  files (Liver, Skin, Spleen, Thymus) checked against that organ's own
  frozen F_developmental_<Organ>.txt gene list.
- Whole-body-style check: D-shared (6 genes) and P-specific (78-gene
  primary tier) checked across all 5 organs' donor-celltype units
  combined (gene x cell type x donor evidence, not organ-pooled).

Per reviewer's explicit interpretive guidance: the priority result is not
the overall detection rate, but which frozen genes show CROSS-DONOR-
CONSISTENT high expression at the adult cell-type level (detected in
ALL eligible donors, with >=2 eligible donors) -- those are the real
candidate false positives worth manual review. Single-donor or low-CPM
sporadic hits are reported but not over-interpreted.

Usage: python3 tier2_validation.py
(run on Argos, argos-codex env)
"""
import os
import zipfile
import tempfile
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

ROOT = "/home/zz950/TWEAKR-OncoPlacental"
TS_DIR = "/home/zz950/DATA/1.Databases/TabulaSapiens/raw"
DFP_DIR = f"{ROOT}/results/04_dfp_signature/dfp_gene_sets"
OUT_DIR = f"{ROOT}/results/05_tier2_validation"
os.makedirs(OUT_DIR, exist_ok=True)

ORGANS = ["Liver", "Skin", "Spleen", "Thymus", "Large_Intestine"]
F_MATCHED_ORGANS = ["Liver", "Skin", "Spleen", "Thymus"]  # organs with an HDMA F-developmental counterpart
MIN_CELLS = 20
NOT_DETECTED_FLOOR = 1.0  # CPM


def load_organ(organ):
    path = f"{TS_DIR}/TS_{organ}.h5ad.zip"
    with zipfile.ZipFile(path) as z:
        h5ad_name = [n for n in z.namelist() if n.endswith(".h5ad")][0]
        with tempfile.TemporaryDirectory() as tmp:
            z.extract(h5ad_name, tmp)
            adata = ad.read_h5ad(os.path.join(tmp, h5ad_name))
    return adata


def build_donor_celltype_pseudobulk(adata):
    """Returns (pb_df: genes x (donor||celltype) columns raw counts,
    meta_df: sample, donor, cell_type, n_cells)."""
    counts = adata.layers["raw_counts"]
    donor = adata.obs["donor"].astype(str).values
    celltype = adata.obs["cell_ontology_class"].astype(str).values
    group = pd.Series(donor).str.cat(pd.Series(celltype), sep="||").values
    uniq_groups = sorted(set(group))
    group_idx = {g: i for i, g in enumerate(uniq_groups)}
    row_idx = np.arange(len(group))
    col_idx = np.array([group_idx[g] for g in group])
    indicator = sparse.csr_matrix(
        (np.ones(len(group)), (row_idx, col_idx)),
        shape=(len(group), len(uniq_groups)),
    )
    pb = indicator.T @ counts  # groups x genes
    pb = np.asarray(pb.todense()) if sparse.issparse(pb) else np.asarray(pb)
    pb_df = pd.DataFrame(pb.T, index=adata.var_names, columns=uniq_groups)
    # Collapse duplicate gene symbols by summing (same discipline as Step 2's
    # canonical_symbol collision handling) -- checked, not assumed absent.
    n_before = pb_df.shape[0]
    pb_df = pb_df.groupby(pb_df.index).sum()
    n_dupes = n_before - pb_df.shape[0]
    if n_dupes > 0:
        print(f"  collapsed {n_dupes} duplicate gene-symbol rows by summing")
    n_cells = pd.Series(group).value_counts().reindex(uniq_groups).values
    meta = pd.DataFrame({
        "sample": uniq_groups,
        "donor": [g.split("||")[0] for g in uniq_groups],
        "cell_type": [g.split("||")[1] for g in uniq_groups],
        "n_cells": n_cells,
    })
    return pb_df, meta


def gene_celltype_report(pb_cpm, meta, genes, min_cells=MIN_CELLS, floor=NOT_DETECTED_FLOOR):
    """For each (gene, cell_type), aggregate across eligible donors
    (n_cells>=min_cells) into: n_eligible_donors, n_donors_detected,
    donor_detection_fraction, median_cpm, max_cpm, cross_donor_consistent
    (all eligible donors detected AND n_eligible_donors>=2)."""
    eligible = meta[meta.n_cells >= min_cells]
    rows = []
    for ct, sub in eligible.groupby("cell_type"):
        samples = sub["sample"].values
        n_donors = len(samples)
        if n_donors == 0:
            continue
        sub_cpm = pb_cpm.loc[[g for g in genes if g in pb_cpm.index], samples]
        for gene in sub_cpm.index:
            vals = sub_cpm.loc[gene]
            detected = vals >= floor
            n_detected = int(detected.sum())
            frac = n_detected / n_donors
            cross_donor_consistent = bool(n_donors >= 2 and n_detected == n_donors)
            rows.append({
                "cell_type": ct, "gene": gene,
                "n_eligible_donors": n_donors, "n_donors_detected": n_detected,
                "donor_detection_fraction": round(frac, 3),
                "median_cpm": round(float(vals.median()), 2),
                "max_cpm": round(float(vals.max()), 2),
                "single_donor_evidence": n_donors == 1,
                "cross_donor_consistent_hit": cross_donor_consistent and n_detected > 0,
            })
    return pd.DataFrame(rows)


def load_gene_list(name):
    with open(f"{DFP_DIR}/{name}.txt") as f:
        return [l.strip() for l in f if l.strip()]


def run_organ_matched_check():
    """Organ-matched check (F-specific)."""
    organ_matched_reports = []
    for organ in F_MATCHED_ORGANS:
        print(f"=== Loading {organ} ===")
        adata = load_organ(organ)
        pb_counts, meta = build_donor_celltype_pseudobulk(adata)
        lib_size = pb_counts.sum(axis=0)
        pb_cpm = pb_counts.div(lib_size, axis=1) * 1e6

        f_genes = load_gene_list(f"F_developmental_{organ}")
        print(f"{organ}: {len(f_genes)} F-specific candidate genes, {pb_counts.shape[1]} donor-celltype samples "
              f"({meta.n_cells.ge(MIN_CELLS).sum()} pass >= {MIN_CELLS} cells)")

        report = gene_celltype_report(pb_cpm, meta, f_genes)
        report.insert(0, "organ", organ)
        organ_matched_reports.append(report)
        del adata

    organ_matched_df = pd.concat(organ_matched_reports, ignore_index=True)
    organ_matched_df.to_csv(f"{OUT_DIR}/organ_matched_F_specific_validation.tsv", sep="\t", index=False)
    n_flagged = int(organ_matched_df.cross_donor_consistent_hit.sum())
    print(f"\nOrgan-matched F-specific validation: {len(organ_matched_df)} (gene, cell_type) pairs tested, "
          f"{n_flagged} cross-donor-consistent hits flagged for review")
    flagged = organ_matched_df[organ_matched_df.cross_donor_consistent_hit]
    flagged.to_csv(f"{OUT_DIR}/organ_matched_F_specific_FLAGGED.tsv", sep="\t", index=False)
    print(flagged.to_string(index=False) if len(flagged) else "(none flagged)")
    return organ_matched_df


def run_whole_body_check():
    """Whole-body-style check (D-shared, P-specific)."""
    d_shared = load_gene_list("D_shared_FINAL")
    p_specific = load_gene_list("P_specific_FINAL")
    whole_body_genes = sorted(set(d_shared) | set(p_specific))
    print(f"\nWhole-body-style check: {len(whole_body_genes)} genes (D-shared {len(d_shared)} + P-specific {len(p_specific)})")

    wb_reports = []
    for organ in ORGANS:
        print(f"=== Loading {organ} (whole-body pass) ===")
        adata = load_organ(organ)
        pb_counts, meta = build_donor_celltype_pseudobulk(adata)
        lib_size = pb_counts.sum(axis=0)
        pb_cpm = pb_counts.div(lib_size, axis=1) * 1e6
        report = gene_celltype_report(pb_cpm, meta, whole_body_genes)
        report.insert(0, "organ", organ)
        wb_reports.append(report)
        del adata

    wb_df = pd.concat(wb_reports, ignore_index=True)
    wb_df["gene_set"] = wb_df["gene"].apply(lambda g: "D_shared" if g in d_shared else "P_specific")
    wb_df.to_csv(f"{OUT_DIR}/whole_body_DP_validation.tsv", sep="\t", index=False)
    n_wb_flagged = int(wb_df.cross_donor_consistent_hit.sum())
    print(f"\nWhole-body-style D-shared/P-specific validation: {len(wb_df)} (organ, cell_type, gene) triples tested, "
          f"{n_wb_flagged} cross-donor-consistent hits flagged for review")
    wb_flagged = wb_df[wb_df.cross_donor_consistent_hit]
    wb_flagged.to_csv(f"{OUT_DIR}/whole_body_DP_FLAGGED.tsv", sep="\t", index=False)
    print(wb_flagged.to_string(index=False) if len(wb_flagged) else "(none flagged)")
    return wb_df


if __name__ == "__main__":
    run_organ_matched_check()
    run_whole_body_check()
    print("\n=== DONE ===")
