#!/usr/bin/env python3
"""
Final D-shared / F-specific / P-specific gene set assembly, per
docs/STEP4_DFP_DESIGN.md. Both halves are now frozen at the per-dataset/
per-organ evidence level (PR #10: P-developmental's replicated_in_placenta;
PR #12: F-developmental's elevated+organ-matched adult_excluded) -- this
script does the two remaining real-data steps the design doc explicitly
deferred, then the final set assembly:

1. P-developmental still needs `adult_excluded(gene, whole_body)` per the
   design doc's compound definition (replicated_in_placenta AND
   adult_excluded(whole_body)) -- never frozen. Computed here with the
   same not-detected-floor + percentile-among-detected + coverage-aware
   provenance design already used for organ-matched (PR #11/#12), across
   GTEx's 68 tissues (primary, per STEP4_STATISTICAL_DESIGN.md) with HPA's
   39 non-placenta tissues reported alongside (fills GTEx's Thymus gap).
2. F-developmental is organ-specific (7 organs), P-developmental is
   global -- STEP4_DFP_DESIGN.md's reviewer note explicitly deferred the
   union-vs-consensus decision to "once real per-organ results exist."
   Computed here: per-organ F-developmental gene sets (frozen
   elev75/adult_excl25/quorum1.0), their pairwise overlap, and both a
   union and an "in >=2 organs" consensus candidate, reported for review.

Usage: python3 build_dfp_gene_sets.py
(run on Argos, argos-codex env)
"""
import os
import pandas as pd

ROOT = "/home/zz950/TWEAKR-OncoPlacental"
EDGE_DIR = f"{ROOT}/results/04_dfp_signature/edgeR"
OUT_DIR = f"{ROOT}/results/04_dfp_signature/dfp_gene_sets"
os.makedirs(OUT_DIR, exist_ok=True)

exec(open(f"{ROOT}/scripts/04_dfp_signature/f_developmental_calibration.py").read().split("# ---- HDMA per-organ elevated calibration ----")[0])

# ============ 1. P-developmental: replicated_in_placenta (frozen, PR #10) ============
a_df = pd.read_csv(f"{EDGE_DIR}/Arutyunyan_edgeR_results.tsv", sep="\t")
n_df = pd.read_csv(f"{EDGE_DIR}/Nature2026_edgeR_results.tsv", sep="\t")
a_pass = set(a_df.loc[(a_df.logFC >= 0.75) & (a_df.FDR < 0.05), "gene"])
n_pass = set(n_df.loc[(n_df.logFC >= 0.75) & (n_df.FDR < 0.05), "gene"])
replicated_in_placenta = a_pass & n_pass
print(f"replicated_in_placenta (frozen, PR #10): {len(replicated_in_placenta)} genes "
      f"(Arutyunyan {len(a_pass)}, Nature2026 {len(n_pass)})")

# ============ 2. P-developmental: adult_excluded(whole_body) -- NOT yet frozen ============
WHOLE_BODY_PCT_CANDIDATES = [10, 25]
WHOLE_BODY_QUORUM = [("all", 0), ("all_but_1", 1), ("all_but_2", 2)]

hpa_non_placenta = [c for c in hpa_pct.columns if c != "placenta"]

wb_rows = []
for pct_cut in WHOLE_BODY_PCT_CANDIDATES:
    fails_gtex = pd.DataFrame(False, index=gtex_not_detected.index, columns=tissue_cols)
    for c in tissue_cols:
        fails_gtex[c] = (~gtex_not_detected[c]) & (gtex_pct[c] > pct_cut)
    fails_hpa = pd.DataFrame(False, index=hpa_not_detected.index, columns=hpa_non_placenta)
    for c in hpa_non_placenta:
        fails_hpa[c] = (~hpa_not_detected[c]) & (hpa_pct[c] > pct_cut)
    for label, allowed in WHOLE_BODY_QUORUM:
        pass_gtex = set(fails_gtex.index[fails_gtex.sum(axis=1) <= allowed])
        pass_hpa = set(fails_hpa.index[fails_hpa.sum(axis=1) <= allowed])
        p_dev_gtex = replicated_in_placenta & pass_gtex
        p_dev_hpa = replicated_in_placenta & pass_hpa
        p_dev_both = replicated_in_placenta & pass_gtex & pass_hpa
        wb_rows.append({
            "pct_cut": pct_cut, "quorum": label,
            "n_pass_gtex_wholebody": len(pass_gtex), "n_P_dev_via_GTEx": len(p_dev_gtex),
            "n_pass_hpa_wholebody": len(pass_hpa), "n_P_dev_via_HPA": len(p_dev_hpa),
            "n_P_dev_GTEx_AND_HPA_agree": len(p_dev_both),
        })
wb_result = pd.DataFrame(wb_rows)
wb_result.to_csv(f"{OUT_DIR}/whole_body_adult_excluded_calibration.tsv", sep="\t", index=False)
print("\n=== Whole-body adult_excluded calibration (P-developmental) ===")
print(wb_result.to_string(index=False))

# Marker check: do known trophoblast markers survive whole-body exclusion?
print("\n=== Marker check: trophoblast markers in whole-body adult_excluded (GTEx, all_but_1, pct=25) ===")
markers = ["ERVFRD-1", "CGA", "CSH1", "CSH2", "PSG1", "PSG3", "GATA3", "KRT7", "HLA-G"]
fails_gtex_25 = pd.DataFrame(False, index=gtex_not_detected.index, columns=tissue_cols)
for c in tissue_cols:
    fails_gtex_25[c] = (~gtex_not_detected[c]) & (gtex_pct[c] > 25)
n_fail_25 = fails_gtex_25.sum(axis=1)
for m in markers:
    if m in n_fail_25.index:
        print(f"{m:10s} in replicated_in_placenta: {m in replicated_in_placenta}, "
              f"GTEx tissues failing whole-body(pct=25): {n_fail_25[m]}/68")

# ============ 3. F-developmental: real per-organ gene sets (frozen, PR #12) ============
ORGANS_LOCAL = ["Adrenal", "Thyroid", "Spleen", "Thymus", "Liver", "Skin", "Stomach"]
f_dev_sets = {}
for organ in ORGANS_LOCAL:
    counts = pd.read_csv(f"{HDMA_DIR}/{organ}_pseudobulk_counts.tsv", sep="\t", index_col=0)
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    nd, pct = detected_percentile(cpm)
    median_pct = pct.median(axis=1, skipna=True)
    detected_frac = (~nd).mean(axis=1)
    elevated = (median_pct >= 75) & (detected_frac >= 1.0)
    mask, prov = adult_excluded_mask(organ, 25)
    elevated_genes = set(elevated.index[elevated])
    excluded_genes = set(mask.index[mask]) if mask is not None else set()
    f_dev_sets[organ] = elevated_genes & excluded_genes
    print(f"\nF-developmental[{organ}]: {len(f_dev_sets[organ])} genes")

union_f_dev = set().union(*f_dev_sets.values())
from collections import Counter
gene_organ_count = Counter()
for organ, genes in f_dev_sets.items():
    for g in genes:
        gene_organ_count[g] += 1
consensus_2plus = {g for g, c in gene_organ_count.items() if c >= 2}
print(f"\nF-developmental UNION (any organ): {len(union_f_dev)} genes")
print(f"F-developmental CONSENSUS (>=2 organs): {len(consensus_2plus)} genes")
print("Organ-count distribution of union genes:", Counter(gene_organ_count.values()))

# ============ 4. Final D/F/P assembly (using GTEx-primary whole-body, all_but_1, pct=25 as the P-developmental adult_excluded choice -- reported, not silently frozen) ============
fails_gtex_final = fails_gtex_25  # pct=25, all_but_1 quorum
p_dev_final = replicated_in_placenta & set(fails_gtex_final.index[fails_gtex_final.sum(axis=1) <= 1])
print(f"\n=== Final P-developmental (replicated_in_placenta AND adult_excluded(whole_body, GTEx, pct=25, all_but_1)) ===")
print(f"P-developmental: {len(p_dev_final)} genes")

d_shared_union = union_f_dev & p_dev_final
f_specific_union = union_f_dev - p_dev_final
p_specific_union = p_dev_final - union_f_dev

d_shared_consensus = consensus_2plus & p_dev_final
f_specific_consensus = consensus_2plus - p_dev_final
p_specific_consensus = p_dev_final - consensus_2plus

print(f"\n=== D/F/P using F-developmental UNION ===")
print(f"D-shared: {len(d_shared_union)}, F-specific: {len(f_specific_union)}, P-specific: {len(p_specific_union)}")
print(f"\n=== D/F/P using F-developmental CONSENSUS (>=2 organs) ===")
print(f"D-shared: {len(d_shared_consensus)}, F-specific: {len(f_specific_consensus)}, P-specific: {len(p_specific_consensus)}")

for m in markers:
    in_p = m in p_dev_final
    in_f_union = m in union_f_dev
    print(f"Marker check -- {m:10s}: P-developmental={in_p}, F-developmental(any organ)={in_f_union}")

# Write out gene lists
for name, geneset in [
    ("P_developmental", p_dev_final),
    ("F_developmental_union", union_f_dev),
    ("F_developmental_consensus2plus", consensus_2plus),
    ("D_shared_union", d_shared_union),
    ("F_specific_union", f_specific_union),
    ("P_specific_union", p_specific_union),
    ("D_shared_consensus", d_shared_consensus),
    ("F_specific_consensus", f_specific_consensus),
    ("P_specific_consensus", p_specific_consensus),
]:
    with open(f"{OUT_DIR}/{name}.txt", "w") as f:
        f.write("\n".join(sorted(geneset)) + "\n")
    print(f"Wrote {OUT_DIR}/{name}.txt ({len(geneset)} genes)")

for organ, genes in f_dev_sets.items():
    with open(f"{OUT_DIR}/F_developmental_{organ}.txt", "w") as f:
        f.write("\n".join(sorted(genes)) + "\n")
