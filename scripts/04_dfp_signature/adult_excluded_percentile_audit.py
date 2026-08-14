#!/usr/bin/env python3
"""
adult_excluded distribution audit, per docs/STEP4_STATISTICAL_DESIGN.md
section 3: computes real within-tissue distributions from GTEx
(organ-matched + whole-body) and HPA (organ-matched, fills GTEx's Thymus
gap) so a concrete cutoff/quorum can be chosen from real numbers -- same
discipline as the P-developmental effect-size calibration (PR #9/#10):
report candidates, don't freeze blind.

Per STEP3_METHOD_CONTRACT.md: GTEx/HPA are used only as within-dataset
rank/percentile "adult-exclusion reference," never compared in raw
magnitude across platforms.

TWO real problems found and fixed before this reached review:

1. Tie-handling: a first draft used pandas `rank(pct=True,
   method="average")`. Bulk RNA-seq tissues routinely have >50% of genes
   at TPM=0 (verified: 56% for GTEx's Adrenal_Gland column alone) --
   "average" tie-breaking dumps that whole zero block at its *mid-rank*
   (~28th percentile for Adrenal_Gland), so 0 genes passed at any
   sensible low-percentile cutoff (5th/10th/25th) for nearly every organ.
2. Even after switching to `method="min"` (zero block -> percentile ~0),
   a deeper problem surfaced: with a zero-tie block covering ~55% of the
   distribution, EVERY candidate cutoff from 5 to 50 returned the exact
   same gene count -- there is no gene between the zero block (percentile
   ~0) and the next distinct expression level (percentile ~55%), so
   "percentile <= 25" and "percentile <= 50" select an identical set.
   Percentile rank is fundamentally degenerate when >50% of a
   distribution is tied at one value -- no tie-breaking convention fixes
   that, the percentile knob is just non-functional in that range.

Fix: split each tissue's genes into "not detected" (TPM/nTPM < 1 -- a
common low-expression convention in bulk RNA-seq; not claimed as any
platform's own official rule, this project's own choice, consistent with
the floor already used for the HPA placenta reference set in
hpa_placenta_enrichment.py) vs "detected" (TPM/nTPM >= 1). Percentile rank
is then computed ONLY within the detected subset, where the distribution
is continuous and percentile is a meaningful, non-degenerate knob. Final
criterion tested at each candidate P:
    adult_excluded(gene, tissue, P) =
        not_detected(gene, tissue) OR (detected AND percentile_among_detected(gene, tissue) <= P)

Usage: python3 adult_excluded_percentile_audit.py
(run on Argos, argos-codex env)
"""
import os
import zipfile
import pandas as pd
import numpy as np

ROOT = "/home/zz950/TWEAKR-OncoPlacental"
GTEX_TSV = "/home/zz950/DATA/1.Databases/GTEx_v11_median_tpm/processed/v0.1/gtex_v11_median_tpm_clean.tsv"
GTEX_MAP = "/home/zz950/DATA/1.Databases/GTEx_v11_median_tpm/processed/v0.1/hdma_organ_to_gtex_tissue_map.tsv"
HPA_ZIP = "/home/zz950/DATA/1.Databases/HPA_RNA_tissue_consensus/raw/rna_tissue_hpa.tsv.zip"
HPA_MAP = "/home/zz950/DATA/1.Databases/HPA_RNA_tissue_consensus/processed/v0.1/hdma_organ_to_hpa_tissue_map.tsv"
OUT_DIR = f"{ROOT}/results/04_dfp_signature/adult_excluded_audit"
os.makedirs(OUT_DIR, exist_ok=True)

ORGANS = ["Adrenal", "Thyroid", "Spleen", "Thymus", "Liver", "Skin", "StomachEsophagus"]
PCT_CANDIDATES = [5, 10, 25, 50]
NOT_DETECTED_FLOOR = 1.0


def detected_percentile(wide_df):
    """Per-column: boolean not_detected mask (value < floor), and percentile
    rank (0-100, min-rank tie-break) computed ONLY among detected (>=floor)
    values -- avoids the degenerate massive-zero-tie-block problem."""
    not_detected = wide_df < NOT_DETECTED_FLOOR
    pct = pd.DataFrame(index=wide_df.index, columns=wide_df.columns, dtype=float)
    for col in wide_df.columns:
        detected_vals = wide_df.loc[~not_detected[col], col]
        if len(detected_vals) > 0:
            pct.loc[detected_vals.index, col] = detected_vals.rank(pct=True, method="min") * 100
    return not_detected, pct


def excluded_mask(not_detected, pct, cols, P):
    """adult_excluded(gene, tissue, P) per the definition above, ANY of
    the given cols failing (i.e. gene must satisfy the criterion in EVERY
    listed tissue -- 'excluded_in_all(cols)')."""
    ok = pd.DataFrame(True, index=not_detected.index, columns=cols)
    for c in cols:
        detected_and_high = (~not_detected[c]) & (pct[c] > P)
        ok[c] = ~detected_and_high  # True = excluded-OK in this tissue (not-detected, or detected-but-low)
    return ok.all(axis=1)


# ---- Load GTEx: gene x tissue TPM, deduped by symbol (max per symbol) ----
gtex = pd.read_csv(GTEX_TSV, sep="\t")
tissue_cols = [c for c in gtex.columns if c not in ("ensembl_id", "ensembl_id_versioned", "symbol")]
print(f"GTEx: {gtex.shape[0]} rows x {len(tissue_cols)} tissues")
gtex_wide = gtex.groupby("symbol")[tissue_cols].max()
print(f"GTEx after symbol dedup (max per symbol): {gtex_wide.shape[0]} genes")
gtex_not_detected, gtex_pct = detected_percentile(gtex_wide)
print(f"GTEx: mean fraction not-detected (TPM<{NOT_DETECTED_FLOOR}) per tissue: {gtex_not_detected.mean().mean():.3f}")

# ---- Load HPA: same pivot as hpa_placenta_enrichment.py ----
with zipfile.ZipFile(HPA_ZIP) as z:
    inner = z.namelist()[0]
    with z.open(inner) as f:
        hpa = pd.read_csv(f, sep="\t")
hpa_wide = hpa.pivot_table(index="Gene name", columns="Tissue", values="nTPM", aggfunc="max")
print(f"HPA: {hpa_wide.shape[0]} genes x {hpa_wide.shape[1]} tissues")
hpa_not_detected, hpa_pct = detected_percentile(hpa_wide)
print(f"HPA: mean fraction not-detected (nTPM<{NOT_DETECTED_FLOOR}) per tissue: {hpa_not_detected.mean().mean():.3f}")

# ---- Load organ-matched mapping tables ----
gtex_map = pd.read_csv(GTEX_MAP, sep="\t")
hpa_map = pd.read_csv(HPA_MAP, sep="\t")


def gtex_cols_for(organ):
    row = gtex_map[gtex_map.hdma_organ_or_crc_relevant == organ]
    if row.empty or row.iloc[0].n_gtex_columns == 0:
        return []
    return row.iloc[0].gtex_tissue_columns.split(",")


def hpa_tissues_for(organ):
    row = hpa_map[hpa_map.hdma_organ_or_crc_relevant == organ]
    if row.empty:
        return []
    return row.iloc[0].hpa_tissue_names.split(",")


rows = []

# ============ Organ-matched (feeds F-developmental) ============
for organ in ORGANS:
    gcols = [c for c in gtex_cols_for(organ) if c in gtex_pct.columns]
    htissues = [t for t in hpa_tissues_for(organ) if t in hpa_pct.columns]
    for pct_cut in PCT_CANDIDATES:
        if gcols:
            excl_gtex = excluded_mask(gtex_not_detected, gtex_pct, gcols, pct_cut)
            n_gtex, n_gtex_universe = int(excl_gtex.sum()), gtex_pct.shape[0]
        else:
            excl_gtex, n_gtex, n_gtex_universe = None, None, None
        if htissues:
            excl_hpa = excluded_mask(hpa_not_detected, hpa_pct, htissues, pct_cut)
            n_hpa, n_hpa_universe = int(excl_hpa.sum()), hpa_pct.shape[0]
        else:
            excl_hpa, n_hpa, n_hpa_universe = None, None, None
        if gcols and htissues:
            shared = gtex_pct.index.intersection(hpa_pct.index)
            both = excl_gtex.reindex(shared).fillna(False) & excl_hpa.reindex(shared).fillna(False)
            n_both, n_both_universe = int(both.sum()), len(shared)
        else:
            n_both, n_both_universe = None, None
        rows.append({
            "check": "organ_matched", "organ": organ, "percentile_cutoff": pct_cut,
            "gtex_cols": ";".join(gcols) if gcols else "NONE",
            "hpa_tissues": ";".join(htissues) if htissues else "NONE",
            "n_pass_gtex": n_gtex, "n_universe_gtex": n_gtex_universe,
            "n_pass_hpa": n_hpa, "n_universe_hpa": n_hpa_universe,
            "n_pass_both_AND": n_both, "n_universe_both": n_both_universe,
        })

# ============ Whole-body (feeds P-developmental) ============
# Per STEP4_STATISTICAL_DESIGN.md section 3: GTEx's 68 tissues is the
# primary whole-body reference. HPA's 39 non-placenta tissues reported
# side by side since GTEx has NO Thymus column (documented gap) -- open
# question for reviewer sign-off, not resolved unilaterally here.
hpa_non_placenta = [c for c in hpa_pct.columns if c != "placenta"]
print(f"Whole-body universes: GTEx {len(tissue_cols)} tissues, HPA (non-placenta) {len(hpa_non_placenta)} tissues")

QUORUM_CANDIDATES = [("all", 0), ("all_but_1", 1), ("all_but_2", 2)]

for pct_cut in PCT_CANDIDATES:
    fails_gtex = pd.DataFrame(False, index=gtex_not_detected.index, columns=tissue_cols)
    for c in tissue_cols:
        fails_gtex[c] = (~gtex_not_detected[c]) & (gtex_pct[c] > pct_cut)
    fails_hpa = pd.DataFrame(False, index=hpa_not_detected.index, columns=hpa_non_placenta)
    for c in hpa_non_placenta:
        fails_hpa[c] = (~hpa_not_detected[c]) & (hpa_pct[c] > pct_cut)
    for label, allowed_fail in QUORUM_CANDIDATES:
        pass_gtex = fails_gtex.sum(axis=1) <= allowed_fail
        pass_hpa = fails_hpa.sum(axis=1) <= allowed_fail
        rows.append({
            "check": "whole_body", "organ": "ALL", "percentile_cutoff": pct_cut,
            "gtex_cols": f"quorum={label} (GTEx {len(tissue_cols)} tissues, no Thymus)",
            "hpa_tissues": f"quorum={label} (HPA {len(hpa_non_placenta)} non-placenta tissues, includes Thymus)",
            "n_pass_gtex": int(pass_gtex.sum()), "n_universe_gtex": gtex_pct.shape[0],
            "n_pass_hpa": int(pass_hpa.sum()), "n_universe_hpa": hpa_pct.shape[0],
            "n_pass_both_AND": None, "n_universe_both": None,
        })

result = pd.DataFrame(rows)
out_path = f"{OUT_DIR}/adult_excluded_percentile_audit.tsv"
result.to_csv(out_path, sep="\t", index=False)
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 20)
print(result.to_string(index=False))
print(f"\nWrote {out_path}")

# ---- Canonical marker cross-check ----
markers = ["ERVFRD-1", "CGA", "CSH1", "CSH2", "PSG1", "PSG3", "GATA3", "KRT7", "HLA-G"]
print("\n=== Canonical marker whole-body GTEx check (sanity) ===")
for m in markers:
    if m in gtex_wide.index:
        nd = gtex_not_detected.loc[m, tissue_cols]
        pc = gtex_pct.loc[m, tissue_cols]
        n_not_detected = int(nd.sum())
        detected_pct_vals = pc[~nd].dropna()
        med = detected_pct_vals.median() if len(detected_pct_vals) else float("nan")
        print(f"{m:10s} not-detected in {n_not_detected}/68 GTEx tissues; "
              f"among the {len(detected_pct_vals)} tissues where detected, median within-tissue percentile={med:.1f}")
    else:
        print(f"{m:10s} not found in GTEx symbol set")
