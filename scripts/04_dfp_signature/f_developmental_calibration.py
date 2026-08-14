#!/usr/bin/env python3
"""
F-developmental "elevated_in_fetal_somatic" calibration, combined with
adult_excluded, per docs/STEP4_STATISTICAL_DESIGN.md section 2 and the
PR #11 reviewer's three design decisions:

1. StomachEsophagus split into Stomach (n=6, used) and Esophagus (n=1,
   insufficient-replication, reported but not used) -- pseudobulk already
   split locally, mapping tables already split
   (hdma_organ_to_{gtex,hpa}_tissue_map.tsv).
2. Thymus's organ-matched adult exclusion uses HPA only (no GTEx column)
   -- provenance tracked explicitly per organ, never silently disguised as
   GTEx-equivalent.
3. Whole-body adult exclusion NOT used here (F-developmental only needs
   organ-matched, per STEP4_STATISTICAL_DESIGN.md section 3) -- no GTEx/
   HPA-as-69th-tissue patching needed at this stage.

Method:
  - HDMA pseudobulk (raw counts) -> CPM per sample -> per organ: same
    not-detected-floor (CPM<1) + percentile-among-detected design already
    used for adult_excluded (avoids the same zero-inflation degeneracy;
    checked directly: 8-10% exact zero, ~32% below CPM=1 in HDMA
    pseudobulk -- less extreme than GTEx's single-tissue TPM but still
    real, so the same fix applies for consistency and correctness).
  - "elevated" candidate: median per-sample percentile-among-detected >=
    cutoff, AND detected (CPM>=1) in >= quorum fraction of the organ's own
    samples (uses the individual-level replicate structure directly).
  - adult_excluded (organ-matched): reused from
    adult_excluded_percentile_audit.py's detected_percentile() logic,
    now against the split Stomach/Esophagus mapping. Provenance
    (GTEx+HPA / GTEx-only / HPA-only) reported explicitly per organ.
  - F-developmental candidate(gene, organ) = elevated AND adult_excluded,
    reported across a small grid of candidate cutoffs -- not frozen here.

Usage: python3 f_developmental_calibration.py
(run on Argos, argos-codex env)
"""
import os
import zipfile
import pandas as pd
import numpy as np

ROOT = "/home/zz950/TWEAKR-OncoPlacental"
HDMA_DIR = f"{ROOT}/results/04_dfp_signature/hdma_pseudobulk"
GTEX_TSV = "/home/zz950/DATA/1.Databases/GTEx_v11_median_tpm/processed/v0.1/gtex_v11_median_tpm_clean.tsv"
GTEX_MAP = "/home/zz950/DATA/1.Databases/GTEx_v11_median_tpm/processed/v0.1/hdma_organ_to_gtex_tissue_map.tsv"
HPA_ZIP = "/home/zz950/DATA/1.Databases/HPA_RNA_tissue_consensus/raw/rna_tissue_hpa.tsv.zip"
HPA_MAP = "/home/zz950/DATA/1.Databases/HPA_RNA_tissue_consensus/processed/v0.1/hdma_organ_to_hpa_tissue_map.tsv"
OUT_DIR = f"{ROOT}/results/04_dfp_signature/f_developmental_audit"
os.makedirs(OUT_DIR, exist_ok=True)

ORGANS = ["Adrenal", "Thyroid", "Spleen", "Thymus", "Liver", "Skin", "Stomach"]  # Esophagus excluded (n=1)
NOT_DETECTED_FLOOR = 1.0
ELEVATED_PCT_CANDIDATES = [75, 90]
QUORUM_CANDIDATES = [0.5, 1.0]  # fraction of organ's own samples gene must be detected in
ADULT_EXCL_PCT_CANDIDATES = [10, 25]


def detected_percentile(wide_df):
    not_detected = wide_df < NOT_DETECTED_FLOOR
    pct = pd.DataFrame(index=wide_df.index, columns=wide_df.columns, dtype=float)
    for col in wide_df.columns:
        vals = wide_df.loc[~not_detected[col], col]
        if len(vals) > 0:
            pct.loc[vals.index, col] = vals.rank(pct=True, method="min") * 100
    return not_detected, pct


# ---- adult reference (GTEx + HPA), same as adult_excluded_percentile_audit.py ----
gtex = pd.read_csv(GTEX_TSV, sep="\t")
tissue_cols = [c for c in gtex.columns if c not in ("ensembl_id", "ensembl_id_versioned", "symbol")]
gtex_wide = gtex.groupby("symbol")[tissue_cols].max()
gtex_not_detected, gtex_pct = detected_percentile(gtex_wide)

with zipfile.ZipFile(HPA_ZIP) as z:
    inner = z.namelist()[0]
    with z.open(inner) as f:
        hpa = pd.read_csv(f, sep="\t")
hpa_wide = hpa.pivot_table(index="Gene name", columns="Tissue", values="nTPM", aggfunc="max")
hpa_not_detected, hpa_pct = detected_percentile(hpa_wide)

gtex_map = pd.read_csv(GTEX_MAP, sep="\t")
hpa_map = pd.read_csv(HPA_MAP, sep="\t")


def gtex_cols_for(organ):
    row = gtex_map[gtex_map.hdma_organ_or_crc_relevant == organ]
    if row.empty or row.iloc[0].n_gtex_columns == 0:
        return []
    return [c for c in row.iloc[0].gtex_tissue_columns.split(",") if c in gtex_pct.columns]


def hpa_tissues_for(organ):
    row = hpa_map[hpa_map.hdma_organ_or_crc_relevant == organ]
    if row.empty:
        return []
    return [t for t in row.iloc[0].hpa_tissue_names.split(",") if t in hpa_pct.columns]


def adult_excluded_mask(organ, pct_cut):
    """Coverage-aware (fixed per PR #12 round 2 review): the original
    version intersected GTEx's and HPA's gene panels first, so a gene
    present in GTEx but simply absent from HPA's gene panel entirely
    (e.g. H19 -- not a detection failure, HPA just doesn't measure it)
    was silently dropped and could never pass -- conflating "platform
    didn't measure this gene" with "gene failed the test". Fixed: each
    gene's provenance is now decided by which platform(s) actually
    include it in their gene panel (gtex_pct.index / hpa_pct.index),
    not by a blanket per-organ label.

    Returns (mask, provenance) -- both pandas Series indexed by gene.
    provenance in {"GTEx+HPA","GTEx-only","HPA-only"}; genes covered by
    neither platform are absent from the index entirely (unresolved,
    can never pass -- not silently marked either True or False)."""
    gcols = gtex_cols_for(organ)
    htissues = hpa_tissues_for(organ)

    def excl(not_det, pct, cols):
        ok = pd.DataFrame(True, index=not_det.index, columns=cols)
        for c in cols:
            fails = (~not_det[c]) & (pct[c] > pct_cut)
            ok[c] = ~fails
        return ok.all(axis=1)

    m_gtex = excl(gtex_not_detected, gtex_pct, gcols) if gcols else None
    m_hpa = excl(hpa_not_detected, hpa_pct, htissues) if htissues else None

    if m_gtex is not None and m_hpa is not None:
        gtex_genes = set(m_gtex.index)
        hpa_genes = set(m_hpa.index)
        both = gtex_genes & hpa_genes
        gtex_only = gtex_genes - hpa_genes
        hpa_only = hpa_genes - gtex_genes
        mask = pd.concat([
            m_gtex.loc[list(both)] & m_hpa.loc[list(both)],
            m_gtex.loc[list(gtex_only)],
            m_hpa.loc[list(hpa_only)],
        ])
        prov = pd.Series(
            ["GTEx+HPA"] * len(both) + ["GTEx-only"] * len(gtex_only) + ["HPA-only"] * len(hpa_only),
            index=list(both) + list(gtex_only) + list(hpa_only),
        )
        return mask, prov
    elif m_gtex is not None:
        return m_gtex, pd.Series("GTEx-only", index=m_gtex.index)
    elif m_hpa is not None:
        return m_hpa, pd.Series("HPA-only", index=m_hpa.index)
    else:
        return None, None


# ---- HDMA per-organ elevated calibration ----
rows = []
for organ in ORGANS:
    counts = pd.read_csv(f"{HDMA_DIR}/{organ}_pseudobulk_counts.tsv", sep="\t", index_col=0)
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    hdma_not_detected, hdma_pct = detected_percentile(cpm)
    n_samples = cpm.shape[1]

    # median percentile-among-detected across samples (NaN where never detected -> treated as not elevated)
    median_pct = hdma_pct.median(axis=1, skipna=True)
    detected_frac = (~hdma_not_detected).mean(axis=1)

    for elev_pct in ELEVATED_PCT_CANDIDATES:
        for quorum in QUORUM_CANDIDATES:
            elevated = (median_pct >= elev_pct) & (detected_frac >= quorum)
            n_elevated = int(elevated.sum())
            for adult_pct in ADULT_EXCL_PCT_CANDIDATES:
                mask, provenance = adult_excluded_mask(organ, adult_pct)
                if mask is None:
                    rows.append({
                        "organ": organ, "n_samples": n_samples, "elevated_pct": elev_pct,
                        "quorum": quorum, "n_elevated_hdma_only": n_elevated,
                        "adult_excl_pct": adult_pct, "adult_ref_provenance": "NONE",
                        "n_F_developmental_candidates": None,
                    })
                    continue
                elevated_genes = set(elevated.index[elevated])
                shared = set(mask.index) & elevated_genes
                excluded_genes = set(mask.index[mask])
                f_dev = shared & excluded_genes
                # provenance breakdown among the actual candidate genes (coverage-aware fix)
                prov_counts = provenance.loc[list(f_dev)].value_counts().to_dict() if f_dev else {}
                prov_summary = ";".join(f"{k}={v}" for k, v in sorted(prov_counts.items())) or "NONE"
                rows.append({
                    "organ": organ, "n_samples": n_samples, "elevated_pct": elev_pct,
                    "quorum": quorum, "n_elevated_hdma_only": n_elevated,
                    "adult_excl_pct": adult_pct, "adult_ref_provenance": prov_summary,
                    "n_F_developmental_candidates": len(f_dev),
                })

result = pd.DataFrame(rows)
out_path = f"{OUT_DIR}/f_developmental_calibration.tsv"
result.to_csv(out_path, sep="\t", index=False)
pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 20)
print(result.to_string(index=False))
print(f"\nWrote {out_path}")

# ---- Sanity check: AFP in Liver -- THE textbook oncofetal gene. Should be
# elevated in fetal liver and adult-excluded (essentially undetectable in
# normal adult liver). ----
print("\n=== Sanity check: AFP in Liver (textbook oncofetal gene) ===")
liver_counts = pd.read_csv(f"{HDMA_DIR}/Liver_pseudobulk_counts.tsv", sep="\t", index_col=0)
liver_cpm = liver_counts.div(liver_counts.sum(axis=0), axis=1) * 1e6
if "AFP" in liver_cpm.index:
    print(f"AFP CPM per fetal Liver sample: {liver_cpm.loc['AFP'].round(1).to_dict()}")
liver_gcols = gtex_cols_for("Liver")
liver_htissues = hpa_tissues_for("Liver")
if "AFP" in gtex_wide.index:
    print(f"AFP GTEx Liver TPM: {gtex_wide.loc['AFP', liver_gcols].to_dict()}")
if "AFP" in hpa_wide.index:
    print(f"AFP HPA liver nTPM: {hpa_wide.loc['AFP', liver_htissues].to_dict()}")
