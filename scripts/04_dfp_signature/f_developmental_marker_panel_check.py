#!/usr/bin/env python3
"""
Multi-organ positive-control marker panel check, per the PR #12 reviewer's
explicit request: AFP alone can't discriminate elevated_pct=75 vs 90 (its
percentile is 99.7, clears both) and can't prove adult_excl_pct=25 is
globally optimal, not just locally sufficient for one gene. Needs a real,
predefined, cross-organ panel checked against all 4 candidate combos.

Marker panel (predefined BEFORE looking at results, not tuned after):
  - Cross-organ pan-fetal/oncofetal reactivation genes, well-established in
    the imprinted-gene and oncofetal-reactivation literature, expected to
    be silenced or markedly reduced in most adult tissues and elevated in
    fetal tissue broadly (not organ-specific developmental biology, which
    is much harder to source per-organ markers for beyond liver):
      DLK1  -- imprinted, classic fetal/preadipocyte marker, silenced postnatally in most tissues
      IGF2  -- imprinted fetal growth factor, silenced postnatally (adult expression retained mainly in choroid plexus/leptomeninges)
      H19   -- imprinted, co-regulated with IGF2, fetal-restricted in most tissues
      LIN28B -- fetal/oncofetal reactivation gene (fetal hematopoiesis, HCC, widely used fetal-vs-adult marker)
      PEG10 -- imprinted, retrotransposon-derived, fetal/placental and broadly reduced postnatally
  - Liver-specific additions (established liver oncofetal markers):
      AFP   -- already used as the primary sanity check
      GPC3  -- oncofetal, hepatoblastoma/HCC marker, fetal liver expression

Confidence note (explicit, not glossed over): the cross-organ panel is
well-supported for broad fetal-vs-adult reactivation biology but is NOT
organ-specific developmental biology for Adrenal/Thyroid/Skin/Stomach --
organ-specific fetal-restricted markers for those organs are much less
established in the literature than liver's AFP/GPC3. Treated as a
pan-fetal biology check, not a per-organ developmental-identity check.

For each of the 4 candidate combos (elev75/adult10, elev75/adult25,
elev90/adult10, elev90/adult25), reports: does each panel gene qualify as
an F-developmental candidate in each organ where it's tested (organ ->
gene -> pass/fail), plus the retention rate.

Usage: python3 f_developmental_marker_panel_check.py
(run on Argos, argos-codex env)
"""
import pandas as pd

exec(open("/home/zz950/TWEAKR-OncoPlacental/scripts/04_dfp_signature/f_developmental_calibration.py").read().split("result = pd.DataFrame")[0])

CROSS_ORGAN_PANEL = ["DLK1", "IGF2", "H19", "LIN28B", "PEG10"]
LIVER_EXTRA = ["AFP", "GPC3"]

COMBOS = [(75, 10), (75, 25), (90, 10), (90, 25)]

rows = []
for organ in ORGANS:
    counts = pd.read_csv(f"{HDMA_DIR}/{organ}_pseudobulk_counts.tsv", sep="\t", index_col=0)
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    nd, pct = detected_percentile(cpm)
    median_pct = pct.median(axis=1, skipna=True)
    detected_frac = (~nd).mean(axis=1)

    panel = CROSS_ORGAN_PANEL + (LIVER_EXTRA if organ == "Liver" else [])

    for elev_pct, adult_pct in COMBOS:
        elevated = (median_pct >= elev_pct) & (detected_frac >= 1.0)
        mask, provenance = adult_excluded_mask(organ, adult_pct)
        n_candidates = int((elevated & mask.reindex(elevated.index).fillna(False)).sum()) if mask is not None else None

        for gene in panel:
            in_hdma = gene in cpm.index
            gene_elevated = bool(elevated.get(gene, False)) if in_hdma else None
            gene_median_pct = float(median_pct.get(gene)) if in_hdma and gene in median_pct.index and pd.notna(median_pct.get(gene)) else None
            gene_excluded = bool(mask.get(gene, False)) if (mask is not None and gene in mask.index) else None
            passes = bool(gene_elevated and gene_excluded) if (gene_elevated is not None and gene_excluded is not None) else False
            rows.append({
                "organ": organ, "provenance": provenance, "elev_pct": elev_pct, "adult_pct": adult_pct,
                "gene": gene, "in_hdma_universe": in_hdma, "hdma_median_pct": gene_median_pct,
                "elevated": gene_elevated, "adult_excluded": gene_excluded, "passes_F_developmental": passes,
                "n_candidates_this_combo": n_candidates,
            })

result = pd.DataFrame(rows)
out_path = "/home/zz950/TWEAKR-OncoPlacental/results/04_dfp_signature/f_developmental_audit/marker_panel_check.tsv"
result.to_csv(out_path, sep="\t", index=False)

pd.set_option("display.width", 220)
pd.set_option("display.max_rows", 500)

# Summary: retention rate per combo (fraction of gene x organ tests that pass)
summary = result.groupby(["elev_pct", "adult_pct"]).agg(
    n_tests=("passes_F_developmental", "size"),
    n_pass=("passes_F_developmental", "sum"),
    mean_n_candidates=("n_candidates_this_combo", "mean"),
).reset_index()
summary["retention_rate"] = (summary["n_pass"] / summary["n_tests"]).round(3)
print("=== Retention summary per combo ===")
print(summary.to_string(index=False))

print("\n=== Full per-gene-per-organ results ===")
print(result[["organ", "provenance", "elev_pct", "adult_pct", "gene", "hdma_median_pct", "elevated", "adult_excluded", "passes_F_developmental"]].to_string(index=False))

print(f"\nWrote {out_path}")
