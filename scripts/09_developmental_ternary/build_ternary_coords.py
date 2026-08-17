#!/usr/bin/env python3
"""
Step 9: Fetal-Placenta-Adult developmental ternary map -- per-gene
ternary coordinate construction for all 3 tracks, per KB note
`2026-GPT-TWEAKR-Oncofetal.md#Fetal-Placenta-Adult ggtern图` and the
user-confirmed decision to reuse this project's own already-reviewed
data (Step 3/4/4a) for Tracks A/B, plus a new independent single-atlas
validation (HCL, GSE134355) as Track C.

Core math object (per the KB note): for each gene, three non-negative
"positive evidence" scores -- S_fetal, S_placenta, S_adult -- each a
ReLU'd, robustly-rescaled effect-size-based statistic (never a raw mean
expression, never a p-value-based statistic -- both explicitly rejected
in the note for the batch-contamination and sample-size-confound reasons
it gives), then closure-normalized (sum to 1) into ternary coordinates.

Track A -- Gut-specific (Fetal gut vs Placenta trophoblast vs Adult gut):
  Fetal/Adult: LargeInt_edgeR_primary.tsv (Gut Cell Atlas, real
  within-atlas fetal-vs-adult colon epithelium DE, Step 4a) -- logFC>0
  is fetal-high (confirmed against AFP=+9.85/+9.89 in that table).
  Placenta: Arutyunyan_edgeR_results.tsv + Nature2026_edgeR_results.tsv
  (Step 4's P-developmental primary DE), combined as the min of the two
  ReLU'd logFCs -- same "2-of-2 concordant replication" spirit already
  used to freeze P-developmental_primary84.

Track B -- Pan-tissue (embryo/fetal somatic vs trophoblast/placenta vs
adult somatic): Fetal/Adult use HDMA's 7-organ per-gene percentile
machinery (Step 4's own f_developmental_calibration.py, reimplemented
here reading the same underlying files) -- NOT a raw cross-platform DE
(explicitly forbidden by this project's own Step 3 method contract: GTEx/
HPA bulk and HDMA scRNA-seq pseudobulk can never be compared by raw
magnitude, only within-dataset percentile rank). diff_organ(g) =
within-organ fetal expression percentile minus organ-matched adult
GTEx/HPA percentile, for each of the 7 organs; S_fetal = ReLU(max over
organs), S_adult = ReLU(-min over organs) -- "best fetal-favoring organ"
and "most adult-favoring organ" respectively, matching F_developmental_
union's own "any organ passes" logic. Placenta: same as Track A.

Track C -- HCL independent validation (GSE134355, Han et al. 2020):
genuine one-vs-rest edgeR DE within ONE uniformly-processed atlas
(hcl_edgeR_onevsrest.R) -- Fetal/Placenta/Adult intestine all same
platform/pipeline, avoiding the cross-dataset batch risk Tracks A/B
cannot fully escape. Real, honest limitation: Placenta n=1 (single
donor, no replication -- effect size only, no meaningful p-value); HCL's
5 fetal-intestine samples are not regionally resolved (no colon/SI
split, unlike Gut Cell Atlas).

Within each track, all 3 raw axis scores are independently rescaled by
their own 99th percentile across that track's gene universe (clipped at
1) before closure -- removes arbitrary absolute-scale differences
between different DE effect-size sources (e.g. Track B mixes a
percentile-differential-based Fetal/Adult score with a logFC-based
Placenta score) without claiming a false statistical equivalence between
them. This is a real methodological simplification for this first
version, stated explicitly here and in the write-up, not hidden.

Usage: python3 build_ternary_coords.py
(runs locally -- reuses already-computed DE/pseudobulk tables, no new
heavy compute)
"""
import zipfile
import numpy as np
import pandas as pd

REPO = "/Volumes/Stelligen_SSD/Stelligen/PR/TWEAKR-OncoPlacental"
DATA = "/Volumes/Stelligen_SSD/Stelligen/DATA"
OUT_DIR = f"{REPO}/results/09_developmental_ternary"
EPS = 0.01

MARKERS = ["AFP", "LIN28B", "CSH1", "CSH2", "CGA", "PSG1", "PSG3",
           "ERVW-1", "ERVFRD-1", "KRT7", "TRIM71"]


def relu(x):
    return np.maximum(x, 0.0)


def robust_rescale(s):
    """Divide by the series' own 99th percentile (min 1e-9 to avoid
    div-by-zero), clip at 1. See module docstring for why this is
    needed (cross-source scale harmonization within one track)."""
    p99 = max(s.quantile(0.99), 1e-9)
    return (s / p99).clip(upper=1.0)


def closure(fetal_raw, placenta_raw, adult_raw):
    df = pd.DataFrame({"Fetal": fetal_raw, "Placenta": placenta_raw, "Adult": adult_raw})
    df = df.add(EPS)
    total = df.sum(axis=1)
    coords = df.div(total, axis=0)
    return coords


def annotate_markers(df):
    # NOTE: consumers reading this TSV back MUST pass keep_default_na=
    # False -- pandas' default read_csv NA-sniffing turns "" (and "NA",
    # "None", "null", etc.) into NaN, which silently breaks any
    # `!= ""` filter on this column downstream. Real bug, caught while
    # reviewing the first rendered plot (see plot_ternary.py).
    df["marker"] = df.index.to_series().apply(lambda g: g if g in MARKERS else "")
    return df


# ============ Track A: Gut-specific ============
def build_track_a():
    large_int = pd.read_csv(f"{REPO}/results/04a_dfp_gut/edgeR/LargeInt_edgeR_primary.tsv", sep="\t").set_index("gene")
    aru = pd.read_csv(f"{REPO}/results/04_dfp_signature/edgeR/Arutyunyan_edgeR_results.tsv", sep="\t").set_index("gene")
    nat = pd.read_csv(f"{REPO}/results/04_dfp_signature/edgeR/Nature2026_edgeR_results.tsv", sep="\t").set_index("gene")

    universe = large_int.index.intersection(aru.index).intersection(nat.index)
    print(f"[Track A] gene universe (tested in all 3 DE tables): {len(universe)}")

    fetal_raw = robust_rescale(relu(large_int.loc[universe, "logFC"]))
    adult_raw = robust_rescale(relu(-large_int.loc[universe, "logFC"]))
    placenta_raw_unscaled = pd.concat([
        relu(aru.loc[universe, "logFC"]), relu(nat.loc[universe, "logFC"])
    ], axis=1).min(axis=1)
    placenta_raw = robust_rescale(placenta_raw_unscaled)

    coords = closure(fetal_raw, placenta_raw, adult_raw)
    coords = annotate_markers(coords)
    return coords


# ============ Track B: Pan-tissue ============
NOT_DETECTED_FLOOR = 1.0
ORGANS = ["Adrenal", "Thyroid", "Spleen", "Thymus", "Liver", "Skin", "Stomach"]


def detected_percentile(wide_df):
    not_detected = wide_df < NOT_DETECTED_FLOOR
    pct = pd.DataFrame(index=wide_df.index, columns=wide_df.columns, dtype=float)
    for col in wide_df.columns:
        vals = wide_df.loc[~not_detected[col], col]
        if len(vals) > 0:
            pct.loc[vals.index, col] = vals.rank(pct=True, method="min") * 100
    return not_detected, pct


def build_track_b():
    gtex = pd.read_csv(f"{DATA}/1.Databases/GTEx_v11_median_tpm/processed/v0.1/gtex_v11_median_tpm_clean.tsv", sep="\t")
    tissue_cols = [c for c in gtex.columns if c not in ("ensembl_id", "ensembl_id_versioned", "symbol")]
    gtex_wide = gtex.groupby("symbol")[tissue_cols].max()
    gtex_not_detected, gtex_pct = detected_percentile(gtex_wide)

    with zipfile.ZipFile(f"{DATA}/1.Databases/HPA_RNA_tissue_consensus/raw/rna_tissue_hpa.tsv.zip") as z:
        inner = z.namelist()[0]
        with z.open(inner) as f:
            hpa = pd.read_csv(f, sep="\t")
    hpa_wide = hpa.pivot_table(index="Gene name", columns="Tissue", values="nTPM", aggfunc="max")
    hpa_not_detected, hpa_pct = detected_percentile(hpa_wide)

    gtex_map = pd.read_csv(f"{DATA}/1.Databases/GTEx_v11_median_tpm/processed/v0.1/hdma_organ_to_gtex_tissue_map.tsv", sep="\t")
    hpa_map = pd.read_csv(f"{DATA}/1.Databases/HPA_RNA_tissue_consensus/processed/v0.1/hdma_organ_to_hpa_tissue_map.tsv", sep="\t")

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

    def adult_pct_for_organ(organ):
        """Max percentile-among-detected across the organ's matched GTEx+HPA
        tissue columns (the 'worst case' adult expression -- matches the
        AND-across-tissues logic the existing adult_excluded_mask gate
        uses: a gene fails whole-organ-matched exclusion if it is highly
        expressed in ANY one matched tissue)."""
        cols = []
        pcts = []
        for c in gtex_cols_for(organ):
            pcts.append(gtex_pct[c].where(~gtex_not_detected[c], 0.0))
        for t in hpa_tissues_for(organ):
            pcts.append(hpa_pct[t].where(~hpa_not_detected[t], 0.0))
        if not pcts:
            return None
        combined = pd.concat(pcts, axis=1)
        return combined.max(axis=1)

    organ_diffs = {}
    for organ in ORGANS:
        counts = pd.read_csv(f"{REPO}/results/04_dfp_signature/hdma_pseudobulk/{organ}_pseudobulk_counts.tsv", sep="\t", index_col=0)
        cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
        hdma_not_detected, hdma_pct = detected_percentile(cpm)
        median_pct = hdma_pct.median(axis=1, skipna=True)  # elevated_pct, per gene

        adult_pct = adult_pct_for_organ(organ)
        if adult_pct is None:
            print(f"[Track B] {organ}: no matched adult reference, skipped")
            continue
        common = median_pct.index.intersection(adult_pct.index)
        diff = median_pct.loc[common] - adult_pct.loc[common]
        organ_diffs[organ] = diff
        print(f"[Track B] {organ}: {len(common)} genes with both fetal-elevation and adult-percentile evidence")

    diff_table = pd.DataFrame(organ_diffs)  # gene x organ
    fetal_diff_max = diff_table.max(axis=1, skipna=True)
    fetal_diff_min = diff_table.min(axis=1, skipna=True)
    # A gene can be in diff_table's index (detected in >=1 organ's own
    # sample set) yet have NO organ where it was ALSO in that organ's
    # matched adult reference's detected set -- max/min over an all-NaN
    # row stays NaN. Real, not a bug: no fetal-vs-adult evidence exists
    # for that gene in any organ, so it cannot get a ternary coordinate.
    has_evidence = fetal_diff_max.notna() & fetal_diff_min.notna()
    n_dropped = (~has_evidence).sum()
    print(f"[Track B] pan-organ gene universe: {len(diff_table)} "
          f"({n_dropped} dropped -- no organ had both fetal AND matched-adult evidence)")
    fetal_diff_max = fetal_diff_max[has_evidence]
    fetal_diff_min = fetal_diff_min[has_evidence]

    aru = pd.read_csv(f"{REPO}/results/04_dfp_signature/edgeR/Arutyunyan_edgeR_results.tsv", sep="\t").set_index("gene")
    nat = pd.read_csv(f"{REPO}/results/04_dfp_signature/edgeR/Nature2026_edgeR_results.tsv", sep="\t").set_index("gene")

    universe = fetal_diff_max.index.intersection(aru.index).intersection(nat.index)
    print(f"[Track B] final gene universe (pan-organ AND placenta-tested): {len(universe)}")

    fetal_raw = robust_rescale(relu(fetal_diff_max.loc[universe]) / 100.0)
    adult_raw = robust_rescale(relu(-fetal_diff_min.loc[universe]) / 100.0)
    placenta_raw_unscaled = pd.concat([
        relu(aru.loc[universe, "logFC"]), relu(nat.loc[universe, "logFC"])
    ], axis=1).min(axis=1)
    placenta_raw = robust_rescale(placenta_raw_unscaled)

    coords = closure(fetal_raw, placenta_raw, adult_raw)
    coords = annotate_markers(coords)
    return coords, diff_table


# ============ Track C: HCL independent ============
def build_track_c():
    f = pd.read_csv(f"{OUT_DIR}/hcl_edgeR_Fetal_vs_rest.tsv", sep="\t").set_index("gene")
    p = pd.read_csv(f"{OUT_DIR}/hcl_edgeR_Placenta_vs_rest.tsv", sep="\t").set_index("gene")
    a = pd.read_csv(f"{OUT_DIR}/hcl_edgeR_Adult_vs_rest.tsv", sep="\t").set_index("gene")

    universe = f.index.intersection(p.index).intersection(a.index)
    print(f"[Track C] gene universe (HCL one-vs-rest, filterByExpr-passed): {len(universe)}")

    fetal_raw = robust_rescale(relu(f.loc[universe, "logFC"]))
    placenta_raw = robust_rescale(relu(p.loc[universe, "logFC"]))
    adult_raw = robust_rescale(relu(a.loc[universe, "logFC"]))

    coords = closure(fetal_raw, placenta_raw, adult_raw)
    coords = annotate_markers(coords)
    return coords


def main():
    print("=== Track A: Gut-specific ===")
    coords_a = build_track_a()
    coords_a.to_csv(f"{OUT_DIR}/track_A_gut_specific_coords.tsv", sep="\t")
    print(f"Wrote track_A_gut_specific_coords.tsv ({len(coords_a)} genes)\n")

    print("=== Track B: Pan-tissue ===")
    coords_b, diff_table = build_track_b()
    coords_b.to_csv(f"{OUT_DIR}/track_B_pantissue_coords.tsv", sep="\t")
    diff_table.to_csv(f"{OUT_DIR}/track_B_organ_diff_table.tsv", sep="\t")
    print(f"Wrote track_B_pantissue_coords.tsv ({len(coords_b)} genes)\n")

    print("=== Track C: HCL independent validation ===")
    coords_c = build_track_c()
    coords_c.to_csv(f"{OUT_DIR}/track_C_hcl_coords.tsv", sep="\t")
    print(f"Wrote track_C_hcl_coords.tsv ({len(coords_c)} genes)\n")

    print("=== Marker sanity check (all 3 tracks) ===")
    for name, coords in [("A (gut)", coords_a), ("B (pan-tissue)", coords_b), ("C (HCL)", coords_c)]:
        print(f"-- Track {name} --")
        for m in MARKERS:
            if m in coords.index:
                row = coords.loc[m]
                print(f"  {m}: Fetal={row.Fetal:.3f} Placenta={row.Placenta:.3f} Adult={row.Adult:.3f}")


if __name__ == "__main__":
    main()
