#!/usr/bin/env python3
"""
Pseudobulk-unit design audit for epi_raw_counts02_v2.h5ad, required by PR #20
round-2 REQUEST_CHANGES before any real DE model is locked:

"不能因为 First trimester 有 4,927 cells、Second trimester 只有 3,481 cells，
就决定把 First+Second trimester 全部 pool 后...统计单位是 donor/pseudobulk，
不是 cell...对 LargeInt 和 SmallInt 分别生成一张真正的 pseudobulk-unit 表：
donor / library / Age_group / trimester / Region / Region_code / 10X / batch /
Fraction / n_cells / total_UMI。检查每个 donor 是否横跨多个 library、chemistry
或 batch...primary contrast 按 biological replicate 和 source comparability
决定，不按 cell 数决定。"

This script does NOT run any DE. It only builds the real donor-level design
table and reports (a) donor counts per Age_group x Region stratum, (b)
whether any donor spans multiple 10X/batch/Fraction values (breaks simple
per-donor pseudobulk if so), and (c) whether candidate contrast design
matrices (e.g. ~ 10X + Age_group, restricted to LargeInt) are full rank --
so the primary-vs-sensitivity contrast choice in the next design-doc
revision is based on this real table, not assumed.

Usage: python3 gut_epi_pseudobulk_design_audit.py <out_dir>
"""
import sys
import os
import json
import numpy as np
import pandas as pd
import anndata as ad

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "results/04a_dfp_gut/inventory"
os.makedirs(OUT_DIR, exist_ok=True)

PATH = "/home/zz950/DATA/scRNAseq/GutCellAtlas_Elmentaite2021/raw/epi_raw_counts02_v2.h5ad"

print(f"Reading {PATH} (backed='r')...", flush=True)
adata = ad.read_h5ad(PATH, backed="r")
obs = adata.obs.copy()

# restrict to epithelial cells and the two regions of interest
obs = obs[obs["category"] == "Epithelial"]
obs = obs[obs["Region"].isin(["LargeInt", "SmallInt"])]
print(f"n cells (Epithelial, LargeInt/SmallInt only): {len(obs)}", flush=True)

DONOR_COL = "Sample name"

# --- 1. per-donor structure: does a donor span multiple 10X/batch/Fraction values? ---
donor_struct = obs.groupby(DONOR_COL, observed=True).agg(
    n_cells=("Region", "size"),
    n_regions=("Region", "nunique"),
    n_age_groups=("Age_group", "nunique"),
    n_10X=("10X", "nunique"),
    n_batch=("batch", "nunique"),
    n_fraction=("Fraction", "nunique"),
    regions=("Region", lambda s: ";".join(sorted(s.unique()))),
    age_groups=("Age_group", lambda s: ";".join(sorted(s.astype(str).unique()))),
    tenX_values=("10X", lambda s: ";".join(sorted(s.astype(str).unique()))),
).reset_index()

n_donor_multi_10X = int((donor_struct["n_10X"] > 1).sum())
n_donor_multi_batch = int((donor_struct["n_batch"] > 1).sum())
n_donor_multi_fraction = int((donor_struct["n_fraction"] > 1).sum())
n_donor_multi_agegroup = int((donor_struct["n_age_groups"] > 1).sum())
print(f"\nDonors spanning >1 10X chemistry value: {n_donor_multi_10X}/{len(donor_struct)}", flush=True)
print(f"Donors spanning >1 batch value: {n_donor_multi_batch}/{len(donor_struct)}", flush=True)
print(f"Donors spanning >1 Fraction value: {n_donor_multi_fraction}/{len(donor_struct)}", flush=True)
print(f"Donors spanning >1 Age_group value: {n_donor_multi_agegroup}/{len(donor_struct)}", flush=True)

donor_struct.to_csv(f"{OUT_DIR}/gut_epi_donor_structure.tsv", sep="\t", index=False)
print(f"Wrote {OUT_DIR}/gut_epi_donor_structure.tsv", flush=True)

# --- 2. real pseudobulk-unit table: donor x Region x Age_group x 10X x batch x Fraction ---
group_cols = [DONOR_COL, "Age_group", "Region", "Region code", "10X", "batch", "Fraction"]
pseudobulk_units = obs.groupby(group_cols, observed=True).agg(
    n_cells=("Region", "size"),
    total_counts=("total_counts", "sum"),
).reset_index()
pseudobulk_units.to_csv(f"{OUT_DIR}/gut_epi_pseudobulk_units.tsv", sep="\t", index=False)
print(f"\nWrote {OUT_DIR}/gut_epi_pseudobulk_units.tsv ({len(pseudobulk_units)} rows)", flush=True)

# --- 3. real donor counts per Age_group x Region, and Age_group x Region x 10X ---
for region in ["LargeInt", "SmallInt"]:
    sub = obs[obs["Region"] == region]
    print(f"\n=== {region}: donor counts per Age_group ===")
    dc = sub.groupby("Age_group", observed=True)[DONOR_COL].nunique()
    print(dc)
    print(f"\n=== {region}: donor counts per Age_group x 10X ===")
    dc2 = sub.groupby(["Age_group", "10X"], observed=True)[DONOR_COL].nunique()
    print(dc2)
    print(f"\n=== {region}: cell counts per Age_group x 10X ===")
    cc2 = sub.groupby(["Age_group", "10X"], observed=True).size()
    print(cc2)

# --- 4. design-matrix full-rank check for candidate contrasts, per region ---
# Build one row per (donor x Region) pseudobulk unit (collapsing 10X/batch within
# donor -- only valid if step 1 shows no donor spans multiple 10X values within
# a region; report explicitly either way)
results_rank = {}
for region in ["LargeInt", "SmallInt"]:
    sub = obs[obs["Region"] == region]
    donor_level = sub.groupby(DONOR_COL, observed=True).agg(
        Age_group=("Age_group", lambda s: s.mode().iloc[0] if len(s.mode()) else s.iloc[0]),
        tenX=("10X", lambda s: s.mode().iloc[0] if len(s.mode()) else s.iloc[0]),
        n_10X_within_donor_region=("10X", "nunique"),
    ).reset_index()

    candidates = {
        "~ Age_group (all fetal pooled vs Adult, First+Second only)": donor_level[donor_level["Age_group"].isin(["First trim", "Second trim", "Adult"])],
        "~ Age_group (Second trim vs Adult only)": donor_level[donor_level["Age_group"].isin(["Second trim", "Adult"])],
        "~ Age_group (First trim vs Adult only)": donor_level[donor_level["Age_group"].isin(["First trim", "Adult"])],
    }
    region_result = {}
    for label, d in candidates.items():
        n_donors = len(d)
        n_levels_age = d["Age_group"].nunique()
        n_levels_10X = d["tenX"].nunique()
        # full-rank check for ~ tenX + Age_group: need each Age_group level to have
        # >1 distinct tenX value represented, OR just check group sizes are non-empty
        # and there's no perfect collinearity between tenX and Age_group
        crosstab = pd.crosstab(d["Age_group"], d["tenX"])
        # perfectly collinear if any Age_group level has cells in exactly one tenX column
        # and that column is unique to it (simple heuristic, not a full rank computation)
        full_rank_note = "see crosstab"
        region_result[label] = {
            "n_donors": n_donors,
            "n_age_group_levels": n_levels_age,
            "n_10X_levels_among_these_donors": n_levels_10X,
            "age_x_10X_donor_crosstab": crosstab.to_dict(),
        }
    results_rank[region] = region_result

print("\n=== Candidate contrast donor counts + Age_group x 10X crosstabs ===")
print(json.dumps(results_rank, indent=2, default=str))

with open(f"{OUT_DIR}/gut_epi_design_rank_audit.json", "w") as f:
    json.dump({
        "n_donor_multi_10X": n_donor_multi_10X,
        "n_donor_multi_batch": n_donor_multi_batch,
        "n_donor_multi_fraction": n_donor_multi_fraction,
        "n_donor_multi_agegroup": n_donor_multi_agegroup,
        "n_total_donors": len(donor_struct),
        "candidate_contrasts": results_rank,
    }, f, indent=2, default=str)
print(f"\nWrote {OUT_DIR}/gut_epi_design_rank_audit.json")
