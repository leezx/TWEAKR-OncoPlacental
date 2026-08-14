#!/usr/bin/env python3
"""
Pseudobulk-unit design audit, v2 -- fixes PR #20 round-3 REQUEST_CHANGES:
the round-2 audit correctly ruled out 10X-chemistry/cohort collinearity but
didn't check a deeper confound the reviewer caught by inspecting the real
batch IDs: within the LargeInt 3' subset, Second-trimester donors F66/F67
are both `Human_colon_16S...` batch family while Adult 3' donors A30/A33
are both `WTDAtest...` -- i.e. Age_group is PERFECTLY collinear with
source/batch family in that stratum, something `~ 10X + Age_group` cannot
separate from real fetal/adult biology.

This version:
  1. Derives `source_family` from `batch` by stripping the trailing numeric
     ID (regex), verified to cleanly partition all 149 real batch values
     into exactly 4 families: 4918STDY, FCA_gut, Human_colon_16S, WTDAtest.
  2. Builds the real Age_group x 10X x source_family donor-level table for
     the primary contrast (Second trim vs Adult) in LargeInt/SmallInt.
  3. Computes REAL design-matrix rank (via patsy dmatrix + np.linalg.matrix_rank),
     not a crosstab heuristic, for both `~ 10X + Age_group` and
     `~ 10X + source_family + Age_group` candidate models.
  4. Builds the mandatory 5'-only donor subset table (LargeInt Second-trim
     5' vs Adult 5' donors) -- the reviewer requires this reported
     unconditionally as a concordance/no-reversal check, not "run
     diagnostics then decide whether to bother."

Usage: python3 gut_epi_pseudobulk_design_audit_v2.py <out_dir>
"""
import sys
import os
import re
import json
import numpy as np
import pandas as pd
import anndata as ad

try:
    from patsy import dmatrix
    HAVE_PATSY = True
except ImportError:
    HAVE_PATSY = False

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "results/04a_dfp_gut/inventory"
os.makedirs(OUT_DIR, exist_ok=True)

PATH = "/home/zz950/DATA/scRNAseq/GutCellAtlas_Elmentaite2021/raw/epi_raw_counts02_v2.h5ad"
DONOR_COL = "Sample name"

print(f"Reading {PATH} (backed='r')...", flush=True)
adata = ad.read_h5ad(PATH, backed="r")
obs = adata.obs.copy()
obs = obs[obs["category"] == "Epithelial"]
obs = obs[obs["Region"].isin(["LargeInt", "SmallInt"])]

# --- source_family extraction, verified against the real batch value list ---
def source_family(batch):
    return re.sub(r"\d+$", "", str(batch))

obs["source_family"] = obs["batch"].apply(source_family)
families = sorted(obs["source_family"].unique())
print(f"\nsource_family values found: {families}", flush=True)
assert len(families) <= 6, "unexpectedly many source families -- regex rule may need revisiting"

# --- donor-level table (one row per donor, mode of each grouping variable) ---
# IMPORTANT bug fix (caught while sanity-checking this script's own rank output
# before reporting it): Age_group/10X are pandas `category` dtype with unused
# levels retained even after filtering (Age_group keeps all 7 original
# categories -- Adult/Adult_MLN/First trim/Pediatric/Pediatric_IBD/Second
# trim/Second trim_MLN -- even in a 2-group subset). patsy's dmatrix builds a
# dummy column per *declared* category, not per *observed* value, so an
# unconverted categorical column silently inflates the design matrix with
# structurally-empty (all-zero) columns and makes a perfectly-identified 2x2
# design look rank-deficient. Cast to plain str explicitly before any
# aggregation/patsy step so only observed levels ever reach the design matrix.
def donor_level_table(sub):
    sub = sub.copy()
    sub["Age_group"] = sub["Age_group"].astype(str)
    sub["10X"] = sub["10X"].astype(str)
    sub["source_family"] = sub["source_family"].astype(str)
    return sub.groupby(DONOR_COL, observed=True).agg(
        Age_group=("Age_group", lambda s: s.mode().iloc[0]),
        tenX=("10X", lambda s: s.mode().iloc[0]),
        source_family=("source_family", lambda s: s.mode().iloc[0]),
        n_source_family_within_donor=("source_family", "nunique"),
        n_cells=("Region", "size"),
    ).reset_index()

def real_rank(df, formula):
    if not HAVE_PATSY:
        return {"error": "patsy not available"}
    try:
        mat = dmatrix(formula, df, return_type="dataframe")
        rank = int(np.linalg.matrix_rank(mat.values))
        return {"formula": formula, "n_rows": len(df), "n_params": mat.shape[1], "rank": rank,
                "full_rank": rank == mat.shape[1]}
    except Exception as exc:
        return {"formula": formula, "error": str(exc)}

results = {}
for region in ["LargeInt", "SmallInt"]:
    sub = obs[obs["Region"] == region]
    donor_level = donor_level_table(sub)

    # primary contrast: Second trim vs Adult
    primary = donor_level[donor_level["Age_group"].isin(["Second trim", "Adult"])].copy()

    print(f"\n=== {region}: Second-trim-vs-Adult donor-level Age_group x 10X x source_family ===")
    ct = primary.groupby(["Age_group", "tenX", "source_family"], observed=True).size()
    print(ct)

    rank_no_family = real_rank(primary, "tenX + Age_group")
    rank_with_family = real_rank(primary, "tenX + source_family + Age_group")
    print(f"\n~ 10X + Age_group : {rank_no_family}")
    print(f"~ 10X + source_family + Age_group : {rank_with_family}")

    # mandatory 5'-only subset (reported unconditionally, per reviewer)
    fiveprime = primary[primary["tenX"] == "5'"]
    fp_counts = fiveprime.groupby("Age_group", observed=True)[DONOR_COL].nunique() if DONOR_COL in fiveprime.columns else fiveprime.groupby("Age_group", observed=True).size()
    print(f"\n5'-only donor counts (mandatory sensitivity subset): \n{fiveprime.groupby('Age_group', observed=True).size()}")

    results[region] = {
        "donor_level_table": primary.to_dict(orient="records"),
        "age_x_10X_x_source_family_donor_counts": {str(k): int(v) for k, v in ct.items()},
        "rank_10X_Age_group": rank_no_family,
        "rank_10X_source_family_Age_group": rank_with_family,
        "five_prime_only_subset": fiveprime[[DONOR_COL, "Age_group", "source_family", "n_cells"]].to_dict(orient="records"),
    }

    # explicit collinearity flag: does every (Age_group, 3') cell have exactly
    # one source_family, and is that family unique to that Age_group?
    threeprime = primary[primary["tenX"] == "3'"]
    tp_family_by_age = threeprime.groupby("Age_group", observed=True)["source_family"].apply(lambda s: sorted(s.unique()))
    print(f"\n3'-only source_family by Age_group (collinearity check): \n{tp_family_by_age}")
    families_by_age = {k: set(v) for k, v in tp_family_by_age.items()}
    if len(families_by_age) == 2:
        ages = list(families_by_age.keys())
        overlap = families_by_age[ages[0]] & families_by_age[ages[1]]
        results[region]["threeprime_source_family_overlap_between_age_groups"] = sorted(overlap)
        results[region]["threeprime_age_source_family_collinear"] = len(overlap) == 0
        print(f"3' source_family overlap between {ages}: {overlap} (collinear={len(overlap)==0})")

with open(f"{OUT_DIR}/gut_epi_design_rank_audit_v2.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nWrote {OUT_DIR}/gut_epi_design_rank_audit_v2.json")
