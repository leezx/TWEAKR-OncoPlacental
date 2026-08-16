#!/usr/bin/env python3
"""
Step 6 gut re-anchor: primary analysis (per docs/STEP6_GUT_SCORING_COMPUTE_DESIGN.md,
PR #26, APPROVE after 2 review rounds). Consumes crc_gut_scoring_full.py's
per-cell percentile output and computes, for each of the 10 locked
revCSC<->D/F/P comparison pairs:

1. Within-study Pearson/Spearman first (36 study_id values), pooled only
   after the donor-aware validation below.
2. Per-donor correlation table (donor_key = (study_id, patient_id)
   composite -- round-2 fix, not bare patient_id, since patient_id is not
   documented as globally unique across this 36-study meta-atlas) with
   n_cells; NOT_ESTIMABLE reported transparently (never filtered by
   outcome) for donors where n_cells is too small or either variable has
   zero variance within that donor.
3. Equal-donor-weighted study summaries: unweighted mean of per-donor
   correlations within each study, reported alongside the cell-weighted
   (pooled-within-study) value, explicitly labeled, never conflated.
4. Leave-one-donor-out and leave-one-study-out pooled sensitivity: the
   full-cohort pooled correlation recomputed once per left-out donor/
   study, reporting the range and which single donor/study shifts it
   most.

Usage: python3 crc_gut_scoring_primary_analysis.py <scoring_full_out_dir> <analysis_out_dir>
(run on Argos, argos-codex env, after crc_gut_scoring_full.py)
"""
import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crc_gut_scoring_core import COMPARISON_PAIRS

MIN_CELLS_FOR_ESTIMATE = 10


def safe_corr(x, y):
    """Returns (pearson_r, spearman_rho) or (None, None) with a reason if
    non-estimable -- never silently coerced to 0 or dropped."""
    n = len(x)
    if n < MIN_CELLS_FOR_ESTIMATE:
        return None, None, f"n_cells={n} < {MIN_CELLS_FOR_ESTIMATE}"
    if np.std(x) == 0 or np.std(y) == 0:
        return None, None, "zero variance in x or y within this group"
    r, _ = pearsonr(x, y)
    rho, _ = spearmanr(x, y)
    return float(r), float(rho), None


def per_donor_table(x, y, donor_key):
    rows = []
    df = pd.DataFrame({"x": x, "y": y, "donor_key": donor_key})
    for dk, sub in df.groupby("donor_key"):
        r, rho, reason = safe_corr(sub["x"].values, sub["y"].values)
        rows.append({
            "donor_key": dk, "n_cells": len(sub),
            "pearson_r": r, "spearman_rho": rho,
            "estimable": r is not None,
            "not_estimable_reason": reason,
        })
    return pd.DataFrame(rows)


def equal_donor_weighted_study_summary(donor_df, study_of_donor):
    donor_df = donor_df.copy()
    donor_df["study_id"] = donor_df["donor_key"].map(study_of_donor)
    rows = []
    for sid, sub in donor_df.groupby("study_id"):
        estimable = sub[sub.estimable]
        rows.append({
            "study_id": sid,
            "n_donors": len(sub),
            "n_donors_estimable": len(estimable),
            "equal_donor_weighted_pearson_r": float(estimable.pearson_r.mean()) if len(estimable) else None,
            "equal_donor_weighted_spearman_rho": float(estimable.spearman_rho.mean()) if len(estimable) else None,
        })
    return pd.DataFrame(rows)


def cell_weighted_within_study(x, y, study):
    rows = []
    df = pd.DataFrame({"x": x, "y": y, "study_id": study})
    for sid, sub in df.groupby("study_id"):
        r, rho, reason = safe_corr(sub["x"].values, sub["y"].values)
        rows.append({
            "study_id": sid, "n_cells": len(sub),
            "cell_weighted_pearson_r": r, "cell_weighted_spearman_rho": rho,
            "estimable": r is not None, "not_estimable_reason": reason,
        })
    return pd.DataFrame(rows)


def leave_one_out_sensitivity(x, y, group, group_name):
    """Full-cohort pooled correlation recomputed once per left-out group
    value. Returns a DataFrame plus the pooled (all-in) value."""
    r_all, rho_all, reason_all = safe_corr(x, y)
    groups = pd.unique(group)
    rows = []
    for g in groups:
        mask = group != g
        r, rho, reason = safe_corr(x[mask], y[mask])
        rows.append({
            f"left_out_{group_name}": g, "n_cells_remaining": int(mask.sum()),
            "pooled_pearson_r_excl": r, "pooled_spearman_rho_excl": rho,
            "delta_pearson_r_vs_full": (r - r_all) if (r is not None and r_all is not None) else None,
        })
    df = pd.DataFrame(rows)
    return df, r_all, rho_all, reason_all


def main():
    scoring_dir = sys.argv[1] if len(sys.argv) > 1 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/gut_scoring_full"
    out_dir = sys.argv[2] if len(sys.argv) > 2 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/gut_scoring_primary_analysis"
    os.makedirs(out_dir, exist_ok=True)

    scores = pd.read_parquet(f"{scoring_dir}/crc_gut_scoring_all_panels.parquet")
    meta = pd.read_parquet(f"{scoring_dir}/crc_gut_scoring_cell_metadata.parquet")
    assert (scores.index == meta.index).all(), "scores/metadata index mismatch"
    print(f"Loaded {len(scores)} cells' scores and metadata", flush=True)

    study_of_donor = meta[["donor_key", "study_id"]].drop_duplicates().set_index("donor_key")["study_id"].to_dict()

    overview_rows = []
    for revcsc_panel, dfp_panel in COMPARISON_PAIRS:
        pair_name = f"{revcsc_panel}__vs__{dfp_panel}"
        print(f"\n=== {pair_name} ===", flush=True)
        x = scores[f"{revcsc_panel}_percentile"].values
        y = scores[f"{dfp_panel}_percentile"].values

        # 1. cell-weighted within-study, pooled
        within_study_df = cell_weighted_within_study(x, y, meta["study_id"].values)
        within_study_path = f"{out_dir}/{pair_name}__within_study_cell_weighted.tsv"
        within_study_df.to_csv(within_study_path, sep="\t", index=False)

        # 2. per-donor table
        donor_df = per_donor_table(x, y, meta["donor_key"].values)
        donor_path = f"{out_dir}/{pair_name}__per_donor.tsv"
        donor_df.to_csv(donor_path, sep="\t", index=False)
        n_donors = len(donor_df)
        n_estimable = int(donor_df.estimable.sum())
        print(f"  {n_donors} donors, {n_estimable} estimable ({n_donors - n_estimable} NOT_ESTIMABLE)", flush=True)

        # 3. equal-donor-weighted study summary
        study_summary_df = equal_donor_weighted_study_summary(donor_df, study_of_donor)
        study_summary_merged = study_summary_df.merge(
            within_study_df[["study_id", "n_cells", "cell_weighted_pearson_r", "cell_weighted_spearman_rho"]],
            on="study_id", how="outer",
        )
        study_summary_path = f"{out_dir}/{pair_name}__study_summary_equal_vs_cell_weighted.tsv"
        study_summary_merged.to_csv(study_summary_path, sep="\t", index=False)

        # 4. leave-one-out sensitivity
        lodo_df, r_all, rho_all, reason_all = leave_one_out_sensitivity(x, y, meta["donor_key"].values, "donor_key")
        losout_df, _, _, _ = leave_one_out_sensitivity(x, y, meta["study_id"].values, "study_id")
        lodo_path = f"{out_dir}/{pair_name}__leave_one_donor_out.tsv"
        losout_path = f"{out_dir}/{pair_name}__leave_one_study_out.tsv"
        lodo_df.to_csv(lodo_path, sep="\t", index=False)
        losout_df.to_csv(losout_path, sep="\t", index=False)

        if r_all is None:
            robust = None
            robust_pearson_only = robust_spearman_only = None
            max_shift_donor = max_shift_study = None
        else:
            max_shift_donor_row = lodo_df.loc[lodo_df.delta_pearson_r_vs_full.abs().idxmax()] \
                if lodo_df.delta_pearson_r_vs_full.notna().any() else None
            max_shift_study_row = losout_df.loc[losout_df.delta_pearson_r_vs_full.abs().idxmax()] \
                if losout_df.delta_pearson_r_vs_full.notna().any() else None
            max_shift_donor = float(max_shift_donor_row.delta_pearson_r_vs_full) if max_shift_donor_row is not None else None
            max_shift_study = float(max_shift_study_row.delta_pearson_r_vs_full) if max_shift_study_row is not None else None
            # BUG FIX (found by reviewer, PR #27 round 2): robustness was
            # determined from Pearson sign stability only, even though both
            # Pearson and Spearman are reported as the primary/secondary
            # metrics throughout this compute. Direct check found a real
            # case (revCSC_extended28_minus_CLU_ASS1 vs F_Gut-specific)
            # where Pearson stays same-sign under leave-one-out but Spearman
            # flips sign when the same influential study (Terekhanova_2023_Nature)
            # is excluded -- a pair reported "robust" that is not actually
            # rank-stable. Fixed: robust now requires same-sign stability
            # for BOTH metrics, both leave-one-donor-out and
            # leave-one-study-out.
            same_sign_donor_r = all(
                (r_all >= 0) == (v >= 0) for v in lodo_df.pooled_pearson_r_excl.dropna()
            )
            same_sign_study_r = all(
                (r_all >= 0) == (v >= 0) for v in losout_df.pooled_pearson_r_excl.dropna()
            )
            same_sign_donor_rho = all(
                (rho_all >= 0) == (v >= 0) for v in lodo_df.pooled_spearman_rho_excl.dropna()
            )
            same_sign_study_rho = all(
                (rho_all >= 0) == (v >= 0) for v in losout_df.pooled_spearman_rho_excl.dropna()
            )
            robust_pearson_only = bool(same_sign_donor_r and same_sign_study_r)
            robust_spearman_only = bool(same_sign_donor_rho and same_sign_study_rho)
            robust = bool(robust_pearson_only and robust_spearman_only)

        # SCHEMA FIX (found by reviewer, PR #27 round 2): these columns
        # are named "max_abs_delta_..." but stored the SIGNED delta at
        # the row of largest absolute shift (so values can be negative),
        # contradicting the "abs" in the name. Renamed to make clear the
        # value is signed, picked by absolute magnitude.
        overview_rows.append({
            "pair": pair_name, "revcsc_panel": revcsc_panel, "dfp_panel": dfp_panel,
            "n_cells_total": len(x), "n_donors": n_donors, "n_donors_estimable": n_estimable,
            "pooled_pearson_r_all_cells": r_all, "pooled_spearman_rho_all_cells": rho_all,
            "pooled_not_estimable_reason": reason_all,
            "signed_delta_pearson_r_at_max_abs_shift_leave_one_donor_out": max_shift_donor,
            "signed_delta_pearson_r_at_max_abs_shift_leave_one_study_out": max_shift_study,
            # Round 2 fix: robust now requires Pearson AND Spearman
            # sign-stability (both leave-one-donor-out and
            # leave-one-study-out); the two component flags are also
            # reported explicitly so a reader can see which metric (if
            # either) drove a non-robust verdict, per the reviewer's
            # suggestion to expose per-metric robustness rather than
            # only a single combined boolean.
            "robust_pearson_sign_stable": robust_pearson_only,
            "robust_spearman_sign_stable": robust_spearman_only,
            "robust_to_no_single_donor_or_study": robust,
        })
        print(f"  pooled r={r_all}, rho={rho_all}, robust={robust}", flush=True)

    overview_df = pd.DataFrame(overview_rows)
    overview_path = f"{out_dir}/primary_analysis_overview.tsv"
    overview_df.to_csv(overview_path, sep="\t", index=False)
    print(f"\nWrote {overview_path}", flush=True)
    print(overview_df.to_string(index=False), flush=True)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
