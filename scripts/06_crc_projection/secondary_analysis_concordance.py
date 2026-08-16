#!/usr/bin/env python3
"""
Step 6 secondary analysis compute (docs/STEP6_SECONDARY_ANALYSIS_DESIGN.md,
PR #28, APPROVE after 4 review rounds), section 5: M11 <-> revCSC
concordance, within the 297,307-cell M11 subset only.

1. Continuous correlation: M11_top50_minus_revCSC_overlap percentile
   (full M11 as sensitivity) vs revCSC_primary27_minus_CLU_ASS1
   percentile -- reuses crc_gut_scoring_primary_analysis.py's
   donor/study-aware machinery directly (same functions, same
   leave-one-donor/study-out sensitivity, same Pearson+Spearman
   sign-stability robustness definition).

2. Enrichment test: M11-high (cross-cell z-score rank of
   M11_top50_minus_revCSC_overlap, ranked WITHIN the 297,307-cell M11
   subset) x revCSC-high (from secondary_analysis_composition.py's
   GLOBAL 665,473-cell z-score rank of revCSC_primary27_minus_CLU_ASS1,
   restricted to the M11 subset -- a fixed, already-computed membership
   set, not re-ranked). Matched cutoff pairs only: 10%x10% primary,
   5%x5%/20%x20% sensitivities -- never one side fixed while the other
   varies (round-2 review fix).

   Effect estimate: Mantel-Haenszel common odds ratio (donor-stratified
   2x2 tables, statsmodels.stats.contingency_tables.StratifiedTable).
   MH is NOT described as solving cell-level pseudoreplication (round-2
   review correction) -- reported alongside: per-donor 2x2/OR table with
   non-estimable donors (zero row/col margin) explicitly flagged, not
   silently dropped from the MH pooling; leave-one-donor-out and
   leave-one-study-out common OR; a donor-cluster bootstrap CI
   (resample donors with replacement, rebuild the stratified table,
   recompute the common OR, 1000 resamples).

Usage: python3 secondary_analysis_concordance.py <m11_scoring_dir> \
    <scoring_full_dir> <composition_dir> <out_dir>
(run on Argos)
"""
import sys
import os
import numpy as np
import pandas as pd
from statsmodels.stats.contingency_tables import StratifiedTable

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crc_gut_scoring_primary_analysis import (
    safe_corr, per_donor_table, equal_donor_weighted_study_summary,
    cell_weighted_within_study, leave_one_out_sensitivity,
)

BOOTSTRAP_N = 1000
BOOTSTRAP_SEED = 20260815
CUTOFF_LABELS = ["top5pct", "top10pct", "top20pct"]
CUTOFF_RANK_THRESH = {"top5pct": 0.95, "top10pct": 0.90, "top20pct": 0.80}


def build_2x2_tables(m11_high, revcsc_high, group):
    """One 2x2 table per group value: [[a,b],[c,d]] with
    a=m11_high&revcsc_high, b=m11_high&~revcsc_high,
    c=~m11_high&revcsc_high, d=~m11_high&~revcsc_high."""
    df = pd.DataFrame({"m11_high": m11_high, "revcsc_high": revcsc_high, "group": group})
    tables = {}
    for g, sub in df.groupby("group"):
        a = int(((sub.m11_high) & (sub.revcsc_high)).sum())
        b = int(((sub.m11_high) & (~sub.revcsc_high)).sum())
        c = int(((~sub.m11_high) & (sub.revcsc_high)).sum())
        d = int(((~sub.m11_high) & (~sub.revcsc_high)).sum())
        tables[g] = np.array([[a, b], [c, d]], dtype=float)
    return tables


def per_donor_or_table(tables_by_donor):
    rows = []
    for donor, t in tables_by_donor.items():
        a, b, c, d = t[0, 0], t[0, 1], t[1, 0], t[1, 1]
        row_ok = (a + b) > 0 and (c + d) > 0
        col_ok = (a + c) > 0 and (b + d) > 0
        estimable = row_ok and col_ok and b > 0 and c > 0
        odds_ratio = (a * d) / (b * c) if (b > 0 and c > 0) else None
        rows.append({
            "donor_key": donor, "n_cells": int(t.sum()),
            "a_m11hi_revcschi": int(a), "b_m11hi_revcsclo": int(b),
            "c_m11lo_revcschi": int(c), "d_m11lo_revcsclo": int(d),
            "estimable": estimable,
            "odds_ratio": odds_ratio,
            "not_estimable_reason": None if estimable else "zero margin or zero off-diagonal cell",
        })
    return pd.DataFrame(rows)


def mh_common_or(tables_by_key):
    tables = list(tables_by_key.values())
    if len(tables) < 2:
        return None, None, None
    st = StratifiedTable(tables)
    or_mh = float(st.oddsratio_pooled)
    ci_lo, ci_hi = st.oddsratio_pooled_confint()
    return or_mh, float(ci_lo), float(ci_hi)


def leave_one_group_out_common_or(tables_by_key, group_of_key):
    """tables_by_key keyed by donor_key; group_of_key maps donor_key ->
    donor_key (leave-one-donor-out) or donor_key -> study_id
    (leave-one-study-out)."""
    or_all, _, _ = mh_common_or(tables_by_key)
    groups = sorted(set(group_of_key.values()))
    rows = []
    for g in groups:
        remaining = {k: t for k, t in tables_by_key.items() if group_of_key[k] != g}
        or_excl, ci_lo, ci_hi = mh_common_or(remaining)
        rows.append({
            "left_out_group": g, "n_donors_remaining": len(remaining),
            "or_mh_excl": or_excl, "delta_or_vs_full": (or_excl - or_all) if (or_excl is not None and or_all is not None) else None,
        })
    return pd.DataFrame(rows), or_all


def donor_cluster_bootstrap_ci(tables_by_donor, n_boot=BOOTSTRAP_N, seed=BOOTSTRAP_SEED):
    donors = list(tables_by_donor.keys())
    rng = np.random.default_rng(seed)
    boot_ors = []
    for _ in range(n_boot):
        sampled = rng.choice(donors, size=len(donors), replace=True)
        tables = [tables_by_donor[d] for d in sampled]
        try:
            st = StratifiedTable(tables)
            or_b = float(st.oddsratio_pooled)
            if np.isfinite(or_b) and or_b > 0:
                boot_ors.append(or_b)
        except Exception:
            continue
    if len(boot_ors) < n_boot // 2:
        return None, None, len(boot_ors)
    lo, hi = np.percentile(boot_ors, [2.5, 97.5])
    return float(lo), float(hi), len(boot_ors)


def run_enrichment_for_cutoff(m11_high, revcsc_high, donor_key, study_id, cutoff_label, out_dir):
    print(f"\n--- Enrichment {cutoff_label} x {cutoff_label} (matched pair) ---", flush=True)
    n_m11_high = int(m11_high.sum())
    n_revcsc_high = int(revcsc_high.sum())
    both = int((m11_high & revcsc_high).sum())
    print(f"  n_m11_high={n_m11_high}, n_revcsc_high={n_revcsc_high}, n_both={both}, "
          f"n_cells={len(m11_high)}", flush=True)

    tables_by_donor = build_2x2_tables(m11_high, revcsc_high, donor_key)
    donor_or_df = per_donor_or_table(tables_by_donor)
    donor_or_path = f"{out_dir}/enrichment_{cutoff_label}_per_donor.tsv"
    donor_or_df.to_csv(donor_or_path, sep="\t", index=False)
    n_estimable = int(donor_or_df.estimable.sum())
    print(f"  {len(donor_or_df)} donors, {n_estimable} estimable ({len(donor_or_df)-n_estimable} NOT_ESTIMABLE)", flush=True)

    or_mh, ci_lo, ci_hi = mh_common_or(tables_by_donor)
    print(f"  OR_MH={or_mh}, 95% CI=({ci_lo}, {ci_hi})", flush=True)

    donor_key_to_donor = {d: d for d in tables_by_donor}
    lodo_df, _ = leave_one_group_out_common_or(tables_by_donor, donor_key_to_donor)
    lodo_path = f"{out_dir}/enrichment_{cutoff_label}_leave_one_donor_out.tsv"
    lodo_df.to_csv(lodo_path, sep="\t", index=False)

    donor_to_study = pd.DataFrame({"donor_key": donor_key, "study_id": study_id}).drop_duplicates().set_index("donor_key")["study_id"].to_dict()
    donor_key_to_study = {d: donor_to_study.get(d) for d in tables_by_donor}
    losout_df, _ = leave_one_group_out_common_or(tables_by_donor, donor_key_to_study)
    losout_path = f"{out_dir}/enrichment_{cutoff_label}_leave_one_study_out.tsv"
    losout_df.to_csv(losout_path, sep="\t", index=False)

    boot_lo, boot_hi, n_valid_boot = donor_cluster_bootstrap_ci(tables_by_donor)
    print(f"  donor-cluster bootstrap 95% CI=({boot_lo}, {boot_hi}) from {n_valid_boot}/{BOOTSTRAP_N} valid resamples", flush=True)

    max_shift_donor = lodo_df.delta_or_vs_full.abs().max() if lodo_df.delta_or_vs_full.notna().any() else None
    max_shift_study = losout_df.delta_or_vs_full.abs().max() if losout_df.delta_or_vs_full.notna().any() else None

    return {
        "cutoff": cutoff_label, "n_m11_high": n_m11_high, "n_revcsc_high": n_revcsc_high,
        "n_both": both, "n_donors": len(tables_by_donor), "n_donors_estimable": n_estimable,
        "or_mh": or_mh, "or_mh_ci_lo_asymptotic": ci_lo, "or_mh_ci_hi_asymptotic": ci_hi,
        "or_mh_ci_lo_donor_cluster_bootstrap": boot_lo, "or_mh_ci_hi_donor_cluster_bootstrap": boot_hi,
        "n_valid_bootstrap_resamples": n_valid_boot,
        "max_abs_delta_or_leave_one_donor_out": max_shift_donor,
        "max_abs_delta_or_leave_one_study_out": max_shift_study,
    }


def main():
    m11_dir = sys.argv[1] if len(sys.argv) > 1 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/m11_scoring_full"
    scoring_full_dir = sys.argv[2] if len(sys.argv) > 2 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/gut_scoring_full"
    composition_dir = sys.argv[3] if len(sys.argv) > 3 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/secondary_analysis_composition"
    out_dir = sys.argv[4] if len(sys.argv) > 4 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/secondary_analysis_concordance"
    os.makedirs(out_dir, exist_ok=True)

    m11_scores = pd.read_parquet(f"{m11_dir}/m11_scores.parquet")
    m11_meta = pd.read_parquet(f"{m11_dir}/m11_cell_metadata.parquet")
    assert (m11_scores.index == m11_meta.index).all(), "M11 scores/metadata index mismatch"
    print(f"Loaded {len(m11_scores)} M11-subset cells' scores and metadata", flush=True)

    full_scores = pd.read_parquet(f"{scoring_full_dir}/crc_gut_scoring_all_panels.parquet")
    revcsc_cols = ["revCSC_primary27_minus_CLU_ASS1_percentile", "revCSC_primary27_minus_CLU_ASS1_zscore"]
    revcsc_sub = full_scores.loc[m11_scores.index, revcsc_cols]
    print(f"Restricted full-atlas revCSC scores to M11 subset: {len(revcsc_sub)} cells "
          f"(must equal {len(m11_scores)})", flush=True)
    assert len(revcsc_sub) == len(m11_scores), "revCSC subset row count mismatch"

    membership = pd.read_parquet(f"{composition_dir}/revcsc_high_cohort_membership.parquet")
    revcsc_high_sub = membership.loc[m11_scores.index]

    # ---- 1. Continuous correlation ----
    print("\n=== Continuous correlation: M11 vs revCSC (within M11 subset) ===", flush=True)
    x_primary = m11_scores["M11_top50_minus_revCSC_overlap_percentile"].values
    x_sensitivity = m11_scores["M11_top50_full_percentile"].values
    y = revcsc_sub["revCSC_primary27_minus_CLU_ASS1_percentile"].values
    donor_key = m11_meta["donor_key"].values
    study_id = m11_meta["study_id"].values

    corr_overview = []
    for label, x in [("M11_top50_minus_revCSC_overlap_PRIMARY", x_primary),
                      ("M11_top50_full_SENSITIVITY", x_sensitivity)]:
        print(f"\n--- {label} vs revCSC_primary27_minus_CLU_ASS1 ---", flush=True)
        within_study_df = cell_weighted_within_study(x, y, study_id)
        within_study_df.to_csv(f"{out_dir}/concordance_{label}__within_study_cell_weighted.tsv", sep="\t", index=False)

        donor_df = per_donor_table(x, y, donor_key)
        donor_df.to_csv(f"{out_dir}/concordance_{label}__per_donor.tsv", sep="\t", index=False)
        n_donors, n_estimable = len(donor_df), int(donor_df.estimable.sum())

        study_of_donor = pd.DataFrame({"donor_key": donor_key, "study_id": study_id}).drop_duplicates().set_index("donor_key")["study_id"].to_dict()
        study_summary_df = equal_donor_weighted_study_summary(donor_df, study_of_donor)
        study_summary_df.to_csv(f"{out_dir}/concordance_{label}__study_summary.tsv", sep="\t", index=False)

        lodo_df, r_all, rho_all, reason_all = leave_one_out_sensitivity(x, y, donor_key, "donor_key")
        losout_df, _, _, _ = leave_one_out_sensitivity(x, y, study_id, "study_id")
        lodo_df.to_csv(f"{out_dir}/concordance_{label}__leave_one_donor_out.tsv", sep="\t", index=False)
        losout_df.to_csv(f"{out_dir}/concordance_{label}__leave_one_study_out.tsv", sep="\t", index=False)

        if r_all is None:
            robust = robust_p = robust_s = None
        else:
            same_sign_donor_r = all((r_all >= 0) == (v >= 0) for v in lodo_df.pooled_pearson_r_excl.dropna())
            same_sign_study_r = all((r_all >= 0) == (v >= 0) for v in losout_df.pooled_pearson_r_excl.dropna())
            same_sign_donor_rho = all((rho_all >= 0) == (v >= 0) for v in lodo_df.pooled_spearman_rho_excl.dropna())
            same_sign_study_rho = all((rho_all >= 0) == (v >= 0) for v in losout_df.pooled_spearman_rho_excl.dropna())
            robust_p = bool(same_sign_donor_r and same_sign_study_r)
            robust_s = bool(same_sign_donor_rho and same_sign_study_rho)
            robust = bool(robust_p and robust_s)

        print(f"  pooled r={r_all}, rho={rho_all}, n_donors={n_donors} ({n_estimable} estimable), robust={robust}", flush=True)
        corr_overview.append({
            "comparison": label, "n_cells": len(x), "n_donors": n_donors, "n_donors_estimable": n_estimable,
            "pooled_pearson_r": r_all, "pooled_spearman_rho": rho_all,
            "robust_pearson_sign_stable": robust_p, "robust_spearman_sign_stable": robust_s,
            "robust_to_no_single_donor_or_study": robust,
        })
    corr_overview_df = pd.DataFrame(corr_overview)
    corr_overview_path = f"{out_dir}/concordance_correlation_overview.tsv"
    corr_overview_df.to_csv(corr_overview_path, sep="\t", index=False)
    print(f"\nWrote {corr_overview_path}", flush=True)
    print(corr_overview_df.to_string(index=False), flush=True)

    # ---- 2. Enrichment test (matched cutoff pairs) ----
    print("\n=== Enrichment test: M11-high x revCSC-high (matched cutoff pairs) ===", flush=True)
    m11_rank = m11_scores["M11_top50_minus_revCSC_overlap_zscore"].rank(pct=True, method="average")

    enrichment_rows = []
    for cutoff_label in CUTOFF_LABELS:
        thresh = CUTOFF_RANK_THRESH[cutoff_label]
        m11_high = (m11_rank >= thresh).values
        revcsc_high_col = f"revCSC_primary27_minus_CLU_ASS1_{cutoff_label}"
        revcsc_high = revcsc_high_sub[revcsc_high_col].values
        row = run_enrichment_for_cutoff(m11_high, revcsc_high, donor_key, study_id, cutoff_label, out_dir)
        enrichment_rows.append(row)

    enrichment_df = pd.DataFrame(enrichment_rows)
    enrichment_path = f"{out_dir}/concordance_enrichment_overview.tsv"
    enrichment_df.to_csv(enrichment_path, sep="\t", index=False)
    print(f"\nWrote {enrichment_path}", flush=True)
    print(enrichment_df.to_string(index=False), flush=True)

    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
