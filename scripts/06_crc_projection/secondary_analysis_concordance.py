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
    """PR #29 round-2 fix (real blocker, confirmed directly): the locked
    design (docs/STEP6_SECONDARY_ANALYSIS_DESIGN.md sec 5) says a donor's
    2x2 table is non-estimable "(a zero cell)" -- literally any of a/b/c/d
    == 0, not just b==0 or c==0 (the round-1 fix's definition, which only
    guards against the individual OR being mathematically undefined via
    division by zero). Confirmed real on real data: donor
    Pelka_2021_Cell.C150 (a=0, b=2, c=2, d=1396) was marked estimable=True
    and kept in MH pooling under the round-1 definition, despite having a
    zero cell -- exactly the literal-noncompliance case the reviewer
    found. There is a genuine statistical argument that a=0/d=0 strata
    (well-defined boundary OR, not undefined) are fine to keep in MH
    pooling -- but this PR explicitly claims "no design changes," so the
    locked text's literal "a zero cell" wording is implemented as written:
    any of a,b,c,d == 0 makes the donor non-estimable and excludes it from
    MH pooling."""
    rows = []
    for donor, t in tables_by_donor.items():
        a, b, c, d = t[0, 0], t[0, 1], t[1, 0], t[1, 1]
        estimable = a > 0 and b > 0 and c > 0 and d > 0
        odds_ratio = (a * d) / (b * c) if (b > 0 and c > 0) else None
        rows.append({
            "donor_key": donor, "n_cells": int(t.sum()),
            "a_m11hi_revcschi": int(a), "b_m11hi_revcsclo": int(b),
            "c_m11lo_revcschi": int(c), "d_m11lo_revcsclo": int(d),
            "estimable": estimable,
            "odds_ratio": odds_ratio,
            "not_estimable_reason": None if estimable else "at least one zero cell (a, b, c, or d == 0)",
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

    # PR #29 round-1 fix (real blocker, confirmed by re-reading the locked
    # design): docs/STEP6_SECONDARY_ANALYSIS_DESIGN.md sec 5 explicitly
    # requires non-estimable (zero-cell) donors to be "excluded from the
    # MH pooling with the exclusion count stated." The original
    # implementation computed the estimable flag for per-donor reporting
    # but then passed ALL donor tables (estimable or not) into
    # mh_common_or/LODO/LOSO/bootstrap -- a real deviation from the
    # approved design, confirmed by the reviewer directly re-reading both
    # the design text and this code. Whether unconditionally retaining
    # zero-cell strata would be more statistically standard is a genuinely
    # separate question (checked directly: it is -- MH's formula handles
    # them without bias, and 11-13 of the ~26/13/6 "non-estimable" donors
    # per cutoff have b=0 XOR c=0, i.e. a real informative margin, not a
    # fully degenerate one; excluding all of them shifts OR_MH by only
    # ~0.3 at top5pct, ~0.08 at top10pct, ~0 at top20pct -- small, no
    # qualitative change). But a compute PR explicitly claiming "no design
    # changes" must execute the design AS APPROVED, not a different
    # (even if arguably better) estimator -- so this is fixed to comply
    # literally: only estimable donors' tables are used for the pooled OR,
    # its CIs, and the leave-one-out/bootstrap sensitivities below. The
    # full per-donor table (all donors, estimable flag included) is still
    # written above for transparency, unchanged.
    tables_estimable = {d: t for d, t in tables_by_donor.items()
                         if donor_or_df.set_index("donor_key").loc[d, "estimable"]}
    n_excluded = len(tables_by_donor) - len(tables_estimable)
    print(f"  MH pooling restricted to {len(tables_estimable)} estimable donors "
          f"({n_excluded} excluded per the locked design's zero-cell exclusion rule)", flush=True)

    or_mh, ci_lo, ci_hi = mh_common_or(tables_estimable)
    print(f"  OR_MH={or_mh}, 95% CI=({ci_lo}, {ci_hi})", flush=True)

    donor_key_to_donor = {d: d for d in tables_estimable}
    lodo_df, _ = leave_one_group_out_common_or(tables_estimable, donor_key_to_donor)
    lodo_path = f"{out_dir}/enrichment_{cutoff_label}_leave_one_donor_out.tsv"
    lodo_df.to_csv(lodo_path, sep="\t", index=False)

    donor_to_study = pd.DataFrame({"donor_key": donor_key, "study_id": study_id}).drop_duplicates().set_index("donor_key")["study_id"].to_dict()
    donor_key_to_study = {d: donor_to_study.get(d) for d in tables_estimable}
    losout_df, _ = leave_one_group_out_common_or(tables_estimable, donor_key_to_study)
    losout_path = f"{out_dir}/enrichment_{cutoff_label}_leave_one_study_out.tsv"
    losout_df.to_csv(losout_path, sep="\t", index=False)

    boot_lo, boot_hi, n_valid_boot = donor_cluster_bootstrap_ci(tables_estimable)
    print(f"  donor-cluster bootstrap 95% CI=({boot_lo}, {boot_hi}) from {n_valid_boot}/{BOOTSTRAP_N} valid resamples", flush=True)

    # PR #29 round-1 fix: the results doc mislabeled which single study
    # produces the largest |delta| at each cutoff (assumed it was always
    # the same study as the one giving the lowest or_mh_excl, which is a
    # different quantity) -- fixed by explicitly identifying and reporting
    # both the max-|delta| study/donor AND the min-or_mh_excl study/donor
    # here, in the machine-generated overview, so the write-up step reads
    # off the correct value instead of re-deriving it by hand.
    max_shift_donor_row = lodo_df.loc[lodo_df.delta_or_vs_full.abs().idxmax()] if lodo_df.delta_or_vs_full.notna().any() else None
    max_shift_study_row = losout_df.loc[losout_df.delta_or_vs_full.abs().idxmax()] if losout_df.delta_or_vs_full.notna().any() else None
    min_or_study_row = losout_df.loc[losout_df.or_mh_excl.idxmin()] if losout_df.or_mh_excl.notna().any() else None

    return {
        "cutoff": cutoff_label, "n_m11_high": n_m11_high, "n_revcsc_high": n_revcsc_high,
        "n_both": both, "n_donors": len(tables_by_donor), "n_donors_estimable": n_estimable,
        "n_donors_excluded_from_MH": n_excluded,
        "or_mh": or_mh, "or_mh_ci_lo_asymptotic": ci_lo, "or_mh_ci_hi_asymptotic": ci_hi,
        "or_mh_ci_lo_donor_cluster_bootstrap": boot_lo, "or_mh_ci_hi_donor_cluster_bootstrap": boot_hi,
        "n_valid_bootstrap_resamples": n_valid_boot,
        # PR #29 round-2 fix (real bug, reappearance of the exact PR #27
        # round-2 naming pattern): these columns are picked by MAXIMUM
        # ABSOLUTE shift but store the SIGNED delta at that row (e.g. a
        # real committed value of -1.762 under a "max_abs_delta_..." name
        # at top10pct) -- contradicting the "abs" in the name, same bug
        # class PR #27 already fixed once. Renamed to make clear the
        # stored value is signed, selected by absolute magnitude.
        "signed_delta_or_at_max_abs_shift_leave_one_donor_out": max_shift_donor_row.delta_or_vs_full if max_shift_donor_row is not None else None,
        "signed_delta_or_at_max_abs_shift_leave_one_donor_out_which": max_shift_donor_row.left_out_group if max_shift_donor_row is not None else None,
        "signed_delta_or_at_max_abs_shift_leave_one_study_out": max_shift_study_row.delta_or_vs_full if max_shift_study_row is not None else None,
        "signed_delta_or_at_max_abs_shift_leave_one_study_out_which": max_shift_study_row.left_out_group if max_shift_study_row is not None else None,
        "min_or_mh_excl_leave_one_study_out": min_or_study_row.or_mh_excl if min_or_study_row is not None else None,
        "min_or_mh_excl_leave_one_study_out_which": min_or_study_row.left_out_group if min_or_study_row is not None else None,
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

    # ---- 1b. Same-population sensitivity: D/F/P<->revCSC correlations,
    # restricted to the SAME 297,307-cell M11 subset (PR #29 round-1 fix,
    # per review) -- PR #27's D/F/P<->revCSC correlations were computed
    # on the full 665,473-cell atlas, so "M11's r=0.318 is stronger than
    # any D/F/P pair's |r|<=0.19" was a cross-population comparison, not
    # a clean same-population effect-size comparison. No new scoring
    # needed (D/F/P panels already scored for all cells in PR #27) --
    # just restricts the existing full-atlas D/F/P percentiles to the
    # M11 subset and reuses the same correlation machinery.
    print("\n=== Same-population sensitivity: D/F/P vs revCSC, M11 subset only ===", flush=True)
    dfp_panels_for_sensitivity = ["D_Gut-shared", "F_Gut-specific", "F_Colon-specific", "F_SI-specific", "P_Gut-specific"]
    dfp_sensitivity_rows = []
    for dfp_panel in dfp_panels_for_sensitivity:
        y_dfp = full_scores.loc[m11_scores.index, f"{dfp_panel}_percentile"].values
        r_all, rho_all, reason_all = safe_corr(y, y_dfp)
        print(f"  revCSC_primary27_minus_CLU_ASS1 vs {dfp_panel} (M11-subset only): "
              f"r={r_all}, rho={rho_all}", flush=True)
        dfp_sensitivity_rows.append({
            "comparison": f"revCSC_primary27_minus_CLU_ASS1__vs__{dfp_panel}__M11_SUBSET_ONLY",
            "n_cells": len(y), "pooled_pearson_r": r_all, "pooled_spearman_rho": rho_all,
        })
    dfp_sensitivity_df = pd.DataFrame(dfp_sensitivity_rows)
    dfp_sensitivity_path = f"{out_dir}/concordance_dfp_vs_revcsc_M11_subset_sensitivity.tsv"
    dfp_sensitivity_df.to_csv(dfp_sensitivity_path, sep="\t", index=False)
    print(f"Wrote {dfp_sensitivity_path}", flush=True)
    print(dfp_sensitivity_df.to_string(index=False), flush=True)

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
