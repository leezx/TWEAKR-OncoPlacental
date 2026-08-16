#!/usr/bin/env python3
"""
Step 6 secondary analysis compute (docs/STEP6_SECONDARY_ANALYSIS_DESIGN.md,
PR #28, APPROVE after 4 review rounds), sections 1-2: revCSC-high cohort
definition + developmental composition. Operates entirely on PR #27's
already-scored, already-verified 665,473-cell parquet -- no new scoring.

Section 1 -- revCSC-high threshold: cross-cell rank (pandas.rank(pct=True))
of revCSC_primary27_minus_CLU_ASS1_zscore across all 665,473 cells. Primary
cutoff = top decile (>=90th percentile of the RANK, i.e. exactly the top
10% of cells by construction); 5%/20% cutoffs reported alongside as
sensitivities, never substituted for the primary. The extended28
equivalent (revCSC_extended28_minus_CLU_ASS1_zscore) is also ranked and
reported as a cohort-definition sensitivity, disclosed at its available
N_PERM=100 precision (vs. the primary cohort's 500).

Section 2 -- developmental composition, computed within each revCSC-high
cohort:
  Step 0 (primary result): per-axis "supported" = null-calibrated
    percentile >=90 (matched-null evidence, NOT a cross-cell rank) for
    each of {D_Gut-shared, F_Gut-specific, P_Gut-specific} -> 8 mutually
    exclusive categories (none/D/F/P/D+F/D+P/F+P/D+F+P).
  Step A (descriptive summary only, explicitly NOT substate evidence):
    3-way argmax among the same 3 percentiles.
  Step A' (predominant axis among Step-0-supported cells only): argmax
    restricted to supported cells, requiring a >=5-percentile-point lead
    over the 2nd-highest axis; else "no_clear_predominance".
  Step B (regional refinement within F-assigned cells: F-supported OR
    F-predominant): F_Colon-specific vs F_SI-specific percentile, same
    Delta5 margin rule; else "no_clear_regional_skew".
All 4 steps reported pooled, and donor/study-aware (unweighted mean
across donor_key / study_id), per cohort cutoff.

Usage: python3 secondary_analysis_composition.py <scoring_full_dir> <out_dir>
(run on Argos; scoring_full_dir = PR #27's gut_scoring_full output dir)
"""
import sys
import os
import numpy as np
import pandas as pd

MARGIN = 5.0  # Delta5 percentile-point margin, locked in the approved design
AXES = ["D_Gut-shared", "F_Gut-specific", "P_Gut-specific"]
SUPPORT_PCT_THRESHOLD = 90.0


def rank_cutoff_membership(zscore, cutoffs=(0.95, 0.90, 0.80)):
    """cross-cell rank(pct=True) of zscore -> {cutoff_label: boolean mask}."""
    rank = pd.Series(zscore).rank(pct=True, method="average")
    out = {}
    for cutoff in cutoffs:
        label = {0.95: "top5pct", 0.90: "top10pct", 0.80: "top20pct"}[cutoff]
        mask = (rank >= cutoff).values
        out[label] = mask
        print(f"    {label} (rank>={cutoff}): n={mask.sum()} "
              f"({mask.sum()/len(mask)*100:.2f}%)", flush=True)
    return out


def axis_supported_status(pct_d, pct_f, pct_p, threshold=SUPPORT_PCT_THRESHOLD):
    d_sup = pct_d >= threshold
    f_sup = pct_f >= threshold
    p_sup = pct_p >= threshold
    n_sup = d_sup.astype(int) + f_sup.astype(int) + p_sup.astype(int)
    status = np.full(len(pct_d), "none", dtype=object)
    status = np.where(n_sup == 0, "none", status)
    status = np.where((n_sup == 1) & d_sup, "D_only", status)
    status = np.where((n_sup == 1) & f_sup, "F_only", status)
    status = np.where((n_sup == 1) & p_sup, "P_only", status)
    status = np.where((n_sup == 2) & d_sup & f_sup, "D+F", status)
    status = np.where((n_sup == 2) & d_sup & p_sup, "D+P", status)
    status = np.where((n_sup == 2) & f_sup & p_sup, "F+P", status)
    status = np.where(n_sup == 3, "D+F+P", status)
    return status, d_sup, f_sup, p_sup


def coarse_argmax(pct_d, pct_f, pct_p):
    stacked = np.vstack([pct_d, pct_f, pct_p])
    idx = np.argmax(stacked, axis=0)
    labels = np.array(["D", "F", "P"])
    return labels[idx]


def predominant_among_supported(pct_d, pct_f, pct_p, d_sup, f_sup, p_sup, margin=MARGIN):
    """Predominant axis among Step-0-supported cells only, requiring a
    >=margin percentile-point lead over the 2nd-highest of the 3 axes.
    Cells with no supported axis get 'not_applicable_none_supported'."""
    stacked = np.vstack([pct_d, pct_f, pct_p])  # 3 x n
    sorted_vals = np.sort(stacked, axis=0)  # ascending
    top = sorted_vals[2]
    second = sorted_vals[1]
    lead = top - second
    argmax_labels = coarse_argmax(pct_d, pct_f, pct_p)

    any_supported = d_sup | f_sup | p_sup
    out = np.full(len(pct_d), "not_applicable_none_supported", dtype=object)
    has_margin = lead >= margin
    out = np.where(any_supported & has_margin, argmax_labels, out)
    out = np.where(any_supported & ~has_margin, "no_clear_predominance", out)
    return out


def regional_refinement(pct_colon, pct_si, f_assigned_mask, margin=MARGIN):
    """Within f_assigned_mask cells only: Colon-biased / SI-biased / no
    clear regional skew (Delta5 margin), else not_applicable."""
    out = np.full(len(pct_colon), "not_applicable_not_F_assigned", dtype=object)
    diff = pct_colon - pct_si
    colon_biased = f_assigned_mask & (diff >= margin)
    si_biased = f_assigned_mask & (diff <= -margin)
    no_skew = f_assigned_mask & (np.abs(diff) < margin)
    out = np.where(colon_biased, "Colon_biased", out)
    out = np.where(si_biased, "SI_biased", out)
    out = np.where(no_skew, "no_clear_regional_skew", out)
    return out


def step0_x_stepA_crosstab(status, argmax_labels, cohort_label):
    """PR #29 round-1 fix (real blocker, confirmed by re-reading the
    approved design): PR #28 explicitly requires the coarse argmax to be
    'restricted to (or at minimum cross-tabulated against) axis-supported
    status' -- the whole point of demoting argmax to descriptive-only is
    to expose exactly how much of its apparent F share comes from cells
    with NO supported axis at all. Reporting Step 0 and Step A as
    separate marginals (the original implementation) does not do this.
    Fixed: an explicit Step0 x StepA cross-tab, pooled counts + fractions
    of the cohort."""
    df = pd.DataFrame({"step0_status": status, "stepA_argmax": argmax_labels})
    ct = pd.crosstab(df["step0_status"], df["stepA_argmax"])
    rows = []
    n_total = len(df)
    for status_cat in ct.index:
        for argmax_cat in ct.columns:
            n = int(ct.loc[status_cat, argmax_cat])
            rows.append({
                "cohort": cohort_label, "step0_status": status_cat, "stepA_argmax": argmax_cat,
                "n_cells": n, "frac_of_cohort": n / n_total,
            })
    return pd.DataFrame(rows)


def categorical_donor_study_summary(labels, donor_key, study_id, label_name):
    """Pooled + unweighted-mean-across-donor + unweighted-mean-across-study
    proportion breakdown for one categorical result."""
    df = pd.DataFrame({"label": labels, "donor_key": donor_key, "study_id": study_id})
    categories = sorted(pd.unique(df["label"]))

    pooled = df["label"].value_counts(normalize=True).reindex(categories, fill_value=0.0)
    pooled_n = df["label"].value_counts().reindex(categories, fill_value=0)

    donor_props = df.groupby("donor_key")["label"].value_counts(normalize=True).unstack(fill_value=0.0)
    donor_props = donor_props.reindex(columns=categories, fill_value=0.0)
    donor_mean = donor_props.mean(axis=0)
    donor_sd = donor_props.std(axis=0, ddof=1)

    study_props = df.groupby("study_id")["label"].value_counts(normalize=True).unstack(fill_value=0.0)
    study_props = study_props.reindex(columns=categories, fill_value=0.0)
    study_mean = study_props.mean(axis=0)
    study_sd = study_props.std(axis=0, ddof=1)

    rows = []
    for c in categories:
        rows.append({
            "result": label_name, "category": c,
            "pooled_n": int(pooled_n[c]), "pooled_frac": float(pooled[c]),
            "n_donors": int(donor_props.shape[0]),
            "donor_unweighted_mean_frac": float(donor_mean[c]),
            "donor_unweighted_sd_frac": float(donor_sd[c]) if donor_props.shape[0] > 1 else None,
            "n_studies": int(study_props.shape[0]),
            "study_unweighted_mean_frac": float(study_mean[c]),
            "study_unweighted_sd_frac": float(study_sd[c]) if study_props.shape[0] > 1 else None,
        })
    return pd.DataFrame(rows)


def run_composition_for_cohort(sub_scores, sub_meta, cohort_label, out_dir):
    print(f"\n=== Composition analysis: {cohort_label} (n={len(sub_scores)}) ===", flush=True)
    pct_d = sub_scores["D_Gut-shared_percentile"].values
    pct_f = sub_scores["F_Gut-specific_percentile"].values
    pct_p = sub_scores["P_Gut-specific_percentile"].values
    pct_colon = sub_scores["F_Colon-specific_percentile"].values
    pct_si = sub_scores["F_SI-specific_percentile"].values
    donor_key = sub_meta["donor_key"].values
    study_id = sub_meta["study_id"].values

    status, d_sup, f_sup, p_sup = axis_supported_status(pct_d, pct_f, pct_p)
    argmax_labels = coarse_argmax(pct_d, pct_f, pct_p)
    predominant = predominant_among_supported(pct_d, pct_f, pct_p, d_sup, f_sup, p_sup)
    f_assigned = f_sup | (predominant == "F")  # F-supported OR F-predominant
    regional = regional_refinement(pct_colon, pct_si, f_assigned)

    all_summaries = []
    all_summaries.append(categorical_donor_study_summary(status, donor_key, study_id, "step0_axis_supported_status"))
    all_summaries.append(categorical_donor_study_summary(argmax_labels, donor_key, study_id, "stepA_coarse_argmax_DESCRIPTIVE_ONLY"))
    all_summaries.append(categorical_donor_study_summary(predominant, donor_key, study_id, "stepA_prime_predominant_among_supported"))
    all_summaries.append(categorical_donor_study_summary(regional, donor_key, study_id, "stepB_regional_refinement_within_F_assigned"))
    summary_df = pd.concat(all_summaries, ignore_index=True)
    summary_df.insert(0, "cohort", cohort_label)
    summary_df.insert(1, "n_cells_in_cohort", len(sub_scores))

    out_path = f"{out_dir}/composition_{cohort_label}.tsv"
    summary_df.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {out_path}", flush=True)
    print(summary_df.to_string(index=False), flush=True)

    crosstab_df = step0_x_stepA_crosstab(status, argmax_labels, cohort_label)
    crosstab_path = f"{out_dir}/composition_{cohort_label}__step0_x_stepA_crosstab.tsv"
    crosstab_df.to_csv(crosstab_path, sep="\t", index=False)
    print(f"Wrote {crosstab_path}", flush=True)
    print(crosstab_df.to_string(index=False), flush=True)

    return summary_df, crosstab_df


def main():
    scoring_dir = sys.argv[1] if len(sys.argv) > 1 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/gut_scoring_full"
    out_dir = sys.argv[2] if len(sys.argv) > 2 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/secondary_analysis_composition"
    os.makedirs(out_dir, exist_ok=True)

    scores = pd.read_parquet(f"{scoring_dir}/crc_gut_scoring_all_panels.parquet")
    meta = pd.read_parquet(f"{scoring_dir}/crc_gut_scoring_cell_metadata.parquet")
    assert (scores.index == meta.index).all(), "scores/metadata index mismatch"
    print(f"Loaded {len(scores)} cells' scores and metadata", flush=True)

    # ---- Section 1: revCSC-high cohort definition(s) ----
    print("\n=== Section 1: revCSC-high threshold (cross-cell z-score rank) ===", flush=True)
    print("  Primary cohort: revCSC_primary27_minus_CLU_ASS1_zscore", flush=True)
    primary_masks = rank_cutoff_membership(scores["revCSC_primary27_minus_CLU_ASS1_zscore"].values)

    print("  Extended sensitivity cohort: revCSC_extended28_minus_CLU_ASS1_zscore "
          "(N_PERM=100, disclosed lower precision)", flush=True)
    extended_masks = rank_cutoff_membership(scores["revCSC_extended28_minus_CLU_ASS1_zscore"].values)

    membership = pd.DataFrame(index=scores.index)
    for label, mask in primary_masks.items():
        membership[f"revCSC_primary27_minus_CLU_ASS1_{label}"] = mask
    for label, mask in extended_masks.items():
        membership[f"revCSC_extended28_minus_CLU_ASS1_{label}"] = mask
    membership_path = f"{out_dir}/revcsc_high_cohort_membership.parquet"
    membership.to_parquet(membership_path)
    print(f"Wrote {membership_path}", flush=True)

    # ---- Section 2: composition, for every primary cutoff + the extended primary cutoff ----
    all_dfs = []
    all_crosstabs = []
    for label, mask in primary_masks.items():
        cohort_name = f"primary27_minus_CLU_ASS1_{label}"
        sub_scores = scores.loc[mask]
        sub_meta = meta.loc[mask]
        summary_df, crosstab_df = run_composition_for_cohort(sub_scores, sub_meta, cohort_name, out_dir)
        all_dfs.append(summary_df)
        all_crosstabs.append(crosstab_df)

    # F-facing overlap-safety sensitivity: extended cohort at its primary (top10pct) cutoff only
    ext_mask = extended_masks["top10pct"]
    summary_df, crosstab_df = run_composition_for_cohort(
        scores.loc[ext_mask], meta.loc[ext_mask],
        "extended28_minus_CLU_ASS1_top10pct_SENSITIVITY", out_dir)
    all_dfs.append(summary_df)
    all_crosstabs.append(crosstab_df)

    combined = pd.concat(all_dfs, ignore_index=True)
    combined_path = f"{out_dir}/composition_all_cohorts_combined.tsv"
    combined.to_csv(combined_path, sep="\t", index=False)
    print(f"\nWrote {combined_path}", flush=True)

    combined_crosstab = pd.concat(all_crosstabs, ignore_index=True)
    combined_crosstab_path = f"{out_dir}/composition_all_cohorts_step0_x_stepA_crosstab.tsv"
    combined_crosstab.to_csv(combined_crosstab_path, sep="\t", index=False)
    print(f"Wrote {combined_crosstab_path}", flush=True)

    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
