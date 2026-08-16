#!/usr/bin/env python3
"""
Step 6 dataset extension: analysis of dataset_extension_scoring.py's
output (per docs/STEP6_DATASET_EXTENSION_DESIGN.md, PR #33, APPROVE
after 3 review rounds). Two analyses, matching the locked contract:

1. Primary-extension correlation (HTAN malignant-only, CRLM
   malignant-only): the same revCSC<->D/F/P continuous-correlation
   construction as crc_gut_scoring_primary_analysis.py (PR #27) --
   percentile-based Pearson+Spearman, donor-aware, robust flag requires
   same-sign stability for BOTH metrics under leave-one-donor-out (the
   PR #27 round-2 bug-fix precedent, reused verbatim, not re-derived).
   No leave-one-STUDY-out here -- HTAN/CRLM are each a single dataset,
   not a multi-study meta-atlas; donor is the only estimable grouping.
   CRLM's LODO is computed identically but reported as a descriptive
   sensitivity diagnostic only (design section 5) -- its `robust_*`
   columns are still written (for transparency) but the results
   write-up must not describe CRLM as "robust"/"not robust", only
   report the LODO range as-is.

2. HTAN patient-matched contrast (joint malignant+normal calibration,
   design section 3): for each of the 6 relevant panels (5 D/F/P + the
   primary revCSC panel), collapse to one malignant summary (median
   percentile) and one normal-epithelial summary (median percentile,
   all subtypes pooled) per patient, for every patient with both cell
   types present, then a paired test (Wilcoxon signed-rank primary,
   paired t-test secondary) across those per-patient pairs -- patients,
   not cells, are the statistical unit (locked contract). Per-named-
   subtype summaries are also written as descriptive/sensitivity-only,
   not a second primary analysis.

Usage: python3 dataset_extension_analysis.py <scoring_out_dir> <analysis_out_dir>
(run on Argos, argos-codex env, after dataset_extension_scoring.py)
"""
import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, wilcoxon, ttest_rel

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crc_gut_scoring_core import COMPARISON_PAIRS, DFP_PANELS

MIN_CELLS_FOR_ESTIMATE = 10
PRIMARY_REVCSC_PANEL = "revCSC_primary27_minus_CLU_ASS1"
HTAN_JOINT_PANELS = DFP_PANELS + [PRIMARY_REVCSC_PANEL]


def safe_corr(x, y):
    """Reused verbatim from crc_gut_scoring_primary_analysis.py (PR #27) --
    never silently coerced to 0 or dropped, reason reported explicitly."""
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


def leave_one_donor_out_sensitivity(x, y, donor_key):
    r_all, rho_all, reason_all = safe_corr(x, y)
    donors = pd.unique(donor_key)
    rows = []
    for d in donors:
        mask = donor_key != d
        r, rho, reason = safe_corr(x[mask], y[mask])
        rows.append({
            "left_out_donor_key": d, "n_cells_remaining": int(mask.sum()),
            "pooled_pearson_r_excl": r, "pooled_spearman_rho_excl": rho,
            "signed_delta_pearson_r_at_this_exclusion": (r - r_all) if (r is not None and r_all is not None) else None,
        })
    return pd.DataFrame(rows), r_all, rho_all, reason_all


def primary_extension_correlation(scores, meta, donor_col, dataset_name, out_dir, report_robustness):
    """donor_col: metadata column giving the patient/donor identity
    (HTAN: 'Patient'; CRLM: 'donor_id'). report_robustness=False for
    CRLM -- LODO is still computed and written, but the robust flag is
    marked N/A_exploratory_per_design in the overview, not a real
    robust/non-robust classification (design section 5)."""
    overview_rows = []
    for revcsc_panel, dfp_panel in COMPARISON_PAIRS:
        pair_name = f"{revcsc_panel}__vs__{dfp_panel}"
        x_col, y_col = f"{revcsc_panel}_percentile", f"{dfp_panel}_percentile"
        if x_col not in scores.columns or y_col not in scores.columns:
            continue
        x = scores[x_col].values
        y = scores[y_col].values

        donor_df = per_donor_table(x, y, meta[donor_col].values)
        n_donors = len(donor_df)
        n_estimable = int(donor_df.estimable.sum())
        donor_path = f"{out_dir}/{dataset_name}__{pair_name}__per_donor.tsv"
        donor_df.to_csv(donor_path, sep="\t", index=False)

        lodo_df, r_all, rho_all, reason_all = leave_one_donor_out_sensitivity(x, y, meta[donor_col].values)
        lodo_path = f"{out_dir}/{dataset_name}__{pair_name}__leave_one_donor_out.tsv"
        lodo_df.to_csv(lodo_path, sep="\t", index=False)

        if r_all is None:
            robust = robust_pearson_only = robust_spearman_only = None
            max_shift_donor = None
        else:
            max_shift_row = lodo_df.loc[lodo_df.signed_delta_pearson_r_at_this_exclusion.abs().idxmax()] \
                if lodo_df.signed_delta_pearson_r_at_this_exclusion.notna().any() else None
            max_shift_donor = float(max_shift_row.signed_delta_pearson_r_at_this_exclusion) if max_shift_row is not None else None
            # Same both-metric-sign-stability requirement as PR #27
            # round 2's bug fix (reused, not re-derived).
            same_sign_r = all((r_all >= 0) == (v >= 0) for v in lodo_df.pooled_pearson_r_excl.dropna())
            same_sign_rho = all((rho_all >= 0) == (v >= 0) for v in lodo_df.pooled_spearman_rho_excl.dropna())
            robust_pearson_only = bool(same_sign_r)
            robust_spearman_only = bool(same_sign_rho)
            robust = bool(same_sign_r and same_sign_rho)

        overview_rows.append({
            "dataset": dataset_name, "pair": pair_name,
            "revcsc_panel": revcsc_panel, "dfp_panel": dfp_panel,
            "n_cells_total": len(x), "n_donors": n_donors, "n_donors_estimable": n_estimable,
            "pooled_pearson_r_all_cells": r_all, "pooled_spearman_rho_all_cells": rho_all,
            "pooled_not_estimable_reason": reason_all,
            "signed_delta_pearson_r_at_max_abs_shift_leave_one_donor_out": max_shift_donor,
            "robust_pearson_sign_stable": robust_pearson_only if report_robustness else None,
            "robust_spearman_sign_stable": robust_spearman_only if report_robustness else None,
            "robust_to_no_single_donor_out": robust if report_robustness else None,
            "robustness_classification_applicable": report_robustness,
        })
    return pd.DataFrame(overview_rows)


def htan_patient_matched_contrast(joint_scores, joint_meta, out_dir):
    """Design section 3's locked contract: joint-calibrated scores
    (already computed with malignant+normal scored together), one
    malignant summary + one normal-epithelial summary per patient per
    panel (median percentile), paired test across patients -- patients,
    not cells, are the statistical unit."""
    is_malignant = joint_meta["cell_type"].astype(str) == "malignant cell"
    patients = joint_meta["Patient"].astype(str)

    overview_rows = []
    subtype_rows = []
    for panel in HTAN_JOINT_PANELS:
        col = f"{panel}_percentile"
        if col not in joint_scores.columns:
            continue
        vals = joint_scores[col].values

        df = pd.DataFrame({
            "patient": patients.values, "cell_type": joint_meta["cell_type"].astype(str).values,
            "is_malignant": is_malignant.values, "score": vals,
        })

        mal_summary = df[df.is_malignant].groupby("patient")["score"].median()
        normal_summary = df[~df.is_malignant].groupby("patient")["score"].median()
        paired_patients = mal_summary.index.intersection(normal_summary.index)

        n_patients_total = patients.nunique()
        n_patients_paired = len(paired_patients)

        if n_patients_paired < 3:
            overview_rows.append({
                "panel": panel, "n_patients_total": n_patients_total,
                "n_patients_paired": n_patients_paired,
                "wilcoxon_stat": None, "wilcoxon_p": None,
                "paired_t_stat": None, "paired_t_p": None,
                "median_malignant_minus_normal_delta": None,
                "not_estimable_reason": f"n_patients_paired={n_patients_paired} < 3",
            })
            continue

        m = mal_summary.loc[paired_patients].values
        n = normal_summary.loc[paired_patients].values
        delta = m - n

        wstat, wp = wilcoxon(m, n)
        tstat, tp = ttest_rel(m, n)

        overview_rows.append({
            "panel": panel, "n_patients_total": n_patients_total,
            "n_patients_paired": n_patients_paired,
            "wilcoxon_stat": float(wstat), "wilcoxon_p": float(wp),
            "paired_t_stat": float(tstat), "paired_t_p": float(tp),
            "median_malignant_minus_normal_delta": float(np.median(delta)),
            "not_estimable_reason": None,
        })

        per_patient_path = f"{out_dir}/htan_patient_matched__{panel}__per_patient.tsv"
        pd.DataFrame({
            "patient": paired_patients, "malignant_median_percentile": m,
            "normal_epithelial_median_percentile": n, "delta": delta,
        }).to_csv(per_patient_path, sep="\t", index=False)

        # Descriptive/sensitivity-only: per-named-normal-subtype summary
        # (not a second primary analysis, per design section 3).
        for subtype, sub in df[~df.is_malignant].groupby("cell_type"):
            subtype_summary = sub.groupby("patient")["score"].median()
            paired_sub = mal_summary.index.intersection(subtype_summary.index)
            if len(paired_sub) < 3:
                continue
            m_sub = mal_summary.loc[paired_sub].values
            n_sub = subtype_summary.loc[paired_sub].values
            wstat_sub, wp_sub = wilcoxon(m_sub, n_sub)
            subtype_rows.append({
                "panel": panel, "normal_subtype": subtype,
                "n_patients_paired": len(paired_sub),
                "wilcoxon_stat": float(wstat_sub), "wilcoxon_p": float(wp_sub),
                "median_malignant_minus_subtype_delta": float(np.median(m_sub - n_sub)),
            })

    overview_df = pd.DataFrame(overview_rows)
    subtype_df = pd.DataFrame(subtype_rows)
    return overview_df, subtype_df


def main():
    scoring_dir = sys.argv[1] if len(sys.argv) > 1 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/dataset_extension"
    out_dir = sys.argv[2] if len(sys.argv) > 2 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/dataset_extension_analysis"
    os.makedirs(out_dir, exist_ok=True)

    print("=== HTAN primary-extension correlation (malignant-only) ===", flush=True)
    htan_scores = pd.read_parquet(f"{scoring_dir}/htan_malignant_only_scores.parquet")
    htan_meta = pd.read_parquet(f"{scoring_dir}/htan_malignant_only_cell_metadata.parquet")
    assert (htan_scores.index == htan_meta.index).all(), "HTAN malignant-only scores/metadata index mismatch"
    htan_overview = primary_extension_correlation(
        htan_scores, htan_meta, donor_col="Patient", dataset_name="HTAN_CRC_progressive_plasticity",
        out_dir=out_dir, report_robustness=True,
    )
    print(htan_overview.to_string(index=False), flush=True)

    print("\n=== CRLM primary-extension correlation (malignant-only, exploratory) ===", flush=True)
    crlm_scores = pd.read_parquet(f"{scoring_dir}/crlm_malignant_only_scores.parquet")
    crlm_meta = pd.read_parquet(f"{scoring_dir}/crlm_malignant_only_cell_metadata.parquet")
    assert (crlm_scores.index == crlm_meta.index).all(), "CRLM scores/metadata index mismatch"
    crlm_overview = primary_extension_correlation(
        crlm_scores, crlm_meta, donor_col="donor_id", dataset_name="CRLM_NMP_ATLAS",
        out_dir=out_dir, report_robustness=False,
    )
    print(crlm_overview.to_string(index=False), flush=True)

    combined_overview = pd.concat([htan_overview, crlm_overview], ignore_index=True)
    overview_path = f"{out_dir}/primary_extension_correlation_overview.tsv"
    combined_overview.to_csv(overview_path, sep="\t", index=False)
    print(f"\nWrote {overview_path}", flush=True)

    print("\n=== HTAN patient-matched contrast (joint malignant+normal calibration) ===", flush=True)
    joint_scores = pd.read_parquet(f"{scoring_dir}/htan_joint_malignant_normal_scores.parquet")
    joint_meta = pd.read_parquet(f"{scoring_dir}/htan_joint_malignant_normal_cell_metadata.parquet")
    assert (joint_scores.index == joint_meta.index).all(), "HTAN joint scores/metadata index mismatch"
    contrast_overview, contrast_subtype = htan_patient_matched_contrast(joint_scores, joint_meta, out_dir)
    print(contrast_overview.to_string(index=False), flush=True)

    contrast_overview_path = f"{out_dir}/htan_patient_matched_contrast_overview.tsv"
    contrast_overview.to_csv(contrast_overview_path, sep="\t", index=False)
    print(f"Wrote {contrast_overview_path}", flush=True)

    contrast_subtype_path = f"{out_dir}/htan_patient_matched_contrast_by_normal_subtype.tsv"
    contrast_subtype.to_csv(contrast_subtype_path, sep="\t", index=False)
    print(f"Wrote {contrast_subtype_path}", flush=True)

    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
