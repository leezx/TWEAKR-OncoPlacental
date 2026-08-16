#!/usr/bin/env python3
"""
Step 6 tertiary analysis compute (docs/STEP6_TERTIARY_ANALYSIS_DESIGN.md,
APPROVE after 4 review rounds): full-atlas, revCSC-independent D/F/P
axis-defined composition/substructure -- also directly answers the
axis-composition component of Q3 of the project's original 6-question
framework.

No new scoring -- reuses PR #27's already-verified full-atlas per-cell
D/F/P + revCSC percentiles. Reuses PR #28/#29's review-tested composition
methodology directly (axis_supported_status, coarse_argmax,
predominant_among_supported, regional_refinement, step0_x_stepA_crosstab,
categorical_donor_study_summary, rank_cutoff_membership -- all imported
from secondary_analysis_composition.py, not reimplemented) but applies it
to EVERY malignant cell in the atlas unconditionally, not the
revCSC-high subset.

Sections implemented (design doc numbering):
  Sec 1: population = all 665,473 cells, no revCSC-based filtering.
  Sec 2: axis-supported status (primary result) + Step0xStepA cross-tab
    (pooled-only, per the design's explicit exemption) + Step A' + Step B.
  Sec 3: Q3 F x P quadrant table (derived from the same d_sup/f_sup/p_sup
    used in sec 2, same percentile>=90 threshold, no new pick).
  Sec 4a: axis-supported category x operational revCSC-high cohort
    (cross-cell z-rank, top 10/5/20%) -- "outside the operational
    revCSC-high cohort" framing, not "revCSC does not capture this cell."
  Sec 4b: axis-supported category x revCSC matched-null support
    (percentile>=90, same semantics as D/F/P) -- "revCSC-not-supported,"
    not "lacks evidence."
  Sec 5: donor/study-aware aggregation for every table except the
    intentionally pooled-only Step0xStepA cross-tab.

Usage: python3 tertiary_analysis_composition.py <scoring_full_dir> <out_dir>
(run on Argos via qsub, for consistency with standing discipline, though
the analysis itself needs no new scoring)
"""
import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from secondary_analysis_composition import (
    axis_supported_status, coarse_argmax, predominant_among_supported,
    regional_refinement, step0_x_stepA_crosstab, categorical_donor_study_summary,
    rank_cutoff_membership, SUPPORT_PCT_THRESHOLD,
)

REVCSC_PANEL = "revCSC_primary27_minus_CLU_ASS1"


def q3_fp_quadrant(f_sup, p_sup):
    """Sec 3: F x P 2D quadrant, collapsed from the same axis-supported
    flags used in sec 2 (same percentile>=90 threshold, no new pick).
    4 categories: Fetal-high/Placenta-low, Fetal-high/Placenta-high,
    Fetal-low/Placenta-high, double-low."""
    quad = np.full(len(f_sup), "double_low", dtype=object)
    quad = np.where(f_sup & ~p_sup, "Fetal_high_Placenta_low", quad)
    quad = np.where(f_sup & p_sup, "Fetal_high_Placenta_high", quad)
    quad = np.where(~f_sup & p_sup, "Fetal_low_Placenta_high", quad)
    return quad


def joint_labels(status, other, other_name):
    """Builds a single composite string label per cell (e.g.
    'P_only|revCSC_high_top10pct') so the joint distribution can be run
    through categorical_donor_study_summary directly -- no new
    statistical procedure, per the design's sec 5."""
    other_str = np.where(other, f"{other_name}", f"not_{other_name}")
    return np.array([f"{s}|{o}" for s, o in zip(status, other_str)], dtype=object)


def main():
    scoring_dir = sys.argv[1] if len(sys.argv) > 1 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/gut_scoring_full"
    out_dir = sys.argv[2] if len(sys.argv) > 2 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/tertiary_analysis_composition"
    os.makedirs(out_dir, exist_ok=True)

    scores = pd.read_parquet(f"{scoring_dir}/crc_gut_scoring_all_panels.parquet")
    meta = pd.read_parquet(f"{scoring_dir}/crc_gut_scoring_cell_metadata.parquet")
    assert (scores.index == meta.index).all(), "scores/metadata index mismatch"
    print(f"Loaded {len(scores)} cells' scores and metadata (sec 1: unconditional population, "
          f"no revCSC-based filtering)", flush=True)

    donor_key = meta["donor_key"].values
    study_id = meta["study_id"].values

    # ---- Section 2: axis-supported status (primary result) ----
    print("\n=== Section 2: axis-supported status (full atlas, unconditional on revCSC) ===", flush=True)
    pct_d = scores["D_Gut-shared_percentile"].values
    pct_f = scores["F_Gut-specific_percentile"].values
    pct_p = scores["P_Gut-specific_percentile"].values
    pct_colon = scores["F_Colon-specific_percentile"].values
    pct_si = scores["F_SI-specific_percentile"].values

    status, d_sup, f_sup, p_sup = axis_supported_status(pct_d, pct_f, pct_p)
    argmax_labels = coarse_argmax(pct_d, pct_f, pct_p)
    predominant = predominant_among_supported(pct_d, pct_f, pct_p, d_sup, f_sup, p_sup)
    f_assigned = f_sup | (predominant == "F")
    regional = regional_refinement(pct_colon, pct_si, f_assigned)

    all_summaries = []
    all_summaries.append(categorical_donor_study_summary(status, donor_key, study_id, "step0_axis_supported_status"))
    all_summaries.append(categorical_donor_study_summary(argmax_labels, donor_key, study_id, "stepA_coarse_argmax_DESCRIPTIVE_ONLY"))
    all_summaries.append(categorical_donor_study_summary(predominant, donor_key, study_id, "stepA_prime_predominant_among_supported"))
    all_summaries.append(categorical_donor_study_summary(regional, donor_key, study_id, "stepB_regional_refinement_within_F_assigned"))

    # Step0xStepA cross-tab: pooled-only by design (sec 5's explicit exemption)
    crosstab_df = step0_x_stepA_crosstab(status, argmax_labels, "full_atlas")
    crosstab_path = f"{out_dir}/tertiary_step0_x_stepA_crosstab.tsv"
    crosstab_df.to_csv(crosstab_path, sep="\t", index=False)
    print(f"Wrote {crosstab_path} (pooled-only, per sec 5's explicit exemption)", flush=True)

    # ---- Section 3: Q3 F x P quadrant table ----
    print("\n=== Section 3: Q3 F x P quadrant table (derived from same sec-2 flags) ===", flush=True)
    quad = q3_fp_quadrant(f_sup, p_sup)
    all_summaries.append(categorical_donor_study_summary(quad, donor_key, study_id, "Q3_FxP_quadrant_occupancy"))

    section23_df = pd.concat(all_summaries, ignore_index=True)
    section23_path = f"{out_dir}/tertiary_composition_sec2_sec3.tsv"
    section23_df.to_csv(section23_path, sep="\t", index=False)
    print(f"Wrote {section23_path}", flush=True)
    print(section23_df.to_string(index=False), flush=True)

    # ---- Section 4a: axis-supported category x operational revCSC-high cohort ----
    print("\n=== Section 4a: axis-supported status x operational revCSC-high cohort (cross-cell rank) ===", flush=True)
    revcsc_zscore = scores[f"{REVCSC_PANEL}_zscore"].values
    rank_masks = rank_cutoff_membership(revcsc_zscore)  # {"top5pct":..., "top10pct":..., "top20pct":...}

    sec4a_summaries = []
    for cutoff_label, mask in rank_masks.items():
        joint = joint_labels(status, mask, f"revCSC_high_{cutoff_label}")
        sec4a_summaries.append(categorical_donor_study_summary(
            joint, donor_key, study_id, f"sec4a_status_x_operational_revCSC_high_{cutoff_label}"))
    sec4a_df = pd.concat(sec4a_summaries, ignore_index=True)
    sec4a_path = f"{out_dir}/tertiary_composition_sec4a_operational_revCSC_cohort.tsv"
    sec4a_df.to_csv(sec4a_path, sep="\t", index=False)
    print(f"Wrote {sec4a_path}", flush=True)

    # ---- Section 4b: axis-supported category x revCSC matched-null support ----
    print("\n=== Section 4b: axis-supported status x revCSC matched-null support (percentile>=90) ===", flush=True)
    revcsc_pct = scores[f"{REVCSC_PANEL}_percentile"].values
    revcsc_supported = revcsc_pct >= SUPPORT_PCT_THRESHOLD
    n_revcsc_supported = int(revcsc_supported.sum())
    print(f"  revCSC-supported (percentile>={SUPPORT_PCT_THRESHOLD}): {n_revcsc_supported}/{len(scores)} "
          f"({n_revcsc_supported/len(scores)*100:.2f}%)", flush=True)

    joint_4b = joint_labels(status, revcsc_supported, "revCSC_supported")
    sec4b_df = categorical_donor_study_summary(joint_4b, donor_key, study_id, "sec4b_status_x_revCSC_matched_null_support")
    sec4b_path = f"{out_dir}/tertiary_composition_sec4b_revCSC_matched_null_support.tsv"
    sec4b_df.to_csv(sec4b_path, sep="\t", index=False)
    print(f"Wrote {sec4b_path}", flush=True)

    # ---- Combined output ----
    combined = pd.concat([section23_df, sec4a_df, sec4b_df], ignore_index=True)
    combined_path = f"{out_dir}/tertiary_composition_all_sections_combined.tsv"
    combined.to_csv(combined_path, sep="\t", index=False)
    print(f"\nWrote {combined_path}", flush=True)

    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
