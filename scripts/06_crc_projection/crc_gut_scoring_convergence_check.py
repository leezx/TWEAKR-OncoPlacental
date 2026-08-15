#!/usr/bin/env python3
"""
Step 6 gut re-anchor: N_PERM=100 vs N_PERM=500 convergence check, run
BEFORE the full 665,473-cell job (per docs/STEP6_GUT_SCORING_COMPUTE_DESIGN.md's
locked contract, PR #26 round 1). Runs the full null-calibration pipeline
on a fixed, stratified-by-study_id 20,000-cell subset at both N_PERM
settings for all 13 panels, then reports:
  1. per-panel correlation between the two N_PERM settings' per-cell
     empirical percentiles (expected ~1.0 if 100 draws is adequate);
  2. for each of the 10 revCSC<->D/F/P comparison pairs, the pooled
     Pearson/Spearman correlation computed at each N_PERM side by side,
     flagging any pair where the two settings' Pearson r differs by
     >0.02 (the locked threshold) -- such a pair is re-scored at
     N_PERM=500 in the full run.

Usage: python3 crc_gut_scoring_convergence_check.py <out_dir>
(run on Argos, argos-codex env, BEFORE crc_gut_scoring_full.py)
"""
import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crc_gut_scoring_core import (
    load_atlas, score_all_panels, ALL_PANELS, COMPARISON_PAIRS, SEED,
)

N_CELLS_SUBSET = 20000
PEARSON_R_THRESHOLD = 0.02


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/gut_scoring_convergence_check"
    os.makedirs(out_dir, exist_ok=True)

    adata = load_atlas(n_cells_subset=N_CELLS_SUBSET, seed=SEED)
    study_counts = adata.obs["study_id"].value_counts()
    print(f"Subset composition: {len(study_counts)} studies, "
          f"min={study_counts.min()}, max={study_counts.max()}", flush=True)

    print("\n=== Scoring at N_PERM=100 ===", flush=True)
    scores_100, n_testable_100 = score_all_panels(adata, n_perm=100, panels=ALL_PANELS, seed=SEED)
    scores_100.to_parquet(f"{out_dir}/scores_nperm100.parquet")

    print("\n=== Scoring at N_PERM=500 ===", flush=True)
    scores_500, n_testable_500 = score_all_panels(adata, n_perm=500, panels=ALL_PANELS, seed=SEED)
    scores_500.to_parquet(f"{out_dir}/scores_nperm500.parquet")

    # ---- 1. per-panel percentile correlation between N_PERM settings ----
    panel_rows = []
    for panel in ALL_PANELS:
        p100 = scores_100[f"{panel}_percentile"].values
        p500 = scores_500[f"{panel}_percentile"].values
        r, _ = pearsonr(p100, p500)
        rho, _ = spearmanr(p100, p500)
        panel_rows.append({
            "panel": panel, "n_testable_100perm": n_testable_100[panel],
            "n_testable_500perm": n_testable_500[panel],
            "percentile_pearson_r_100_vs_500": round(float(r), 4),
            "percentile_spearman_rho_100_vs_500": round(float(rho), 4),
        })
    panel_df = pd.DataFrame(panel_rows)
    panel_path = f"{out_dir}/convergence_per_panel_percentile_correlation.tsv"
    panel_df.to_csv(panel_path, sep="\t", index=False)
    print(f"\nWrote {panel_path}")
    print(panel_df.to_string(index=False), flush=True)

    # ---- 2. pooled revCSC<->D/F/P correlation at each N_PERM ----
    pair_rows = []
    for revcsc_panel, dfp_panel in COMPARISON_PAIRS:
        x100 = scores_100[f"{revcsc_panel}_percentile"].values
        y100 = scores_100[f"{dfp_panel}_percentile"].values
        x500 = scores_500[f"{revcsc_panel}_percentile"].values
        y500 = scores_500[f"{dfp_panel}_percentile"].values

        r100, _ = pearsonr(x100, y100)
        rho100, _ = spearmanr(x100, y100)
        r500, _ = pearsonr(x500, y500)
        rho500, _ = spearmanr(x500, y500)
        delta_r = abs(r100 - r500)
        needs_500 = bool(delta_r > PEARSON_R_THRESHOLD)

        pair_rows.append({
            "revcsc_panel": revcsc_panel, "dfp_panel": dfp_panel,
            "pearson_r_nperm100": round(float(r100), 4),
            "pearson_r_nperm500": round(float(r500), 4),
            "abs_delta_pearson_r": round(float(delta_r), 4),
            "spearman_rho_nperm100": round(float(rho100), 4),
            "spearman_rho_nperm500": round(float(rho500), 4),
            "exceeds_0.02_threshold": needs_500,
        })
    pair_df = pd.DataFrame(pair_rows)
    pair_path = f"{out_dir}/convergence_pairwise_correlation_comparison.tsv"
    pair_df.to_csv(pair_path, sep="\t", index=False)
    print(f"\nWrote {pair_path}")
    print(pair_df.to_string(index=False), flush=True)

    flag_col = pair_df["exceeds_0.02_threshold"]
    n_flagged = int(flag_col.sum())
    flagged_revcsc = pair_df.loc[flag_col, "revcsc_panel"].unique().tolist() if n_flagged else []
    flagged_dfp = pair_df.loc[flag_col, "dfp_panel"].unique().tolist() if n_flagged else []
    gate_path = f"{out_dir}/nperm500_required_panels.txt"
    flagged_panels = sorted(set(flagged_revcsc) | set(flagged_dfp))
    with open(gate_path, "w") as f:
        f.write("\n".join(flagged_panels) + ("\n" if flagged_panels else ""))
    print(f"\n{n_flagged}/10 comparison pairs exceed the 0.02 Pearson-r threshold.")
    print(f"Panels requiring N_PERM=500 in the full run (gate file: {gate_path}): {flagged_panels}")
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
