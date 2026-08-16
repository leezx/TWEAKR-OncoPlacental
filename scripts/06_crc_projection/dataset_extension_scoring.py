#!/usr/bin/env python3
"""
Step 6 dataset extension: real scoring compute for
HTAN_CRC_progressive_plasticity + CRLM_NMP_ATLAS (per
docs/STEP6_DATASET_EXTENSION_DESIGN.md, PR #33, APPROVE after 3 review
rounds). Produces per-cell null-calibrated percentile + z-score for all
13 panels, using the locked, corrected contract:

- N_PERM=500 for every panel on both datasets (round-1 review fix --
  PR #27's own convergence gate did NOT certify N_PERM=100 for every
  panel, and that certification is population-specific to the primary
  atlas's detectability strata anyway).
- HTAN gets TWO separate scoring passes on two different populations
  (round-2 review fix -- these are not interchangeable calibration
  bases): (1) malignant-cells-only, all 13 panels -- the true
  PR #27-equivalent primary-extension analysis; (2) malignant +
  normal-epithelial cells scored JOINTLY, the 6 relevant panels (5 D/F/P
  + the primary revCSC panel) -- feeds the patient-matched contrast
  analysis (a separate script, dataset_extension_analysis.py, run after
  this compute).
- CRLM gets one scoring pass, malignant-cells-only, all 13 panels,
  explicitly exploratory (docs section 5) -- no composition re-run.

Also runs the coverage check, canonical-marker sentinel check, and HTAN
provenance-overlap audit (all required pre-compute/pre-writeup steps per
the design) and writes their results alongside the scores.

Usage: python3 dataset_extension_scoring.py <out_dir>
(run on Argos, argos-codex env, qsub only -- never locally)
"""
import sys
import os
import time
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crc_gut_scoring_core import (
    score_all_panels, load_panel_ensembl_ids, ALL_PANELS, DFP_PANELS, SEED,
)
from dataset_extension_core import (
    EXTENSION_DATASETS, load_extension_dataset, coverage_check,
    canonical_marker_sentinel_check, htan_provenance_overlap_audit,
)

N_PERM = 500  # locked for ALL panels on both extension datasets, round-1 fix
PRIMARY_REVCSC_PANEL = "revCSC_primary27_minus_CLU_ASS1"
HTAN_JOINT_PANELS = DFP_PANELS + [PRIMARY_REVCSC_PANEL]  # 6 panels, patient-matched-contrast pass


def score_and_write(adata, panels, n_perm, out_dir, tag, meta_cols):
    checkpoint_dir = f"{out_dir}/_checkpoints_{tag}"
    t0 = time.time()
    scores, n_testable = score_all_panels(
        adata, n_perm=n_perm, panels=panels, seed=SEED,
        checkpoint_dir=checkpoint_dir, fast=True,
    )
    print(f"[{tag}] scored {len(panels)} panels x {adata.n_obs} cells in "
          f"{time.time()-t0:.1f}s", flush=True)

    scores_path = f"{out_dir}/{tag}_scores.parquet"
    scores.to_parquet(scores_path)
    print(f"[{tag}] wrote {scores_path} ({scores.shape[0]} cells x {scores.shape[1]} cols)", flush=True)

    meta_cols_present = [c for c in meta_cols if c in adata.obs.columns]
    meta = adata.obs[meta_cols_present].copy()
    meta_path = f"{out_dir}/{tag}_cell_metadata.parquet"
    meta.to_parquet(meta_path)
    print(f"[{tag}] wrote {meta_path}", flush=True)

    n_testable_path = f"{out_dir}/{tag}_n_testable_genes_per_panel.tsv"
    pd.DataFrame([
        {"panel": p, "n_testable": n_testable[p], "n_perm_used": n_perm}
        for p in panels
    ]).to_csv(n_testable_path, sep="\t", index=False)
    print(f"[{tag}] wrote {n_testable_path}", flush=True)
    return scores_path, meta_path, n_testable_path


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/dataset_extension"
    os.makedirs(out_dir, exist_ok=True)

    # ---- Pre-compute gates ----
    print("\n=== Loading datasets ===", flush=True)
    htan = load_extension_dataset(EXTENSION_DATASETS["HTAN_CRC_progressive_plasticity"],
                                   "HTAN_CRC_progressive_plasticity")
    crlm = load_extension_dataset(EXTENSION_DATASETS["CRLM_NMP_ATLAS"], "CRLM_NMP_ATLAS")

    print("\n=== Coverage check (all 13 panels, both datasets) ===", flush=True)
    cov_htan = coverage_check(htan, "HTAN_CRC_progressive_plasticity", ALL_PANELS, load_panel_ensembl_ids)
    cov_crlm = coverage_check(crlm, "CRLM_NMP_ATLAS", ALL_PANELS, load_panel_ensembl_ids)
    cov_all = pd.concat([cov_htan, cov_crlm], ignore_index=True)
    cov_path = f"{out_dir}/coverage_check.tsv"
    cov_all.to_csv(cov_path, sep="\t", index=False)
    print(f"Wrote {cov_path}", flush=True)
    # Round-1 review correction: a flat 50% floor missed the real case
    # this project actually hit -- CRLM's P_Gut-specific at 40/76=52.6%
    # is well above 50% but far below every other panel/dataset pair
    # (next lowest: 87.5%), i.e. "unexpectedly low" relative to its
    # peers per the design's actual wording, not relative to an
    # arbitrary absolute floor. Flagged here as a per-dataset relative
    # deviation instead (>15 points below that dataset's own median
    # coverage). The design's contract is a REQUIRED INVESTIGATION gate,
    # not an automatic compute-blocking gate or an automatically-waived
    # one -- this script does not itself decide the panel remains usable;
    # it only surfaces the deviation loudly so the investigation (done
    # once, by hand, and reported in docs/STEP6_DATASET_EXTENSION_RESULTS.md)
    # is not silently skipped.
    cov_all["coverage_frac"] = cov_all["n_testable"] / cov_all["n_panel_genes"]
    cov_all["dataset_median_coverage_frac"] = cov_all.groupby("dataset")["coverage_frac"].transform("median")
    flagged = cov_all[cov_all["coverage_frac"] < cov_all["dataset_median_coverage_frac"] - 0.15]
    if len(flagged) > 0:
        print(f"REQUIRES INVESTIGATION (per design's pre-compute gate, not "
              f"auto-waived): {len(flagged)} panel/dataset pairs have coverage "
              f">15 points below their own dataset's median -- see {cov_path}. "
              f"Investigation result MUST be documented in the results write-up "
              f"before proceeding to claims based on that panel:\n{flagged}", flush=True)

    print("\n=== Canonical-marker sentinel check ===", flush=True)
    htan_group_col = "cell_type"
    crlm_group_col = "cell_type"
    marker_htan = canonical_marker_sentinel_check(htan, "HTAN_CRC_progressive_plasticity", htan_group_col)
    marker_crlm = canonical_marker_sentinel_check(crlm, "CRLM_NMP_ATLAS", crlm_group_col)
    marker_all = pd.concat([marker_htan, marker_crlm], ignore_index=True)
    marker_path = f"{out_dir}/canonical_marker_sentinel_check.tsv"
    marker_all.to_csv(marker_path, sep="\t", index=False)
    print(f"Wrote {marker_path}", flush=True)

    print("\n=== HTAN study-provenance overlap audit ===", flush=True)
    prov_summary, prov_conclusion, prov_studies = htan_provenance_overlap_audit()
    prov_path = f"{out_dir}/htan_provenance_overlap_audit.tsv"
    prov_summary.to_csv(prov_path, sep="\t", index=False)
    with open(f"{out_dir}/htan_provenance_overlap_audit_conclusion.txt", "w") as f:
        f.write(f"conclusion: {prov_conclusion}\n")
        f.write(f"studies_with_id_prefix_hits: {prov_studies}\n")
    print(f"Wrote {prov_path}; conclusion={prov_conclusion}, studies={prov_studies}", flush=True)

    # ---- HTAN pass 1: malignant-only, all 13 panels (primary-extension) ----
    print("\n=== HTAN pass 1: malignant-only, all 13 panels ===", flush=True)
    htan_malignant_mask = htan.obs["cell_type"].astype(str) == "malignant cell"
    htan_malignant = htan[htan_malignant_mask].copy()
    print(f"HTAN malignant-only: {htan_malignant.n_obs} cells", flush=True)
    meta_cols_htan = ["cell_type", "Patient", "donor_id", "Sample ID",
                      "Tumor Status", "Sample Type", "Primary Site", "Site"]
    score_and_write(htan_malignant, ALL_PANELS, N_PERM, out_dir,
                     "htan_malignant_only", meta_cols_htan)

    # ---- HTAN pass 2: malignant + normal, jointly, 6 relevant panels ----
    print("\n=== HTAN pass 2: malignant+normal jointly, 6 panels (patient-matched contrast) ===", flush=True)
    print(f"HTAN joint (all epithelial): {htan.n_obs} cells", flush=True)
    score_and_write(htan, HTAN_JOINT_PANELS, N_PERM, out_dir,
                     "htan_joint_malignant_normal", meta_cols_htan)

    # ---- CRLM: malignant-only, all 13 panels, exploratory ----
    print("\n=== CRLM: malignant-only, all 13 panels (exploratory) ===", flush=True)
    crlm_malignant_mask = crlm.obs["cell_type"].astype(str) == "malignant cell"
    crlm_malignant = crlm[crlm_malignant_mask].copy()
    print(f"CRLM malignant-only: {crlm_malignant.n_obs} cells", flush=True)
    meta_cols_crlm = ["cell_type", "donor_id", "timepoint"]
    score_and_write(crlm_malignant, ALL_PANELS, N_PERM, out_dir,
                     "crlm_malignant_only", meta_cols_crlm)

    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
