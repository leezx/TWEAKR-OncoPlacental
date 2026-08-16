#!/usr/bin/env python3
"""
Step 6 gut re-anchor: full 665,473-cell primary scoring compute (per
docs/STEP6_GUT_SCORING_COMPUTE_DESIGN.md, PR #26, APPROVE after 2 review
rounds). MUST be run after crc_gut_scoring_convergence_check.py -- reads
that check's `nperm500_required_panels.txt` gate file and scores any
flagged panel at N_PERM=500 instead of the default N_PERM=100 (the
locked, checkable convergence contract, not a discretionary choice made
here).

Writes per-cell empirical percentile + null-calibrated z-score for all 13
panels (665,473 cells x 13 panels x 2 metrics), plus the obs metadata
(study_id, patient_id, platform, atlas_cell_type_middle) needed for the
donor/study-aware primary analysis (crc_gut_scoring_primary_analysis.py).
Checkpointed per-panel (see crc_gut_scoring_core.score_all_panels) so a
job that fails partway can be resumed without re-scoring completed panels.

Uses score_genes_fast (crc_gut_scoring_core.py), not repeated real
scanpy.tl.score_genes calls: profiled empirically (timing probe, not
guessed) that the naive per-call approach costs ~45-55s/call at this
scale (dominated by score_genes' own internal full-genome average-
expression recomputation, done from scratch every call regardless of
gene-list size -- confirmed by reading scanpy 1.11's _score_genes_bins
source directly), which extrapolates to ~18h+ for the full 13-panel job
-- infeasible. score_genes_fast precomputes that binning once and reuses
it, empirically ~12s/call for the largest panel at full scale (~4-5x
faster for large panels, more for small ones), and is numerically
validated against real scanpy.tl.score_genes on the EXACT production
code path -- CSC-converted object, bins precomputed in production's own
call order -- before use here
(validate_score_genes_fast_csc_production_path.py: byte-identical
control-gene-set selection AND exact score agreement,
max_abs_diff=0.000e+00, after also fixing _nan_means_dense_or_sparse to
match scanpy's float64 summation exactly rather than relying on sparse
.mean()'s unspecified accumulation dtype). The convergence check (crc_gut_scoring_convergence_check.py)
deliberately uses the real, non-fast implementation throughout, since
its whole purpose is validating the statistical design against the
reference method.

Usage: python3 crc_gut_scoring_full.py <out_dir> [<convergence_check_dir>]
(run on Argos, argos-codex env)
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crc_gut_scoring_core import load_atlas, score_all_panels, ALL_PANELS, SEED

DEFAULT_N_PERM = 100
GATED_N_PERM = 500


def load_gate_panels(convergence_check_dir):
    """BUG FIX (found by reviewer, PR #27 round 2): a missing gate file
    previously only printed a WARNING and silently fell back to
    N_PERM=100 for every panel -- exactly the failure mode the locked
    contract ("crc_gut_scoring_convergence_check.py MUST be run first")
    exists to prevent, and the warning could easily scroll past unnoticed
    in a long qsub log. Fixed to raise instead: a missing gate file is
    always a real problem (either the convergence check hasn't run, or a
    path is wrong), never a legitimate reason to proceed ungated."""
    gate_path = f"{convergence_check_dir}/nperm500_required_panels.txt"
    if not os.path.exists(gate_path):
        raise FileNotFoundError(
            f"Convergence-check gate file not found at {gate_path}. "
            f"crc_gut_scoring_convergence_check.py MUST be run first, per the "
            f"locked contract -- refusing to silently proceed with N_PERM={DEFAULT_N_PERM} "
            f"for every panel ungated."
        )
    with open(gate_path) as f:
        panels = [l.strip() for l in f if l.strip()]
    return panels


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/gut_scoring_full"
    convergence_check_dir = sys.argv[2] if len(sys.argv) > 2 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/06_crc_projection/gut_scoring_convergence_check"
    os.makedirs(out_dir, exist_ok=True)
    checkpoint_dir = f"{out_dir}/_checkpoints"

    gated_panels = load_gate_panels(convergence_check_dir)
    print(f"Panels requiring N_PERM={GATED_N_PERM} (per convergence-check gate): "
          f"{gated_panels if gated_panels else 'none'}", flush=True)
    default_panels = [p for p in ALL_PANELS if p not in gated_panels]

    adata = load_atlas(n_cells_subset=None, to_csc=True)

    all_scores = None
    all_n_testable = {}

    if default_panels:
        print(f"\n=== Scoring {len(default_panels)} panels at N_PERM={DEFAULT_N_PERM} (fast path) ===", flush=True)
        scores, n_testable = score_all_panels(
            adata, n_perm=DEFAULT_N_PERM, panels=default_panels, seed=SEED,
            checkpoint_dir=checkpoint_dir, fast=True,
        )
        all_scores = scores
        all_n_testable.update(n_testable)

    if gated_panels:
        print(f"\n=== Scoring {len(gated_panels)} gated panels at N_PERM={GATED_N_PERM} (fast path) ===", flush=True)
        scores_gated, n_testable_gated = score_all_panels(
            adata, n_perm=GATED_N_PERM, panels=gated_panels, seed=SEED,
            checkpoint_dir=checkpoint_dir, fast=True,
        )
        all_scores = scores_gated if all_scores is None else all_scores.join(scores_gated)
        all_n_testable.update(n_testable_gated)

    scores_path = f"{out_dir}/crc_gut_scoring_all_panels.parquet"
    all_scores.to_parquet(scores_path)
    print(f"\nWrote {scores_path} ({all_scores.shape[0]} cells x {all_scores.shape[1]} columns)", flush=True)

    # obs metadata needed for the primary analysis
    meta_cols = ["study_id", "patient_id", "platform", "atlas_cell_type_middle",
                 "donor_id", "sample_id"]
    meta_cols = [c for c in meta_cols if c in adata.obs.columns]
    meta = adata.obs[meta_cols].copy()
    meta["donor_key"] = meta["study_id"].astype(str) + "||" + meta["patient_id"].astype(str)
    meta_path = f"{out_dir}/crc_gut_scoring_cell_metadata.parquet"
    meta.to_parquet(meta_path)
    print(f"Wrote {meta_path}", flush=True)

    n_testable_path = f"{out_dir}/n_testable_genes_per_panel.tsv"
    import pandas as pd
    n_perm_used = {p: (GATED_N_PERM if p in gated_panels else DEFAULT_N_PERM) for p in ALL_PANELS}
    pd.DataFrame([
        {"panel": p, "n_testable": all_n_testable[p], "n_perm_used": n_perm_used[p]}
        for p in ALL_PANELS
    ]).to_csv(n_testable_path, sep="\t", index=False)
    print(f"Wrote {n_testable_path}", flush=True)

    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
