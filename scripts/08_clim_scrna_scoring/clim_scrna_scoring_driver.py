#!/usr/bin/env python3
"""
Step 8 compute: D/F/P/revCSC scoring for the Step 7 scRNA-seq cohorts,
per docs/STEP8_CLIM_SCRNA_SCORING_DESIGN.md (PR #37, APPROVE after 3
review rounds). Runs the locked coverage check + required empirical
gene investigation, then scores all 4 locked populations (GSE231559
primary [15 samples], GSE231559 normal [11 samples], GSE285990 [10
samples], GSE225857 non-immune [GSM7058755 only]) against all 13 frozen
D/F/P/revCSC panels, reusing crc_gut_scoring_core.py unchanged
(N_PERM=500, same null-calibration method).

Usage: python3 clim_scrna_scoring_driver.py <out_dir>
(run on Argos, argos-codex env, qsub only -- never locally)
"""
import sys
import os
import time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "06_crc_projection"))

from crc_gut_scoring_core import score_all_panels, load_panel_ensembl_ids, ALL_PANELS, SEED
from clim_scrna_scoring_core import (
    load_gse231559_population, load_gse285990_population,
    load_gse225857_nonimmune_population, coverage_check,
    canonicalize_ensembl_ids,
    GSE231559_PRIMARY_POPULATION_GSMS, GSE231559_NORMAL_GSMS,
)

N_PERM = 500  # locked for all panels x all populations, per design

# Small panels the design requires extending the coverage investigation
# to (beyond P_Gut-specific), since proportionally large losses on a
# small panel don't trigger the aggregate ">15pts below dataset median"
# rule but are still worth checking per-gene.
SMALL_PANELS = ["D_Gut-shared"] + [p for p in ALL_PANELS if p.startswith("revCSC_")]


def investigate_missing_genes(population_name, missing_by_panel, reference_adatas):
    """Required empirical (not mechanism-assumed) investigation for
    GSE225857's coverage deviations, per the design's round-2-corrected
    contract: check each missing gene's real detection/expression status
    in the OTHER populations that DO have it (GSE231559 primary,
    GSE285990 -- both ~96%+ coverage, confirmed in the design), rather
    than assigning absence to any specific named GSE225857 preprocessing
    mechanism this project could not confirm applies to the object being
    scored (round-2 review finding)."""
    rows = []
    for panel, missing_genes in missing_by_panel.items():
        for gene in missing_genes:
            row = {"population": population_name, "panel": panel, "gene": gene}
            for ref_name, adata in reference_adatas.items():
                if gene in adata.var_names:
                    counts = adata.layers["counts"][:, adata.var_names.get_loc(gene)]
                    counts = np.asarray(counts.todense()).ravel() if hasattr(counts, "todense") else np.asarray(counts).ravel()
                    detect_frac = float((counts > 0).mean())
                    mean_count = float(counts.mean())
                    row[f"{ref_name}_present"] = True
                    row[f"{ref_name}_detect_frac"] = round(detect_frac, 4)
                    row[f"{ref_name}_mean_count"] = round(mean_count, 4)
                else:
                    row[f"{ref_name}_present"] = False
                    row[f"{ref_name}_detect_frac"] = None
                    row[f"{ref_name}_mean_count"] = None
            rows.append(row)
    return pd.DataFrame(rows)


def score_and_write(adata, panels, n_perm, out_dir, tag, seed=SEED):
    checkpoint_dir = f"{out_dir}/_checkpoints_{tag}"
    t0 = time.time()
    scores, n_testable = score_all_panels(
        adata, n_perm=n_perm, panels=panels, seed=seed,
        checkpoint_dir=checkpoint_dir, fast=True,
    )
    print(f"[{tag}] scored {len(panels)} panels x {adata.n_obs} cells in "
          f"{time.time()-t0:.1f}s", flush=True)

    scores_path = f"{out_dir}/{tag}_scores.parquet"
    scores.to_parquet(scores_path)
    print(f"[{tag}] wrote {scores_path} ({scores.shape[0]} cells x {scores.shape[1]} cols)", flush=True)

    n_testable_path = f"{out_dir}/{tag}_n_testable_genes_per_panel.tsv"
    pd.DataFrame([
        {"panel": p, "n_testable": n_testable[p], "n_perm_used": n_perm}
        for p in panels
    ]).to_csv(n_testable_path, sep="\t", index=False)
    print(f"[{tag}] wrote {n_testable_path}", flush=True)
    return scores_path, n_testable_path


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/08_clim_scrna_scoring"
    os.makedirs(out_dir, exist_ok=True)

    print("\n=== Loading all 4 locked scoring populations ===", flush=True)
    gse231559_primary = load_gse231559_population(
        GSE231559_PRIMARY_POPULATION_GSMS, "GSE231559_primary")
    gse231559_normal = load_gse231559_population(
        GSE231559_NORMAL_GSMS, "GSE231559_normal")
    gse285990 = load_gse285990_population()
    gse225857_nonimmune, gse225857_mapping_stats = load_gse225857_nonimmune_population()

    populations = {
        "gse231559_primary": gse231559_primary,
        "gse231559_normal": gse231559_normal,
        "gse285990": gse285990,
        "gse225857_nonimmune": gse225857_nonimmune,
    }

    print("\n=== Coverage check (all 13 panels x all 4 populations) ===", flush=True)
    cov_all = pd.concat([
        coverage_check(adata, name, ALL_PANELS, load_panel_ensembl_ids)
        for name, adata in populations.items()
    ], ignore_index=True)
    cov_all["coverage_frac"] = cov_all["n_testable"] / cov_all["n_panel_genes"]
    cov_all["population_median_coverage_frac"] = cov_all.groupby("population")["coverage_frac"].transform("median")
    cov_all["deviation_flagged"] = (
        cov_all["coverage_frac"] < cov_all["population_median_coverage_frac"] - 0.15
    )
    cov_path = f"{out_dir}/coverage_check.tsv"
    cov_all.to_csv(cov_path, sep="\t", index=False)
    print(f"Wrote {cov_path}", flush=True)
    flagged = cov_all[cov_all["deviation_flagged"]]
    print(f"Flagged panel/population pairs (>15pts below population median): "
          f"{len(flagged)}", flush=True)
    for _, row in flagged.iterrows():
        print(f"  {row['population']}/{row['panel']}: "
              f"{row['n_testable']}/{row['n_panel_genes']} "
              f"({row['coverage_frac']:.1%}, median {row['population_median_coverage_frac']:.1%})",
              flush=True)

    print("\n=== Required empirical gene investigation (GSE225857 non-immune, "
          "small panels + any flagged deviation) ===", flush=True)
    # Per the design's round-2-corrected contract: investigate D_Gut-shared,
    # all revCSC panels, and any flagged deviation for GSE225857
    # non-immune specifically -- empirically, using the reference
    # populations that DO have near-complete coverage (GSE231559 primary,
    # GSE285990), not assigned to an assumed preprocessing mechanism.
    investigate_panels = set(SMALL_PANELS) | set(
        flagged[flagged["population"] == "gse225857_nonimmune"]["panel"]
    )
    gse225857_genes = set(gse225857_nonimmune.var_names)
    missing_by_panel = {}
    for panel in investigate_panels:
        panel_genes = canonicalize_ensembl_ids(load_panel_ensembl_ids(panel))
        missing = [g for g in panel_genes if g not in gse225857_genes]
        if missing:
            missing_by_panel[panel] = missing
    investigation_df = investigate_missing_genes(
        "gse225857_nonimmune", missing_by_panel,
        {"gse231559_primary": gse231559_primary, "gse285990": gse285990},
    )
    investigation_path = f"{out_dir}/gse225857_missing_gene_investigation.tsv"
    investigation_df.to_csv(investigation_path, sep="\t", index=False)
    print(f"Wrote {investigation_path} ({len(investigation_df)} missing "
          f"gene x panel rows investigated)", flush=True)
    if len(investigation_df):
        both_present = investigation_df[
            investigation_df["gse231559_primary_present"] & investigation_df["gse285990_present"]
        ]
        if len(both_present):
            print(f"  Missing genes present in BOTH reference populations: "
                  f"{len(both_present)}/{len(investigation_df)}, median detect_frac "
                  f"gse231559_primary={both_present['gse231559_primary_detect_frac'].median():.3f}, "
                  f"gse285990={both_present['gse285990_detect_frac'].median():.3f}", flush=True)

    print("\n=== Scoring: GSE231559 primary (15 samples, primary "
          "population) ===", flush=True)
    score_and_write(gse231559_primary, ALL_PANELS, N_PERM, out_dir, "gse231559_primary")

    print("\n=== Scoring: GSE231559 normal (11 samples, separate "
          "calibration pass) ===", flush=True)
    score_and_write(gse231559_normal, ALL_PANELS, N_PERM, out_dir, "gse231559_normal")

    print("\n=== Scoring: GSE285990 (10 samples) ===", flush=True)
    score_and_write(gse285990, ALL_PANELS, N_PERM, out_dir, "gse285990")

    print("\n=== Scoring: GSE225857 non-immune (GSM7058755 only) ===", flush=True)
    score_and_write(gse225857_nonimmune, ALL_PANELS, N_PERM, out_dir, "gse225857_nonimmune")

    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
