#!/usr/bin/env python3
"""
Three-layer statistical validation, layers 1+3: hypergeometric enrichment
(with fold-enrichment/odds-ratio/95% CI/BH-FDR) and size-and-expression-
matched permutation nulls, per the user's explicit statistical-rigor
directive when redirecting the project onto real fetal gut data
(docs/STEP4A_GUT_FDEV_DESIGN.md's "What this design does NOT do yet" ->
now being done).

Tests each of the 5 independent mike_verzi mouse/human fetal-development
signatures against F_Colon-developmental (LargeInt) and F_SI-developmental
(SmallInt) -- the real, frozen gut-specific F-developmental gene sets from
PR #21, never used to construct these signatures or vice versa (non-
circularity requirement from the design doc).

Universe (locked before any test, per the user's requirement that the
background be "genes actually measurable/orthology-resolved in the
relevant data, not whole genome"):
  U = {human genes with >=1 Compara-verified one2one mouse ortholog}
      ^ {genes passing filterByExpr in that region's primary edgeR fit}
  (same "eligible to appear in either the signature or the target" logic
  already used to fix F_Gut-core's universe in PR #21 round 1.)

Layer 1 -- hypergeometric enrichment:
  fold_enrichment = (k/n) / (K/N)
  odds_ratio + 95% CI (Haldane-Anscombe continuity-corrected log-OR CI)
  p_value = hypergeom.sf(k-1, N, K, n) (one-sided over-representation)
  BH-FDR across all 10 tests (5 signatures x 2 regions) together

Layer 3 -- size-and-expression-matched permutation null (10,000 draws):
  For each real signature gene, its logCPM decile bin (from the real
  edgeR fit, a real measured expression-level covariate, not an assumed
  one) is recorded; each permutation draws one random gene per bin
  (without replacement within a draw) from the same universe, preserving
  BOTH exact gene-set size and expression-decile composition. Reports
  observed overlap, null-expected overlap, fold-enrichment vs. null,
  empirical p-value, and null Z-score -- the strongest defense against
  the "gene-set size explains the overlap" critique, per the user's
  explicit request.

Usage: python3 mike_verzi_gut_enrichment_permutation.py <out_dir>
(run on Argos, argos-codex env; real permutation compute, not gene-list
arithmetic, so run via qsub not locally)
"""
import sys
import os
import csv
import json
import numpy as np
from scipy.stats import hypergeom, fisher_exact
import pandas as pd

REPO = "/home/zz950/TWEAKR-OncoPlacental"
OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "results/04a_dfp_gut/mike_verzi_validation"
os.makedirs(OUT_DIR, exist_ok=True)

MIKE_VERZI_DIR = os.path.join(REPO, "results/06a_normal_context")
GUT_DIR = os.path.join(REPO, "results/04a_dfp_gut")
SIGNATURES = ["YAP_SIGNALING_GENES", "REVIVAL_STEM_CELL_GENES", "FETAL_SPHEROID_EPITHELIUM_GENES",
              "REGENERATIVE_EPITHELIUM", "FETAL_INTESTINE_GENES"]
REGIONS = {"Colon": ("LargeInt", "F_Colon-developmental.txt"), "SI": ("SmallInt", "F_SI-developmental.txt")}
N_PERM = 10000
N_BINS = 20
RNG_SEED = 20260815  # fixed seed, stated explicitly for reproducibility

rng = np.random.default_rng(RNG_SEED)


def load_gene_set(path):
    with open(path) as f:
        return {l.strip() for l in f if l.strip()}


def load_edger(path):
    df = pd.read_csv(path, sep="\t")
    return df.set_index("gene")


def bh_fdr(pvals):
    pvals = np.asarray(pvals)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    fdr = ranked * n / (np.arange(n) + 1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr = np.clip(fdr, 0, 1)
    out = np.empty(n)
    out[order] = fdr
    return out


def odds_ratio_ci(a, b, c, d):
    """2x2 table: a=k, b=n-k, c=K-k, d=N-K-n+k. Haldane-Anscombe continuity
    correction (+0.5 to all cells) applied uniformly, not just when a cell
    is zero, for a consistent CI computation across all tests."""
    a2, b2, c2, d2 = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orr = (a2 * d2) / (b2 * c2)
    se_log_or = np.sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
    log_or = np.log(orr)
    ci_low = np.exp(log_or - 1.96 * se_log_or)
    ci_high = np.exp(log_or + 1.96 * se_log_or)
    return orr, ci_low, ci_high


def one2one_human_targets(mouse_biomart_path):
    rows = list(csv.DictReader(open(mouse_biomart_path), delimiter="\t"))
    return {r["Human gene name"] for r in rows
            if r["Human gene name"] and r["Human homology type"] == "ortholog_one2one"}


def permutation_null(signature_genes_in_U, target_in_U, U_expr_bins, k_observed, n_perm, rng):
    """signature_genes_in_U: list of genes (in U) belonging to the real signature.
    U_expr_bins: dict gene -> bin index, covers all of U.
    target_in_U: set of target genes (already intersected with U)."""
    bin_to_genes = {}
    for g in U_expr_bins:
        bin_to_genes.setdefault(U_expr_bins[g], []).append(g)
    bin_arrays = {b: np.array(genes) for b, genes in bin_to_genes.items()}

    sig_bins = [U_expr_bins[g] for g in signature_genes_in_U]  # multiset of bins for the real signature
    bin_counts = pd.Series(sig_bins).value_counts().to_dict()

    null_overlaps = np.empty(n_perm, dtype=int)
    for p in range(n_perm):
        drawn = []
        for b, count in bin_counts.items():
            pool = bin_arrays[b]
            idx = rng.choice(len(pool), size=count, replace=False)
            drawn.extend(pool[idx].tolist())
        null_overlaps[p] = len(set(drawn) & target_in_U)

    expected_null = float(np.mean(null_overlaps))
    std_null = float(np.std(null_overlaps))
    fold_vs_null = k_observed / expected_null if expected_null > 0 else float("inf")
    empirical_p = (np.sum(null_overlaps >= k_observed) + 1) / (n_perm + 1)
    null_z = (k_observed - expected_null) / std_null if std_null > 0 else float("nan")
    return {
        "observed_overlap": int(k_observed),
        "null_expected_overlap": round(expected_null, 3),
        "null_std_overlap": round(std_null, 3),
        "fold_enrichment_vs_null": round(fold_vs_null, 3),
        "empirical_p": float(empirical_p),
        "null_z_score": round(null_z, 3) if not np.isnan(null_z) else None,
    }


def main():
    one2one = one2one_human_targets(f"{MIKE_VERZI_DIR}/mouse_biomart_full.tsv")
    print(f"Human genes with >=1 Compara one2one mouse ortholog: {len(one2one)}", flush=True)

    sig_sets = {sig: load_gene_set(f"{MIKE_VERZI_DIR}/{sig}_human_primary.txt") for sig in SIGNATURES}

    enrichment_rows = []
    permutation_rows = []

    for label, (region, target_file) in REGIONS.items():
        edger = load_edger(f"{GUT_DIR}/edgeR/{region}_edgeR_primary.tsv")
        region_testable = set(edger.index)
        U = one2one & region_testable
        N = len(U)
        print(f"\n=== {label} ({region}): U = one2one ^ filterByExpr-tested = {N} ===", flush=True)

        target = load_gene_set(f"{GUT_DIR}/dfp_gut_gene_sets/{target_file}")
        target_in_U = target & U
        K = len(target_in_U)
        print(f"  F_{label}-developmental in U: {K} (of {len(target)} total)", flush=True)

        # expression bins from real logCPM, computed on U only (region-specific, real data)
        logcpm_U = edger.loc[list(U), "logCPM"]
        bin_edges = np.quantile(logcpm_U, np.linspace(0, 1, N_BINS + 1))
        bin_edges[0] -= 1e-9
        bin_edges[-1] += 1e-9
        bin_idx = np.digitize(logcpm_U.values, bin_edges) - 1
        U_expr_bins = dict(zip(logcpm_U.index, bin_idx))

        for sig in SIGNATURES:
            sig_in_U = sig_sets[sig] & U
            n = len(sig_in_U)
            k = len(sig_in_U & target_in_U)
            expected = n * K / N if N else 0.0
            fold = (k / n) / (K / N) if (k > 0 and n > 0 and K > 0) else 0.0
            p = hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0

            a, b, c, d = k, n - k, K - k, N - K - n + k
            orr, ci_low, ci_high = odds_ratio_ci(a, b, c, d)

            enrichment_rows.append({
                "region": label, "signature": sig, "n_signature_in_U": n, "K_target_in_U": K,
                "k_overlap": k, "N_universe": N, "expected": round(expected, 3),
                "fold_enrichment": round(fold, 3), "odds_ratio": round(orr, 3),
                "or_ci_low": round(ci_low, 3), "or_ci_high": round(ci_high, 3), "p_value": p,
            })

            print(f"  {sig}: n={n}, k={k}, expected={expected:.2f}, fold={fold:.2f}, "
                  f"OR={orr:.2f} [{ci_low:.2f}-{ci_high:.2f}], p={p:.3e}", flush=True)

            perm_result = permutation_null(list(sig_in_U), target_in_U, U_expr_bins, k, N_PERM, rng)
            perm_result["region"] = label
            perm_result["signature"] = sig
            perm_result["n_signature_in_U"] = n
            permutation_rows.append(perm_result)
            print(f"    permutation null (n={N_PERM}): expected={perm_result['null_expected_overlap']}, "
                  f"fold_vs_null={perm_result['fold_enrichment_vs_null']}, "
                  f"empirical_p={perm_result['empirical_p']:.4f}, Z={perm_result['null_z_score']}", flush=True)

    pvals = [r["p_value"] for r in enrichment_rows]
    fdrs = bh_fdr(pvals)
    for r, fdr in zip(enrichment_rows, fdrs):
        r["fdr"] = fdr

    enrichment_path = f"{OUT_DIR}/mike_verzi_gut_hypergeometric_enrichment.tsv"
    with open(enrichment_path, "w", newline="") as f:
        cols = ["region", "signature", "n_signature_in_U", "K_target_in_U", "k_overlap", "N_universe",
                "expected", "fold_enrichment", "odds_ratio", "or_ci_low", "or_ci_high", "p_value", "fdr"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in enrichment_rows:
            w.writerow(r)
    print(f"\nWrote {enrichment_path}")

    permutation_path = f"{OUT_DIR}/mike_verzi_gut_permutation_null.tsv"
    with open(permutation_path, "w", newline="") as f:
        cols = ["region", "signature", "n_signature_in_U", "observed_overlap", "null_expected_overlap",
                "null_std_overlap", "fold_enrichment_vs_null", "empirical_p", "null_z_score"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in permutation_rows:
            w.writerow(r)
    print(f"Wrote {permutation_path}")

    with open(f"{OUT_DIR}/validation_run_metadata.json", "w") as f:
        json.dump({"n_permutations": N_PERM, "n_expression_bins": N_BINS, "rng_seed": RNG_SEED,
                    "n_tests": len(enrichment_rows)}, f, indent=2)


if __name__ == "__main__":
    main()
