#!/usr/bin/env python3
"""
Gene-set enrichment (not just raw overlap) of the 5 mike_verzi normal-context
fetal/revival signatures against D-shared / F-specific (global + 7 lineage
modules) / P-specific, per user request: raw overlap counts aren't
comparable across D/F/P targets that differ ~400x in size (D-shared=6 vs
F-specific=2,504) -- a hypergeometric enrichment score (fold-enrichment +
BH-FDR-corrected p-value) is the properly size-normalized version of the
same question, and is the standard statistic for this kind of gene-set
comparison (same category of test as GSEA/Fisher's-exact overlap analysis).

Background gene universe: 23,272 human protein-coding genes, queried
directly from Ensembl BioMart (gene_biotype=protein_coding, GRCh38,
jun2026 archive) -- not assumed from a textbook figure.

For each (signature, target) pair:
  k = observed overlap (primary one2one human genes only, same set used
      in the raw-overlap audit)
  n = signature's primary human gene-set size
  K = target's gene-set size
  N = background (23,272)
  fold_enrichment = (k/n) / (K/N)
  p = hypergeometric survival function P(X >= k), one-sided over-representation
  fdr = BH-corrected p across all signature x target tests (55 tests:
        5 signatures x 11 targets [D-shared, F-specific-global, 7 lineage, P-specific])

Usage: python3 mike_verzi_dfp_enrichment.py
"""
import csv
import json
import os
from scipy.stats import hypergeom
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT_DIR = os.path.join(REPO, "results/06a_normal_context")
DFP_DIR = os.path.join(REPO, "results/04_dfp_signature/dfp_gene_sets")
ORGANS = ["Adrenal", "Liver", "Skin", "Spleen", "Stomach", "Thymus", "Thyroid"]

BACKGROUND_N = 23272  # human protein-coding genes, Ensembl BioMart, verified this session


def load_gene_set(path):
    with open(path) as f:
        return {l.strip() for l in f if l.strip()}


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


def main():
    sets = json.load(open(f"{OUT_DIR}/mike_verzi_sets_raw.json"))

    # reload primary human sets (written by build_mike_verzi_human_final.py)
    primary_human_sets = {}
    for sig in sets:
        with open(f"{OUT_DIR}/{sig}_human_primary.txt") as f:
            primary_human_sets[sig] = {l.strip() for l in f if l.strip()}

    d_shared = load_gene_set(f"{DFP_DIR}/D_shared_FINAL.txt")
    f_specific = load_gene_set(f"{DFP_DIR}/F_specific_FINAL.txt")
    p_specific = load_gene_set(f"{DFP_DIR}/P_specific_FINAL.txt")
    f_lineage = {}
    for organ in ORGANS:
        organ_set = load_gene_set(f"{DFP_DIR}/F_developmental_{organ}.txt")
        f_lineage[organ] = organ_set & f_specific

    targets = {"D-shared": d_shared, "F-specific-global": f_specific}
    for organ in ORGANS:
        targets[f"F-lineage-{organ}"] = f_lineage[organ]
    targets["P-specific"] = p_specific

    results = []
    for sig, human_set in primary_human_sets.items():
        n = len(human_set)
        for tname, tset in targets.items():
            K = len(tset)
            k = len(human_set & tset)
            expected = n * K / BACKGROUND_N
            fold = (k / n) / (K / BACKGROUND_N) if k > 0 else 0.0
            # hypergeom.sf(k-1, N, K, n) = P(X >= k)
            p = hypergeom.sf(k - 1, BACKGROUND_N, K, n) if k > 0 else 1.0
            results.append({
                "signature": sig, "target": tname, "n_signature": n,
                "K_target": K, "k_overlap": k, "expected": round(expected, 3),
                "fold_enrichment": round(fold, 3), "p_value": p,
            })

    pvals = [r["p_value"] for r in results]
    fdrs = bh_fdr(pvals)
    for r, fdr in zip(results, fdrs):
        r["fdr"] = fdr
        r["neg_log10_fdr"] = round(-np.log10(max(fdr, 1e-300)), 3)

    out_path = f"{OUT_DIR}/mike_verzi_dfp_enrichment.tsv"
    with open(out_path, "w", newline="") as f:
        cols = ["signature", "target", "n_signature", "K_target", "k_overlap",
                "expected", "fold_enrichment", "p_value", "fdr", "neg_log10_fdr"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in results:
            w.writerow(r)

    print(f"Background: {BACKGROUND_N} human protein-coding genes (Ensembl BioMart)")
    print(f"Total tests: {len(results)} (5 signatures x 11 targets)")
    print(f"Wrote {out_path}")

    # console summary: per signature, which target has max fold-enrichment
    print("\nPer-signature top enrichment target (by fold-enrichment among FDR<0.05):")
    from collections import defaultdict
    by_sig = defaultdict(list)
    for r in results:
        by_sig[r["signature"]].append(r)
    for sig, rs in by_sig.items():
        sig_results = [r for r in rs if r["fdr"] < 0.05]
        sig_results.sort(key=lambda r: -r["fold_enrichment"])
        top = sig_results[:3]
        print(f"  {sig}:")
        for r in top:
            print(f"    {r['target']}: fold={r['fold_enrichment']}, k={r['k_overlap']}/{r['n_signature']}, FDR={r['fdr']:.2e}")


if __name__ == "__main__":
    main()
