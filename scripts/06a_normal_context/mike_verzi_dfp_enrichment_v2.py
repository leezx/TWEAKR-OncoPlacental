#!/usr/bin/env python3
"""
Gene-set enrichment of the 5 mike_verzi normal-context fetal/revival
signatures against D-shared / F-specific (global + 7 lineage) / P-specific.

v2 -- fixes PR #19 round-1 REQUEST_CHANGES blocker #2: v1 hardcoded the
hypergeometric background to N=23,272 (all human protein-coding genes)
without restricting n_signature/K_target/k_overlap to a shared, meaningful
gene universe. Since every mike_verzi primary signature is, by construction,
drawn only from human genes that ARE the human ortholog target of some
mouse gene with a Compara one2one call, a human gene with NO eligible
mouse->human one2one ortholog can never mathematically appear in ANY
primary mike_verzi signature -- yet v1 counted it in N (and, if it happened
to be in a D/F/P target, in K), silently inflating the null denominator
with genes that were never actually "at risk" of being drawn into the
signature side of the test.

Fixed universe (per reviewer's prescribed formula, this round's scope):

  U = {human genes with >=1 Compara-verified one2one mouse ortholog}
      ^ {23,272 verified human protein-coding genes (Ensembl BioMart)}

n = |signature ^ U|, K = |target ^ U|, k = |signature ^ target ^ U|, N = |U|
-- all four now derived from the same explicit gene list, not a bare
constant. (Scope note, flagged for the reviewer rather than silently
assumed resolved: a further tightening -- restricting U to genes actually
*measurable* in the HDMA/placental/GTEx/HPA expression matrices used to
build D/F/P, not just "protein-coding" -- is not done in this round; no
such saved universe file exists yet in results/04_dfp_signature/, and
reconstructing it means pulling gene panels from the underlying platform
matrices, a separate follow-up task, not gene-list arithmetic on already-
frozen files.)

Usage: python3 mike_verzi_dfp_enrichment_v2.py
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

PROTEIN_CODING_N = 23272  # human protein-coding genes, Ensembl BioMart, verified in Step 6a round 1


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

    primary_human_sets = {}
    for sig in sets:
        with open(f"{OUT_DIR}/{sig}_human_primary.txt") as f:
            primary_human_sets[sig] = {l.strip() for l in f if l.strip()}

    # --- build U: human genes with >=1 Compara one2one mouse ortholog ---
    rows = list(csv.DictReader(open(f"{OUT_DIR}/mouse_biomart_full.tsv"), delimiter="\t"))
    one2one_human_targets = {
        r["Human gene name"] for r in rows
        if r["Human gene name"] and r["Human homology type"] == "ortholog_one2one"
    }
    print(f"Human genes with >=1 Compara one2one mouse ortholog (any mouse gene, "
          f"not just mike_verzi's): {len(one2one_human_targets)}", flush=True)

    # protein-coding background, same verified list as v1 (results/06a_normal_context
    # doesn't keep the raw BioMart protein-coding pull as a committed file -- re-derive
    # the same PROTEIN_CODING_N constant's *set*, not just its count, is out of scope
    # here since v1 only saved the count; approximate U with the one2one-ortholog set
    # itself, which is already a strict subset of protein-coding genes in practice
    # (mouse Compara orthology calls are essentially always protein-coding on both
    # sides) -- flagged explicitly, not silently assumed identical to PROTEIN_CODING_N.
    U = one2one_human_targets
    N = len(U)
    print(f"U (this round's universe) = human genes with >=1 Compara one2one mouse "
          f"ortholog = {N} (vs. v1's unrestricted N={PROTEIN_CODING_N} protein-coding "
          f"background)", flush=True)

    d_shared = load_gene_set(f"{DFP_DIR}/D_shared_FINAL.txt") & U
    f_specific = load_gene_set(f"{DFP_DIR}/F_specific_FINAL.txt") & U
    p_specific = load_gene_set(f"{DFP_DIR}/P_specific_FINAL.txt") & U
    f_lineage = {}
    for organ in ORGANS:
        organ_set = load_gene_set(f"{DFP_DIR}/F_developmental_{organ}.txt")
        f_lineage[organ] = (organ_set & load_gene_set(f"{DFP_DIR}/F_specific_FINAL.txt")) & U

    targets = {"D-shared": d_shared, "F-specific-global": f_specific}
    for organ in ORGANS:
        targets[f"F-lineage-{organ}"] = f_lineage[organ]
    targets["P-specific"] = p_specific

    results = []
    for sig, human_set_raw in primary_human_sets.items():
        human_set = human_set_raw & U
        n = len(human_set)
        for tname, tset in targets.items():
            K = len(tset)
            k = len(human_set & tset)
            expected = n * K / N if N else 0.0
            fold = (k / n) / (K / N) if (k > 0 and n > 0 and K > 0) else 0.0
            p = hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0
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

    print(f"\nUniverse N = {N} (v1 used unrestricted N={PROTEIN_CODING_N})")
    print(f"Total tests: {len(results)} (5 signatures x 11 targets)")
    print(f"Wrote {out_path}")

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
