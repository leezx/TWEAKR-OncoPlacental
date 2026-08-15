#!/usr/bin/env python3
"""
Build F_Colon-developmental / F_SI-developmental from the real edgeR DE
results (Second-trimester-fetal vs. Adult epithelium, LargeInt/SmallInt,
docs/STEP4A_GUT_FDEV_DESIGN.md PR #20 APPROVE), and redefine
D_Colon-shared/F_Colon-specific/P_Colon-specific (+ SI parallel) against
the unchanged, already-frozen P_developmental_primary84.txt.

v2 -- fixes 3 real blockers caught in PR #21 round-1 REQUEST_CHANGES:

1. Threshold calibration claim was mathematically false. The real
   calibration table is now built by build_threshold_calibration.py
   (results/04a_dfp_gut/dfp_gut_gene_sets/threshold_calibration.tsv) --
   every detectable marker's logFC is well above 2, so FDR<0.05 & logFC>2
   retains the full panel too, not just logFC>1. Explicit, pre-stated
   tie-break rule (decided before looking at which candidate "wins",
   avoiding the outcome-dependent-selection mistake flagged in PR #17
   round 3): among candidates retaining 100% of the marker panel, prefer
   the larger gene set -- a genuine sensitivity-oriented recall choice,
   not a false "only this one works" claim. FDR<0.05 & logFC>1 wins on
   that stated basis (1,456 genes both regions), not because logFC>2
   "fails" (it doesn't).

2. Marker panel was missing IGF2 due to a real gene-ID bug (see
   gene_id_audit.py): the h5ad has both "IGF2" (-> ENSG00000284779, a
   BioMart-confirmed "novel protein", NOT HGNC-curated IGF2) and "IGF2-1"
   (-> ENSG00000167244, HGNC:5466, the real canonical IGF2) as distinct
   var_names -- anndata's duplicate-symbol uniquification. IGF2-1 is a
   massive fetal-up signal (logFC=12.53 LargeInt, 6.59 SmallInt) that the
   naive "IGF2" string lookup missed entirely. Verified via BioMart
   (jun2026.archive.ensembl.org), not assumed. Checked: 0/84
   P_developmental genes collide with any of the 104 duplicate-symbol
   base names in this h5ad, so D_Colon-shared/D_SI-shared/P_*-specific
   set arithmetic itself was NOT affected by this bug -- only the marker
   calibration check was.

3. F_Gut-core's hypergeometric enrichment used N=33,538 (every gene in
   the h5ad var) as the background, but only ~15,149/~14,982 genes ever
   entered LargeInt/SmallInt's filterByExpr-passing set and were
   therefore even eligible to appear in either F-developmental set. Using
   the full unfiltered gene count as background artificially deflates the
   chance-expected overlap and inflates the apparent enrichment. Fixed:
   N = |common testable universe| = genes passing filterByExpr in BOTH
   regions. Also reports bidirectional overlap fraction and logFC
   concordance (Pearson/Spearman) on the shared filterByExpr genes, not
   just Jaccard + a single enrichment p-value.

Usage: python3 build_dfp_gut_gene_sets.py
(pure gene-list arithmetic on already-computed frozen edgeR results --
run locally per the project's narrow exception, no qsub needed)
"""
import csv
import os

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
EDGER_DIR = os.path.join(REPO, "results/04a_dfp_gut/edgeR")
DFP_DIR = os.path.join(REPO, "results/04_dfp_signature/dfp_gene_sets")
OUT_DIR = os.path.join(REPO, "results/04a_dfp_gut/dfp_gut_gene_sets")
os.makedirs(OUT_DIR, exist_ok=True)

# Chosen per the real calibration table (build_threshold_calibration.py) and
# the explicit tie-break rule stated in the module docstring above.
FDR_CUTOFF = 0.05
LOGFC_CUTOFF = 1.0
# real symbol as it appears in this dataset's edgeR output (IGF2 -> IGF2-1,
# per the gene-ID audit)
MARKER_PANEL = {"AFP": "AFP", "DLK1": "DLK1", "PEG10": "PEG10", "LGR5": "LGR5", "IGF2": "IGF2-1"}


def load_edger(path):
    with open(path) as f:
        return {r["gene"]: (float(r["logFC"]), float(r["FDR"])) for r in csv.DictReader(f, delimiter="\t")}


def load_gene_set(path):
    with open(path) as f:
        return {l.strip() for l in f if l.strip()}


def build_f_developmental(edger_results, label):
    genes = {g for g, (lfc, fdr) in edger_results.items() if fdr < FDR_CUTOFF and lfc > LOGFC_CUTOFF}
    with open(f"{OUT_DIR}/F_{label}-developmental.txt", "w") as f:
        f.write("\n".join(sorted(genes)) + "\n")
    markers_present = {d: s for d, s in MARKER_PANEL.items() if s in edger_results}
    markers_retained = [d for d, s in markers_present.items() if s in genes]
    print(f"F_{label}-developmental: {len(genes)} genes (FDR<{FDR_CUTOFF} & logFC>{LOGFC_CUTOFF})")
    print(f"  Marker panel: {markers_retained} retained of {list(markers_present)} detectable "
          f"(of full panel {list(MARKER_PANEL)})")
    return genes


def main():
    p_dev = load_gene_set(f"{DFP_DIR}/P_developmental_primary84.txt")
    print(f"P_developmental (unchanged, frozen): {len(p_dev)} genes\n")

    edger_by_region = {}
    for label, region in [("Colon", "LargeInt"), ("SI", "SmallInt")]:
        edger = load_edger(f"{EDGER_DIR}/{region}_edgeR_primary.tsv")
        edger_by_region[label] = edger
        f_dev = build_f_developmental(edger, label)

        d_shared = f_dev & p_dev
        f_specific = f_dev - p_dev
        p_specific = p_dev - f_dev

        for name, geneset in [(f"D_{label}-shared", d_shared),
                               (f"F_{label}-specific", f_specific),
                               (f"P_{label}-specific", p_specific)]:
            with open(f"{OUT_DIR}/{name}.txt", "w") as fh:
                fh.write("\n".join(sorted(geneset)) + "\n")

        print(f"  D_{label}-shared = F_{label}-developmental ^ P_developmental = {len(d_shared)} genes: {sorted(d_shared)}")
        print(f"  F_{label}-specific = F_{label}-developmental \\ P_developmental = {len(f_specific)} genes")
        print(f"  P_{label}-specific = P_developmental \\ F_{label}-developmental = {len(p_specific)} genes")
        print()

    print(f"Wrote all gene sets to {OUT_DIR}/")

    # --- tertiary summary: Colon/SI concordance + Gut-core (intersection) ---
    # Per design doc: "A pooled F_Gut-developmental ... is reported only as a
    # tertiary summary, and only if Colon and SI turn out concordant enough
    # for pooling to be meaningful -- not assumed up front." Check that now,
    # against real data, with the CORRECT background universe (fix #3 above).
    colon_dev = load_gene_set(f"{OUT_DIR}/F_Colon-developmental.txt")
    si_dev = load_gene_set(f"{OUT_DIR}/F_SI-developmental.txt")
    intersection = colon_dev & si_dev
    union = colon_dev | si_dev

    common_testable = set(edger_by_region["Colon"]) & set(edger_by_region["SI"])
    N = len(common_testable)
    K = len(colon_dev & common_testable)  # should equal len(colon_dev) since F sets are subsets of their own filterByExpr set, but colon_dev genes not in SI's testable set can't count toward SI-side chance calc -- use the common universe consistently on both sides
    n = len(si_dev & common_testable)
    k = len(intersection & common_testable)

    from scipy.stats import hypergeom
    p_overlap = hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0
    expected = K * n / N if N else 0.0
    fold = (k / n) / (K / N) if (k > 0 and n > 0 and K > 0) else 0.0

    with open(f"{OUT_DIR}/F_Gut-core.txt", "w") as f:
        f.write("\n".join(sorted(intersection)) + "\n")

    # bidirectional overlap fraction + logFC concordance on the common testable genes
    frac_of_colon = len(intersection) / len(colon_dev)
    frac_of_si = len(intersection) / len(si_dev)
    common_list = sorted(common_testable)
    colon_lfc = [edger_by_region["Colon"][g][0] for g in common_list]
    si_lfc = [edger_by_region["SI"][g][0] for g in common_list]
    from scipy.stats import pearsonr, spearmanr
    pear_r, _ = pearsonr(colon_lfc, si_lfc)
    spear_r, _ = spearmanr(colon_lfc, si_lfc)

    print(f"\n=== Colon/SI concordance (tertiary summary, corrected universe) ===")
    print(f"F_Colon-developmental ({len(colon_dev)}) vs F_SI-developmental ({len(si_dev)}):")
    print(f"  intersection (F_Gut-core) = {len(intersection)} genes, Jaccard = {len(intersection)/len(union):.3f}")
    print(f"  bidirectional overlap: {frac_of_colon:.3f} of Colon set, {frac_of_si:.3f} of SI set")
    print(f"  common testable universe (filterByExpr-passing in BOTH regions): N={N}")
    print(f"  hypergeometric on common universe: K(Colon-in-universe)={K}, n(SI-in-universe)={n}, "
          f"k(overlap-in-universe)={k}, expected={expected:.1f}, fold={fold:.2f}, p={p_overlap:.3e}")
    print(f"  logFC concordance on the {len(common_list)} common testable genes: "
          f"Pearson r={pear_r:.3f}, Spearman rho={spear_r:.3f}")
    print(f"Wrote {OUT_DIR}/F_Gut-core.txt")


if __name__ == "__main__":
    main()
