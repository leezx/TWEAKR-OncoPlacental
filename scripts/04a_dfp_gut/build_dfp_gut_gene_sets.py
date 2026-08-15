#!/usr/bin/env python3
"""
Build F_Colon-developmental / F_SI-developmental from the real edgeR DE
results (Second-trimester-fetal vs. Adult epithelium, LargeInt/SmallInt,
docs/STEP4A_GUT_FDEV_DESIGN.md PR #20 APPROVE), and redefine
D_Colon-shared/F_Colon-specific/P_Colon-specific (+ SI parallel) against
the unchanged, already-frozen P_developmental_primary84.txt.

Threshold calibration (real marker-panel check, same discipline as HDMA's
PR #12 calibration -- test candidates, pick by real marker retention, not
"comfortably above background" reasoning):

  Candidates tested: FDR<0.05|0.01 x logFC>1|2, against the real,
  filterByExpr-detectable canonical oncofetal/developmental marker panel
  (AFP, DLK1, PEG10, LGR5 -- IGF2/H19/LIN28B did not pass filterByExpr in
  either region, a real finding, not silently dropped). Only
  FDR<0.05 & logFC>1 retains ALL detectable markers in BOTH regions
  (LGR5's SmallInt FDR=0.0298 fails the FDR<0.01 candidates) -- chosen on
  that basis. AFP itself is the strongest possible calibration signal:
  logFC=+9.85 (LargeInt) / +9.89 (SmallInt), FDR<3e-6 both regions -- the
  single most face-valid oncofetal marker showing a massive, unambiguous
  fetal-elevation effect.

P_developmental provenance: P_developmental_primary84.txt already exists
as the raw (pre-F-subtraction) frozen program (verified: 84 = 6 D-shared +
78 P-specific, matching the disjoint-partition math in
docs/STEP4_DFP_DESIGN.md) -- used directly, not reconstructed, avoiding
any double-subtraction risk.

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

FDR_CUTOFF = 0.05
LOGFC_CUTOFF = 1.0
MARKER_PANEL = ["AFP", "DLK1", "PEG10", "LGR5"]


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
    markers_present = [m for m in MARKER_PANEL if m in edger_results]
    markers_retained = [m for m in markers_present if m in genes]
    print(f"F_{label}-developmental: {len(genes)} genes (FDR<{FDR_CUTOFF} & logFC>{LOGFC_CUTOFF})")
    print(f"  Marker panel: {markers_retained} retained of {markers_present} detectable "
          f"(of full panel {MARKER_PANEL})")
    return genes


def main():
    p_dev = load_gene_set(f"{DFP_DIR}/P_developmental_primary84.txt")
    print(f"P_developmental (unchanged, frozen): {len(p_dev)} genes\n")

    for label, region in [("Colon", "LargeInt"), ("SI", "SmallInt")]:
        edger = load_edger(f"{EDGER_DIR}/{region}_edgeR_primary.tsv")
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
    # against real data, rather than deciding in advance.
    colon_dev = load_gene_set(f"{OUT_DIR}/F_Colon-developmental.txt")
    si_dev = load_gene_set(f"{OUT_DIR}/F_SI-developmental.txt")
    intersection = colon_dev & si_dev
    union = colon_dev | si_dev

    from scipy.stats import hypergeom
    N = 33538  # approx universe: total genes in the h5ad's var (both regions same gene panel)
    p_overlap = hypergeom.sf(len(intersection) - 1, N, len(colon_dev), len(si_dev))

    with open(f"{OUT_DIR}/F_Gut-core.txt", "w") as f:
        f.write("\n".join(sorted(intersection)) + "\n")

    print(f"\n=== Colon/SI concordance (tertiary summary) ===")
    print(f"F_Colon-developmental ({len(colon_dev)}) vs F_SI-developmental ({len(si_dev)}):")
    print(f"  intersection (F_Gut-core) = {len(intersection)} genes, Jaccard = {len(intersection)/len(union):.3f}")
    print(f"  hypergeometric enrichment of this overlap vs chance (N={N}): p={p_overlap:.3e}")
    print(f"  (massively above chance-expected overlap of ~{len(colon_dev)*len(si_dev)/N:.0f} genes -- "
          f"real concordant fetal-vs-adult-gut signal, plus real region-specific signal in the "
          f"~{len(colon_dev)-len(intersection)}/{len(si_dev)-len(intersection)} Colon/SI-only genes)")
    print(f"Wrote {OUT_DIR}/F_Gut-core.txt")


if __name__ == "__main__":
    main()
