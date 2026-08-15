#!/usr/bin/env python3
"""
Real threshold-calibration table, per PR #21 REQUEST_CHANGES blocker #1:
the previous claim "Only FDR<0.05 & logFC>1 retains every detectable
marker" was never actually computed as a saved artifact -- the code
hardcoded FDR_CUTOFF=0.05/LOGFC_CUTOFF=1.0 directly -- and the claim itself
is false: every detected marker has logFC well above 2, so
FDR<0.05 & logFC>2 retains them all too. This script actually builds the
4-candidate table (now with 5 markers, IGF2 added via IGF2-1 -- see
gene_id_audit.py / PR #21 review: IGF2-1 is the real, HGNC-confirmed IGF2,
BioMart-verified ENSG00000167244; the unsuffixed "IGF2" var_name is a
different, non-canonically-named "novel protein" locus that happens to
share the symbol string in this dataset's annotation).

Marker panel (real, filterByExpr-detectable, gene-ID-verified):
  AFP, DLK1 (LargeInt only), PEG10, LGR5, IGF2 (via IGF2-1 symbol)

Tie-break rule (stated explicitly, decided BEFORE looking at which
candidate "wins" anything beyond marker retention -- avoiding the same
outcome-dependent-selection mistake flagged in PR #17 round 3):
  1. Prefer candidates that retain 100% of the detectable marker panel
     (a real, pre-registered developmental-biology sanity check).
  2. Among those, prefer the LARGER gene set (logFC>1 over logFC>2) --
     this is a genuine sensitivity-oriented choice for maximizing
     developmental-program recall, stated explicitly as such, not
     justified by a false "logFC>2 loses markers" claim.

Usage: python3 build_threshold_calibration.py <out_dir>
"""
import csv
import os

REPO = "/home/zz950/TWEAKR-OncoPlacental"
EDGER_DIR = os.path.join(REPO, "results/04a_dfp_gut/edgeR")
OUT_DIR = os.path.join(REPO, "results/04a_dfp_gut/dfp_gut_gene_sets")
os.makedirs(OUT_DIR, exist_ok=True)

# real symbol as it appears in this dataset's edgeR output (IGF2 -> IGF2-1,
# per the gene-ID audit -- see docstring)
MARKER_PANEL = {"AFP": "AFP", "DLK1": "DLK1", "PEG10": "PEG10", "LGR5": "LGR5", "IGF2": "IGF2-1"}

CANDIDATES = [
    ("FDR<0.05 & logFC>1", 0.05, 1.0),
    ("FDR<0.05 & logFC>2", 0.05, 2.0),
    ("FDR<0.01 & logFC>1", 0.01, 1.0),
    ("FDR<0.01 & logFC>2", 0.01, 2.0),
]


def load_edger(path):
    with open(path) as f:
        return {r["gene"]: (float(r["logFC"]), float(r["FDR"])) for r in csv.DictReader(f, delimiter="\t")}


def main():
    rows = []
    for region, file in [("LargeInt", "LargeInt_edgeR_primary.tsv"), ("SmallInt", "SmallInt_edgeR_primary.tsv")]:
        edger = load_edger(f"{EDGER_DIR}/{file}")
        markers_present = {display: symbol for display, symbol in MARKER_PANEL.items() if symbol in edger}
        print(f"\n=== {region}: detectable markers = {list(markers_present.keys())} ===")
        for display, symbol in markers_present.items():
            lfc, fdr = edger[symbol]
            print(f"  {display} ({symbol}): logFC={lfc:.2f} FDR={fdr:.2e}")

        for label, fdr_cut, lfc_cut in CANDIDATES:
            n_pass = sum(1 for lfc, fdr in edger.values() if fdr < fdr_cut and lfc > lfc_cut)
            retained = [d for d, s in markers_present.items() if edger[s][1] < fdr_cut and edger[s][0] > lfc_cut]
            rows.append({
                "region": region, "candidate": label, "fdr_cutoff": fdr_cut, "logfc_cutoff": lfc_cut,
                "n_genes": n_pass, "n_markers_detectable": len(markers_present),
                "n_markers_retained": len(retained), "markers_retained": ";".join(retained),
                "all_markers_retained": len(retained) == len(markers_present),
            })

    out_path = f"{OUT_DIR}/threshold_calibration.tsv"
    with open(out_path, "w", newline="") as f:
        cols = ["region", "candidate", "fdr_cutoff", "logfc_cutoff", "n_genes",
                "n_markers_detectable", "n_markers_retained", "markers_retained", "all_markers_retained"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nWrote {out_path}")

    print("\n=== Candidates retaining ALL detectable markers in BOTH regions ===")
    by_candidate = {}
    for r in rows:
        by_candidate.setdefault(r["candidate"], []).append(r)
    for label, _, _ in CANDIDATES:
        rs = by_candidate[label]
        all_ok = all(r["all_markers_retained"] for r in rs)
        sizes = {r["region"]: r["n_genes"] for r in rs}
        print(f"  {label}: all_markers_retained_both_regions={all_ok}, gene set sizes={sizes}")


if __name__ == "__main__":
    main()
