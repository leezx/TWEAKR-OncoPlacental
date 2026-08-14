#!/usr/bin/env python3
"""
Final D-shared / F-specific / P-specific gene set assembly, per
docs/STEP4_DFP_DESIGN.md. Single, re-runnable, single-source-of-truth
pipeline -- consolidates what were two separate scripts (an earlier
build_dfp_gene_sets.py using an unjustified GTEx-only-primary whole-body
choice, and a follow-up whole_body_disagreement_audit.py that derived the
real, reviewer-approved 84-gene primary tier) into one, per PR #13 round
2 review: having two scripts produce two different, inconsistent answers
was a real reproducibility gap -- rerunning "the" assembly script must
always reproduce the frozen numbers, not silently regress to a
superseded definition.

Pipeline:
1. P-developmental = replicated_in_placenta (frozen, PR #10) AND
   adult_excluded(whole_body). Whole-body evidence is coverage-aware
   (same fix as organ-matched adult_excluded_mask, PR #12): of the
   1,007 replicated_in_placenta genes, only those measured by BOTH GTEx
   and HPA are a valid disagreement population (815 genes) -- genes
   covered by only one platform are not silently penalized for the
   other platform's absence.
   - P_primary (84 genes): both GTEx and HPA agree the gene is
     whole-body adult-excluded (pct_cut=25, quorum=all_but_1). Highest
     confidence -- both bulk references agree.
   - P_extended (174 genes): P_primary plus the 90 genuine
     single-platform-disagreement genes (15 GTEx-pass/HPA-fail + 75
     HPA-pass/GTEx-fail) -- lower confidence, reported separately, never
     silently merged into the primary tier (PR #13 round 1 fix: neither
     disagreement direction was confidently explainable by a bespoke
     combination rule, so the reviewer's own fallback was used).
2. F-developmental = union of the 7 organs' frozen per-organ sets
   (elevated_pct=75, adult_excl_pct=25, quorum=1.0 -- PR #12). Union used
   as the primary set-logic input (consensus reported as a stricter
   sensitivity alternative, approved as-is by the PR #13 reviewer).
3. Final D/F/P = simple set operations between P_primary and F_union
   (STEP4_DFP_DESIGN.md): D-shared = P AND F; F-specific = F AND NOT P;
   P-specific = P AND NOT F.

Usage: python3 build_dfp_gene_sets.py
(run on Argos, argos-codex env)
"""
import os
import pandas as pd
from collections import Counter

ROOT = "/home/zz950/TWEAKR-OncoPlacental"
EDGE_DIR = f"{ROOT}/results/04_dfp_signature/edgeR"
OFFICIAL_65 = f"{ROOT}/results/04_dfp_signature/hpa_placenta_official_tissue_enriched_65genes.tsv"
OUT_DIR = f"{ROOT}/results/04_dfp_signature/dfp_gene_sets"
os.makedirs(OUT_DIR, exist_ok=True)

exec(open(f"{ROOT}/scripts/04_dfp_signature/f_developmental_calibration.py").read().split("# ---- HDMA per-organ elevated calibration ----")[0])

# ============ 1. P-developmental: replicated_in_placenta (frozen, PR #10) ============
a_df = pd.read_csv(f"{EDGE_DIR}/Arutyunyan_edgeR_results.tsv", sep="\t")
n_df = pd.read_csv(f"{EDGE_DIR}/Nature2026_edgeR_results.tsv", sep="\t")
a_pass = set(a_df.loc[(a_df.logFC >= 0.75) & (a_df.FDR < 0.05), "gene"])
n_pass = set(n_df.loc[(n_df.logFC >= 0.75) & (n_df.FDR < 0.05), "gene"])
replicated_in_placenta = a_pass & n_pass
print(f"replicated_in_placenta (frozen, PR #10): {len(replicated_in_placenta)} genes "
      f"(Arutyunyan {len(a_pass)}, Nature2026 {len(n_pass)})")
assert len(replicated_in_placenta) == 1007, f"replicated_in_placenta count drifted: {len(replicated_in_placenta)} != 1007 (frozen PR #10 input changed?)"

# ============ 2. Whole-body adult_excluded, coverage-aware (PR #13 round 1/2 fix) ============
PCT_CUT = 25
ALLOWED_FAIL = 1  # all_but_1
hpa_non_placenta = [c for c in hpa_pct.columns if c != "placenta"]

fails_gtex = pd.DataFrame(False, index=gtex_not_detected.index, columns=tissue_cols)
for c in tissue_cols:
    fails_gtex[c] = (~gtex_not_detected[c]) & (gtex_pct[c] > PCT_CUT)
fails_hpa = pd.DataFrame(False, index=hpa_not_detected.index, columns=hpa_non_placenta)
for c in hpa_non_placenta:
    fails_hpa[c] = (~hpa_not_detected[c]) & (hpa_pct[c] > PCT_CUT)

pass_gtex = set(fails_gtex.index[fails_gtex.sum(axis=1) <= ALLOWED_FAIL])
pass_hpa = set(fails_hpa.index[fails_hpa.sum(axis=1) <= ALLOWED_FAIL])

gtex_covered = set(fails_gtex.index)
hpa_covered = set(fails_hpa.index)
rip = replicated_in_placenta

both_covered = rip & gtex_covered & hpa_covered
both_pass = both_covered & pass_gtex & pass_hpa
gtex_pass_hpa_fail = (both_covered & pass_gtex) - pass_hpa
hpa_pass_gtex_fail = (both_covered & pass_hpa) - pass_gtex

P_primary = both_pass
P_extended = both_pass | gtex_pass_hpa_fail | hpa_pass_gtex_fail

print(f"\nWhole-body coverage-aware disagreement audit (pct={PCT_CUT}, quorum=all_but_1):")
print(f"  both-platform-covered: {len(both_covered)} of {len(rip)}")
print(f"  P_primary (both-pass): {len(P_primary)}")
print(f"  GTEx-pass/HPA-fail (real disagreement): {len(gtex_pass_hpa_fail)}")
print(f"  HPA-pass/GTEx-fail (real disagreement): {len(hpa_pass_gtex_fail)}")
print(f"  P_extended (primary + both disagreement directions): {len(P_extended)}")
assert len(P_primary) == 84, f"P_primary count drifted: {len(P_primary)} != 84"
assert len(P_extended) == 174, f"P_extended count drifted: {len(P_extended)} != 174"

# Marker validation
markers = ["ERVFRD-1", "CGA", "CSH1", "CSH2", "PSG1", "PSG3", "GATA3", "KRT7", "HLA-G"]
print("\nMarker check against P_primary (84 genes):")
for m in markers:
    print(f"  {m:10s} in P_primary: {m in P_primary}, in P_extended: {m in P_extended}")

official_65 = set(pd.read_csv(OFFICIAL_65, sep="\t")["Gene"])
overlap = P_primary & official_65
print(f"\nP_primary overlap with HPA official 65-gene placenta set: {len(overlap)}/{len(P_primary)} "
      f"({100*len(overlap)/len(P_primary):.1f}%, vs ~0.3% background)")

# ============ 3. F-developmental: real per-organ gene sets (frozen, PR #12) ============
ORGANS_LOCAL = ["Adrenal", "Thyroid", "Spleen", "Thymus", "Liver", "Skin", "Stomach"]
f_dev_sets = {}
for organ in ORGANS_LOCAL:
    counts = pd.read_csv(f"{HDMA_DIR}/{organ}_pseudobulk_counts.tsv", sep="\t", index_col=0)
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    nd, pct = detected_percentile(cpm)
    median_pct = pct.median(axis=1, skipna=True)
    detected_frac = (~nd).mean(axis=1)
    elevated = (median_pct >= 75) & (detected_frac >= 1.0)
    mask, prov = adult_excluded_mask(organ, 25)
    elevated_genes = set(elevated.index[elevated])
    excluded_genes = set(mask.index[mask]) if mask is not None else set()
    f_dev_sets[organ] = elevated_genes & excluded_genes
    print(f"\nF-developmental[{organ}]: {len(f_dev_sets[organ])} genes")

F_union = set().union(*f_dev_sets.values())
gene_organ_count = Counter()
for organ, genes in f_dev_sets.items():
    for g in genes:
        gene_organ_count[g] += 1
F_consensus = {g for g, c in gene_organ_count.items() if c >= 2}
print(f"\nF_union: {len(F_union)}, F_consensus(>=2 organs): {len(F_consensus)}")
assert len(F_union) == 2510, f"F_union count drifted: {len(F_union)} != 2510"

# ============ 4. Final D/F/P assembly (P_primary x F_union, per PR #13 review) ============
D_shared = P_primary & F_union
F_specific = F_union - P_primary
P_specific = P_primary - F_union
print(f"\n=== FINAL (P_primary x F_union) ===")
print(f"D-shared: {len(D_shared)}, F-specific: {len(F_specific)}, P-specific: {len(P_specific)}")
assert len(D_shared) == 6, f"D_shared count drifted: {len(D_shared)} != 6"
assert len(F_specific) == 2504, f"F_specific count drifted: {len(F_specific)} != 2504"
assert len(P_specific) == 78, f"P_specific count drifted: {len(P_specific)} != 78"
print(f"D-shared genes: {sorted(D_shared)}")

for m in markers:
    print(f"Marker check -- {m:10s}: P_primary={m in P_primary}, F_union={m in F_union}")

# ============ Write outputs (single source of truth -- old ambiguous
# GTEx-only-primary files from the superseded first draft removed) ============
for stale in ["P_developmental.txt", "D_shared_union.txt", "F_specific_union.txt", "P_specific_union.txt"]:
    stale_path = f"{OUT_DIR}/{stale}"
    if os.path.exists(stale_path):
        os.remove(stale_path)
        print(f"Removed stale/superseded file: {stale_path}")

outputs = {
    "P_developmental_primary84": P_primary,
    "P_developmental_extended174": P_extended,
    "F_developmental_union": F_union,
    "F_developmental_consensus2plus": F_consensus,
    "D_shared_FINAL": D_shared,
    "F_specific_FINAL": F_specific,
    "P_specific_FINAL": P_specific,
}
for name, geneset in outputs.items():
    path = f"{OUT_DIR}/{name}.txt"
    with open(path, "w") as f:
        f.write("\n".join(sorted(geneset)) + "\n")
    print(f"Wrote {path} ({len(geneset)} genes)")

for organ, genes in f_dev_sets.items():
    with open(f"{OUT_DIR}/F_developmental_{organ}.txt", "w") as f:
        f.write("\n".join(sorted(genes)) + "\n")

print("\n=== DONE: single-source-of-truth D/F/P assembly, all cardinalities asserted ===")
