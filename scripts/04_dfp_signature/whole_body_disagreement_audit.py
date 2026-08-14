#!/usr/bin/env python3
"""
GTEx vs HPA whole-body adult_excluded disagreement audit, per the PR #13
reviewer's explicit request: the final P-developmental assembly used GTEx
as an arbitrary "primary" reference and only reported HPA alongside, but
GTEx-pass (158) and HPA-pass (160) whole-body sets only overlap ~84 genes
-- real platform sensitivity that must be explained before freezing, not
glossed over by picking one platform as primary without justification.

Categorizes the 1,007 replicated_in_placenta genes into 4 buckets at
pct_cut=25, quorum=all_but_1:
  - both_pass:  passes whole-body exclusion in both GTEx and HPA
  - gtex_only:  passes in GTEx, fails in HPA
  - hpa_only:   passes in HPA, fails in GTEx
  - both_fail:  fails in both

Cross-checks each bucket against the real official HPA 65-gene placenta
"Tissue enriched" classification (results/04_dfp_signature/
hpa_placenta_official_tissue_enriched_65genes.tsv, from PR #10) and the 9
canonical trophoblast/placental markers, to understand what's actually
driving the gtex_only/hpa_only discordant genes before deciding whether a
combination rule is justified or the safest floor is the both-pass set.

Usage: python3 whole_body_disagreement_audit.py
(run on Argos, argos-codex env)
"""
import pandas as pd

ROOT = "/home/zz950/TWEAKR-OncoPlacental"
EDGE_DIR = f"{ROOT}/results/04_dfp_signature/edgeR"
OFFICIAL_65 = f"{ROOT}/results/04_dfp_signature/hpa_placenta_official_tissue_enriched_65genes.tsv"
OUT_DIR = f"{ROOT}/results/04_dfp_signature/dfp_gene_sets"

exec(open(f"{ROOT}/scripts/04_dfp_signature/build_dfp_gene_sets.py").read().split("# ============ 3. F-developmental")[0])

PCT_CUT = 25
ALLOWED_FAIL = 1  # all_but_1

fails_gtex = pd.DataFrame(False, index=gtex_not_detected.index, columns=tissue_cols)
for c in tissue_cols:
    fails_gtex[c] = (~gtex_not_detected[c]) & (gtex_pct[c] > PCT_CUT)
fails_hpa = pd.DataFrame(False, index=hpa_not_detected.index, columns=hpa_non_placenta)
for c in hpa_non_placenta:
    fails_hpa[c] = (~hpa_not_detected[c]) & (hpa_pct[c] > PCT_CUT)

pass_gtex = set(fails_gtex.index[fails_gtex.sum(axis=1) <= ALLOWED_FAIL])
pass_hpa = set(fails_hpa.index[fails_hpa.sum(axis=1) <= ALLOWED_FAIL])

# Coverage-aware bucketing (same fix as organ-matched adult_excluded_mask,
# PR #12 round 2): a gene only in fails_gtex.index/fails_hpa.index if that
# platform's gene panel actually includes it. Genes covered by only one
# platform are NOT a real cross-platform disagreement -- they're a
# coverage gap, same category as the H19 case, must not be conflated with
# genes both platforms measure but disagree on.
rip = replicated_in_placenta  # 1,007 genes, from the exec'd script
gtex_covered = set(fails_gtex.index)
hpa_covered = set(fails_hpa.index)

both_covered = rip & gtex_covered & hpa_covered
gtex_only_covered = (rip & gtex_covered) - hpa_covered  # HPA doesn't measure this gene at all
hpa_only_covered = (rip & hpa_covered) - gtex_covered   # GTEx doesn't measure this gene at all
neither_covered = rip - gtex_covered - hpa_covered

print(f"replicated_in_placenta total: {len(rip)}")
print(f"Coverage: both platforms measure {len(both_covered)}, GTEx-only-covered (HPA doesn't measure) {len(gtex_only_covered)}, "
      f"HPA-only-covered (GTEx doesn't measure) {len(hpa_only_covered)}, neither {len(neither_covered)}")

# Within the genes BOTH platforms actually measure, the real disagreement:
both_pass = both_covered & pass_gtex & pass_hpa
gtex_pass_hpa_fail = both_covered & pass_gtex - pass_hpa   # true disagreement
hpa_pass_gtex_fail = both_covered & pass_hpa - pass_gtex   # true disagreement
both_fail = both_covered - pass_gtex - pass_hpa

print(f"\nOf the {len(both_covered)} genes BOTH platforms measure (the only valid disagreement population):")
print(f"both_pass: {len(both_pass)}, GTEx-pass/HPA-fail (real disagreement): {len(gtex_pass_hpa_fail)}, "
      f"HPA-pass/GTEx-fail (real disagreement): {len(hpa_pass_gtex_fail)}, both_fail: {len(both_fail)}")
assert len(both_pass) + len(gtex_pass_hpa_fail) + len(hpa_pass_gtex_fail) + len(both_fail) == len(both_covered)

# Coverage-gap-only genes (not real disagreement, just measured by one
# platform): pass/fail per that platform's own evidence.
gtex_only_pass = gtex_only_covered & pass_gtex
gtex_only_fail = gtex_only_covered - pass_gtex
hpa_only_pass = hpa_only_covered & pass_hpa
hpa_only_fail = hpa_only_covered - pass_hpa
print(f"\nOf the {len(gtex_only_covered)} genes only GTEx measures: {len(gtex_only_pass)} pass, {len(gtex_only_fail)} fail (on GTEx's own evidence)")
print(f"Of the {len(hpa_only_covered)} genes only HPA measures: {len(hpa_only_pass)} pass, {len(hpa_only_fail)} fail (on HPA's own evidence)")

# Keep the old bucket names for the rest of this script (cross-checks
# below), now correctly scoped to genuine cross-platform disagreement only.
gtex_only = gtex_pass_hpa_fail
hpa_only = hpa_pass_gtex_fail

official_65 = set(pd.read_csv(OFFICIAL_65, sep="\t")["Gene"])
markers = ["ERVFRD-1", "CGA", "CSH1", "CSH2", "PSG1", "PSG3", "GATA3", "KRT7", "HLA-G"]

print("\n=== Cross-check against HPA official 65-gene placenta set ===")
for name, bucket in [("both_pass", both_pass), ("gtex_only", gtex_only), ("hpa_only", hpa_only), ("both_fail", both_fail)]:
    overlap = bucket & official_65
    print(f"{name}: {len(bucket)} genes, {len(overlap)} in HPA-official-65 ({100*len(overlap)/max(len(bucket),1):.1f}%): {sorted(overlap)}")

print("\n=== Marker bucket membership ===")
for m in markers:
    for name, bucket in [("both_pass", both_pass), ("gtex_only", gtex_only), ("hpa_only", hpa_only), ("both_fail", both_fail)]:
        if m in bucket:
            print(f"{m}: {name}")
            break
    else:
        print(f"{m}: not in replicated_in_placenta at all")

print(f"\n=== gtex_only genes ({len(gtex_only)}) ===")
print(sorted(gtex_only))
print(f"\n=== hpa_only genes ({len(hpa_only)}) ===")
print(sorted(hpa_only))

# For the discordant genes, show WHY: which specific tissues are failing
print("\n=== Why gtex_only genes fail HPA: which HPA tissues, and how many ===")
for g in sorted(gtex_only)[:15]:  # sample first 15 for inspection
    failing_hpa_tissues = [c for c in hpa_non_placenta if fails_hpa.loc[g, c]] if g in fails_hpa.index else ["NOT_IN_HPA_PANEL"]
    print(f"{g}: fails in {len(failing_hpa_tissues)} HPA tissues: {failing_hpa_tissues[:5]}{'...' if len(failing_hpa_tissues)>5 else ''}")

print("\n=== Why hpa_only genes fail GTEx: which GTEx tissues, and how many ===")
for g in sorted(hpa_only)[:15]:
    failing_gtex_tissues = [c for c in tissue_cols if fails_gtex.loc[g, c]] if g in fails_gtex.index else ["NOT_IN_GTEX_PANEL"]
    print(f"{g}: fails in {len(failing_gtex_tissues)} GTEx tissues: {failing_gtex_tissues[:5]}{'...' if len(failing_gtex_tissues)>5 else ''}")

# Write outputs
for name, bucket in [("both_pass_84", both_pass), ("gtex_only", gtex_only), ("hpa_only", hpa_only), ("both_fail", both_fail)]:
    with open(f"{OUT_DIR}/whole_body_{name}.txt", "w") as f:
        f.write("\n".join(sorted(bucket)) + "\n")
print(f"\nWrote bucket files to {OUT_DIR}/whole_body_{{both_pass_84,gtex_only,hpa_only,both_fail}}.txt")
