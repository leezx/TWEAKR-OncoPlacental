# HPA trophoblast/placenta enrichment calibration — freezing the P-developmental effect-size cutoff

Per the PR #9 reviewer's explicit next step: "先把 0.75 作为 leading candidate、
再用独立的 HPA trophoblast/placenta enrichment 决定是否冻结". Independent
check against HPA's own bulk RNA consensus data (not derived from the
scRNA-seq datasets the edgeR DE itself used), to calibrate which of the
directional `logFC ≥ cutoff` candidates (0.5 / 0.75 / 1.0, all already
`FDR<0.05`) to freeze for `replicated_in_placenta`'s per-dataset pass
criterion.

## Method

- **HPA "placenta-enriched" reference set**: computed directly from
  `rna_tissue_hpa.tsv` (Tier-1 adult reference, `STEP3_METHOD_CONTRACT.md`;
  placenta used here strictly as a **positive cross-check**, per the `role`
  column fix in PR #5 — never as adult-negative reference). Uses HPA's own
  published "Tissue enriched" definition: nTPM in placenta ≥ 4× every other
  of HPA's 39 tissues, with a floor (placenta nTPM ≥ 1) to exclude
  near-zero-everywhere noise. Result: **82 genes** (of 20,151 total HPA
  genes) meet this bar — a small, high-confidence reference set, not the
  thing being validated (it's built from HPA's independent bulk data, not
  from Arutyunyan/Nature2026's scRNA-seq).
- **Enrichment test**: for each dataset's edgeR-tested gene universe
  (`filterByExpr`-passing genes) and each candidate cutoff, checked what
  fraction of the directional pass-set (`logFC≥cutoff & FDR<0.05`) falls in
  the 82-gene HPA-placenta-enriched set, vs. the background rate in the full
  universe. Reported as fold-enrichment and a one-sided hypergeometric
  p-value. Also computed for the 2-of-2 overlap sets (genes passing in both
  datasets).
- Run on Argos (`scripts/04_dfp_signature/hpa_placenta_enrichment.py`),
  results pulled back and verified byte-exact (md5) against Argos output.

## Results

| Cutoff | Arutyunyan fold (p) | Nature2026 fold (p) | 2-of-2 overlap fold (p) | 2-of-2 overlap gene count | Canonical marker retention (from `trophoblast_edgeR_audit.md`) |
|---|---|---|---|---|---|
| 0.5 | 3.44 (7.6e-22) | 4.27 (1.3e-26) | 6.44 (2.4e-32) | 1,742 | 9/9 both datasets |
| 0.75 | 4.42 (1.6e-26) | 7.23 (1.2e-38) | 11.15 (5.1e-44) | 1,007 | 9/9 both datasets |
| 1.0 | 5.51 (1.6e-30) | 10.87 (1.6e-40) | 18.85 (6.3e-49) | 536 | 7/9 Nature2026 (loses ERVFRD-1, KRT7) |

All three cutoffs show massive, highly significant enrichment for the
independent HPA placenta reference (background rate ~0.3-0.5% vs pass-set
rates of 1.5-8.4%) — confirming the pseudobulk DE signal is real placental
biology at every tested cutoff, not just at the top of the range.

Enrichment purity increases monotonically with cutoff (as expected — a
stricter logFC bar should enrich further for strong, specific placental
genes). But the purity gain from 0.75→1.0 has a real cost: **1.0 loses two
textbook trophoblast markers in Nature2026** (ERVFRD-1/syncytin-2, the
fusogen defining syncytiotrophoblast formation; KRT7, the defining
trophoblast cytokeratin) — both are core, uncontroversial trophoblast
biology, not edge cases. Losing them at the discovery stage (before any
independent replication beyond 2 datasets) is a worse trade than accepting
somewhat lower purity.

## Proposal (for reviewer sign-off, not yet frozen)

**`logFC ≥ 0.75` and `FDR < 0.05` per dataset, 2-of-2 quorum**, as the
`replicated_in_placenta(gene)` pass criterion:

- Full canonical-marker retention (9/9) in both datasets — 1.0 fails this.
- Substantially higher HPA-independent enrichment purity than 0.5 (11.15×
  vs 6.44× fold at the 2-of-2 overlap level; p=5.1e-44 vs 2.4e-32) — 0.5 is
  the weaker of the three candidates on this axis.
- 1,007 genes at the 2-of-2 overlap level — a workable candidate pool size
  for the downstream D/F/P assembly, not so large it looks unfiltered, not
  so small it risks losing real biology the way 1.0 demonstrably does.

This is a proposal, not a unilateral freeze — flagged here for the same
review-before-lock discipline used for every other threshold in this
project.
