# HPA placenta tissue-level enrichment calibration — freezing the P-developmental effect-size cutoff

Per the PR #9 reviewer's explicit next step: "先把 0.75 作为 leading candidate、
再用独立的 HPA trophoblast/placenta enrichment 决定是否冻结". Independent
check against HPA's own **bulk tissue-level** RNA consensus data (not
derived from the scRNA-seq datasets the edgeR DE itself used), to calibrate
which of the directional `logFC ≥ cutoff` candidates (0.5 / 0.75 / 1.0, all
already `FDR<0.05`) to freeze for `replicated_in_placenta`'s per-dataset
pass criterion.

**Scope note (PR #10 round 1 review)**: this is an HPA **placenta
tissue-level** calibration (bulk RNA consensus across HPA's 40 tissues),
not an HPA single-cell trophoblast cell-type calibration. The two should
not be conflated in naming — this doc covers only the tissue-level check
unless a separate HPA single-cell trophoblast reference is incorporated
later.

## Method

- **HPA reference sets** — computed directly from `rna_tissue_hpa.tsv`
  (Tier-1 adult reference, `STEP3_METHOD_CONTRACT.md`; placenta used here
  strictly as a **positive cross-check**, per the `role` column fix in
  PR #5 — never as adult-negative reference). Two sets, kept explicitly
  separate per PR #10 round 1 review (the first draft mislabeled the
  floored set as "HPA's own/official definition", which it is not):
  - **`no_floor`**: HPA's own published "Tissue enriched" rule, unmodified
    — nTPM in placenta ≥ 4× every other of HPA's 39 tissues. **604 genes.**
  - **`floored`**: the same rule *plus* a project-defined nTPM≥1 noise
    floor, to exclude near-zero-everywhere genes that can trivially clear
    a 4× ratio on noise. **82 genes.** This floor is a project addition,
    not part of HPA's definition — labeled as such everywhere in this repo
    going forward.
- **Enrichment test**: for each dataset's edgeR-tested gene universe
  (`filterByExpr`-passing genes) and each candidate cutoff, checked what
  fraction of the directional pass-set (`logFC≥cutoff & FDR<0.05`) falls in
  each HPA reference set, vs. the background rate in the full universe.
  Reported as fold-enrichment and a one-sided hypergeometric p-value, run
  **against both reference sets** as a sensitivity check on whether the
  floor is doing the work of picking the cutoff. Also computed for the
  2-of-2 overlap sets (genes passing in both datasets).
- Run on Argos (`scripts/04_dfp_signature/hpa_placenta_enrichment.py`),
  results pulled back and verified byte-exact (md5) against Argos output.

## Results

| Cutoff | Arutyunyan fold (floored / no_floor) | Nature2026 fold (floored / no_floor) | 2-of-2 overlap fold (floored / no_floor) | 2-of-2 gene count | Canonical marker retention |
|---|---|---|---|---|---|
| 0.5 | 3.44 / 3.28 | 4.27 / 3.92 | 6.44 / 6.14 | 1,742 | 9/9 both datasets |
| 0.75 | 4.42 / 4.21 | 7.23 / 6.65 | 11.15 / 10.62 | 1,007 | 9/9 both datasets |
| 1.0 | 5.51 / 5.26 | 10.87 / 10.11 | 18.85 / 18.00 | 536 | 7/9 Nature2026 (loses ERVFRD-1, KRT7) |

All p-values remain ≪1e-19 for both reference sets at every cutoff (full
table in `hpa_placenta_enrichment.tsv`). **The floored and un-floored
reference sets give the same monotonic trend and the same relative gap
between cutoffs** (floored fold is 5-8% higher than no_floor throughout,
never enough to change which cutoff looks best) — confirms the nTPM≥1
project floor is not driving the calibration result; HPA's own unmodified
fourfold rule alone supports the same conclusion.

Enrichment purity increases monotonically with cutoff under both reference
sets (as expected — a stricter logFC bar should enrich further for strong,
specific placental genes). But the purity gain from 0.75→1.0 has a real
cost: **1.0 loses two textbook trophoblast markers in Nature2026**
(ERVFRD-1/syncytin-2, the fusogen defining syncytiotrophoblast formation;
KRT7, the defining trophoblast cytokeratin) — both are core,
uncontroversial trophoblast biology, not edge cases. Losing them at the
discovery stage (before any independent replication beyond 2 datasets) is
a worse trade than accepting somewhat lower purity.

## Proposal (for reviewer sign-off, not yet frozen)

**`logFC ≥ 0.75` and `FDR < 0.05` per dataset, 2-of-2 quorum**, as the
`replicated_in_placenta(gene)` pass criterion:

- Full canonical-marker retention (9/9) in both datasets — 1.0 fails this.
- Substantially higher HPA-independent enrichment purity than 0.5 (11.15×
  vs 6.44× fold at the 2-of-2 overlap level, floored reference; 10.62× vs
  6.14× under the un-floored sensitivity check) — 0.5 is the weaker of the
  three candidates on this axis, under both reference definitions.
- 1,007 genes at the 2-of-2 overlap level — a workable candidate pool size
  for the downstream D/F/P assembly, not so large it looks unfiltered, not
  so small it risks losing real biology the way 1.0 demonstrably does.
- Conclusion is robust to whether the nTPM≥1 project floor is applied or
  not — not an artifact of that choice.

This is a proposal, not a unilateral freeze — flagged here for the same
review-before-lock discipline used for every other threshold in this
project.
