# HPA placenta tissue-level enrichment calibration — freezing the P-developmental effect-size cutoff

Per the PR #9 reviewer's explicit next step: "先把 0.75 作为 leading candidate、
再用独立的 HPA trophoblast/placenta enrichment 决定是否冻结". Independent
check against HPA's own **bulk tissue-level** RNA data (not derived from
the scRNA-seq datasets the edgeR DE itself used), to calibrate which of the
directional `logFC ≥ cutoff` candidates (0.5 / 0.75 / 1.0, all already
`FDR<0.05`) to freeze for `replicated_in_placenta`'s per-dataset pass
criterion.

**Scope note**: this is an HPA **placenta tissue-level** calibration (bulk
RNA consensus across HPA's 40 tissues), not an HPA single-cell trophoblast
cell-type calibration — the two should not be conflated in naming.

## Method — three reference sets, kept explicitly separate

Two earlier drafts of this doc conflated a project recomputation with
HPA's own official classification (PR #10 round 1 and round 2 review both
caught real attribution/methodology problems here — see "Review history"
below). Final version uses three clearly-labeled reference sets:

1. **`official_65`** — the real HPA classification. Fetched directly from
   proteinatlas.org's search API (`tissue_category_rna:placenta;Tissue
   enriched`), not recomputed from the raw matrix. **65 genes** — matches
   HPA's own reported count exactly. Saved as
   `results/04_dfp_signature/hpa_placenta_official_tissue_enriched_65genes.tsv`.
   This is the ground truth for this calibration.
2. **`floored`** — a project recomputation from `rna_tissue_hpa.tsv`
   (Tier-1 adult reference, placenta used strictly as a positive
   cross-check per the PR #5 `role` fix): nTPM in placenta ≥ 4× every other
   of HPA's 39 tissues, plus a project-defined nTPM≥1 floor to exclude
   zero/zero artifacts (see below). **82 genes.** Not HPA's official
   classification — a sensitivity check.
3. **`no_floor`** — the same recomputation without the floor. **604
   genes.** Also not HPA's official classification; kept only as a second
   sensitivity check, with a known artifact (see below).
- **Enrichment test**: for each dataset's edgeR-tested gene universe
  (`filterByExpr`-passing genes) and each candidate cutoff, checked what
  fraction of the directional pass-set (`logFC≥cutoff & FDR<0.05`) falls in
  each of the 3 reference sets, vs. the background rate in the full
  universe. Reported as fold-enrichment and a one-sided hypergeometric
  p-value. Also computed for the 2-of-2 overlap set (genes passing in both
  datasets).
- Run on Argos (`scripts/04_dfp_signature/hpa_placenta_enrichment.py`),
  results pulled back and verified byte-exact (md5) against Argos output.

## Results — all three references agree

| Cutoff | official_65 fold (p) | floored fold (p) | no_floor fold (p) | 2-of-2 gene count | Marker retention |
|---|---|---|---|---|---|
| 0.5 | 6.63 (3.6e-28) | 6.44 (2.4e-32) | 6.14 (2.4e-31) | 1,742 | 9/9 both datasets |
| 0.75 | 11.46 (5.6e-38) | 11.15 (5.1e-44) | 10.62 (3.4e-43) | 1,007 | 9/9 both datasets |
| 1.0 | 19.49 (2.8e-42) | 18.85 (6.3e-49) | 18.00 (1.3e-48) | 536 | 7/9 Nature2026 (loses ERVFRD-1, KRT7) |

(2-of-2-overlap fold-enrichment shown; full per-dataset breakdown in
`hpa_placenta_enrichment.tsv`.) The real HPA official 65-gene set gives
numbers essentially identical to both recomputed sensitivity sets — same
monotonic trend, same relative gap between cutoffs, all p≪1e-19 throughout.
Notably, HPA's real official list includes several of the canonical
markers already used in the DE sanity check (CSH1, CSH2, PSG1, PSG3,
ERVFRD-1), which is an independent confirmation that the marker panel and
the HPA reference are measuring the same underlying biology.

Enrichment purity increases monotonically with cutoff under all three
references (as expected — a stricter logFC bar should enrich further for
strong, specific placental genes). But the purity gain from 0.75→1.0 has a
real cost: **1.0 loses two textbook trophoblast markers in Nature2026**
(ERVFRD-1/syncytin-2, the fusogen defining syncytiotrophoblast formation;
KRT7, the defining trophoblast cytokeratin) — both are core,
uncontroversial trophoblast biology, not edge cases. Losing them at the
discovery stage (before any independent replication beyond 2 datasets) is
a worse trade than accepting somewhat lower purity.

## Proposal (for reviewer sign-off, not yet frozen)

**`logFC ≥ 0.75` and `FDR < 0.05` per dataset, 2-of-2 quorum**, as the
`replicated_in_placenta(gene)` pass criterion:

- Full canonical-marker retention (9/9) in both datasets — 1.0 fails this.
- Substantially higher enrichment purity than 0.5 across all three
  independent reference sets, including the real HPA official
  classification (11.46× vs 6.63× fold at the 2-of-2 overlap level) — 0.5
  is the weaker of the three candidates on this axis under every
  reference tested.
- 1,007 genes at the 2-of-2 overlap level — a workable candidate pool size
  for the downstream D/F/P assembly, not so large it looks unfiltered, not
  so small it risks losing real biology the way 1.0 demonstrably does.
- Conclusion is robust across the real HPA official set and two
  independent project recomputations (with and without a noise floor) —
  not an artifact of any one reference-set construction choice.

This is a proposal, not a unilateral freeze — flagged here for the same
review-before-lock discipline used for every other threshold in this
project.

## Review history (PR #10)

- **Round 1**: reviewer caught that a project-defined `nTPM≥1` noise floor
  had been mislabeled as part of "HPA's own/official definition." Fixed by
  separating `floored` from `no_floor` and relabeling both as project
  recomputations.
- **Round 2**: reviewer looked up HPA's actual published placenta
  "Tissue enriched" classification directly (65 genes) and caught that the
  round-1 `no_floor` set (604 genes) still didn't reproduce it and was
  still being described as "HPA's rule exactly as published" — the real
  bug: when placenta and every other tissue are both 0, `0 >= 4*0` is
  trivially true, so zero-everywhere genes flood into the naive
  recomputation. Fixed by fetching HPA's real classification from
  proteinatlas.org's search API as a third, authoritative reference
  (`official_65`) and relabeling both recomputed sets explicitly as
  project recomputations, never "HPA's rule as published." All three
  references now agree on the same conclusion.
