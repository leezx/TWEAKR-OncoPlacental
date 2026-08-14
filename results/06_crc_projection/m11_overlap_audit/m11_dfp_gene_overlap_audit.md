# M11 × D/F/P gene-overlap audit

Per PR #16 review round 4 requirement: before using M11/revCSC as an independent
Oncofetal anchor whose score is correlated against the frozen D/F/P scores, check
directly whether M11's own gene list (top50/top100/top200, mapped Ensembl→symbol
via `adata_nmf.h5ad`'s own `var['GeneSymbol']` column — 50/100/200 of 50/100/200
genes mapped cleanly, no losses) shares any genes with D-shared, F-specific
(global and each of the 7 per-organ lineage modules, intersected with the frozen
`F_specific_FINAL.txt`), or P-specific. If shared genes exist, the same gene
would contribute to both the M11 score and the target D/F/P score, producing
mechanical (not biological) correlation.

## Result

| M11 version (n genes) | D-shared (6) | F-specific global (2,504) | P-specific (78) | F-lineage modules (7) |
|---|---|---|---|---|
| top50 | 0 | 0 | 0 | 0 in all 7 |
| top100 | 0 | 0 | 0 | 0 in all 7 |
| top200 | 0 | 1 (`MALAT1`) | 0 | `MALAT1` in 6 of 7 (all except Thymus) |

**Overlap is essentially negligible.** Across all three M11 gene-list versions
and every D/F/P signature (84 target genes total, plus 7 lineage modules),
the only shared gene anywhere is `MALAT1` — and only at the most permissive
top200 M11 cutoff. `MALAT1` is a near-ubiquitously highly-expressed nuclear
lncRNA well known in the single-cell literature to appear in the top-loading
gene list of almost any NMF/PCA factor regardless of biological specificity
(a common "usual suspect," not a developmental marker) — its presence in 6 of
7 F-lineage modules (all sharing a single gene, not a coordinated block) is
consistent with this generic technical pattern rather than real shared biology.

## Contract adopted for compute (per review)

- **Primary correlation uses the overlap-excluded M11 score**: for each pairwise
  comparison (M11 vs. D-shared, vs. F-specific global, vs. each F-lineage
  module, vs. P-specific), remove any gene present in both M11's scoring set
  and the target signature before computing M11's `score_genes` value for that
  specific comparison. Given the near-zero overlap found, this in practice
  only affects the M11-top200-vs-F-specific(-related) comparisons (drop
  `MALAT1` from the M11 gene set for those comparisons only).
- **Full (non-excluded) M11 score retained as a sensitivity analysis** —
  reported alongside the overlap-excluded primary result; given the overlap is
  ≤1 gene, the two are expected to be nearly identical, and any material
  difference between them would itself be a signal worth investigating.
- **Cutoff fallback rule**: if overlap-exclusion at a given cutoff leaves too
  few genes for a stable score (not expected here, since the maximum overlap
  found is 1 of 200 genes), fall back to a larger top-N version rather than
  forcing the smaller one.
- Only M11↔F/P/D associations that survive overlap exclusion may be
  interpreted as a genuine program-level developmental relationship, not an
  artifact of shared gene content.

## Files

- `M11_top50_symbols.txt` / `M11_top100_symbols.txt` / `M11_top200_symbols.txt` —
  M11's own gene lists, Ensembl IDs mapped to gene symbols via the source
  atlas's own `var['GeneSymbol']` column (pulled from Argos, verified byte-exact).
