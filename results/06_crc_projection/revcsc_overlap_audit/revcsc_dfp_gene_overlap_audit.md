# revCSC x D/F/P gene-overlap audit

Per Phase I design pivot: revCSC (not M11) is the primary Oncofetal anchor. Checked directly whether revCSC's own mapped human gene set shares genes with D-shared, F-specific (global and each of the 7 per-organ lineage modules, intersected with the frozen F_specific_FINAL.txt), or P-specific, before any null-calibrated revCSC score is correlated against D/F/P scores.

## revCSC signature extraction

- Raw rows in `CSC_subtype_signatures.ensembl_mapping.tsv` for cluster=revCSC: 42 (includes duplicate rows -- the source table repeats some genes).
- Distinct gene symbols (raw, before mapping filter): 32.
- Successfully mapped to a human Ensembl ID: 31.
- Failed mouse->human mapping (dropped from scoring set): 1 (CTLA2A).

## Overlap result

| Target signature (n genes) | revCSC overlap (n) | Overlapping genes |
|---|---|---|
| D-shared (6) | 0 | - |
| F-specific global (2504) | 2 | ACTA1, ANKRD1 |
| P-specific (78) | 0 | - |
| F-lineage Adrenal (680) | 0 | - |
| F-lineage Liver (797) | 0 | - |
| F-lineage Skin (1123) | 0 | - |
| F-lineage Spleen (1087) | 1 | ANKRD1 |
| F-lineage Stomach (749) | 0 | - |
| F-lineage Thymus (1189) | 2 | ACTA1, ANKRD1 |
| F-lineage Thyroid (682) | 1 | ACTA1 |

revCSC scoring-set size used going forward: 31 genes (this was the pre-ortholog-audit set, listed here for traceability only — ACTA1, ANKRD1, ANXA1, AREG, ASS1, BASP1, CCN1, CCN2, CD44, CLDN4, CLU, CTSE, CTSL, ECM1, F3, FN1, ITGA2, KRT18, LY6A, MARCKSL1, PMEPA1, PRDX2, PYY, SFN, SOX4, SOX9, SPRR1A, TMSB4X, TNFRSF12A, TPM1, TUBA1A).

**Superseded (PR #17 round 2)**: this 31-gene set proved symbol-matching success only, not verified mouse→human orthology. See `revcsc_ortholog_provenance_audit.md` for the Ensembl-Compara-verified correction — the frozen scoring set is now 28 genes (27 confirmed `one2one` + `Ly6a` flagged `one2many`), with `Cldn4`/`Ctsl`/`Sprr1a` dropped (no Compara-confirmed ortholog despite a valid-looking symbol match) and `Ccn1`'s Ensembl ID corrected (the old ID pointed to `CCNA2`, a different gene). The D/F/P overlap numbers above (D=0, F-specific=2 [`ACTA1`,`ANKRD1`], P=0) are unaffected by this correction and do not need to be re-run — all 4 excluded/corrected genes are absent from every D/F/P set. All future steps should use `revCSC_human_FINAL.tsv` / `revCSC_symbols.primary28.txt`, not this file's 31-gene list.
