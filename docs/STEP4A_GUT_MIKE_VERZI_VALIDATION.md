# Step 4a external validation: F_Colon/SI-developmental vs. 5 independent mike_verzi signatures

Executes the three-layer statistical validation planned in `docs/STEP4A_GUT_FDEV_DESIGN.md`'s "What this design does NOT do yet" and `docs/STEP4A_GUT_COMPUTE_RESULTS.md`'s "what's next": tests whether the real, same-atlas fetal-vs-adult gut DE signal (`F_Colon-developmental` / `F_SI-developmental`, PR #21) is independently corroborated by 5 published mouse intestinal fetal/regenerative signatures (`mike_verzi`, Step 6a, PR #19) that were never used to build the gut D/F/P and never touched the Gut Cell Atlas data.

**Non-circularity**: `mike_verzi` signatures come from an independent mouse study; `F_{Colon,SI}-developmental` come from independent human fetal-vs-adult edgeR DE. Neither informed the other. This is genuine external validation, not a consistency check on the same data.

**Revised after PR #22 round-1 REQUEST_CHANGES** — the reviewer caught 3 real, substantive issues in the first submission (all statistical/interpretive, none touching the underlying gut D/F/P or edgeR results). See "Round-1 corrections" below; the numbers in this doc are the corrected ones.

## Design

- **Universe** (locked before testing, region-specific, and now applied consistently across **all three layers** — see round-1 correction #1): `U = {human genes with ≥1 Compara one2one mouse ortholog} ∩ {genes passing filterByExpr in that region's primary edgeR fit}`. Same "eligible to appear in either set" logic used to fix `F_Gut-core`'s universe in PR #21. `N = 11,646` (`LargeInt`/Colon), `N = 11,661` (`SmallInt`/SI).
- **Layer 1 — hypergeometric enrichment**: fold-enrichment, odds ratio with 95% CI (Haldane-Anscombe continuity-corrected), one-sided hypergeometric p-value, BH-FDR across all 10 tests (5 signatures × 2 regions) together.
- **Layer 2 — preranked GSEA (primary evidence, per the user's explicit directive)**: `fgsea`, ranked by the real edgeR `logFC` from each region's primary fit, **restricted to the same universe `U`** as layers 1/3 — the continuous fetal-vs-adult differential statistic, not a binarized cutoff. `minSize=5, maxSize=5000, eps=0`, `set.seed(20260815)`.
- **Layer 3 — size-and-expression-matched permutation null**: 10,000 draws per (signature × region), universe `U`. Each real signature gene's `logCPM` 20-quantile bin (from the real edgeR fit) is recorded; each permutation draws one random gene per bin from `U`, preserving exact gene-set size and expression-strata composition. Reports observed overlap, null-expected overlap, fold-vs-null, empirical p (`(#perm≥observed + 1)/(n_perm+1)`), **BH-FDR across the same 10 tests** (round-1 correction #2), and null Z-score.
- Fixed seed `20260815` throughout, stated for reproducibility.

Code: `scripts/04a_dfp_gut/mike_verzi_gut_enrichment_permutation.py` (layers 1+3), `scripts/04a_dfp_gut/mike_verzi_gut_gsea.R` (layer 2), run via `scripts/04a_dfp_gut/run_mike_verzi_gut_validation.sh` (qsub job 3620772, completed cleanly, stderr empty). Outputs in `results/04a_dfp_gut/mike_verzi_validation/`.

## Results

All numbers below are the real, unedited script output from the corrected re-run, pulled back from Argos and md5-verified byte-exact.

### Layer 1 — hypergeometric enrichment

| Region | Signature | n (sig∩U) | k (overlap) | Fold | OR [95% CI] | p | FDR |
|---|---|---|---|---|---|---|---|
| Colon | FETAL_INTESTINE_GENES | 998 | 169 | **1.71** | 2.00 [1.68–2.39] | 4.0e-13 | **4.0e-12** |
| Colon | REGENERATIVE_EPITHELIUM | 119 | 33 | **2.80** | 3.59 [2.40–5.38] | 2.8e-08 | **9.2e-08** |
| Colon | YAP_SIGNALING_GENES | 288 | 40 | 1.40 | 1.50 [1.07–2.10] | 0.0175 | 0.0350 |
| Colon | FETAL_SPHEROID_EPITHELIUM_GENES | 205 | 25 | 1.23 | 1.29 [0.85–1.96] | 0.161 | 0.230 (ns) |
| Colon | REVIVAL_STEM_CELL_GENES | 206 | 12 | 0.59 | 0.58 [0.33–1.03] | 0.987 (depleted) | 0.994 (ns) |
| SI | FETAL_INTESTINE_GENES | 993 | 151 | **1.54** | 1.73 [1.44–2.08] | 2.2e-08 | **9.2e-08** |
| SI | REGENERATIVE_EPITHELIUM | 114 | 19 | 1.68 | 1.87 [1.15–3.06] | 0.0161 | 0.0350 |
| SI | YAP_SIGNALING_GENES | 284 | 39 | 1.39 | 1.48 [1.05–2.08] | 0.0218 | 0.0363 |
| SI | FETAL_SPHEROID_EPITHELIUM_GENES | 203 | 24 | 1.20 | 1.25 [0.81–1.91] | 0.206 | 0.258 (ns) |
| SI | REVIVAL_STEM_CELL_GENES | 206 | 11 | 0.54 | 0.53 [0.29–0.97] | 0.994 (depleted) | 0.994 (ns) |

### Layer 3 — permutation null (10,000 draws, size + logCPM-20-quantile-bin matched, universe `U`)

BH-FDR now applied across the same 10 tests as layer 1 (round-1 correction #2) — several nominally-under-0.05 empirical p's do **not** survive correction.

| Region | Signature | Observed | Null expected | Fold vs. null | Empirical p | Empirical FDR | Z |
|---|---|---|---|---|---|---|---|
| Colon | FETAL_INTESTINE_GENES | 169 | 107.7 | **1.57** | 1e-4 | **3.3e-4** | **6.61** |
| Colon | REGENERATIVE_EPITHELIUM | 33 | 13.8 | **2.39** | 1e-4 | **3.3e-4** | **5.65** |
| Colon | YAP_SIGNALING_GENES | 40 | 30.1 | 1.33 | 0.0349 | 0.070 (ns) | 1.97 |
| Colon | FETAL_SPHEROID_EPITHELIUM_GENES | 25 | 21.6 | 1.16 | 0.247 | 0.353 (ns) | 0.80 |
| Colon | REVIVAL_STEM_CELL_GENES | 12 | 14.4 | 0.83 | 0.786 | 0.873 (ns) | -0.67 |
| SI | FETAL_INTESTINE_GENES | 151 | 102.8 | **1.47** | 1e-4 | **3.3e-4** | **5.30** |
| SI | REGENERATIVE_EPITHELIUM | 19 | 12.4 | 1.54 | 0.0350 | 0.070 (ns) | 2.02 |
| SI | YAP_SIGNALING_GENES | 39 | 29.6 | 1.32 | 0.0427 | 0.071 (ns) | 1.89 |
| SI | FETAL_SPHEROID_EPITHELIUM_GENES | 24 | 21.3 | 1.13 | 0.289 | 0.361 (ns) | 0.65 |
| SI | REVIVAL_STEM_CELL_GENES | 11 | 19.6 | 0.56 (depleted) | 0.991 | 0.991 (ns) | -2.10 |

Only the 3 strongest hits (Colon FETAL_INTESTINE, Colon REGENERATIVE_EPITHELIUM, SI FETAL_INTESTINE) survive 10-test FDR at the permutation layer; the 3 borderline nominal-p<0.05 results (Colon YAP, SI REGENERATIVE_EPITHELIUM, SI YAP) all land at FDR≈0.07 and are **not** significant once corrected for the same multiple-testing burden layer 1 was already held to.

### Layer 2 — preranked GSEA (primary evidence, ranked by real logFC, universe `U`)

| Region | Signature | Size | NES | p | padj | Direction |
|---|---|---|---|---|---|---|
| Colon | FETAL_INTESTINE_GENES | 998 | **+1.42** | 1.5e-05 | **7.6e-05** | fetal-up ✓ |
| Colon | REGENERATIVE_EPITHELIUM | 119 | **+1.63** | 1.6e-03 | **2.0e-03** | fetal-up ✓ |
| Colon | YAP_SIGNALING_GENES | 288 | **-1.47** | 3.7e-04 | **9.2e-04** | **adult-up** |
| Colon | FETAL_SPHEROID_EPITHELIUM_GENES | 205 | **-1.54** | 5.8e-04 | **9.6e-04** | **adult-up** |
| Colon | REVIVAL_STEM_CELL_GENES | 206 | -1.14 | 0.125 | 0.125 (ns) | — |
| SI | FETAL_INTESTINE_GENES | 993 | **+1.44** | 1.1e-05 | **2.7e-05** | fetal-up ✓ |
| SI | REGENERATIVE_EPITHELIUM | 114 | **+1.40** | 0.0224 | 0.0280 | fetal-up ✓ |
| SI | YAP_SIGNALING_GENES | 284 | **-1.45** | 1.7e-03 | **2.8e-03** | **adult-up** |
| SI | FETAL_SPHEROID_EPITHELIUM_GENES | 203 | **-1.85** | 1.1e-06 | **5.3e-06** | **adult-up** |
| SI | REVIVAL_STEM_CELL_GENES | 206 | -1.18 | 0.109 | 0.109 (ns) | — |

## Interpretation — honest, not cherry-picked, FDR-consistent across all three layers

**`FETAL_INTESTINE_GENES` fully triangulates as genuinely fetal-enriched across all three FDR-corrected layers, in both regions; `REGENERATIVE_EPITHELIUM` fully triangulates in Colon, with 2-of-3-layers support in SI:**

- **`FETAL_INTESTINE_GENES`** — the strongest and most consistent result. Significant hypergeometric enrichment (FDR<1e-7 both regions), permutation FDR=3.3e-4 both regions (Z≈5.3–6.6), and GSEA gives significant positive NES (padj<8e-5 both regions). All three layers agree, FDR-corrected, in both regions.
- **`REGENERATIVE_EPITHELIUM`** — full triangulation in **Colon** (hypergeometric FDR=9.2e-8, permutation FDR=3.3e-4, GSEA padj=2.0e-3). In **SI**, the pattern is weaker once permutation is FDR-corrected: hypergeometric FDR=0.035 and GSEA padj=0.028 are both significant, but permutation is only nominal (p=0.035, FDR=0.070, ns). So SI REGENERATIVE_EPITHELIUM is 2-of-3-layers significant after correction, directionally consistent on the third, not full triangulation.

**One signature shows a real, unresolved discordance, reported plainly: `YAP_SIGNALING_GENES`.** Hypergeometric shows nominally significant positive overlap enrichment (FDR=0.035 Colon, 0.036 SI). Permutation shows the same direction but does **not** survive FDR correction in either region (FDR=0.070/0.071, ns) — so the overlap-based signal is weak and not robust to the same multiple-testing standard used elsewhere. GSEA — the primary evidence layer, using the full 284–288-gene set against the complete continuous ranking — gives significant **negative** NES in both regions (padj=9.2e-4 Colon, 2.8e-3 SI): the bulk of the signature's genes skew toward the adult-up end of the ranking. Taken together, `YAP_SIGNALING_GENES` does **not** support "YAP signaling as a whole is part of this fetal gut program" — the only FDR-robust signal for this set (GSEA) points the other way.

**`FETAL_SPHEROID_EPITHELIUM_GENES` shows no positive support on the overlap-based layers, and significant negative GSEA in *both* regions.** Hypergeometric and permutation are both non-significant in both regions (hypergeometric FDR=0.230/0.258, permutation FDR=0.353/0.361). GSEA is significant negative in **both** Colon (NES=-1.54, padj=9.6e-4) and SI (NES=-1.85, padj=5.3e-6) — the strongest and most consistent negative finding in the whole validation.

**`REVIVAL_STEM_CELL_GENES` shows no signal anywhere** — non-significant (numerically depleted) on hypergeometric, non-significant on permutation, non-significant on GSEA, in both regions. A clean, honest null.

**Bottom line**: the re-anchored `F_Colon-developmental`/`F_SI-developmental` gene sets are robustly, FDR-consistently externally validated by `FETAL_INTESTINE_GENES` in both regions and by `REGENERATIVE_EPITHELIUM` in Colon (full triangulation) with weaker, partial support in SI. Of the remaining three signatures: `REVIVAL_STEM_CELL_GENES` shows no signal anywhere; `YAP_SIGNALING_GENES` and `FETAL_SPHEROID_EPITHELIUM_GENES` both show a significant *adult*-skewed signal on GSEA (the primary evidence layer) in both regions, with no FDR-robust positive support on the overlap-based layers. This is reported as-is; no signature was dropped, reweighted, or re-run to improve the pattern.

## Round-1 corrections (PR #22 REQUEST_CHANGES, 3 real issues)

1. **Blocker — GSEA's ranked universe didn't match layers 1/3.** The original GSEA script ranked all ~15k `filterByExpr`-tested genes per region, not the `one2one-ortholog-eligible ∩ filterByExpr-tested` universe (`U`, N=11,646/11,661) locked for layers 1/3. Genes without a one-to-one mouse ortholog can structurally never be a `mike_verzi` signature member (the signature files are themselves built only from one2one orthologs), so including them as guaranteed misses shifted GSEA's ES/NES/p reference space away from layers 1/3's. **Fixed**: `mike_verzi_gut_gsea.R` now restricts the ranked list to `U` before running `fgsea`, using the same `mouse_biomart_full.tsv` one2one filter as the Python script. Re-run confirms the conclusions are unchanged in direction, and in fact strengthened (e.g. Colon `FETAL_SPHEROID_EPITHELIUM_GENES` padj tightened from 0.019 to 9.6e-4 — see correction #3).
2. **Blocker — the permutation layer's 10 tests weren't BH-FDR corrected, while layer 1's identical 10 tests were**, so the two overlap-based layers were being held to inconsistent statistical standards. Three nominal-p<0.05 permutation results (Colon YAP p=0.035, SI REGENERATIVE_EPITHELIUM p=0.035, SI YAP p=0.043) do **not** survive BH-FDR across the same 10 tests (all land at FDR≈0.07). **Fixed**: `mike_verzi_gut_enrichment_permutation.py` now computes and reports `empirical_fdr` (BH across all 10 tests) alongside `empirical_p`; this doc's interpretation now treats permutation FDR, not nominal p, as the significance criterion, matching layer 1.
3. **A factual contradiction in the original write-up**: the results table itself already showed Colon `FETAL_SPHEROID_EPITHELIUM_GENES` GSEA as significant (`padj=0.019 < 0.05`), but the interpretation text described it as "non-significant in Colon." **Fixed**: after the universe correction (#1), Colon is now unambiguously significant (`padj=9.6e-4`), consistent with SI — the interpretation above states this correctly, and a related miscount in the original "bottom line" (which said "two show no signal, and two show adult skew" for what were actually three remaining signatures) is also corrected.

Also fixed alongside (non-blocking, reviewer-suggested): `20-quantile bin` terminology used consistently instead of the imprecise `decile bin` (the code always used 20 bins, not 10).

## Outputs

`results/04a_dfp_gut/mike_verzi_validation/`:
- `mike_verzi_gut_hypergeometric_enrichment.tsv`
- `mike_verzi_gut_permutation_null.tsv` (now includes `empirical_fdr`)
- `mike_verzi_gut_gsea_results.tsv` (universe-restricted; includes leading-edge gene lists per signature/region)
- `validation_run_metadata.json` (n_permutations=10000, n_expression_bins=20, rng_seed=20260815, n_tests=10)
