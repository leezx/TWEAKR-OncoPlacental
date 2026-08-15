# Step 4a external validation: F_Colon/SI-developmental vs. 5 independent mike_verzi signatures

Executes the three-layer statistical validation planned in `docs/STEP4A_GUT_FDEV_DESIGN.md`'s "What this design does NOT do yet" and `docs/STEP4A_GUT_COMPUTE_RESULTS.md`'s "what's next": tests whether the real, same-atlas fetal-vs-adult gut DE signal (`F_Colon-developmental` / `F_SI-developmental`, PR #21) is independently corroborated by 5 published mouse intestinal fetal/regenerative signatures (`mike_verzi`, Step 6a, PR #19) that were never used to build the gut D/F/P and never touched the Gut Cell Atlas data.

**Non-circularity**: `mike_verzi` signatures come from an independent mouse study; `F_{Colon,SI}-developmental` come from independent human fetal-vs-adult edgeR DE. Neither informed the other. This is genuine external validation, not a consistency check on the same data.

## Design

- **Universe** (locked before testing, region-specific): `U = {human genes with ≥1 Compara one2one mouse ortholog} ∩ {genes passing filterByExpr in that region's primary edgeR fit}`. Same "eligible to appear in either set" logic used to fix `F_Gut-core`'s universe in PR #21. `N = 11,646` (`LargeInt`/Colon), `N = 11,661` (`SmallInt`/SI).
- **Layer 1 — hypergeometric enrichment**: fold-enrichment, odds ratio with 95% CI (Haldane-Anscombe continuity-corrected), one-sided hypergeometric p-value, BH-FDR across all 10 tests (5 signatures × 2 regions) together.
- **Layer 2 — preranked GSEA (primary evidence, per the user's explicit directive)**: `fgsea`, ranked by the real edgeR `logFC` from each region's primary fit — the continuous fetal-vs-adult differential statistic, not a binarized cutoff. `minSize=5, maxSize=5000, eps=0`, `set.seed(20260815)`.
- **Layer 3 — size-and-expression-matched permutation null**: 10,000 draws per (signature × region). Each real signature gene's `logCPM` decile bin (20 bins, from the real edgeR fit) is recorded; each permutation draws one random gene per bin from `U`, preserving exact gene-set size and expression-strata composition. Reports observed overlap, null-expected overlap, fold-vs-null, empirical p (`(#perm≥observed + 1)/(n_perm+1)`), null Z-score. This is the strongest available defense against "the overlap is just gene-set size/detectability."
- Fixed seed `20260815` throughout, stated for reproducibility.

Code: `scripts/04a_dfp_gut/mike_verzi_gut_enrichment_permutation.py` (layers 1+3), `scripts/04a_dfp_gut/mike_verzi_gut_gsea.R` (layer 2), run via `scripts/04a_dfp_gut/run_mike_verzi_gut_validation.sh` (qsub job 3620768, completed cleanly, stderr empty). Outputs in `results/04a_dfp_gut/mike_verzi_validation/`.

## Results

All numbers below are the real, unedited script output, pulled back from Argos and md5-verified byte-exact.

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

### Layer 3 — permutation null (10,000 draws, size + logCPM-decile matched)

| Region | Signature | Observed | Null expected | Fold vs. null | Empirical p | Z |
|---|---|---|---|---|---|---|
| Colon | FETAL_INTESTINE_GENES | 169 | 107.5 | **1.57** | **1e-4** | **6.76** |
| Colon | REGENERATIVE_EPITHELIUM | 33 | 13.8 | **2.39** | **1e-4** | **5.64** |
| Colon | YAP_SIGNALING_GENES | 40 | 30.2 | 1.33 | 0.0352 | 1.96 |
| Colon | FETAL_SPHEROID_EPITHELIUM_GENES | 25 | 21.6 | 1.16 | 0.243 (ns) | 0.79 |
| Colon | REVIVAL_STEM_CELL_GENES | 12 | 14.3 | 0.84 | 0.779 (ns) | -0.65 |
| SI | FETAL_INTESTINE_GENES | 151 | 102.8 | **1.47** | **1e-4** | **5.30** |
| SI | REGENERATIVE_EPITHELIUM | 19 | 12.4 | 1.53 | 0.0363 | 2.02 |
| SI | YAP_SIGNALING_GENES | 39 | 29.4 | 1.33 | 0.0396 | 1.90 |
| SI | FETAL_SPHEROID_EPITHELIUM_GENES | 24 | 21.1 | 1.14 | 0.286 (ns) | 0.66 |
| SI | REVIVAL_STEM_CELL_GENES | 11 | 19.7 | 0.56 (depleted) | 0.990 (ns) | -2.08 |

### Layer 2 — preranked GSEA (primary evidence, ranked by real logFC)

| Region | Signature | Size | NES | p | padj | Direction |
|---|---|---|---|---|---|---|
| Colon | FETAL_INTESTINE_GENES | 998 | **+1.50** | 1.5e-07 | **7.7e-07** | fetal-up ✓ |
| Colon | REGENERATIVE_EPITHELIUM | 119 | **+1.70** | 3.6e-04 | **9.0e-04** | fetal-up ✓ |
| Colon | YAP_SIGNALING_GENES | 288 | **-1.30** | 0.0199 | 0.0249 | **adult-up** |
| Colon | FETAL_SPHEROID_EPITHELIUM_GENES | 205 | **-1.38** | 0.0114 | 0.0190 | **adult-up** |
| Colon | REVIVAL_STEM_CELL_GENES | 206 | -1.01 | 0.444 | 0.444 (ns) | — |
| SI | FETAL_INTESTINE_GENES | 993 | **+1.46** | 1.5e-06 | **7.4e-06** | fetal-up ✓ |
| SI | REGENERATIVE_EPITHELIUM | 114 | **+1.46** | 0.0131 | 0.0164 | fetal-up ✓ |
| SI | YAP_SIGNALING_GENES | 284 | **-1.32** | 0.0091 | 0.0152 | **adult-up** |
| SI | FETAL_SPHEROID_EPITHELIUM_GENES | 203 | **-1.72** | 2.8e-05 | **7.0e-05** | **adult-up** |
| SI | REVIVAL_STEM_CELL_GENES | 206 | -1.06 | 0.301 | 0.301 (ns) | — |

## Interpretation — honest, not cherry-picked

**Two of five signatures triangulate as genuinely fetal-enriched across all three independent layers, in both regions:**

- **`FETAL_INTESTINE_GENES`** — the strongest and most consistent result. Significant positive hypergeometric enrichment (FDR<1e-7 both regions), the permutation null rejects at its floor (empirical p=1e-4, the minimum measurable at 10,000 draws) with Z≈5.3–6.8, and GSEA gives strongly significant positive NES (padj<1e-5 both regions) — genes are ranked toward the fetal-up end of the real, continuous edgeR effect, not just present in the significant-gene overlap. This is the cleanest possible triangulation.
- **`REGENERATIVE_EPITHELIUM`** — same pattern, strong in Colon (fold=2.80, permutation Z=5.64, GSEA padj=9.0e-4) and directionally consistent but weaker in SI (fold=1.68, permutation Z=2.02, GSEA padj=0.016) — still nominally significant on all three layers in both regions.

**One signature shows a real, unresolved discordance that is being reported plainly rather than smoothed over: `YAP_SIGNALING_GENES`.** Its hypergeometric and permutation tests show weak but nominally significant *positive* overlap enrichment with the fetal-up significant-gene set (fold≈1.4, permutation Z≈1.9–2.0, both FDR/empirical-p just under 0.04). But GSEA — using the full 284–288-gene set against the complete continuous ranking, not just the binary significant-overlap subset — gives significant **negative** NES in both regions (padj=0.025 Colon, 0.015 SI): the bulk of the YAP signature's genes actually skew toward the adult-up end of the real fetal-vs-adult ranking. These are not contradictory bugs; they are two different, both-real statistics describing the same gene set: a minority of YAP-signaling genes are among the significant fetal-up hits (driving the overlap test), while the set as a whole is, on average, shifted toward the adult side of the ranking (driving GSEA). Read together, `YAP_SIGNALING_GENES` does **not** support "YAP signaling as a whole is part of this fetal gut program" — if anything the GSEA result points the other way.

**`FETAL_SPHEROID_EPITHELIUM_GENES` shows no positive support, and a significant negative GSEA signal in SI.** Hypergeometric and permutation are both non-significant in both regions (fold≈1.2, permutation p≈0.24–0.29). GSEA is non-significant in Colon (p=0.011, but that's *negative* direction, i.e. this set skews adult-up, weakly) and clearly significant negative in SI (NES=-1.72, padj=7.0e-5) — a real, reportable adult-skewed signal, the strongest negative finding in the whole validation.

**`REVIVAL_STEM_CELL_GENES` shows no signal anywhere** — non-significant (in fact numerically depleted) on hypergeometric, non-significant on permutation, non-significant on GSEA, in both regions. A clean, honest null.

**Bottom line**: the re-anchored `F_Colon-developmental`/`F_SI-developmental` gene sets are genuinely, robustly externally validated by 2 of 5 independent published mouse fetal/regenerative intestinal signatures (`FETAL_INTESTINE_GENES`, `REGENERATIVE_EPITHELIUM`) — a real, triangulated, non-circular positive result. The other 3 signatures (`YAP_SIGNALING_GENES`, `FETAL_SPHEROID_EPITHELIUM_GENES`, `REVIVAL_STEM_CELL_GENES`) do not support a naive "all fetal/revival mouse signatures should be fetal-up in this human DE" expectation — two show no signal, and two (YAP, spheroid) show a genuine skew toward the *adult* side on the primary (GSEA) evidence layer. This is reported as-is; no signature was dropped, reweighted, or re-run to "improve" the pattern.

## Outputs

`results/04a_dfp_gut/mike_verzi_validation/`:
- `mike_verzi_gut_hypergeometric_enrichment.tsv`
- `mike_verzi_gut_permutation_null.tsv`
- `mike_verzi_gut_gsea_results.tsv` (includes leading-edge gene lists per signature/region)
- `validation_run_metadata.json` (n_permutations=10000, n_expression_bins=20, rng_seed=20260815, n_tests=10)
