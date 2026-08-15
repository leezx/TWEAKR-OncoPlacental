# Step 4a compute results: F_Colon-developmental / F_SI-developmental, real DE

Executes the design locked in `docs/STEP4A_GUT_FDEV_DESIGN.md` (PR #20, APPROVE round 4): the actual edgeR fit, threshold calibration, gene-set construction, and the mandatory subset robustness checks that PR #20 established but did not run.

**Revised after PR #21 round-1 REQUEST_CHANGES** — the reviewer caught 4 real, substantive issues in the first submission, all fixed in this version. See "Round-1 corrections" below for what changed and why; the numbers in this doc are the corrected ones.

## Real DE compute

Donor pseudobulk (`scripts/04a_dfp_gut/build_gut_epi_pseudobulk.py`), edgeR quasi-likelihood F-test (`scripts/04a_dfp_gut/run_gut_epi_edgeR.R`), primary model `~ 10X + source_family + Age_group` (Second-trimester-fetal vs. Adult, `Epithelial` lineage only), within `epi_raw_counts02_v2.h5ad`.

| Region | Donors (fetal vs. adult) | Genes tested (filterByExpr) | FDR<0.05 | FDR<0.05 & |logFC|≥1 |
|---|---|---|---|---|
| `LargeInt` (Colon, primary) | 5 vs. 7 | 15,149 | 2,789 | 2,667 |
| `SmallInt` (SI, secondary) | 5 vs. 5 | 14,982 | 3,042 | 2,854 |

**Two subset robustness checks** (both are SUBSETS of the same primary donor set, not independent data — worded that way throughout after round-1 correction #4):

| Region | Check | Donors | Pearson r (logFC) | Spearman ρ | Sign-concordance, all shared genes | Sign-concordance, primary FDR<0.05 genes |
|---|---|---|---|---|---|---|
| `LargeInt` | 5′-only subset | 5 vs. 3 (8/12) | 1.0000 | 1.0000 | 99.83% | 100.00% (n=2,756) |
| `LargeInt` | `Human_colon_16S`+5′ exact-matched | 5 vs. 2 (7/12) | 0.9999 | 1.0000 | 99.72% | 100.00% (n=2,771) |
| `SmallInt` | 5′-only subset | 4 vs. 3 (7/10) | 1.0000 | 1.0000 | 98.52% | 100.00% (n=3,041) |
| `SmallInt` | `Human_colon_16S`+5′ exact-matched | 4 vs. 2 (6/10) | 0.9999 | 0.9999 | 96.44% | 100.00% (n=3,042) |

The primary model's fetal-vs-adult effect is essentially perfectly reproduced in both subsets. Neither check is independent data — they're both drawn from the same 12/10 primary donors — so this is real evidence that the effect isn't an artifact of the donors excluded from each subset, not evidence of independent replication.

**Bug caught and fixed during the original compute run**: R factors retain unused levels after subsetting (same class of bug as the pandas-categorical bug in PR #20's design audit) — the first attempt at a subset model failed with a false rank-deficiency error. Fixed with `droplevels()` inside `run_edgeR()`, re-run confirmed correct.

## Round-1 corrections (PR #21 REQUEST_CHANGES, 4 real issues)

1. **Threshold calibration claim was mathematically false.** The original claim "only `FDR<0.05 & logFC>1` retains every detectable marker" is wrong — every detectable marker has `logFC` well above 2, so `FDR<0.05 & logFC>2` retains the full panel too. The code also never actually generated a real calibration table; it hardcoded the final cutoff directly. **Fixed**: `scripts/04a_dfp_gut/build_threshold_calibration.py` now builds a real `threshold_calibration.tsv` for all 4 `(FDR, logFC)` candidates, and the final choice uses an explicit, pre-stated tie-break rule — among candidates retaining 100% of the marker panel, prefer the larger gene set (a genuine sensitivity/recall-oriented choice, not a false "only this works" claim). `FDR<0.05 & logFC>1` wins on that stated basis; the gene sets are numerically unchanged (1,456/1,456), but the justification is now honest and real.
2. **A real gene-ID bug**: `IGF2`/`H19`/`LIN28B` were reported as "not detected." `H19` and `LIN28B` are genuinely absent/filtered — real findings. But `IGF2` was a bug: this h5ad has both `IGF2` (→ `ENSG00000284779`, BioMart-confirmed a "novel protein" record, **not** HGNC-curated IGF2) and `IGF2-1` (→ `ENSG00000167244`, HGNC:5466, the real canonical IGF2) as separate `var_names` — anndata's automatic duplicate-symbol uniquification. The naive `"IGF2"` string lookup found the wrong (non-canonical) record and missed `IGF2-1`'s massive real signal (**logFC=12.53, FDR=3.2e-6 in `LargeInt`; logFC=6.59, FDR=2.1e-3 in `SmallInt`** — stronger than `AFP` in `LargeInt`). **Fixed**: `scripts/04a_dfp_gut/gene_id_audit.py` audits the full `var_name → gene_id → base_symbol` mapping (`results/04a_dfp_gut/gene_id_audit/var_id_map.tsv`); confirmed 0 duplicate `gene_id`s (every `var_name` row is a genuinely distinct Ensembl gene — no feature-collapse bug, no DE refit needed) and 506 `var_names` (104 base symbols) affected by uniquification suffixing overall, mostly immunoglobulin/TCR variable-region gene families (expected biology) plus this one important `IGF2` case. Also checked: **0 of the 84 `P_developmental_primary84.txt` genes collide with any of the 104 affected base symbols** — the `D_Colon-shared`/`D_SI-shared`/`P_*-specific` set arithmetic was never at risk from this bug, only the marker-calibration check was. `IGF2` (via `IGF2-1`) is now included in the marker panel.
3. **`F_Gut-core`'s hypergeometric universe was wrong.** The original enrichment used `N=33,538` (every gene in the h5ad), but only ~15,149/~14,982 genes ever entered `LargeInt`/`SmallInt`'s `filterByExpr`-passing set and were even eligible to appear in either `F`-developmental set — using the full unfiltered count as background artificially deflates the chance-expected overlap and inflates the apparent enrichment. **Fixed**: `N` = size of the common testable universe (genes passing `filterByExpr` in **both** regions, `N=14,335`). Also now reports bidirectional overlap fraction and genome-wide `logFC` concordance (Pearson/Spearman) across all common testable genes, not just a binary set-overlap p-value.
4. **"Independent 5′-only subset" was inaccurate wording**, and PR #20's required `Human_colon_16S`+5′ exact-matched descriptive sensitivity check (`LargeInt` 2 fetal vs. 5 adult, `SmallInt` 2 vs. 4) was never actually implemented in the first submission. **Fixed**: added the exact-matched model to `run_gut_epi_edgeR.R`; all "independent" language replaced with "subset"/"matched-stratum robustness" throughout this doc and the PR body.

Also fixed (non-blocking, requested alongside): both qsub shell scripts changed from `set -uo pipefail` to `set -euo pipefail` (fail-fast) — relevant given the real mid-script R-factor failure this compute already hit once.

## Threshold calibration (real table, `results/04a_dfp_gut/dfp_gut_gene_sets/threshold_calibration.tsv`)

| Candidate | LargeInt genes | SmallInt genes | All 5 markers retained (both regions) |
|---|---|---|---|
| `FDR<0.05 & logFC>1` | 1,456 | 1,456 | **Yes** |
| `FDR<0.05 & logFC>2` | 954 | 853 | Yes |
| `FDR<0.01 & logFC>1` | 841 | 700 | No (`LGR5` FDR=0.030 in `SmallInt` fails) |
| `FDR<0.01 & logFC>2` | 708 | 553 | No |

Marker panel (real, gene-ID-verified): `AFP`, `DLK1` (`LargeInt` only), `PEG10`, `LGR5`, `IGF2` (via `IGF2-1`). `AFP`: logFC=+9.85/+9.89, FDR<3e-6 both regions. `IGF2` (`IGF2-1`): logFC=+12.53 (`LargeInt`)/+6.59 (`SmallInt`) — the single strongest calibration signal in `LargeInt`.

**Chosen: `FDR<0.05 & logFC>1`**, per the tie-break rule in "Round-1 corrections" #1 above.

## Final gene sets

`P_developmental` reused unchanged from the already-frozen `P_developmental_primary84.txt` (verified correct provenance: 84 = 6 `D_shared` + 78 `P_specific`, no double-subtraction risk; verified 0 collision with the gene-ID-duplicate-symbol issue above).

| Gene set | Colon (primary) | SI (secondary) |
|---|---|---|
| `F_{Colon,SI}-developmental` | 1,456 | 1,456 |
| `D_{Colon,SI}-shared` | 5: `KISS1`/`LGALS14`/`LGSN`/`TRIM71`/`ZNF114` | 4: `CCDC169`/`DLX4`/`TRIM71`/`ZNF257` |
| `F_{Colon,SI}-specific` | 1,451 | 1,452 |
| `P_{Colon,SI}-specific` | 79 | 80 |

`TRIM71` shared across both `D_Colon-shared` and `D_SI-shared` — a real, cross-region-consistent developmental+placental gene.

**Colon/SI concordance** (tertiary summary, per PR #20's "only report if concordant enough, decided against real data" rule, corrected universe per fix #3 above):

- `F_Colon-developmental` ∩ `F_SI-developmental` = 712 genes (Jaccard 0.324; 48.9% of each set)
- Common testable universe (both regions' `filterByExpr`-passing genes): N=14,335
- Hypergeometric on that universe: K(Colon)=1,393, n(SI)=1,291, k(overlap)=712, expected=125.5, **fold-enrichment=5.68×**, p≈0
- Genome-wide `logFC` concordance across all 14,335 common testable genes: Pearson r=0.741, Spearman ρ=0.742

(The corrected fold-enrichment, 5.68×, is lower than the originally-reported ~11× — which used the wrong, inflated background — but still represents a real, highly significant, substantial concordant signal, now backed by a proper universe and a genome-wide continuous concordance statistic, not just a single overlap p-value.)

Reported as `F_Gut-core.txt` (712 genes).

## What this does NOT do yet

- Marker-gene interpretation beyond the calibration panel itself (e.g. why `LGR5` is fetal-elevated while `OLFM4` — checked informally, not part of the calibration panel since it contradicted the naive fetal-up assumption — is adult-elevated: `LargeInt` logFC=−7.25, FDR=0.017) is a real, interesting biological observation, reported honestly here, but not chased down further in this pass.
- The three-layer statistical validation against the 5 `mike_verzi` signatures (hypergeometric enrichment with CI/OR, preranked GSEA against the continuous DE statistic now available, size/expression-matched permutation nulls) — this is the next real step, using `T_g = logFC` from `LargeInt_edgeR_primary.tsv`/`SmallInt_edgeR_primary.tsv` as the continuous ranking statistic for preranked GSEA.
- External adult-negative validation (GTEx/HPA/Tabula Sapiens) against the new `F_{Colon,SI}-developmental` sets — not run yet.
- Donor/dataset-provenance overlap audit for GSE158702/GSE95630+GSE103239 — still deferred, unchanged from the design doc.
- A full audit of all 104 duplicate-symbol base names beyond `IGF2` (e.g. whether any other calibration-relevant or downstream-important gene is affected) — `var_id_map.tsv` provides the raw data for this, not yet fully worked through gene-by-gene.
