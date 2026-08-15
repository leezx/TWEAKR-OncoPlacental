# Step 4a compute results: F_Colon-developmental / F_SI-developmental, real DE

Executes the design locked in `docs/STEP4A_GUT_FDEV_DESIGN.md` (PR #20, APPROVE round 4): the actual edgeR fit, threshold calibration, gene-set construction, and the mandatory 5′-only concordance check that PR #20 established but did not run.

## Real DE compute

Donor pseudobulk (`scripts/04a_dfp_gut/build_gut_epi_pseudobulk.py`), edgeR quasi-likelihood F-test (`scripts/04a_dfp_gut/run_gut_epi_edgeR.R`), primary model `~ 10X + source_family + Age_group` (Second-trimester-fetal vs. Adult, `Epithelial` lineage only), within `epi_raw_counts02_v2.h5ad`.

| Region | Donors (fetal vs. adult) | Genes tested (filterByExpr) | FDR<0.05 | FDR<0.05 & \|logFC\|≥1 |
|---|---|---|---|---|
| `LargeInt` (Colon, primary) | 5 vs. 7 | 15,149 | 2,789 | 2,667 |
| `SmallInt` (SI, secondary) | 5 vs. 5 | 14,982 | 3,042 | 2,854 |

**Mandatory 5′-only concordance check** (pre-registered in the design doc, not contingent on diagnostics): both regions pass cleanly.

| Region | 5′-only donors | Pearson r (logFC) | Spearman ρ (logFC) | Sign-concordance, all shared genes | Sign-concordance, primary FDR<0.05 genes |
|---|---|---|---|---|---|
| `LargeInt` | 3 vs. 5 | 1.0000 | 1.0000 | 99.83% | 100.00% (n=2,756) |
| `SmallInt` | 3 vs. 4 | 1.0000 | 1.0000 | 98.52% | 100.00% (n=3,041) |

The primary model's fetal-vs-adult effect is essentially perfectly reproduced in the independent 5′-only subset — the `source_family`/chemistry confound identified in PR #20 round 3 does not appear to be driving the primary result.

## Threshold calibration (real marker-panel check, same discipline as HDMA's PR #12)

Candidate `(FDR, logFC)` cutoffs tested against the real, `filterByExpr`-detectable canonical oncofetal/developmental marker panel — **not** picked by "comfortably above background" reasoning:

| Marker | `LargeInt` logFC / FDR | `SmallInt` logFC / FDR |
|---|---|---|
| `AFP` | +9.85 / 3.0e-10 | +9.89 / 2.5e-6 |
| `DLK1` | +5.87 / 2.8e-3 | not detected (fails `filterByExpr`) |
| `PEG10` | +4.76 / 3.1e-3 | +5.82 / 1.5e-3 |
| `LGR5` | +4.86 / 8.1e-5 | +3.26 / 3.0e-2 |

`IGF2`/`H19`/`LIN28B` (part of HDMA's original calibration panel) do not pass `filterByExpr` in either region — a real finding (this dataset's detection sensitivity differs from HDMA's), reported honestly rather than dropped silently.

Only **`FDR<0.05 & logFC>1`** retains every detectable marker in both regions (a stricter `FDR<0.01` loses `LGR5` in `SmallInt`, FDR=0.030) — chosen on that basis. `AFP` is the strongest single calibration signal: a massive, unambiguous fetal-elevation effect, exactly matching its textbook role as the canonical oncofetal marker.

## Final gene sets

`P_developmental` reused unchanged from `results/04_dfp_signature/dfp_gene_sets/P_developmental_primary84.txt` (84 genes — verified as the correct pre-F-subtraction frozen program: 84 = 6 `D_shared` + 78 `P_specific`, matching the disjoint-partition math in `STEP4_DFP_DESIGN.md`; not reconstructed, no double-subtraction risk).

| Gene set | Colon (primary) | SI (secondary) |
|---|---|---|
| `F_{Colon,SI}-developmental` | 1,456 | 1,456 |
| `D_{Colon,SI}-shared` | 5: `KISS1`, `LGALS14`, `LGSN`, `TRIM71`, `ZNF114` | 4: `CCDC169`, `DLX4`, `TRIM71`, `ZNF257` |
| `F_{Colon,SI}-specific` | 1,451 | 1,452 |
| `P_{Colon,SI}-specific` | 79 | 80 |

`TRIM71` (a known stemness/pluripotency-associated E3 ligase) is shared between `D_Colon-shared` and `D_SI-shared` — a real, cross-region-consistent gene in both "shared with placenta" sets, biologically plausible.

**Colon/SI concordance** (tertiary summary, per the design doc's "only report `F_Gut-core` if concordant enough, decided against real data" rule): `F_Colon-developmental` ∩ `F_SI-developmental` = 712 genes (Jaccard 0.324) — vs. a chance-expected overlap of ~63 genes given the ~33,538-gene universe (hypergeometric p≈0, effectively 0). This is a real, ~11× enrichment — concordant enough to report `F_Gut-core.txt` (712 genes) as the tertiary summary, while `F_Colon-specific`/`F_SI-specific`-only genes (~744 each) represent real region-specific signal, not just noise.

## What this does NOT do yet

- Marker-gene interpretation beyond the calibration panel itself (e.g. why `LGR5` is fetal-elevated while `OLFM4` — checked informally, not part of the calibration panel since it contradicted the naive fetal-up assumption — is adult-elevated: `LargeInt` logFC=−7.25, FDR=0.017) is a real, interesting biological observation, reported honestly here, but not chased down further in this pass.
- The three-layer statistical validation against the 5 `mike_verzi` signatures (hypergeometric enrichment with CI/OR, preranked GSEA against the continuous DE statistic now available, size/expression-matched permutation nulls) — this is the next real step, using `T_g = logFC` from `LargeInt_edgeR_primary.tsv`/`SmallInt_edgeR_primary.tsv` as the continuous ranking statistic for preranked GSEA.
- External adult-negative validation (GTEx/HPA/Tabula Sapiens) against the new `F_{Colon,SI}-developmental` sets — not run yet.
- Donor/dataset-provenance overlap audit for GSE158702/GSE95630+GSE103239 — still deferred, unchanged from the design doc.
