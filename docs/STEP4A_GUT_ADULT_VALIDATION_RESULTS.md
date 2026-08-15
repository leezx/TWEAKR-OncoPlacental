# Step 4a external adult-expression audit: results

Executes the design approved in `docs/STEP4A_GUT_ADULT_VALIDATION.md`
(PR #23, APPROVE after 3 real review rounds). Job provenance, per file
(all stderr empty, all pulled back and md5-verified byte-exact from
whichever job actually produced them — PR #24 went through 2 real
review rounds that each required a real re-run, not just a doc edit):
GTEx outputs (Check 2) — **job 3620899** (initial, unaffected by either
fix below); Check 1's permutation output
(`check1_tabula_sapiens_permutation_null.tsv`) — **job 3620901**
(round-1 fix: a single-donor degenerate-test bug in the original
3620899 run, see Check 1 below); Check 3's direct-flags output
(`check3_tabula_sapiens_direct_flags.tsv`) — **job 3620903** (round-2
fix: a mutually-exclusive gene-set label had silently dropped dual
membership for `TRIM71`, see Check 3 below; the 54-flag headline number
is unchanged by this fix). Two different scientific questions
throughout, per the design — never conflated:

- **F-arm** (`F_Colon-developmental`, `F_SI-developmental`): **adult-expression / adult-specificity audit**. A flag is a red flag worth manual scrutiny, not proof the gene is a false positive (the frozen definition never required near-zero adult expression). A non-flag is inconclusive, not proof of adult-negativity.
- **D/P-arm** (`D_Colon-shared`, `D_SI-shared`, `P_Colon-specific` ∪ `P_SI-specific`): GTEx = **construction-consistency audit** (these genes already inherited whole-body GTEx+HPA adult-exclusion evidence from `P_developmental`'s own construction); Tabula Sapiens = genuinely independent held-out check.

## Check 2 — GTEx bulk

**Gene-ID mapping** (frozen contract, Ensembl-primary for GTEx): `F_Colon-developmental` 1,435/1,456 testable (1,434 via Ensembl, 1 via symbol fallback, 21 `NOT_TESTABLE`); `F_SI-developmental` 1,441/1,456 testable (15 `NOT_TESTABLE`); D/P union (92 gene-set-rows, 84 unique genes) 100% testable.

**F-arm adult-expression audit** — most fetal-up genes ARE detected (median bulk TPM≥1) in adult colon/SI, which is expected and not itself concerning (see reframing above); the informative number is *how high* within the adult tissue's own distribution:

| Set | Tissue | Below floor | Detected | Median percentile (of detected) | >90th percentile (% of detected) | >75th percentile (% of detected) |
|---|---|---|---|---|---|---|
| `F_Colon-developmental` | Colon_Sigmoid | 343 (23.9%) | 1,092 (76.1%) | 52.8 | 143 (13.1%) | 302 (27.7%) |
| `F_Colon-developmental` | Colon_Transverse | 332 (23.1%) | 1,103 (76.9%) | 47.5 | 93 (8.4%) | 233 (21.1%) |
| `F_SI-developmental` | Small_Intestine_Terminal_Ileum | 283 (19.6%) | 1,158 (80.4%) | 53.4 | 128 (11.1%) | 305 (26.3%) |

Honest reading: the typical detected fetal-up gene sits near the *middle* of the adult tissue's own expression distribution (median percentile ≈48–53), not at an extreme — consistent with "fetal significantly higher than adult" without adult expression being trivial. A real minority — 8–13% of *detected* genes (the table's `>90th percentile` column) — sit above the 90th percentile of the adult tissue's own detected distribution; as a fraction of all *testable* genes (the more conservative denominator) that's ~6.5–10% (143/1,435=10.0% Colon_Sigmoid, 93/1,435=6.5% Colon_Transverse, 128/1,441=8.9% SI). These are the genuine adult-expression-burden candidates worth flagging for manual review if a stricter "F-specific" tier is ever wanted later; full gene-level list in `gtex_F_adult_expression_audit.tsv`.

**D/P construction-consistency audit** — reused the exact frozen `P_developmental` whole-body criterion (`pct_cut=25, quorum=all_but_1`) unchanged: **92/92 tested gene-set-rows (84 unique genes) still pass**. This is the expected, reassuring result — nothing about `P_developmental`'s own construction or GTEx has changed since Step 4, so the frozen gut D/P membership (which is just a re-partition of the same fixed 84-gene `P_developmental_primary84` by the new gut-specific F sets) should still be internally consistent, and it is.

## Check 1 — Tabula Sapiens, `Large_Intestine`, epithelial-restricted (`F_Colon-developmental` only)

1,344/1,456 genes testable (112 `NOT_TESTABLE`). 13 of 15 donor×epithelial-cell-type samples pass the ≥20-cell floor (2 donors, `TSP14`/`TSP2`), spanning 11 epithelial cell types.

**Bug found and fixed before these numbers were trustworthy**: the first run treated all 11 cell types as one inferential family, but 9 of them are backed by only **1** eligible donor — the cross-donor-consistent flag criterion requires ≥2 donors by definition, so for those 9 a flag is *structurally impossible*, not evidence of anything. Both the observed data and every permutation draw are identically `False` for those columns, producing a degenerate `p=1.0` that isn't a real test — including them inflated the reported "9/11 clean" story with fake, uninformative tests and understated the true (smaller) FDR family. **Fixed**: the permutation/FDR family is now restricted to the 2 cell types where a hit is even structurally possible.

**Primary evidence (direct flags)**: 14,784 (gene, cell type) pairs tested across all 11 cell types (11 including the 9 single-donor ones, reported honestly as `single_donor_evidence=True` in the raw output — not used for inference), 1,009 cross-donor-consistent adult-expression flags, concentrated entirely in the 2 genuinely-testable cell types:

| Cell type | Eligible donors | Observed flag rate | Permutation null median | Empirical p | FDR (2-test family) |
|---|---|---|---|---|---|
| `paneth cell of epithelium of large intestine` | 2 | 38.8% | 38.7% | 0.445 | 0.717 (ns) |
| `transit amplifying cell of large intestine` | 2 | 36.2% | 36.5% | 0.717 | 0.717 (ns) |
| *9 other epithelial cell types* (enterocyte, immature/mature enterocyte, goblet cell, intestinal crypt stem cell ×2, enteroendocrine cell, tuft cell, large intestine goblet cell) | *1 each* | 0.0% | — | *structurally untestable, excluded from FDR family* | — |

**Honest reading**: the only two epithelial cell types with enough donor replication (`TSP14` + `TSP2` both eligible) to run a real cross-donor-consistent test — Paneth cells and transit-amplifying cells — show **no statistically-supported adult-expression anomaly**: observed rates are indistinguishable from the size-and-detectability-matched permutation null (enrichment ≈1.0×, p=0.445/0.717, not significant even nominally). The other 9 epithelial cell types, including the most CRC-relevant ones (enterocytes, goblet cells, intestinal stem cells), show zero flags in the direct evidence — but this is now correctly reported as **single-donor evidence, not a statistically tested negative** (Tabula Sapiens' `Large_Intestine` file has only 2 donors total, a real data limitation, not a design choice). The genuinely-testable subset is smaller than originally reported, but its conclusion — no adult-expression anomaly detected — is unchanged.

## Check 3 — Tabula Sapiens, all 5 organs, D/P construction-consistency + independent audit

84 unique genes tested (`D_Colon-shared` 5, `D_SI-shared` 4, `P_Colon-specific` ∪ `P_SI-specific` union-deduplicated 83 — these are always exactly the same 84 genes as `P_developmental_primary84`, since `D_{region} ∪ P_{region}-specific = P_developmental` by construction for any region). 8,568 (organ, cell_type, gene) triples tested, **54 cross-donor-consistent flags** (`D_Colon-shared` 2, `D_SI-shared` 14, `P_union` 38) — concentrated in a small set of genes already known from Step 5 (`RPA4`, `ZNF257`, `ZNF850`, `ZNF695`, `ZNF114`, `GJB7`, `ERVFRD-1`, `TMEM191C`, `GCM1`, `HTRA4`, `LAIR2`, `KISS1`, `TSKS`, `ZBTB8B`).

**Verified directly, not assumed**: this result numerically reproduces Step 5's original pan-organ whole-body validation exactly (same 8,568 triples, same 54 flags, same flagged genes — confirmed `P_developmental_primary84.txt` is literally identical to Step 5's `D_shared_FINAL ∪ P_specific_FINAL`, `diff` clean). This is expected, not a bug: Check 3 tests the same fixed 84-gene `P_developmental` set against the same Tabula Sapiens data with the same methodology as Step 5 — the gut-specific D/F/P re-anchoring only changes how these 84 genes are *labeled* (Colon-D vs Colon-P-specific vs SI-D vs SI-P-specific), not which 84 genes they are. A useful internal-consistency confirmation that this compute pipeline reproduces a known-good prior result correctly, though it adds no new information beyond what Step 5 already established for these particular 84 genes.

**Bug found and fixed** (round-2 review): the output's `gene_set` label was originally computed as mutually-exclusive ("`D_Colon-shared` else `D_SI-shared` else `P_union`"), which silently drops dual membership for any gene in more than one set — `TRIM71` is the one gene shared by `D_Colon-shared` and `D_SI-shared` and was only ever labeled `D_Colon-shared`. Doesn't change the 54-flag headline (`TRIM71` isn't among them, reconfirmed after the fix), but a future region-specific summary sliced by that label alone would have silently misrepresented `TRIM71`'s membership. **Fixed**: added explicit, non-exclusive `in_D_Colon_shared`/`in_D_SI_shared` boolean columns (matching the pattern already used for `in_P_Colon_specific`/`in_P_SI_specific`); `TRIM71` now correctly shows both `True`. Re-run (job 3620903), pulled back and md5-verified.

## Bottom line

- **`F_Colon-developmental`** shows no statistically-supported adult-expression anomaly in either of the 2 adult epithelial cell types with enough donor replication to genuinely test (Tabula Sapiens, genuinely independent) — the most direct, single-cell-resolution external check available, and it comes back clean. The other 9 epithelial cell types (including the most CRC-relevant ones) have only single-donor coverage in this dataset and are reported as such, not claimed as tested negatives. GTEx bulk shows the expected middle-of-distribution detection pattern for most genes, with a real ~10% minority worth flagging for future scrutiny, not evidence against the panel as a whole.
- **`F_SI-developmental`** has no single-cell check available (Tabula Sapiens has no Small_Intestine organ, a stated design limitation) — GTEx terminal-ileum bulk shows a similar pattern to Colon, with the same partial-anatomical-coverage caveat.
- **D/P sets** remain internally consistent with the whole-body adult-exclusion evidence that helped define them in the first place (GTEx, expected) and reproduce Step 5's already-established, small set of genuine adult-expression flags exactly (Tabula Sapiens, genuinely independent but not new information for this particular 84-gene set).

No gene's membership in any frozen set was changed by this audit — purely additive reporting, per the design's own "does not gate anything" contract.

## Outputs

`results/04a_dfp_gut/adult_validation/`:
- `gtex_F_adult_expression_audit.tsv`, `gtex_DP_construction_consistency_audit.tsv`, `gtex_gene_id_mapping_summary.tsv`
- `check1_tabula_sapiens_direct_flags.tsv`, `check1_tabula_sapiens_permutation_null.tsv`
- `check3_tabula_sapiens_direct_flags.tsv`, `check3_gene_id_mapping_summary.tsv`, `tabula_sapiens_gene_id_mapping_summary.tsv`
