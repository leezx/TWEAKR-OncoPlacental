# Dataset: HGNC symbol ↔ Ensembl gene ID mapping

- Source: HGNC (genenames.org) custom download, user-supplied URL
- Role in project: resolves the gene-ID convention gap found in Step 1 ([[TWEAKR-Worklog#Step 1 results]] / `results/01_inventory/SUMMARY.md` Finding #1) — blocks the D-shared/F-specific/P-specific pseudobulk merge until placental (symbol) and HDMA (mixed symbol/Ensembl) datasets share one gene-ID space.

## What's on disk

`DATA/1.Databases/HGNC_gene_id_mapping/`:
- `raw/hgnc_custom_download.tsv` (8.1MB, 50,321 rows, all HGNC statuses)
- `processed/v0.1/hgnc_symbol_ensembl_map.tsv` (45,032 `Approved`-status rows; 42,348 with a resolvable Ensembl ID) — columns: `hgnc_id`, `symbol`, `ensembl_id`, `previous_symbols`, `alias_symbols`
- Also pushed to Argos (`~/DATA/1.Databases/HGNC_gene_id_mapping/`) for the Step 2 pipeline

## Correction to Step 1's coarse gene-ID characterization

Step 1's inventory only checked whether `rownames(obj)[1]` started with `ENSG` — a single-gene spot check, not a real characterization. Verified against the actual full gene lists for **all 7 HDMA organs** (Argos job 3620286, `scripts/02_gene_id_mapping/extract_all_organ_gene_lists.sh`), not just Adrenal as in the first pass: HDMA gene names are **not purely Ensembl**, they're a mix in every organ — native-symbol fraction ~68–81%, ENSG fallback fraction ~19–32% (e.g. Adrenal 23.5%, Skin 31.6%). The placental h5ad datasets remain 100% symbol, as Step 1 found. The union of distinct ENSG IDs across all 7 organs is 11,015.

## Coverage check (real, not assumed — revised after PR #4 review, now full-organ + real biotype query)

Of the 11,015 union `ENSG...` fallback IDs:
- **1,023 (9.3%) resolve** via this HGNC table (any status, both Ensembl-ID columns)
- **9,992 (90.7%) don't** — genuinely absent from HGNC entirely, not a gap in this table.

**Biotype of the 9,992 unresolved IDs**, queried against the Ensembl REST API (`scripts/02_gene_id_mapping/query_biotype.py`, `results/02_gene_id_mapping/union_ensg_biotypes.tsv`) instead of assumed:

| biotype | count | % |
|---|---|---|
| lncRNA | 8,808 | 88.2% |
| NOT_FOUND_IN_ENSEMBL | 724 | 7.2% |
| protein_coding | **459** | **4.6%** |
| TR_V_gene | 1 | 0.0% |

459 real protein-coding genes is **not negligible** — the first-pass "almost certainly" language overclaimed by not measuring this. Checked whether these 459 have an Ensembl-native `display_name` bypassing HGNC entirely: only 7/459 do; the other 452 are genuinely unnamed loci in Ensembl itself (novel/predicted genes, readthrough transcripts) — so no symbol source has them, confirming "not a coverage bug" but with the actual count now stated rather than assumed.

**Handling**: none of these are dropped. `scripts/02_gene_id_mapping/build_canonical_feature_map.py` builds a per-organ `canonical_feature_map.tsv` that keeps every original feature — unresolved ENSG IDs keep their own ID as `canonical_symbol` (`mapping_status = unmapped_kept_as_ensembl_id`), so nothing is silently lost. Also checked for **symbol collisions** created by the mapping (two original features landing on the same `canonical_symbol` within one organ): found **29 across all 7 organs** (3–5 per organ, e.g. `PDE8B` present as both a native symbol and `ENSG00000284762`) — see `collision_report.tsv`. Small but real; Step 3 pseudobulk needs an explicit aggregation rule for these before summing counts by symbol (not yet decided).

## Full record

`DATA/1.Databases/HGNC_gene_id_mapping/link.md`
