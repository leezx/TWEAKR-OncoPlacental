# Dataset: HGNC symbol ↔ Ensembl gene ID mapping

- Source: HGNC (genenames.org) custom download, user-supplied URL
- Role in project: resolves the gene-ID convention gap found in Step 1 ([[TWEAKR-Worklog#Step 1 results]] / `results/01_inventory/SUMMARY.md` Finding #1) — blocks the D-shared/F-specific/P-specific pseudobulk merge until placental (symbol) and HDMA (mixed symbol/Ensembl) datasets share one gene-ID space.

## What's on disk

`DATA/1.Databases/HGNC_gene_id_mapping/`:
- `raw/hgnc_custom_download.tsv` (8.1MB, 50,321 rows, all HGNC statuses)
- `processed/v0.1/hgnc_symbol_ensembl_map.tsv` (45,032 `Approved`-status rows; 42,348 with a resolvable Ensembl ID) — columns: `hgnc_id`, `symbol`, `ensembl_id`, `previous_symbols`, `alias_symbols`
- Also pushed to Argos (`~/DATA/1.Databases/HGNC_gene_id_mapping/`) for the Step 2 pipeline

## Correction to Step 1's coarse gene-ID characterization

Step 1's inventory only checked whether `rownames(obj)[1]` started with `ENSG` — a single-gene spot check, not a real characterization. Verified against the actual full gene list for `HDMA_Adrenal` (loaded on Argos): HDMA gene names are **not purely Ensembl**, they're a mix — **76.5% already gene symbols** (19,368/25,314), **23.5% Ensembl IDs** (5,946/25,314) used only where no resolved symbol was available at the source pipeline's processing time. The placental h5ad datasets remain 100% symbol, as Step 1 found.

## Coverage check (real, not assumed)

Of the 5,946 `ENSG...` fallback IDs in `HDMA_Adrenal`:
- **651 (~11%) resolve** via this HGNC table
- **5,295 (~89%) don't** — checked against the full raw dump (any HGNC status, both Ensembl-ID columns), same result. These aren't a gap in this specific table; they're **genuinely absent from HGNC entirely** — predominantly lncRNA/pseudogene/novel-locus Ensembl entries that Ensembl's gene model includes but HGNC hasn't curated an approved symbol for. This is expected: HGNC deliberately curates narrower than Ensembl's full annotation, especially for non-protein-coding loci.

**Practical read for Step 2**: don't chase the unresolvable ~89% with a different source — HGNC not having them is the ground truth. The protein-coding markers that actually matter for the trophoblast/placental signature work (`CGA`, `CGB`, `ERVFRD-1`, etc.) are almost certainly already in the 76.5% that's already a plain symbol, not in the ENSG-only tail. Genes that stay unresolved after this mapping should just be dropped from any symbol-keyed merge (or kept as their Ensembl ID if the downstream step doesn't require symbols) rather than treated as a blocker.

## Full record

`DATA/1.Databases/HGNC_gene_id_mapping/link.md`
