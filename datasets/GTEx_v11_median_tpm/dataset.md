# Dataset: GTEx v11 adult gene-level median TPM by tissue

- Source: GTEx Portal, Google Cloud Storage `adult-gtex` bucket
- Role in project: Tier-1 adult normal reference — the prerequisite the ChatGPT reviewer flagged repeatedly as missing before formal D-shared/F-specific/P-specific (D/F/P) signature construction can start. Without an adult baseline there's no way to separate "developmentally re-expressed" genes from "just normally expressed in adults too."

## Two-tier adult reference design (per project discussion, 2026-08-13)

- **Tier 1 (this dataset, mandatory, used now to define signatures)**: GTEx v11 bulk median TPM, 68 adult tissues. Used two ways, deliberately different designs:
  - **F-specific** = organ-matched: fetal HDMA organ vs its matched adult GTEx tissue(s).
  - **P-specific** = whole-body: trophoblast vs **all 68** GTEx tissues, not just the 7 HDMA-matched organs — otherwise a gene absent from the 7 matched organs but high in, say, testis would be wrongly called placenta-specific.
  - **CRC sanity check**: normal adult colon columns (`Colon_Sigmoid`, `Colon_Transverse`+Mucosa/Muscularis) — the eventual target is CRC malignant epithelial cells, so a "placenta-specific" gene that's already high in normal adult colon epithelium isn't a useful discriminator.
- **Tier 2 (deferred, independent validation only)**: Tabula Sapiens (adult single-cell) + a dedicated adult colon epithelial single-cell atlas, held out of signature *definition* entirely so "P-specific genes are largely absent across adult cell types" is an independent post-hoc validation, not circular. Not downloaded yet — revisit after Step 3's first D/F/P version is frozen.

**Cross-platform comparison rule (fixed 2026-08-13, PR #5 round-2)**: GTEx (bulk TPM) is never directly compared in raw magnitude against HDMA/trophoblast (single-cell counts). It's used only as an **adult-exclusion reference** — is a candidate gene still meaningfully expressed in the matched tissue(s), judged by GTEx's own internal rank/percentile/threshold — combined as a separate evidence axis alongside developmental evidence computed entirely within the scRNA-seq data. Full rule: `docs/STEP3_METHOD_CONTRACT.md`. `hdma_organ_to_gtex_tissue_map.tsv` now has an explicit `role` column (`adult_negative_reference` for every row — GTEx has no placenta column, so unlike HPA there's no exclusion case needed here).

## What's on disk

`DATA/1.Databases/GTEx_v11_median_tpm/`:
- `raw/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz` (10,129,906 bytes, byte-exact vs the GCS bucket's declared object size)
- `processed/v0.1/gtex_v11_median_tpm_clean.tsv` — 74,628 genes × 68 tissues, `ensembl_id` (version stripped) + `ensembl_id_versioned` + `symbol` + all tissue columns
- `processed/v0.1/hdma_organ_to_gtex_tissue_map.tsv` — the HDMA-organ ↔ GTEx-tissue-column mapping, decided by inspecting the real 68 column names, not assumed
- Also pushed to Argos (`~/DATA/1.Databases/GTEx_v11_median_tpm/`)

## Verification

- Byte-exact: downloaded size matches the GCS bucket's declared object size (checked via the JSON listing API before downloading, not just curl's own exit code).
- Structural: `gzip -t` passes; GCT header declares 74,628 rows × 68 columns — both the raw file's actual row count and the cleanup script's independently-counted written-row count match exactly.
- Gene ID format confirmed by direct inspection: `Name` = versioned Ensembl ID, `Description` = symbol — same convention as HDMA, plugs directly into the `[[HGNC_gene_id_mapping]]` / `canonical_feature_map` pipeline from Step 2.

## Real gap found (not glossed over)

**No `Thymus` column in GTEx at all** — checked directly against all 68 tissue names, absent. Genuine collection gap: GTEx's adult donors skew older and the thymus involutes with age, so it isn't represented. There is no organ-matched adult baseline for HDMA's Thymus from this source; Step 3 needs an explicit fallback (e.g. Tabula Sapiens, if it covers thymus) for that one organ's F-specific comparison. Whole-body P-specific comparisons are unaffected.

## Full record

`DATA/1.Databases/GTEx_v11_median_tpm/link.md`
