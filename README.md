# TWEAKR-OncoPlacental

Code and documentation for building a CNS-grade human developmental
reference (fetal-somatic-specific / placental-specific / shared-developmental
modules) and testing whether the CRC "Oncofetal" cell state is actually a
mixture of an embryonic/fetal program and a separate extraembryonic/
placental program.

**Start here:**
- [`docs/PROJECT_SUMMARY.md`](docs/PROJECT_SUMMARY.md) — structured
  "what was done, how, why, what's been answered" summary. Read this
  first for the current state of the project.
- [`Worklog.md`](Worklog.md) — full chronological log (every review
  round, every bug caught and fixed, the progress tracker). Read this for
  the detailed history or to resume active work.

## Current status (see `docs/PROJECT_SUMMARY.md` for details)

**Step 4 (D-shared / F-specific / P-specific signature construction) is
CLOSED** — frozen, reviewer-approved gene sets exist:
`results/04_dfp_signature/dfp_gene_sets/{D_shared,F_specific,P_specific}_FINAL.txt`.
Next: Tier-2 validation (Tabula Sapiens) and applying the frozen
signature to real CRC Oncofetal single-cell/spatial data.

## Scope

This repo covers building and freezing the normal-development reference
(placenta/trophoblast + fetal-somatic + adult) entirely from normal
tissue — no cancer data touches signature *definition*. Applying the
frozen signature to CRC Oncofetal cells, and the remaining Q2–Q6
questions, come after.

## Data

**No large data files live in this repo.** All raw/processed/result data
lives under `/Volumes/Stelligen_SSD/Stelligen/DATA/<data_type>/<dataset_id>/`,
following the workspace-wide convention in `DATA/README.md`. This repo
only holds:

- `datasets/<dataset_id>/dataset.md` — manifest/metadata for each dataset
  (source URLs, sizes, what was downloaded and why, what wasn't)
- `scripts/<step>/` — analysis code, one subdirectory per pipeline step
  (`01_inventory`, `02_gene_id_mapping`, `03_pseudobulk_prep`,
  `04_dfp_signature`) — every script that ran on real data lives here,
  committed alongside the results it produced
- `results/<step>/` — every real compute output (gene lists, calibration
  tables, audit docs), mirroring the `scripts/` step structure
- `docs/` — design docs, method contracts, and this summary

See `datasets/` for the datasets acquired so far.

## Related repos

- `REPOS/TWEAKR_Oncofetal-workflow` — separate, pre-existing repo; mostly
  `GSE131696` placenta scRNA analysis. Not migrated here; left as-is.
