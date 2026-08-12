# TWEAKR-OncoPlacental

Code and documentation for building a CNS-grade human developmental reference (fetal-specific / placental-specific / shared-developmental modules) and testing whether the CRC "Oncofetal" cell state is actually a mixture of an embryonic/fetal program and a separate extraembryonic/placental program.

**Start here:** [`Worklog.md`](Worklog.md) — session-by-session log of what was done, why, and what's left. Read it before doing anything else in this repo.

## Scope

This repo covers **Aim 1** of the project: acquiring and freezing the normal-development reference datasets (placenta/trophoblast + fetal somatic + adult) and building the P1 (placenta-core) / P2 (early-placenta) / P3 (placenta-specific) / D (shared-developmental) signature construction pipeline — entirely from normal tissue, no cancer data.

Aim 2 (decomposing existing Oncofetal/revCSC signatures against this reference) and Aim 3 (SPP1-TAM → TWEAK → Fn14 → YAP mechanism) come after the reference is frozen.

Full scientific framing lives in the KB, not here:
- `Zhixins-KB/3.Distill/1.Projects/TWEAKR/TWEAKR-Worklog.md`
- `Zhixins-KB/3.Distill/1.Projects/TWEAKR/2026-GPT-TWEAKR-Oncofetal.md`

## Data

**No large data files live in this repo.** All raw/processed/result data lives under `/Volumes/Stelligen_SSD/Stelligen/DATA/<data_type>/<dataset_id>/`, following the workspace-wide convention in `DATA/README.md`. This repo only holds:

- `datasets/<dataset_id>/dataset.md` — manifest/metadata for each dataset (source URLs, sizes, what was downloaded and why, what wasn't)
- `scripts/`, `notebooks/` — analysis code
- `docs/` — design docs

See `datasets/` for the datasets acquired so far.

## Related repos

- `REPOS/TWEAKR_Oncofetal-workflow` — separate, pre-existing repo; mostly `GSE131696` placenta scRNA analysis. Not migrated here; left as-is.
