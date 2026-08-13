# Step 1 — Inventory pass results

Ran via SGE job `3620272` on Argos (`argos-codex` env). All 15/15 files succeeded — verified directly against each JSON's `error` field (not the job's own exit status; see `Worklog.md` for why that distinction matters here). Raw JSON per dataset in this directory.

## Dataset shapes and gene-ID convention

| Dataset | Cells | Genes | Gene IDs |
|---|---|---|---|
| Arutyunyan primary_tissue | 325,665 | 30,800 | **symbol** |
| Arutyunyan organoid PTO | 26,853 | 23,281 | symbol |
| Arutyunyan organoid TSC | 9,957 | 22,523 | symbol |
| Arutyunyan organoid Fig3 | 37,480 | 23,281 | symbol |
| Nature2026 scPlacenta_host | 191,735 | 32,981 | symbol |
| Nature2026 snRNA_raw_counts | 243,103 | 36,601 | symbol |
| VentoTormo decidua-v3 | 70,325 | 31,764 | symbol |
| Greenbaum RNA (mtx) | 36,456 cells | 36,546 | symbol (from earlier compute-feasibility check) |
| HDMA Adrenal | 2,883 | 25,314 | **Ensembl** |
| HDMA Thyroid | 9,299 | 26,163 | Ensembl |
| HDMA Spleen | 33,177 | 32,513 | Ensembl |
| HDMA Thymus | 41,702 | 33,648 | Ensembl |
| HDMA Skin | 57,801 | 34,185 | Ensembl |
| HDMA Liver | 71,807 | 32,748 | Ensembl |
| HDMA StomachEsophagus | 83,663 | 33,035 | Ensembl |

**Finding #1 — gene-ID mismatch, must be resolved before any cross-dataset comparison**: every placental/trophoblast dataset (h5ad, all 7 + Greenbaum) uses gene **symbols**; every HDMA fetal-somatic dataset (RDS) uses **Ensembl IDs**. Building the D-shared/F-specific/P-specific comparison means either mapping HDMA to symbols or the placental side to Ensembl — needs a decision + a mapping table (e.g. via `biomaRt`/a static Ensembl-to-symbol table) before Step 2.

## Cell-type / annotation columns (the actually useful part)

| Dataset | Annotation column(s) | Trophoblast-relevant labels present |
|---|---|---|
| Arutyunyan primary_tissue | `coarse_annot` (8 classes) / `cell_type` (41 fine types) | `coarse_annot == "Trophoblast"` (75,042 cells) is a clean top-level filter; fine `cell_type` has SCT/VCT/VCT_p/VCT_CCC/VCT_fusing/EVT_1/EVT_2/iEVT/eEVT/GC |
| Arutyunyan organoid PTO/TSC/Fig3 | `cell_annotation` / `final_annot_v9` | Pure trophoblast by construction (organoid) — all values are VCT/EVT/SCT-lineage subtypes, no filtering needed |
| Nature2026 scPlacenta_host | `major_class` (15 classes) / `celltype_fullname` (33 fine) / `origin` (Fetal/Maternal/Unknown) | `major_class` has SCT/VCT/EVT directly. **Caution**: `origin=="Fetal"` here means fetal *side of the placenta* (includes fetal endothelium `fEC`, Hofbauer `HB`) — not the same thing as "fetal somatic organ tissue" from HDMA. Don't conflate the two when building F-specific vs. P-specific. |
| Nature2026 snRNA_raw_counts | only `ID`/`dataset`/`BC` | **No usable annotation in this file** — likely needs joining to `scPlacenta_host`'s `celltype_fullname` via barcode (`BC`), or is a redundant raw companion. Flag for Step 2: check `BC` overlap between the two 2026-Nature files before assuming this one is independently usable. |
| VentoTormo decidua-v3 | `CellType` (32 types) / `Location` (Blood/Decidua/Placenta) | Has VCT/EVT/SCT directly in `CellType`, confirming (from the earlier data-audit PR) this dataset is not decidua-only |
| Greenbaum (cluster.csv, not obs) | `cell_type` | vCTB/STB/EVT/EVT-progenitor/STB-progenitor + Hofbauer/fibroblast/endothelial/erythroblast — clean, already checked in an earlier run |
| HDMA (all 7 organs) | `annotv1` (coarse: `_end`/`_epi`/`_imm`/`_str` suffix = endothelial/epithelial/immune/stromal) / `annotv2` (fine, organ-specific) | No trophoblast labels present, as expected — confirms these are clean fetal-somatic-organ references with no placental contamination |

**Finding #2 — nomenclature is consistent enough to build a shared trophoblast filter across 5 independent studies without inventing a new vocabulary**: `VCT`/`EVT`/`SCT` (and close variants: `iEVT`/`eEVT`/`proEVT` for EVT subtypes, `VCT_p`/`VCT_CCC`/`VCT_fusing`/`cycling_VCT`/`avVCT` for VCT subtypes, `SCT_A`/`SCT_B`/`proSCT` for SCT subtypes) appear across Arutyunyan, Nature2026, VentoTormo, and Greenbaum independently. This is a strong, literature-consistent trophoblast vocabulary — Step 2 can define one canonical trophoblast filter and apply near-identical logic to all 5 placental datasets rather than 5 bespoke ones.

**Finding #3 — `Nature2026_snRNA_raw_counts` needs a follow-up check**, not full inventory data as-is (see table above).

## What this unblocks for Step 2

The three-way D-shared/F-specific/P-specific comparison (Q1) can now be scoped concretely:
- **P (placental/trophoblast) pseudobulk**: filter each of the 5 placental datasets to trophoblast-lineage cells using the columns above, aggregate to per-sample/per-donor pseudobulk.
- **F (fetal-somatic) pseudobulk**: HDMA's 7 organs, all cells (no filtering needed — these are already whole-organ, non-placental).
- **Gene-ID harmonization**: resolve symbol-vs-Ensembl before merging (Finding #1) — first concrete task for Step 2.
- **Adult reference**: still the open gap from earlier in this project — not resolved by this inventory pass, still needed before the full 3-way comparison.
