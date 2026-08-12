# Dataset: Human Development Multiomic Atlas (HDMA) — fetal 12-organ reference

- Paper: Multiomics and deep learning dissect regulatory syntax in human development
- Authors: Liu*, Jessa*, Kim*, Ng*, …, Kundaje+, Farh+, Greenleaf+ (Nature, 2026)
- DOI: `10.1038/s41586-026-10326-9`
- Role in project: "第四梯队" fetal negative/reference universe ([[TWEAKR-Worklog#Placenta数据集]]) — the `Fetal-specific` side of the `Fetal-specific` vs `Placenta-specific` orthogonal axis (Q1, [[2026-GPT-TWEAKR-Oncofetal#定义清楚Placenta的问题]]), so the project isn't just comparing placenta to adult tissue.

**Not the same paper as** `2026_human_maternal_fetal_Nature` (Wang et al., placenta-specific, DOI `10.1038/s41586-026-10316-x`) — easy to confuse since both are "Nature 2026", do not conflate.

## Status history

1. User originally dropped just this repo's own GitHub README into `DATA/scRNAseq/HumanDevelopmentMultiomicAtlas/` — no actual data, a placeholder.
2. This session: the full collection is **~2.2TB** across ArchR projects, bigwigs, ChromBPNet models/contribution-scores, motif compendium, BPCells objects, fragments+count matrices, and RNA Seurat objects — way too much to pull wholesale, and most of it (chromatin models, tracks, motifs) is irrelevant to RNA-level signature construction.
3. Resolved the actual file-level structure via the Zenodo API: the "RNA Seurat Objects Part1-4" Zenodo records each bundle a few organs, but **files within each record are individually downloadable** — no need to pull an entire ~46-49GB "Part" to get one organ.

## Full per-organ RNA Seurat object listing (resolved via Zenodo API, 2026-08-12)

| Organ | Size | Germ layer | Zenodo record | Downloaded? |
|---|---|---|---|---|
| Adrenal | 0.65GB | mesoderm | 15014813 | ✅ |
| Thyroid | 1.58GB | endoderm | 15014809 | ✅ |
| Spleen | 7.41GB | mesoderm | 15014813 | ✅ |
| Thymus | 10.62GB | endoderm | 15014809 | ✅ |
| Skin | 14.36GB | ectoderm | 15014809 | ✅ |
| Eye | 14.24GB | ectoderm | 15014815 | — |
| Liver | 17.31GB | endoderm | 15014815 | ✅ |
| Brain | 17.71GB | ectoderm | 15014815 | — |
| StomachEsophagus | 19.92GB | endoderm (GI-adjacent) | 15014809 | ✅ |
| Muscle | 19.96GB | mesoderm | 15014977 | — |
| Heart | 27.79GB | mesoderm | 15014977 | — |
| Lung | 38.09GB | endoderm | 15014813 | — |

**Note: this atlas has no intestine/colon organ.** It's a broad-spectrum "is this gene just generically fetal, in any organ" negative reference, not a gut-specific one — the project's existing fetal-intestine references (e.g. via the CRC atlas work) cover that role separately.

## Decision (user-selected, 2026-08-12)

Downloaded 7 of 12 organs (≈57GB): Adrenal, Thyroid, Spleen, Thymus, Liver, Skin, StomachEsophagus — covers all 3 germ layers plus a GI-adjacent organ, while skipping the 5 largest (Lung/Heart/Muscle/Brain/Eye, ≈133GB combined). Can add any of the skipped 5 later via the same per-file Zenodo URL pattern if a specific organ becomes analytically relevant.

**Deliberately skipped entirely** (out of scope for RNA-level signature construction): ArchR projects (166GB), Bigwigs (392GB), ChromBPNet contribution scores (~473GB) + models (89GB), BPCells objects (40GB), TF-MoDISco motifs (22GB), Fragments+count matrices (~334GB — raw multiome inputs, redundant with the already-processed Seurat objects), cell metadata/regulatory bundle (13.7GB). Revisit only if chromatin-level corroboration becomes necessary (parallel to the same open question flagged for `Arutyunyan2023_MFI_multiome`).

## How to get more organs later

Each Zenodo record's files are individually addressable:

```
https://zenodo.org/api/records/<record_id>/files/<Organ>_RNA_obj_clustered_final.rds/content
```

Record IDs: `15014813` (Adrenal/Spleen/Lung), `15014815` (Brain/Eye/Liver), `15014977` (Heart/Muscle), `15014809` (Skin/StomachEsophagus/Thymus/Thyroid). Full index of all HDMA data types (not just RNA Seurat objects): `table_s14.tsv` from https://greenleaflab.github.io/HDMA/tables/table_s14.tsv.

## Full record

`DATA/scRNAseq/HumanDevelopmentMultiomicAtlas/link.md`
