# Dataset: Human Development Multiomic Atlas (HDMA) — fetal 12-organ reference

- Paper: Multiomics and deep learning dissect regulatory syntax in human development
- Authors: Liu*, Jessa*, Kim*, Ng*, …, Kundaje+, Farh+, Greenleaf+ (Nature, 2026)
- DOI: `10.1038/s41586-026-10326-9`
- Role in project: "第四梯队" fetal negative/reference universe ([[TWEAKR-Worklog#Placenta数据集]]) — needed to build the `Fetal-specific` vs `Placenta-specific` orthogonal axis (Q1/D-shared/F-specific/P-specific framework, [[2026-GPT-TWEAKR-Oncofetal#定义清楚Placenta的问题]]) instead of just comparing placenta to adult tissue.

**⚠️ Placeholder only — no data downloaded yet.** What exists in `DATA/scRNAseq/HumanDevelopmentMultiomicAtlas/raw/` right now is just the project's own GitHub README (data-availability section + code map), saved as `source_notes.md`. This is not the same paper as `2026_human_maternal_fetal_Nature` (Wang et al., a different Nature 2026 placenta-specific paper, DOI `10.1038/s41586-026-10316-x`) — easy to confuse since both are "Nature 2026" — do not conflate them.

## Source / how to actually get data

- Docs: https://greenleaflab.github.io/HDMA/
- GitHub: https://github.com/GreenleafLab/HDMA
- Data: https://zenodo.org/communities/hdma (fragment files, count matrices, cell annotations, caCRE annotations, ChromBPNet models, motif lexicon, genomic tracks — all separate Zenodo records)
- Full data-type → Zenodo-record index: `tables/table_s14.tsv` in the repo (Table S14 of the manuscript)
- Raw-genomic-data metadata: https://doi.org/10.5281/zenodo.17259745

This atlas is much bigger in scope than what's needed here — it includes trained ChromBPNet models, motif compendia, and WashU genome-browser tracks alongside the actual expression/accessibility data. **Do not bulk-download the whole Zenodo community.** Next step: pull `table_s14.tsv`, find the Seurat object / count matrix records for the specific fetal somatic organs relevant to the project (likely a handful, not all 12), and download just those.

## Full record

`DATA/scRNAseq/HumanDevelopmentMultiomicAtlas/link.md`
