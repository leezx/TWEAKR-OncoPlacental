# Dataset: Vento-Tormo et al. 2018, Nature — Single-cell reconstruction of the early maternal-fetal interface

- Journal: Nature, 2018
- Role in project: "第二梯队" independent replication ([[TWEAKR-Worklog#Placenta数据集]]) — the classic, fully independent maternal-fetal interface atlas; strengthens any gene that's also trophoblast-specific in the 2023/2024/2026 atlases already in this project

**Downloaded manually by the user**, not by this session's automated workflow. This file documents what's there after the fact.

## Source

- Interactive browser: http://data.teichlab.org
- ArrayExpress accessions: `E-MTAB-6701` (droplet-based), `E-MTAB-6678` (Smart-seq2), `E-MTAB-7304` (whole-genome sequencing, not relevant here)
- CellPhoneDB: https://www.cellphonedb.org (companion resource from the same paper, not pulled)

## What's on disk

`DATA/scRNAseq/VentoTormo_Nature_2018/raw/`:

- `decidua-v3.h5ad` (800MB) — processed AnnData, downloaded from the data.teichlab.org portal
- `source_notes.md` — the paper's data-availability paragraph, saved verbatim

## Open question

The filename `decidua-v3.h5ad` suggests this might be **decidua (maternal) cells only**, not trophoblast. This paper's title is about the whole maternal-fetal interface (which does include trophoblast), but the specific portal file downloaded needs to be checked before treating it as a trophoblast reference. **TODO next session:** load once, inspect `obs` cell-type column, document findings in `link.md` (same pattern as `DATA/scRNAseq/2026_human_maternal_fetal_Nature/raw/README.md`). If it turns out to be decidua-only, it's still useful as a maternal-negative-reference, just not as a trophoblast-positive reference.

## Full record

`DATA/scRNAseq/VentoTormo_Nature_2018/link.md`
