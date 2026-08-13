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

## Cell-type coverage — resolved 2026-08-13 (data-audit PR)

Despite the `decidua-v3.h5ad` filename, this is **not** decidua-only. Inspected via `h5py` directly (`obs['CellType']`, 32 categories):

- Trophoblast: `EVT` 3,626 + `SCT` 1,261 + `VCT` 9,479 ≈ **14,366 cells**
- Decidual stromal/immune: `dS1-3`, `dP1-2`, `dM1-3`, `dNK1-3`/`dNK p`, etc.
- Fetal fibroblasts (`fFB1-2`), endothelium (`Endo (f)`/`Endo (m)`/`Endo L`), and more (32 types total)
- `obs['Location']`: `Blood` 11,266 / `Decidua` 40,512 / `Placenta` 18,547

Safe to use as an independent trophoblast reference, not just a maternal/decidual-negative one.

## Full record

`DATA/scRNAseq/VentoTormo_Nature_2018/link.md`
