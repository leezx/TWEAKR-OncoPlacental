# Dataset: Han et al. 2020, Nature — Human Cell Landscape (HCL)

- Paper: *Construction of a human cell landscape at single-cell level*
- Journal: Nature, 2020
- DOI: `10.1038/s41586-020-2157-4`
- GEO: `GSE134355`
- Portal: http://bis.zju.edu.cn/HCL/ (mirror: https://db.cngb.org/HCL/)
- Role in project: Step 9's independent, single-atlas validation track
  (Track C) for the Fetal-Placenta-Adult developmental ternary map — the
  one candidate atlas containing Fetal intestine, Placenta, and Adult
  intestine on one uniformly-processed platform (Microwell-seq),
  avoiding the cross-dataset/cross-platform batch risk this project's
  other developmental-reference tracks (built from separate studies)
  cannot fully escape. See `docs/STEP9_DEVELOPMENTAL_TERNARY_MAP.md`.

## Access

Public, no login/authentication required. Single supplementary file:
`GSE134355_RAW.tar`, exact `Content-Length` 927,109,120 bytes
(byte-verified on download, not GEO's rounded MB display). 141 GSM
samples total (human tissues; the series also nominally covers a mouse
cell atlas companion dataset, not used here).

Downloaded to `DATA/scRNAseq/GSE134355/raw/GSE134355_RAW.tar`, extracted
to `DATA/scRNAseq/GSE134355/raw/extracted/` (141 `GSM*_dge.txt.gz`
files, one per sample). Archive integrity: `tar -tf` + `gzip -t` on
every extracted member — 141/141 clean, 0 failures.

## Format

Each `_dge.txt.gz` is a plain tab-separated digital gene expression
matrix: `GENE` row header + gene-symbol rows, one column per cell
barcode. Confirmed by direct inspection (not filename/format guess):
real integer raw UMI counts (values like 0, 1, 2 — genuine sparse
counts, not normalized/log-transformed).

**Real, honest data characteristic found during inventory**: every one
of the 141 files has exactly 10,000 barcode columns — HCL's own
processing pads/caps each sample's DGE export at a fixed 10,000 columns
rather than a variable barcode-rank-knee cutoff. A real (unquantified)
fraction of those 10,000 columns per sample are likely low/near-empty
droplets, not all genuine cells. Not treated as "10,000 real cells" in
any downstream use.

## Cell-type annotation (added round-1 review fix)

The GEO `RAW.tar` contains only the 141 `_dge.txt.gz` whole-tissue
matrices, no cell-type annotation. The real per-cell annotation lives
on Figshare, not GEO: article
[7235471](https://figshare.com/articles/dataset/HCL_DGE_Data/7235471),
file `HCL_Fig1_cell_Info.xlsx` (19,772,723 bytes, md5
`fe73a9b7129abb10d09dfcd355c19f12`, both verified on download) —
599,926 cells with `cellnames`/`sample`/`cluster`/`stage`/`batch`/
`donor`/`celltype` columns. Barcode suffixes in `cellnames` (after the
`.`) confirmed to be an exact subset of the corresponding `_dge.txt.gz`
file's column headers (direct membership test, not assumed).

**Real, honest finding this annotation revealed**: `Placenta1`'s
dominant cell type is **Fibroblast (72.7%)**, not trophoblast — HCL's
own annotation has **no "Trophoblast" label at all**; "Epithelial cell"
(11.4% of the sample) is the best available proxy. Downloaded and used
after a round-1 ChatGPT PR review correctly flagged that scoring the
whole-tissue pseudobulk (this project's first version) confounds
cell-composition with developmental state — see
`docs/STEP9_DEVELOPMENTAL_TERNARY_MAP.md` "Round 1 review" for the full
account, including 2 further self-caught bugs (a trailing-whitespace
bug in HCL's own `celltype` labels; a real `"AdultJeJunum"` vs
`"AdultJejunum"` naming quirk) found while implementing the fix.

Also real: `AdultTransverseColon2` in this annotation's own curated
batch structure combines what this project's GSM-level download
treated as two separate samples (`Adult-Transverse-Colon2-1`,
`Adult-Transverse-Colon2-2`) — confirmed by direct barcode set-overlap
testing against both underlying files, not assumed.

## Samples used in this project (of 141 total)

Only 14 of the 141 GSMs are used (this project's Step 9 needs Fetal
intestine, Placenta, and Adult intestine only — not a general-purpose
HCL ingestion):

| Group | Samples | n |
|---|---|---|
| Placenta | `Placenta1` | 1 (single donor — real limitation, treated as directional/exploratory support, not primary replicated evidence, matching this project's precedent for other single-donor sources) |
| Fetal intestine | `Fetal-Intestine1` through `Fetal-Intestine5` | 5 (not regionally resolved — no colon/small-intestine split, unlike this project's Gut Cell Atlas reference) |
| Adult intestine | `Adult-Ascending-Colon1`, `Adult-Sigmoid-Colon1`, `Adult-Transverse-Colon1`, `Adult-Transverse-Colon2-1`, `Adult-Transverse-Colon2-2`, `Adult-Duodenum1`, `Adult-Ileum2`, `Adult-Jejunum2` | 8 (colon + small-intestine segments pooled as one "Adult" group, to match Fetal-Intestine's own lack of regional split) |

The other 127 GSMs (other adult organs, mouse atlas) are not downloaded
into `processed/` or used anywhere in this project.
