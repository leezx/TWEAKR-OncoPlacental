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
