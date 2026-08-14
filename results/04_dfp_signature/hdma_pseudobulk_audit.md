# HDMA per-organ, per-sample pseudobulk: build + verification

Per `docs/STEP4_STATISTICAL_DESIGN.md` section 2 (F-developmental's
positive evidence): per-organ, per-individual-`Sample` pseudobulk raw-count
matrices for all 7 downloaded HDMA organs, built via
`scripts/04_dfp_signature/build_hdma_pseudobulk.R`
(`run_build_hdma_pseudobulk.sh`, Argos qsub). This required first finding
and fixing a real gene-mapping gap — see
`hdma_rna_assay_gene_mapping_gap.md`.

## Build result

All 7 organs succeeded, clean stderr, raw-counts confirmed integer-valued
(fraction of sampled nonzero values integer = 1.0000) for every organ —
true raw counts, not normalized:

| Organ | Canonical genes | Samples | Sample cell counts |
|---|---|---|---|
| Adrenal | 28,371 | 4 | 612 / 1,381 / 112 / 778 |
| Thyroid | 30,153 | 3 | 2,109 / 5,293 / 1,897 |
| Spleen | 34,346 | 3 | 17,874 / 12,366 / 2,937 |
| Thymus | 34,809 | 3 | 24,826 / 15,839 / 1,037 |
| Liver | 33,934 | 7 | 7,702 / 12,271 / 6,963 / 6,279 / 16,280 / 14,781 / 7,531 |
| Skin | 34,978 | 3 | 22,386 / 13,860 / 21,555 |
| StomachEsophagus | 34,418 | 7 | 4,130 / 16,310 / 17,481 / 8,950 / 18,106 / 16,464 / 2,222 |

Sample counts per organ **match exactly** the replicate-structure numbers
already established (`results/04_dfp_signature/replicate_structure_audit.md`
context / Worklog's progress-tracker table: Adrenal 4, Thyroid/Spleen/
Thymus/Skin 3 each, Liver/StomachEsophagus 7 each).

## Independent verification

- All 14 output files (7 organs × counts + meta) pulled back from Argos
  and verified **byte-exact (md5)**.
- Gene/sample row and column counts independently recomputed from the raw
  TSVs (not trusted from the R script's own printed summary) — match
  exactly.
- Column order in each `*_pseudobulk_counts.tsv` verified to exactly match
  the corresponding `*_pseudobulk_meta.tsv`'s `sample` column order.

## Biological sanity check (organ-identity markers, independent of the build's own claims)

CPM-normalized within each organ, checked 2-3 canonical identity markers
per organ against their within-organ expression percentile:

| Organ | Markers checked | Result |
|---|---|---|
| Adrenal | STAR, CYP11B1, CYP17A1 (steroidogenesis) | 93.9–99.9th percentile in every sample |
| Thyroid | TG, TPO, PAX8 | TG at the **100.0th percentile in all 3 samples** (thyroglobulin is one of the most abundant transcripts in thyroid tissue — textbook result); TPO/PAX8 90.7–100.0th |
| Spleen | PTPRC, CD19, MS4A1 | PTPRC (pan-immune/CD45) ~98.5th in all samples; B-cell markers moderate-high (62–92nd), consistent with spleen's mixed immune-cell composition |
| Thymus | CD3D, CD3E, PTCRA | 67.7–93.3rd percentile, consistent with thymocyte-lineage identity |
| Liver | ALB, APOA1, TTR | ALB (albumin) at the **100.0th percentile in all 7 samples** — textbook liver identity marker |
| Skin | KRT14, KRT5, COL17A1 | 77.8–88.2nd percentile, consistent with a heterogeneous tissue where keratinocytes are one of several cell types present |
| StomachEsophagus | PGC, MUC5AC, KRT13 | see below — real signal, and it directly surfaces the known tissue-mixing issue |

**StomachEsophagus result is worth calling out specifically**: MUC5AC
(gastric mucin) is high (95.4–98.9th percentile) in all 6 Stomach-labeled
samples but markedly lower (84.8th) in the single Esophagus-labeled
sample; **KRT13** (esophageal squamous marker) shows the **opposite**
pattern — highest (99.5th percentile) in the Esophagus sample, lower
(86.3–97.6th) in the Stomach samples. This is real data evidence, not
speculation, that the previously-flagged
(`docs/STEP4_STATISTICAL_DESIGN.md` section 2) StomachEsophagus
organ-mixing concern is genuine: the one Esophagus sample has a measurably
different expression profile from the six Stomach samples. Sample-label
breakdown (from the pseudobulk build's own sample names): 6 Stomach + 1
Esophagus, not an even split — confirms this needs an explicit decision
(split into two separate organs vs. keep pooled with this asymmetry
acknowledged) before F-developmental's per-organ "elevated" computation
runs for this organ, flagged for reviewer input.

## What this doesn't do yet

No "elevated_in_fetal_somatic" criterion applied — this is the raw
pseudobulk build + verification step. Computing the actual within-organ
percentile + majority-of-samples detection criterion
(`STEP4_STATISTICAL_DESIGN.md` section 2) is the next sub-task, likely
needing its own calibration round (same discipline as P-developmental's
effect-size cutoff and `adult_excluded`'s percentile cutoff).
