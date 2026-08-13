# Trophoblast-vs-rest edgeR QLF results: distribution audit

First real compute for Step 4's primary P-developmental DE model, per
`docs/STEP4_STATISTICAL_DESIGN.md` §1 (edgeR QLF, paired `~donor+status`
design, fit independently per dataset). Run on the pseudobulk matrices in
`results/04_dfp_signature/pseudobulk/` (built by
`scripts/04_dfp_signature/build_trophoblast_pseudobulk.py`, verified
byte-exact against Argos output and cross-checked against the known
replicate-structure-audit donor counts before this DE was run).

Only 2 of the 3 datasets named in the original design have real raw counts —
see `raw_counts_availability_audit.md` for the VentoTormo finding and the
resulting 2-of-2 interim quorum revision.

## Run summary

| Dataset | Paired donors used | Genes passing `filterByExpr` | Common dispersion | FDR<0.05 | FDR<0.05 & \|logFC\|≥1 |
|---|---|---|---|---|---|
| Arutyunyan | 17/18 | 16,919 / 30,800 | 0.2505 | 10,227 | 6,806 |
| Nature2026 | 23/23 | 21,958 / 36,601 | 0.1227 | 14,244 | 4,848 |

Both runs: 0 NaN in logFC/PValue/FDR; clean stderr (no warnings/errors from
edgeR); output row counts and FDR-significant counts independently
reverified locally against the qsub job's own printed summary (byte-exact
md5 match between Argos and local copies).

## Biological sanity check (independent of the design's own claims)

Canonical trophoblast markers vs. canonical non-trophoblast (immune,
endothelial) markers, checked directly against the result tables — not
assumed:

| Marker | Expected direction | Arutyunyan logFC (FDR) | Nature2026 logFC (FDR) |
|---|---|---|---|
| ERVFRD-1 (syncytin-2, fusogenic STB gene) | up in troph | +5.44 (6.8e-08) | +0.93 (9.9e-06) |
| CGA (hCG alpha) | up | +5.88 (3.8e-07) | +1.59 (9.1e-05) |
| CSH1 / CSH2 (placental lactogen) | up | +5.81 / +4.33 | +1.44 / +1.36 |
| PSG1 / PSG3 (pregnancy-specific glycoproteins) | up | +4.63 / +5.37 | +1.61 / +1.70 |
| GATA3 (trophoblast TF) | up | +4.13 (4.3e-08) | +1.40 (1.7e-05) |
| KRT7 (trophoblast cytokeratin) | up | +5.74 (2.4e-07) | +0.89 (6.1e-06) |
| HLA-G (EVT marker) | up | +6.05 (1.4e-06) | +1.48 (3.3e-06) |
| PTPRC (CD45, immune) | down in troph | -5.59 (1.8e-07) | -3.12 (9.2e-06) |
| PECAM1 (CD31, endothelial) | down | -3.99 (1.4e-09) | -2.56 (6.8e-06) |

All 9 markers land in the expected direction with FDR well below 0.05 in
both datasets. Effect sizes are systematically larger in Arutyunyan than
Nature2026 (consistent with Arutyunyan's non-trophoblast comparator likely
being a broader/less trophoblast-adjacent cell mix — worth keeping in mind
when picking a single cross-dataset effect-size cutoff below, since a fixed
logFC threshold will behave asymmetrically across the two datasets).

## Correction (PR #9 round 1 review): effect-size criterion must be directional

The first draft of this doc proposed `abs(logFC)≥1` as the pass criterion
and claimed "all 9 markers clear it several-fold." Both statements were
wrong, caught by review:

1. **`abs(logFC)` is the wrong statistic.** The edgeR contrast is
   `troph vs nontroph`; a large *negative* logFC means strongly
   depleted-in-trophoblast. PTPRC and PECAM1 above are exactly this case —
   they clear `abs(logFC)≥1` by a wide margin precisely because they go the
   *wrong* direction. An absolute-value threshold would have let
   immune/endothelial-depletion signal masquerade as placenta-developmental
   evidence. The correct criterion is **`logFC ≥ cutoff`** (positive,
   trophoblast-enriched only). Negative-logFC genes can be tracked
   separately as a `trophoblast_depleted` QC set but never enter the
   P-developmental positive program this way. Fixed in
   `docs/STEP4_STATISTICAL_DESIGN.md` §1.
2. **The "all 9 markers clear +1" claim was checked wrong.** Re-checked
   directly against the per-dataset logFC values: in Nature2026, ERVFRD-1
   (+0.93) and KRT7 (+0.89) do **not** clear a +1 cutoff — only 7/9 canonical
   markers survive at that threshold in that dataset (Arutyunyan clears all
   9/9 at every cutoff tested below, since its effect sizes run much larger
   — see the asymmetry noted above).

## Calibration curve (directional `logFC ≥ cutoff` & FDR<0.05), computed before freezing anything

| Cutoff | Arutyunyan pass | Nature2026 pass | 2-of-2 overlap | Canonical markers retained (Arutyunyan / Nature2026 / both) |
|---|---|---|---|---|
| 0.5 | 3,684 | 3,775 | 1,742 | 9/9 / 9/9 / 9/9 |
| 0.75 | 2,819 | 2,227 | 1,007 | 9/9 / 9/9 / 9/9 |
| 1.0 | 2,216 | 1,293 | 536 | 9/9 / 7/9 / 7/9 (Nature2026 loses ERVFRD-1, KRT7) |

Computed directly from the two edgeR result tables (no new DE run needed —
same discipline the reviewer confirmed: "本轮只需要修改 threshold 逻辑和 audit
的错误表述，不需要重跑 edgeR").

## What this unlocks for the still-placeholder cutoffs

- **Effect-size minimum** (§1 "TBD"): **not frozen yet.** The calibration
  curve shows 0.5 and 0.75 both retain all 9 canonical markers in both
  datasets, while 1.0 loses 2 markers in Nature2026 — consistent with the
  cross-dataset effect-size scale asymmetry already noted (Arutyunyan's
  non-trophoblast comparator likely being a broader/less trophoblast-adjacent
  mix than Nature2026's). **0.75 is the current leading candidate**
  (full marker retention, meaningfully tighter than 0.5's larger candidate
  set), but per reviewer guidance this should not be locked until an
  independent check is also run: HPA trophoblast/placenta known-gene
  enrichment against the pass sets at each cutoff (not yet done — needs the
  HPA processed mapping table from Step 3, `datasets/HPA_RNA_tissue_consensus/`,
  which hasn't been joined against these results yet). That HPA-enrichment
  check is the next concrete sub-task before this number freezes.
- **Quorum**: with only 2 usable datasets (see raw-counts audit), "≥2 of 3"
  from the original design collapses to a **2-of-2 interim discovery rule**
  (both must pass) — explicitly *not* CNS-grade replication with only 2
  independent primary datasets left. VentoTormo/Greenbaum/HPA/organoid
  follow-up validation become correspondingly more important as secondary
  support, per reviewer note.

No gene list is finalized here — this is the distribution/threshold audit
step the reviewer asked for after PR #8's approval, reporting what the real
data look like before locking exact numbers into `STEP4_STATISTICAL_DESIGN.md`.
