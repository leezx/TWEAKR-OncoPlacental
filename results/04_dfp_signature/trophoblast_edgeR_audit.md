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

## What this unlocks for the still-placeholder cutoffs

- **Effect-size minimum** (§1 "TBD"): the marker panel above suggests
  |logFC|≥1 is comfortably conservative for real trophoblast biology in both
  datasets (all 9 markers clear it several-fold), while still leaving a
  sizeable but not overwhelming candidate set (6,806 / 4,848 genes at
  FDR<0.05). Proposing **|logFC|≥1 and FDR<0.05 per dataset** as the
  per-dataset pass criterion feeding the 2-of-2 quorum vote — flagged here
  for reviewer sign-off before it's locked into the design doc.
- **Quorum**: with only 2 usable datasets (see raw-counts audit), "≥2 of 3"
  from the original design collapses to a **2-of-2** requirement (both must
  pass) as the strict interim rule, consistent with the recommendation in
  `raw_counts_availability_audit.md`.

No gene list is finalized here — this is the distribution/threshold audit
step the reviewer asked for after PR #8's approval, reporting what the real
data look like before locking exact numbers into `STEP4_STATISTICAL_DESIGN.md`.
