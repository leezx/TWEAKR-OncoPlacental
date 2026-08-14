# adult_excluded distribution audit

Per `docs/STEP4_STATISTICAL_DESIGN.md` section 3, and the PR #10 reviewer's
explicit next step after freezing P-developmental's effect-size threshold:
"下一步进入 F-developmental 和 adult-excluded threshold calibration". This
covers `adult_excluded` — real within-tissue distributions from GTEx and
HPA, both organ-matched (feeds F-developmental) and whole-body (feeds
P-developmental), so a concrete percentile cutoff and whole-body quorum
can be chosen from real numbers rather than picked blind.

## Two real problems found and fixed before this reached review

**1. Tie-handling bug.** A first draft used `pandas.rank(pct=True,
method="average")`. Bulk RNA-seq tissues routinely have >50% of genes at
TPM=0 (verified directly: 56% for GTEx's `Adrenal_Gland` column alone) —
"average" tie-breaking dumps that entire zero block at its *mid-rank*
(~28th percentile for Adrenal_Gland), not the bottom of the distribution.
Result: 0 genes passed at any sensible low-percentile cutoff (5th/10th/
25th) for nearly every organ — caught by direct inspection of the output,
not assumed from a clean-looking summary table.

**2. Percentile is fundamentally degenerate when >50% of a distribution is
tied.** Switching to `method="min"` (zero block → percentile ≈0) fixed
problem 1, but exposed a deeper one: with a zero-tie block covering ~55%
of the distribution, *every* candidate cutoff from 5 to 50 returned the
exact same gene count. There is no gene between the zero block (percentile
≈0) and the next distinct expression level (percentile ≈55%) — no
tie-breaking convention fixes this, the percentile knob is simply
non-functional across that whole range.

**Fix**: split each tissue's genes into **not-detected** (TPM/nTPM < 1 — a
common low-expression convention in bulk RNA-seq, not attributed to any
platform's own official rule; the same floor value already used earlier
for the HPA placenta reference set in `hpa_placenta_enrichment.py`, kept
consistent) vs **detected** (≥1). Percentile rank is then computed *only*
within the detected subset, where the distribution is continuous and
percentile is a meaningful, non-degenerate knob:

```
adult_excluded(gene, tissue, P) =
    not_detected(gene, tissue) OR (detected AND percentile_among_detected(gene, tissue) <= P)
```

GTEx: 77.7% of gene×tissue pairs are not-detected (TPM<1) on average.
HPA: 34.6% (nTPM<1). Both platforms show the same zero-inflation pattern,
at different magnitudes — plausibly reflecting GTEx's per-sample-median
aggregation vs HPA's own consensus computation, not investigated further
here since it doesn't change the fix.

## Results

Full table: `adult_excluded_percentile_audit.tsv`. With the fix, all
candidate cutoffs (5/10/25/50) now give distinct, monotonically increasing
pass counts — the percentile knob works as intended.

**Organ-matched** (GTEx column(s) + HPA tissue(s) per
`hdma_organ_to_{gtex,hpa}_tissue_map.tsv`), pass counts at cutoff=10:

| Organ | GTEx pass / 73,321 | HPA pass / 20,151 | Both (AND) / 19,569 shared |
|---|---|---|---|
| Adrenal | 59,604 | 8,411 | 6,876 |
| Thyroid | 56,381 | 8,778 | 6,140 |
| Spleen | 57,098 | 8,322 | 6,237 |
| Thymus | N/A (no GTEx column) | 8,501 | N/A |
| Liver | 55,970 | 8,818 | 5,445 |
| Skin | 57,857 | 8,296 | 6,250 |
| StomachEsophagus | 56,031 | 6,830 | 4,894 |

**Whole-body** (GTEx's 68 tissues primary per the design doc; HPA's 39
non-placenta tissues reported alongside since GTEx has no Thymus column),
pass counts at cutoff=10:

| Quorum | GTEx pass / 73,321 | HPA pass / 20,151 |
|---|---|---|
| all 68/39 | 39,104 | 1,972 |
| all but 1 | 47,668 | 3,452 |
| all but 2 | 49,706 | 4,242 |

## Canonical marker sanity check

The 9 trophoblast markers from the P-developmental DE sanity check behave
exactly as expected — genuinely placenta-restricted genes should be
adult-excluded almost everywhere:

| Marker | Not-detected in GTEx tissues | Median percentile where detected |
|---|---|---|
| ERVFRD-1 | 65/68 | 1.9 |
| CGA | 67/68 | 99.9 (1 tissue, likely pituitary — biologically real cross-reactivity, not a bug) |
| CSH1 | 68/68 | n/a (never detected) |
| CSH2 | 66/68 | 10.9 |
| PSG1 | 68/68 | n/a (never detected) |
| PSG3 | 67/68 | 7.9 |
| GATA3 | 40/68 | 38.7 |
| KRT7 | 29/68 | 45.6 |
| HLA-G | 17/68 | 22.3 |

The strongly placenta-specific genes (CSH1/2, PSG1/3, ERVFRD-1) are
almost universally not-detected in GTEx, while the broader-function genes
(GATA3, KRT7, HLA-G — general epithelial/TF/immune-modulatory roles with
known adult expression elsewhere) are detected more widely — matches known
biology, not an artifact.

## Open questions for reviewer sign-off (not resolved here)

1. **Exact percentile cutoff** among 5/10/25/50 — not chosen here, same
   discipline as P-developmental's calibration (real numbers reported,
   final pick deferred to review, likely needs its own marker/enrichment
   cross-check once F-developmental's HDMA pseudobulk is available too).
2. **Whole-body quorum**: "all 68" vs "all but 1/2" — the design doc's
   default ("every one of GTEx's 68 tissues... or a high fraction") isn't
   pinned to a specific number.
3. **GTEx's Thymus gap**: whole-body currently uses GTEx alone per the
   design doc, which has no Thymus column — HPA's 39-tissue whole-body
   numbers are reported side by side as a candidate supplement/replacement
   but not adopted here.
4. **StomachEsophagus organ/tissue mixing** (flagged previously, still
   open): HDMA's 7 StomachEsophagus samples mix Stomach- and
   Esophagus-labeled tissue while GTEx/HPA treat them as separate adult
   tissues — this audit used the combined organ-level mapping as-is;
   splitting may be warranted before F-developmental's final computation
   for this organ.
