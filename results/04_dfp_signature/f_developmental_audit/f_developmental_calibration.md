# F-developmental calibration: elevated + adult_excluded combined

Per `docs/STEP4_STATISTICAL_DESIGN.md` section 2, and the three design
decisions from the PR #11 reviewer. This is the first real compute for
F-developmental's positive evidence, combining the HDMA pseudobulk (built
+ verified in PR #11) with `adult_excluded` (real distributions computed
in PR #11) into an actual candidate gene count per organ.

## The three PR #11 decisions, implemented

1. **StomachEsophagus split.** Pseudobulk split locally by sample name
   (`T*_Stomach_*` vs `T*_Esophagus_*`, confirmed unambiguous): Stomach
   n=6 (used below), Esophagus n=1 (excluded — insufficient-replication,
   not run through F-developmental). `hdma_organ_to_{gtex,hpa}_tissue_map.tsv`
   updated with separate `Stomach`/`Esophagus` rows.
2. **Thymus provenance**: no GTEx column (documented gap), so Thymus's
   `adult_ref_provenance` is explicitly `HPA-only` in the output table
   below — never silently presented as GTEx-equivalent.
3. **Whole-body GTEx/HPA independence**: not relevant to this step —
   F-developmental only uses organ-matched exclusion per the design doc;
   whole-body is P-developmental's concern and untouched here.

## Method

- HDMA pseudobulk → CPM per sample → same not-detected-floor (CPM<1) +
  percentile-among-detected design as `adult_excluded` (checked directly:
  HDMA pseudobulk has 8-10% exact zeros and ~32% below CPM=1 per sample —
  less extreme than GTEx's single-tissue TPM but real, so the same fix
  applies for consistency).
- **Elevated candidate**: median per-sample percentile-among-detected ≥
  cutoff (75 or 90 tested) AND detected in ≥ quorum fraction of the
  organ's own samples (0.5 or 1.0 tested — uses the individual-level
  replicate structure directly, per `STEP4_STATISTICAL_DESIGN.md`).
- **F-developmental candidate** = elevated AND `adult_excluded` (organ-
  matched, cutoff 10 or 25 tested, provenance-tagged per organ).

## Results

Full grid: `f_developmental_calibration.tsv`. Summary at the two elevated
cutoffs (quorum=1.0, i.e. detected in every sample):

| Organ | Provenance | elev=75, adult=10 | elev=75, adult=25 | elev=90, adult=10 | elev=90, adult=25 |
|---|---|---|---|---|---|
| Adrenal | GTEx+HPA | 259 | 466 | 75 | 153 |
| Thyroid | GTEx+HPA | 287 | 461 | 104 | 178 |
| Spleen | GTEx+HPA | 462 | 754 | 142 | 263 |
| Thymus | **HPA-only** | 671 | 1,191 | 203 | 379 |
| Liver | GTEx+HPA | 315 | 516 | 89 | 156 |
| Skin | GTEx+HPA | 476 | 786 | 176 | 310 |
| Stomach | GTEx+HPA | 320 | 545 | 127 | 209 |

**Quorum (0.5 vs 1.0) makes almost no difference** at either elevated-
percentile level — a gene that clears the 75th/90th within-organ
percentile is almost always detected across *all* the organ's samples,
not just a majority. The quorum knob isn't doing discriminating work in
this data; the elevated-percentile cutoff is what actually matters.

**Thymus (HPA-only) shows systematically higher candidate counts** than
GTEx-covered organs at the same cutoffs (671 vs 259-476 at elev=75/adult=10)
— plausibly because HPA's within-tissue percentile behaves differently
than GTEx's (different assay/aggregation, and HPA generally has fewer
near-zero genes per tissue than GTEx, per `adult_excluded_percentile_audit.md`'s
34.6% vs 77.7% not-detected rates) — flagged as a real asymmetry to keep
in mind if Thymus's numbers end up looking out of line with GTEx-covered
organs downstream, not something to correct here.

## Calibration signal: AFP in Liver (textbook oncofetal gene)

**AFP (alpha-fetoprotein)** — the single most famous oncofetal gene,
massively expressed in fetal liver and clinically used as an HCC
biomarker precisely because it's normally silent in adult liver — is a
built-in, unplanned test of whether this whole pipeline is doing the
right thing:

- **Fetal Liver CPM**: 612.7–1,484.0 across all 7 samples (median
  percentile-among-detected = 99.7, detected in 7/7 samples) — massively
  elevated, exactly as expected.
- **GTEx Liver TPM**: 1.33 (main `Liver` column, within-tissue percentile
  7.4 — passes adult_excl_pct≤10) / 0.81, 0.96, 0.26 (Hepatocyte/
  Mixed_Cell/Portal_Tract columns, all below the CPM<1 not-detected
  floor).
- **HPA liver nTPM**: 3.4 (within-tissue percentile 20.5 — **fails**
  adult_excl_pct=10, **passes** adult_excl_pct=25).

**Direct check**: AFP is **excluded** from the Liver F-developmental
candidate set at `adult_excl_pct=10`, but **included** at
`adult_excl_pct=25` — verified directly, not inferred from the summary
counts. This is the same calibration signal pattern as P-developmental's
ERVFRD-1/KRT7 case (PR #9/#10): a too-strict adult-exclusion cutoff loses
a textbook true positive. **Recommends adult_excl_pct=25 over 10** on this
evidence, consistent with the earlier P-developmental calibration's
overall lesson that the stricter end of these percentile ranges tends to
cut real biology.

## What this doesn't do yet

No cutoffs frozen — this is the calibration report, same discipline as
every prior threshold in this project. Proposal for reviewer sign-off:
**elevated_pct=90** (comfortably above background, AFP clears it at 99.7),
**quorum irrelevant at this stage** (doesn't discriminate; could be
dropped or kept at 1.0 for simplicity), **adult_excl_pct=25** (retains
AFP, matches the P-developmental calibration's lesson). Under this combo:
Adrenal 153, Thyroid 178, Spleen 263, Thymus 379 (HPA-only), Liver 156,
Skin 310, Stomach 209 F-developmental candidates.
