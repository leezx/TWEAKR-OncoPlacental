# Q1/Q2/Q4 cross-reference: mapping this project onto the 6-question framework

Source: `2026-GPT-TWEAKR-Oncofetal.md#定义清楚Placenta的问题` (external,
`Zhixins-KB/3.Distill/1.Projects/TWEAKR/`), read in full for the first
time 2026-08-16, previously only referenced by filename in `Worklog.md`.
The framework poses six questions in two stages: Q1/Q2/Q3/Q4 build and
apply a normal-development reference; Q5/Q6 test a causal mechanism
(macrophage/TWEAK signaling) and functional consequences. This document
maps this project's already-frozen, already-reviewed work onto Q1, Q2,
and Q4 — the three questions that turn out to already be substantively
answered under different labels. Q3 is answered directly by the tertiary
analysis and is cross-referenced, not re-derived, here. **Q5 and Q6 are
explicitly out of scope** (see final section) — this doc does not attempt
them.

No new compute or gene-set construction happens in this document. Every
number below is quoted from an already-committed, already-reviewed result
(design + compute PRs listed inline) and re-verified against the
underlying frozen files before being written here.

## Q1 — 正常人类发育中，placental trophoblast program 和 fetal somatic program 到底共享什么，又分别拥有什么？

*("In normal human development, what does the placental trophoblast
program share with the fetal somatic program, and what does each hold
independently?" — answered as three modules: D-shared, F-specific,
P-specific.)*

This is exactly what Steps 4/4a build, frozen at two successive
resolutions. **Round-1 review correction**: the first draft called these
"three orthogonal modules" — D/F/P are a mutually exclusive gene-set
*partition* by construction (each gene is assigned to at most one set by
the pipeline's own logic), not three axes whose statistical independence
was tested or demonstrated. Restated below as "three mutually exclusive
gene-set components."

**Resolution 1 — pan-organ (Step 4, `results/04_dfp_signature/dfp_gene_sets/`,
frozen)**: built from HDMA (7 fetal organs: Adrenal, Thyroid, Spleen,
Thymus, Liver, Skin, StomachEsophagus) vs. Arutyunyan/VentoTormo/Nature2026
placental trophoblast data. **Round-1 review correction**: the two arms
do not use the same evidence type. The placental (P-developmental) arm
does use real DE — trophoblast-vs-other-placental-cell edgeR, per donor
— against HPA/GTEx adult-exclusion calibration. The fetal-somatic
(F-developmental) arm does **not** use internal fetal-vs-adult DE on the
HDMA side (HDMA has no matched adult comparator built in) — it uses a
threshold-based construction instead: within-organ expression percentile
among detected genes (calibrated at `elevated_pct=75`), combined with
organ-matched `adult_excluded` calibration against GTEx/HPA (`adult_excl_pct=25`),
per `docs/STEP4_STATISTICAL_DESIGN.md`. Frozen counts (verified against
the committed files): `D_shared_FINAL.txt` = **6 genes**,
`F_specific_FINAL.txt` = **2,504 genes**, `P_specific_FINAL.txt` = **78
genes**. Validated Tier-2 against Tabula Sapiens (PR #15): 69/84
whole-body D/P genes clean, 11/82 organ-matched F-specific cell types
survive matched-null+BH-FDR correction.

**Resolution 2 — gut-re-anchored (Step 4a, `results/04a_dfp_gut/dfp_gut_gene_sets/`
+ `results/06_crc_projection/revcsc_gut_overlap_audit/`, frozen)**: the
pan-organ reference has zero gut/colon organ representation — a real
mismatch for a CRC project, surfaced by an honest null result against 5
independent published fetal-mouse signatures (see below). Here
F-developmental genuinely is DE-based: rebuilt from real human fetal gut
data (Gut Cell Atlas, Elmentaite et al. *Nature* 2021) via same-atlas
donor-pseudobulk fetal-vs-adult epithelial edgeR DE, in two regions — a
different, and for this project's CRC application more directly
relevant, evidence type than Resolution 1's threshold-based F. Frozen
per-region counts: `F_Colon-specific` = 1,451, `F_SI-specific`
= 1,452, `D_Colon-shared` = 5, `D_SI-shared` = 4. Combined into the
global gut-coordinate D/F/P actually used for CRC scoring (Step 6's Layer
2 substitution, PR #25): `D_Gut-shared` = **8 genes** (union of regional
D sets), `F_Gut-specific` = **2,192 genes** (coarse/secondary; regional
F stays two sets — `F_Colon-specific` primary, `F_SI-specific`
secondary), `P_Gut-specific` = **76 genes** (intersection of regional P
sets). This resolution is the project's **primary CRC coordinate
system** — the one used in every Step 6 primary/secondary/tertiary
scoring result below.

Both resolutions are externally validated, not just internally
consistent, though the validation outcome for gut F-developmental is
mixed, not uniformly positive — see `docs/STEP4A_GUT_MIKE_VERZI_VALIDATION.md`
(PR #22) for the full picture: benchmarked against 5 independent
`mike_verzi` mouse fetal-intestine signatures across three FDR-corrected
evidence layers (hypergeometric, permutation, GSEA), `FETAL_INTESTINE_GENES`
fully triangulates (all three layers, both regions) and `REGENERATIVE_EPITHELIUM`
fully triangulates in Colon with weaker, partial (2-of-3-layers) support
in SI — genuine positive external validation for these two. The other
three are not positive validations: `REVIVAL_STEM_CELL_GENES` is a clean
null throughout; `YAP_SIGNALING_GENES` and `FETAL_SPHEROID_EPITHELIUM_GENES`
both show a significant *adult*-skewed signal on GSEA (the primary
evidence layer) in both regions, with no FDR-robust positive support on
the overlap-based layers — a real, reported-as-is discordance, not
smoothed over. Separately, an adult-expression audit against GTEx (68
tissues) + Tabula Sapiens (PR #23/#24) found no statistically-supported
adult-expression anomaly in the 2 epithelial cell types with real donor
replication (other coverage limitations explicitly retained, not
resolved).

**Answer to Q1**: this project delivers the three-component decomposition
the question asks for, at two resolutions (general whole-body, and
gut-specific for the CRC application) — not a single "fetal" and a
single "placenta" list, which the framework document explicitly warns
against building. The two resolutions differ in evidence type per
component (Resolution 1's F is threshold-based, its P and Resolution 2's
F are DE-based), and external validation of gut F is genuinely mixed
(2 of 5 benchmark signatures positively triangulate, 1 is null, 2 show a
discordant adult-skewed signal) — reported honestly above rather than
summarized as uniform success.

## Q2 — 目前所谓的 Oncofetal 到底是什么？

*("What actually is the so-called 'Oncofetal' state?" — decompose
existing Oncofetal/revCSC/revival signatures against the Q1 D/F/P
reference, rather than assuming they equal it.)*

This project does not decompose "15 signatures" as the framework's
illustrative example suggests — it decomposes the two Oncofetal proxies
actually in scope here: the published mouse **`revCSC`** signature
(*An oncogenic phenoscape of colonic stem cell polarization*,
ortholog-mapped to 27 human genes, `docs/STEP6_CRC_PROJECTION_DESIGN.md`)
and **M11**, an independently-discovered NMF meta-program from this
project's own CRC atlas decomposition that is *annotated* (not defined)
as revCSC-like via a Jaccard best-match. **Round-2 review correction**:
the first draft said both are scored against gut D/F/P "on the same
665,473-cell population" — inaccurate for M11. `revCSC` and D/F/P are
scored across all 665,473 `CRC_single_cell_atlas_2025` malignant cells;
`M11` is scored only on its own 297,307-cell subset (the cells with a
precomputed M11 score), and the M11↔revCSC same-population concordance
comparison below is performed within that 297,307-cell subset, not the
full atlas. Both use the same null-calibrated empirical-percentile
machinery (`score_genes_fast`, `crc_gut_scoring_core.py`) throughout.

**Primary analysis** (PR #27, `docs/STEP6_GUT_SCORING_COMPUTE_RESULTS.md`):
all 10 revCSC↔D/F/P pooled correlations are weak, |r|≤0.19. The 4 D/P
pairs are weakest and the only ones robust to leave-one-donor/study-out;
5 of 6 F pairs are not robust (one pair's entire pooled signal is driven
by a single study). **No clean single-axis "Oncofetal = F" or "= P"
result at this resolution.**

**Secondary analysis** (PR #28/#29, `docs/STEP6_SECONDARY_ANALYSIS_RESULTS.md`,
revCSC-high cohort only): developmental composition of revCSC-high cells
is **bimodal, not F-dominant** — 41.1% show no gut-developmental axis
evidence at all, 42.6% show `F_Gut-specific` alone, and only small
minorities show `P_Gut-specific` or multi-axis support. **M11 shows real,
robust concordance with revCSC substantially stronger than any single
D/F/P axis**: continuous r=0.318 in the same population where every D/F/P
panel is |r|≤0.074 against the identical revCSC score; enrichment
OR≈5–6 at every matched cutoff, every CI excluding 1, robust to any
single donor/study exclusion.

**Tertiary analysis** (PR #30/#31, `docs/STEP6_TERTIARY_ANALYSIS_RESULTS.md`,
full 665,473-cell atlas, unconditional on revCSC status): the same
bimodal pattern holds population-wide — not an artifact of restricting to
the revCSC-high subset. F-supported biology is, if anything, marginally
*less* concentrated in revCSC-high cells (52.9%) than in the general
population (56.7%); P-supported biology is only modestly more
concentrated there (8.7% vs 8.0%, though the *continuous* revCSC↔P
correlation is weak and slightly negative, r=−0.025/−0.030 — a discrete
thresholded-cohort comparison and a continuous correlation, not
contradictory but not sharing the same direction either).

**Answer to Q2**: `revCSC` does not decompose cleanly into a single
D/F/P axis — it is weakly and heterogeneously related to the normal
gut-developmental reference, with `F_Gut-specific` showing the strongest
(if still modest) association and `P_Gut-specific` showing a small
P-supported component. This is itself the kind of finding the framework
document anticipated: an existing "Oncofetal" proxy turns out not to be
well-explained by any single normal-developmental module in this
reference. **Round-1 review correction**: the first draft went further
and stated that `revCSC`/`M11`'s shared structure is therefore
"substantially non-developmental-axis structure" — that inference was
not actually tested. No multivariate D+F+P variance decomposition,
residualization, or nonlinear/multi-axis model was run; only marginal
(single-axis) correlations were computed. `M11`, discovered completely
independently (unsupervised NMF on this project's own atlas, not derived
from or fit to D/F/P), does concord with `revCSC` far more strongly
(r=0.318) than any single D/F/P axis does against the same revCSC score
in the same population (every D/F/P marginal correlation |r|≤0.074) —
but the correct, supported statement is narrower: **this concordance is
not captured by any single D/F/P axis in the tested marginal
comparisons; what explains the remaining concordance is open.** The
secondary-results document itself is explicit on this point (`docs/STEP6_SECONDARY_ANALYSIS_RESULTS.md`,
"What this shows, honestly"): it "does not establish that M11 *is* the
'Oncofetal' fetal-gut program." This document follows that same
restraint rather than the stronger, untested framing (see "What this
does not answer" below).

## Q3 — pointer only (already answered elsewhere)

*("Does Oncofetal split into separable Fetal×Placental quadrant
states?")* — directly answered by the tertiary analysis's Q3 F×P
quadrant table (`docs/STEP6_TERTIARY_ANALYSIS_RESULTS.md` §3, full atlas,
n=665,473): Fetal-high/Placenta-low 51.9%, double-low 40.1%,
Fetal-low/Placenta-high 3.2%, Fetal-high/Placenta-high 4.8%. Per that
document's explicit, review-established boundary: this is quadrant
occupancy / axis-defined composition, **not** on its own evidence of
formal biological separability (a continuous gradient thresholded at the
same bar produces an identical-looking table; a genuine separability
claim needs a pre-specified clustering/mixture model, out of scope for
Step 6). Not re-derived here — see the source document for the full
table and discussion.

## Q4 — 我们到底需要输出哪些可以直接用于 scRNA/spatial 的 signatures？

*("Which signatures should we actually output for direct use on
scRNA/spatial data?")* The framework asks for (at least) 4 primary
signatures plus 2 secondary (early- vs. term-placenta) modules. Mapping
onto what this project has actually frozen:

| Framework's ask | This project's deliverable | Status |
|---|---|---|
| Signature 1 — Developmental Shared | `D_Gut-shared` (8 genes) | **Frozen**, reviewed (PR #25) |
| Signature 2 — Fetal Somatic | `F_Gut-specific` (2,192, coarse) + `F_Colon-specific`/`F_SI-specific` (regional, primary/secondary) | **Frozen**, reviewed (PR #25), externally validated (PR #22) |
| Signature 3 — Placental/Trophoblast | `P_Gut-specific` (76 genes) | **Frozen**, reviewed (PR #25) |
| Signature 4 — Consensus Oncofetal | Not built as a literature-consensus list. This project instead has a reviewed alternative: the published `revCSC` 27-gene ortholog-mapped signature as the primary empirical Oncofetal comparator, plus the independently-discovered `M11` NMF module as a secondary concordance anchor — see Q2. | **Gap relative to framework; existing reviewed alternative in place** |
| Secondary — Early Placenta (1st trimester > term) | Not built. | **Gap** |
| Secondary — Term/Late Placenta (term > 1st trimester) | Not built. | **Gap** |

All frozen gene sets are directly usable on new scRNA/spatial data as-is
(`scanpy.tl.score_genes`/`score_genes_fast` + null calibration is the
project's own scoring recipe, already applied at atlas scale) — this
part of Q4 is satisfied. **Round-1 review correction**: the first draft
described Signature 4's gap as a "deliberate, reviewed design decision"
made "relative to" this framework's ask — that overstates what PR
#16/#17 actually decided. Those PRs chose `revCSC` as the primary
Layer-1 Oncofetal anchor (explicitly *not* asserted as "a field-consensus
'the' Oncofetal definition," per `docs/STEP6_CRC_PROJECTION_DESIGN.md`)
and `M11` as a secondary check — but that choice was made before this
external 6-question framework had ever been read in this project (it was
only located and read in full 2026-08-16, this session). It cannot have
been a decision made *in light of* Q4's specific "Consensus Oncofetal"
ask, because that ask wasn't known yet. The accurate framing: Signature 4
is a genuine gap relative to the framework, for which this project
already has a reviewed alternative analysis strategy (`revCSC` primary
comparator + `M11` secondary anchor) that happens to serve a similar
purpose, not a substitution chosen in response to the framework itself.
The two secondary modules (early- vs. term-placenta) were never in scope
for any Step 4/4a/6 design document — a genuine gap relative to the
framework's full ask, not previously flagged. Closing it would require
re-analyzing the placental reference data (Arutyunyan/VentoTormo/Nature2026)
with an explicit gestational-age split, which none of this project's
frozen P-developmental work currently does.

## What this document does not answer

- Does not build the two missing early/term-placenta secondary modules
  (flagged above as a genuine gap, not attempted here).
- Does not build a new literature-consensus Oncofetal signature (Q4's
  Signature 4 as literally specified) — a genuine gap, for which
  `revCSC`+`M11` is an existing reviewed alternative, not the same
  deliverable (see Q4 above for why this is a gap, not a substitution
  made in response to the framework).
- Does not investigate what explains the `revCSC`↔`M11` concordance
  beyond the tested single-axis D/F/P marginal comparisons (flagged
  under Q2, not tested — **round-2 review correction**: the first draft
  of this bullet said "non-developmental structure," reintroducing the
  exact untested inference round 1 had already removed from Q2's main
  text; this bullet must match that same restraint).
- Does not extend any of the above to the 2 additional CRC datasets
  (`HTAN_CRC_progressive_plasticity`, `CRLM_NMP_ATLAS`) — separate,
  already-scoped future work.
- Does not touch Q5 or Q6 (see below).

## Q5/Q6 — explicitly out of scope

**Q5** (does SPP1+ TAM/TWEAK/Fn14/YAP macrophage signaling drive a
*specific* developmental component — generic reversion, F-specific,
P-specific, or a mixture in some ratio?) and **Q6** (does the placental
component have independent functional consequences — invasion, immune
exclusion, metastasis, antigen presentation, TAM proximity, TWEAK/Fn14/YAP
dependence — versus F-specific cells' regenerative/plastic/stem-like
phenotype?) both require data this project does not have: spatial or
perturbation data linking macrophage/TWEAK signaling to developmental
state (Q5), and functional/phenotypic assay data — invasion, metastasis,
immune-exclusion readouts (Q6). **User-confirmed 2026-08-16**: both are
explicitly out of scope for project completion, documented here as a
distinct future aim rather than a blocker. Any future work on Q5/Q6
should reuse this project's frozen `D_Gut-shared`/`F_Gut-specific`/
`P_Gut-specific`/`revCSC`/`M11` scores directly as the developmental-state
axis — no new signature construction should be needed to start that
work, only new data (macrophage colocalization/perturbation, functional
assays) layered against the existing per-cell scores.

## Review history

- **Round 1 (REQUEST_CHANGES, documentation-only — no compute,
  frozen-gene-set, or numerical provenance problem found; 4 wording/
  provenance issues, all independently verified against the real
  committed docs before fixing)**: (1) Q1 called D/F/P "three orthogonal
  modules" (not demonstrated statistically independent — they are a
  mutually exclusive gene-set partition by construction) and glossed over
  the fact that pan-organ F-developmental is threshold-based, not
  DE-based like pan-organ P-developmental and gut-re-anchored F — fixed
  to "three mutually exclusive gene-set components" with the evidence
  types distinguished per resolution. (2) Q2 claimed the revCSC/M11
  shared structure is "substantially non-developmental-axis structure" —
  an inference not actually tested (no multivariate/residualization
  analysis was run) — fixed to the narrower, supported statement: not
  captured by any single D/F/P axis in the tested marginal comparisons,
  with what explains the rest left open; also softened "independent
  signal" to "P-supported component" for the small P-enrichment finding.
  (3) Q4's Signature 4 gap was framed as a "deliberate, reviewed design
  decision" made "relative to" this framework — but PR #16/#17's
  `revCSC`+`M11` choice predates this project ever reading the external
  framework, so it cannot have been a decision made in response to Q4's
  specific ask — fixed to "gap relative to framework; existing reviewed
  alternative in place." (4) The mike_verzi triangulation summary
  ("triangulates against 5 independent signatures") compressed a
  genuinely mixed result into an apparent uniform success — fixed to
  spell out that only 2 of 5 signatures positively triangulate
  (`FETAL_INTESTINE_GENES` fully, `REGENERATIVE_EPITHELIUM` fully in
  Colon/partially in SI), 1 is a clean null, and 2 show a discordant
  adult-skewed GSEA signal.

- **Round 2 (REQUEST_CHANGES, 2 small residual documentation issues,
  both independently verified against `docs/STEP6_SECONDARY_ANALYSIS_RESULTS.md`
  before fixing)**: (1) Q2's setup said `revCSC` and `M11` are "scored
  against gut D/F/P on the same 665,473-cell population" — inaccurate for
  M11, which is scored only on its own 297,307-cell subset (confirmed
  against the secondary-results doc's job-provenance table); fixed to
  state the two population scopes separately. (2) The "What this
  document does not answer" section's Q2-related bullet still said
  "non-developmental structure shared between revCSC and M11," silently
  reintroducing the exact untested inference round 1 had already removed
  from Q2's main text — fixed to match round 1's restraint ("what
  explains the concordance beyond the tested single-axis comparisons").

Submitting for round-3 review before merge, same discipline as every
prior step this session.
