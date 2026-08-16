# Project summary: what was done, how, why, and what's been answered

Living document — updated after every reviewed step. For the full
blow-by-blow (every review round, every bug caught and fixed), see
`Worklog.md`. This doc is the structured "what/how/why/answered" entry
point; `Worklog.md` is the chronological record.

## The question this project answers

Is CRC's "Oncofetal" malignant epithelial cell state actually one
program, or a mixture of two distinct normal-development programs —
a **fetal-somatic** (embryonic organ development) program and a
separate **placental/trophoblast** program — plus a possible third,
genuinely **shared** program between the two? Before that question can
be asked of cancer data, a clean, adult-corrected reference for all
three has to exist from normal tissue alone. That reference is what
Steps 1–4 build. Applying it to real CRC Oncofetal cells (Step 6) is
what actually answers the question — the primary analysis found a weak,
largely donor/study-heterogeneous association (not a clean single-axis
"= F" or "= P" result); the secondary analysis (revCSC-high cells'
developmental composition + M11 concordance) is closed, real, and
honest: composition is bimodal (not F-dominant), and M11 — an
independent NMF-derived module — shows real, robust concordance with
revCSC substantially stronger than any single developmental axis; and
the tertiary analysis (the same composition question asked of the
*entire* 665,473-cell atlas, unconditional on revCSC status) is now also
closed and confirms this pattern is not an artifact of subsetting to
revCSC-high cells — it holds population-wide, with the SI-over-Colon
regional skew getting *stronger*, not weaker, at full scale. A
cross-reference against the project's original external 6-question
framework (Q1/Q2/Q3/Q4) confirms Q1/Q2/Q4 are already substantively
answered by this work and Q3 by the tertiary analysis's F×P quadrant
table — see "Step 6 Q1/Q2/Q4 cross-reference" below. See "What's next"
for the one remaining item that closes out the project: extending this
same pipeline to the 2 additional CRC datasets.

## Step 1 — Data inventory

**What**: catalogued every dataset's actual structure (assay names, gene
ID convention, cell/donor counts, annotation columns) by opening the
real files, not by trusting prior notes.
**Why**: earlier documentation about several datasets (e.g. Nature2026's
`snRNA_raw_counts.h5ad` annotation status) turned out to be stale or
wrong once checked directly — repeated pattern this whole project,
worth stating as a standing lesson.
**Result**: `results/01_inventory/*.json` + `SUMMARY.md`.

## Step 2 — Gene-ID mapping

**What**: built a canonical gene-symbol mapping (`canonical_feature_map`)
for all 7 HDMA organs, resolving ENSG IDs to HGNC symbols where possible.
**Why**: HDMA's gene identifiers mix native symbols and raw ENSG IDs
inconsistently across organs; every downstream comparison needs one
consistent gene key.
**Found later (Step 4) and fixed**: the original mapping was built from
each Seurat object's *default* assay (a QC-reduced gene set), not the
`RNA` assay's raw-counts feature space needed for pseudobulk — extended
with `canonical_feature_map_rna_assay/`, verified as a clean superset
(217,606 shared mappings, 0 mismatches).
**Result**: `results/02_gene_id_mapping/`.

## Step 3 — Cross-platform method contract + adult reference

**What**: locked the rule that GTEx (bulk TPM), HPA (bulk nTPM), and the
placental/HDMA scRNA-seq datasets can never be compared by raw magnitude
across platforms — only used as within-dataset rank/percentile
"adult-exclusion reference." Downloaded and processed GTEx v11 (68
tissues) and HPA RNA tissue consensus (40 tissues) as the Tier-1 adult
reference; deferred Tabula Sapiens to Tier-2 (post-freeze validation).
**Why**: without an adult baseline there's no way to tell "developmentally
re-expressed" genes from "just normally expressed in adults too" — and
without the platform-comparison rule, a spurious fold-change between
different assay types could masquerade as biology.
**Result**: `docs/STEP3_METHOD_CONTRACT.md`, `datasets/GTEx_v11_median_tpm/`,
`datasets/HPA_RNA_tissue_consensus/`.

## Step 4 — D-shared / F-specific / P-specific signature construction (CLOSED)

This is the project's core deliverable so far. Two independent
"developmental programs" are each built, adult-corrected, then combined
by simple set logic (`docs/STEP4_DFP_DESIGN.md`):

- **F-developmental(gene, organ)** = elevated in that HDMA fetal organ's
  own expression distribution AND excluded from that organ's matched
  adult GTEx/HPA tissue(s).
- **P-developmental(gene)** = significantly enriched in trophoblast vs.
  non-trophoblast within placental scRNA-seq datasets (replicated across
  ≥2 independent datasets) AND excluded from adult expression body-wide
  (GTEx + HPA).
- **D-shared** = F-developmental AND P-developmental. **F-specific** =
  F-developmental AND NOT P-developmental. **P-specific** = P-developmental
  AND NOT F-developmental.

### How each half was actually calibrated (not assumed — computed, checked against real markers, corrected when review caught problems)

**P-developmental**: paired pseudobulk edgeR QLF DE (`~donor+status`) on
2 datasets with confirmed raw counts (Arutyunyan, Nature2026 — a third,
VentoTormo, was found to have no raw counts anywhere and was demoted to
secondary support). Frozen at `logFC≥0.75 & FDR<0.05` in both datasets
(2-of-2) — the cutoff was originally proposed higher, then corrected
down after the review caught that `|logFC|` was picking up
trophoblast-*depleted* genes by mistake, and that the first "±1" cutoff
choice would have failed two textbook trophoblast markers (ERVFRD-1,
KRT7). Whole-body adult exclusion: computed real GTEx/HPA percentile
distributions, caught and fixed a tie-handling bug where >50% zero-TPM
values in bulk tissue made naive percentile ranking degenerate, then
found that GTEx and HPA disagreed on ~50% of candidate genes — traced
this down to (a) most of the disagreement being a coverage-gap artifact
(genes only one platform measures at all, not real disagreement) and
(b) the true 11% residual disagreement being explainable (testis-driven
noise on one side, genuinely broadly-expressed genes GTEx's larger panel
catches on the other). Final: 84-gene high-confidence primary tier
(both platforms agree), 174-gene lower-confidence extended tier kept
separate.

**F-developmental**: per-organ, per-sample pseudobulk from HDMA (7
organs), with a real gene-mapping gap found and fixed along the way (see
Step 2 above). Frozen at `elevated_pct=75` (within-organ percentile) and
`adult_excl_pct=25` (organ-matched) — the cutoff was corrected from an
original higher proposal after a predefined cross-organ marker panel
(DLK1/IGF2/H19/LIN28B/PEG10 + liver-specific AFP/GPC3) showed the
stricter cutoff lost 3 of 4 real marker hits. AFP — the textbook
oncofetal gene — was used as a decisive unplanned validation: correctly
excluded from Liver's candidate set at the stricter adult-exclusion
cutoff, correctly included at the looser one that was ultimately frozen.
StomachEsophagus was split into separate Stomach/Esophagus organs after
real marker data (MUC5AC vs. KRT13) confirmed they're genuinely different
tissue identities, not just a labeling quirk.

### What questions Step 4 has answered so far

- **Is the pipeline capturing real biology, not just running without
  errors?** Yes — repeatedly validated against markers with independently
  known biology (AFP in fetal vs. adult liver; ERVFRD-1/CGA/CSH1/CSH2/
  PSG1/PSG3 in trophoblast vs. adult tissues; organ-identity markers
  like TG in thyroid, ALB in liver hitting literally the 100th
  percentile) — never assumed, always checked directly against the
  actual output.
- **Do the frozen gene sets partition cleanly, with real placental genes
  staying placental and not leaking into the "shared" or "fetal-organ"
  buckets?** Yes, directly verified: all 6 core placental hormone/fusion
  markers land as P-developmental=True, F-developmental=False, never in
  D-shared or F-specific.
- **What does the final signature look like?** D-shared: 6 genes
  (`IL1RAPL2`, `KCNH5`, `PCDH11X`, `TMC1`, `TRIM71`, `ZNF730`).
  F-specific: 2,504 genes. P-specific: 78 genes (includes `CGA`, `CGB3`,
  `CGB5`, `CSH1`, `CSH2`, `CSHL1`, `DLX4` — textbook placental hormone/TF
  genes).

**Result**: `results/04_dfp_signature/dfp_gene_sets/` —
`D_shared_FINAL.txt`, `F_specific_FINAL.txt`, `P_specific_FINAL.txt`,
`P_developmental_primary84.txt` / `_extended174.txt`,
`F_developmental_union.txt` / `_consensus2plus.txt`, per-organ
`F_developmental_<Organ>.txt`. Single reproducible pipeline:
`scripts/04_dfp_signature/build_dfp_gene_sets.py` (hard cardinality
assertions on every output — rerunning it cannot silently drift to a
stale definition). Full design/audit trail: `docs/STEP4_DFP_DESIGN.md`,
`docs/STEP4_STATISTICAL_DESIGN.md`, and the `*_audit.md` files under
`results/04_dfp_signature/`.

**Status update (2026-08-15)**: everything below this point (Tier-2
validation, M11/revCSC CRC decomposition against this pan-organ D/F/P) is
now real and CLOSED/merged — see `Worklog.md` for the full detail. More
importantly, the pan-organ D/F/P built here is **no longer the primary
CRC coordinate system** — see "Step 4a — Gut re-anchoring" below.

## Step 4a — Re-anchoring F-developmental on real fetal gut data (CLOSED; NEW primary CRC coordinate system)

The pan-organ HDMA F-lineage above (Adrenal/Liver/Skin/Spleen/
StomachEsophagus/Thymus/Thyroid) has **zero gut/colon organ
representation** — a real scientific mismatch for a CRC project, surfaced
when a hypergeometric enrichment check of 5 independent published
fetal/revival mouse signatures (`mike_verzi_fetal_signature.gmt`) against
this F-lineage came back honestly null (1/50 tests FDR<0.05). User
directed a full re-anchoring rather than a small patch: `F_Colon-
developmental` (large intestine, primary) and `F_SI-developmental` (small
intestine, secondary parallel), built from the **Gut Cell Atlas**
(Elmentaite et al., *Nature* 2021), using the same-atlas's own combined
fetal+pediatric+adult epithelium object for a real fetal-vs-adult
pseudobulk DE (edgeR, `~10X + source_family + Age_group`, Second-
trimester-fetal vs. Adult) — not an HDMA-style percentile port against an
external adult reference.

Went through 3 PRs (#19 close-out of the enrichment finding, #20 design,
#21 real compute) and ~10 real review rounds total, catching genuine bugs
each time — not wording nitpicks: a case-sensitivity bug that had wrongly
excluded 382 mouse genes from an ortholog mapping; a pseudobulk-donor
design audit that found and fixed a real chemistry/cohort confound; a
mathematically false threshold-calibration claim; a real gene-ID mix-up
(`IGF2` vs. `IGF2-1`, BioMart-verified) that had hidden the single
strongest calibration marker; a wrong hypergeometric background inflating
an enrichment claim from a real 5.68× to a false ~11×; and an unreliable
regex-based duplicate-symbol heuristic replaced with a real, authoritative
BioMart lookup. Every one of these was caught by the reviewer or by
directly re-checking the data, not assumed away.

**Frozen result** (`results/04a_dfp_gut/dfp_gut_gene_sets/`):
`F_Colon-developmental`=1,456 genes, `F_SI-developmental`=1,456 genes,
`D_Colon-shared`=5 (`KISS1`/`LGALS14`/`LGSN`/`TRIM71`/`ZNF114`),
`D_SI-shared`=4 (`CCDC169`/`DLX4`/`TRIM71`/`ZNF257`), `F_Gut-core`=712
genes (5.68× enrichment over the correctly-computed chance overlap).
Calibrated against real oncofetal markers — `AFP` (logFC=+9.85/+9.89,
FDR<3e-6 both regions) is the single strongest signal, matching its
textbook role. Full design/audit trail: `docs/STEP4A_GUT_FDEV_DESIGN.md`,
`docs/STEP4A_GUT_COMPUTE_RESULTS.md`.

The original pan-organ HDMA D/F/P above is **not deleted** — demoted to a
secondary "pan-fetal / non-gut developmental validation framework."

### External validation — three-layer statistical test vs. 5 independent `mike_verzi` signatures (CLOSED)

Tested `F_Colon-developmental`/`F_SI-developmental` against the 5
independent `mike_verzi` mouse fetal/regenerative intestinal signatures
(never used in construction) with three layers sharing one locked
universe (`one2one-ortholog-eligible ∩ filterByExpr-tested`, Colon
N=11,646, SI N=11,661): hypergeometric enrichment (fold/OR/95% CI/
BH-FDR), preranked GSEA (`fgsea`, ranked by the real edgeR `logFC`,
primary evidence), and a 10,000-draw size-and-expression-matched
permutation null (BH-FDR across the same 10 tests as layer 1). Went
through 2 real review rounds (PR #22): round 1 caught GSEA's ranked
universe not matching layers 1/3, the permutation layer missing
multiple-testing correction (which flips 3 borderline results from
nominally-significant to not, once corrected), and a factual
contradiction between the results table and the interpretation text;
round 2 caught a remaining overclaim in the interpretation's own
section header. **APPROVED and merged** (`f1a59a0`).

**Honest, FDR-consistent result** (`docs/STEP4A_GUT_MIKE_VERZI_VALIDATION.md`):
`FETAL_INTESTINE_GENES` fully triangulates (all 3 layers significant,
both regions) — the strongest validation of the re-anchoring.
`REGENERATIVE_EPITHELIUM` fully triangulates in Colon, partial (2-of-3
layers) in SI. `YAP_SIGNALING_GENES` and `FETAL_SPHEROID_EPITHELIUM_GENES`
both show significant *adult*-skewed GSEA (the primary evidence layer)
in both regions, with no FDR-robust positive signal on the overlap-based
layers — a real discordance, reported as-is. `REVIVAL_STEM_CELL_GENES`
shows no signal on any layer, either region.

### External adult-expression audit — CLOSED (GTEx + Tabula Sapiens, PR #23/#24)

Checked the frozen gut D/F/P against adult references beyond the Gut Cell
Atlas's own internal contrast. Went through 5 review rounds total (3 design,
2 compute) — see `Worklog.md`'s PR #23/#24 entry for the full list of real
bugs caught (estimand mismatches, a false "reused contract" claim, a gene-ID
mapping contract, a single-donor degenerate-permutation-test bug). Result
(`docs/STEP4A_GUT_ADULT_VALIDATION_RESULTS.md`): `F_Colon-developmental`
shows no statistically-supported adult-expression anomaly in the 2 epithelial
cell types with real Tabula Sapiens donor replication; GTEx bulk shows the
expected mid-distribution detection pattern with a real ~6.5–10% minority
worth future scrutiny; D/P sets remain internally consistent and reproduce
Step 5's known-good result exactly.

### Step 6 gut re-anchor — Layer 2 substitution contract locked (PR #25)

The already-approved Step 6 CRC decomposition design (`docs/STEP6_CRC_PROJECTION_DESIGN.md`,
PR #16/#17) is being re-run against gut-specific D/F/P instead of
pan-organ — the method itself is unchanged, only which gene sets are
Layer 2. Went through 4 review rounds to lock: global D = `D_Gut-shared`
(8 genes), global F = `F_Gut-specific` (2,192, coarse/secondary),
regional F = `F_Colon-specific` (primary) + `F_SI-specific` (secondary),
global P = `P_Gut-specific` (76 genes). Gene-overlap audit against
revCSC done — negligible overlap (`CLU`/`ASS1` only), matching the
original pan-organ finding. See `docs/STEP6_GUT_REANCHOR_DELTA.md`.

### Step 6 primary scoring compute — CLOSED (PR #27, real result, answers Q1)

`docs/STEP6_GUT_SCORING_COMPUTE_DESIGN.md` (PR #26, 2 review rounds)
locked the implementation contract: 13 scored gene sets (8 revCSC panels
— 27-primary + 28-extended/sensitivity, each up to 4 overlap-exclusion
forms — + 5 D/F/P panels), null-calibrated empirical percentile as the
primary common-scale metric (null-calibrated z-score secondary),
`N_PERM=100` gated by a 100-vs-500 convergence check, and donor/study-aware
validation (composite `donor_key=(study_id,patient_id)`).

The real compute (PR #27, 2 review rounds) ran on all 665,473
`CRC_single_cell_atlas_2025` malignant cells. Two real compute-feasibility
problems were found and fixed mid-run (naive `scanpy.tl.score_genes`'s
~18h+ full-genome-recomputation-per-call cost, replaced with a
numerically-**exactly**-validated `score_genes_fast` + CSR→CSC
conversion, cutting the run to ~50 min), and one real reproducibility bug
(per-panel RNG seeds used Python's process-randomized `hash()`, not
actually reproducible from the nominal fixed seed — fixed with
`hashlib.sha256`-based seeding, both prior runs discarded and re-run from
scratch to confirm).

**Real, honest result** (`docs/STEP6_GUT_SCORING_COMPUTE_RESULTS.md`):
all 10 revCSC↔D/F/P pooled correlations are weak (|r| ≤ 0.19). The 4 D/P
pairs are weakest and the only ones robust to leave-one-donor/study-out;
5 of 6 F pairs are not robust (for one, a single study's exclusion shifts
the pooled r by >11× its own magnitude). **This does not find a strong
single-axis "Oncofetal = fetal-gut program" or "= placental program"
result** — the honest answer to Q1 at this primary-analysis resolution is
a weak, largely donor/study-heterogeneous association, not a clean
single-axis one.

### Step 6 secondary analysis — design + compute CLOSED (PR #28/#29)

`docs/STEP6_SECONDARY_ANALYSIS_DESIGN.md` (4 review rounds) locks the
design for the analysis explicitly deferred out of PR #27: whether
revCSC-high cells show a distinguishable developmental composition
(D/F/P axis-supported status, not a forced single-label argmax), and
whether M11 (an independent NMF-derived Oncofetal-annotated module)
concords with revCSC. Key locks: revCSC-high defined by cross-cell rank
of the null-calibrated z-score (10% primary, 5%/20% sensitivities, using
`revCSC_primary27_minus_CLU_ASS1` uniformly to stay overlap-safe against
every gut F axis); composition reported as axis-supported status
(none/D/F/P/D+F/D+P/F+P/D+F+P at a pre-specified per-axis threshold) with
argmax kept only as a disclaimed descriptive summary; M11 scored using
`M11_minus_revCSC_overlap` (M11 shares 5 of its top50 genes with revCSC
itself — the same genes that originally identified M11 as revCSC-like);
concordance tested both as a continuous correlation and as a
donor-stratified (Mantel-Haenszel) enrichment test with the same
per-donor/leave-one-out robustness discipline PR #27 established.

Two rounds required real new compute to close, not just wording fixes: a
qsub job (3621108) independently certified `N_PERM=500` is adequate for
the actual 45-gene M11 panel (500-vs-1000 z-rank Jaccard 0.94–0.98,
Pearson r=0.9997, byte-exact/md5-verified), and a later round caught that
the overlap-audit script itself was circular (hardcoded the expected
overlap genes instead of deriving them independently) — fixed and
re-verified to produce identical output for the right reason. **PR #28
merged.**

**PR #29 (the actual compute) then went through 3 more real review
rounds**, catching real implementation deviations from the just-approved
design that the design text alone hadn't guaranteed correct execution of:
a missing Step0×StepA cross-tab (needed to show precisely how much of
the unconditional argmax's apparent F share comes from cells with *no*
supported axis — turned out to be 21.4 of 71.7 percentage points, ~30%),
an MH zero-cell-donor exclusion rule implemented too loosely twice in a
row (fixed to the design's literal wording only on the second pass), a
reappeared signed/max-abs column-naming bug (same class PR #27 already
fixed once), plus several write-up errors (an arithmetic mistake, an
unsupported causal claim, a misattributed study, and a cross-population
r-comparison overstatement that — once fixed by direct same-population
recomputation — made the original claim stronger, not weaker). **PR #29
merged** (`docs/STEP6_SECONDARY_ANALYSIS_RESULTS.md`).

**Real, honest result**: revCSC-high cells' developmental composition is
bimodal, not F-dominant (41.1% show no gut-developmental axis evidence at
all, 42.6% show `F_Gut-specific` alone). **M11 shows real, robust
concordance with revCSC that is substantially stronger than any single
gut-developmental axis** — continuous r=0.318 (gene-disjoint,
same-population comparison against all 5 D/F/P panels: |r|≤0.074),
enrichment OR≈5–6 at every matched cutoff, every CI (asymptotic and
donor-cluster bootstrap) excluding 1, robust to any single donor/study
exclusion.

### Step 6 tertiary analysis — design + compute CLOSED (PR #30/#31)

`docs/STEP6_TERTIARY_ANALYSIS_DESIGN.md` (PR #30, 4 review rounds) locks
the last item in Step 6's original 3-analysis plan: apply the same
axis-supported-status composition methodology from the secondary
analysis (PR #28/#29) to **all 665,473 malignant cells**, unconditional
on revCSC status, rather than only the revCSC-high cohort — directly
answering the axis-composition component of Q3 of the original 6-question
framework. Reuses PR #27's already-verified per-cell percentiles and
PR #28/#29's already-reviewed composition machinery (imported, not
reimplemented). Rounds 2–3 were narrow wording/internal-consistency
cleanups (dominance-vs-support language, absence-of-evidence-vs-threshold
language, a self-introduced inconsistency between two sections) — no new
conceptual issues past round 1, which fixed a "formally separable states"
overclaim and required splitting one conflated revCSC cross-tab into two
correctly-labeled ones (operational cohort membership vs. matched-null
support — the same distinction PR #28 round 1 first established).

The real compute (PR #31, 2 review rounds; real qsub job 3621204, <1 min,
byte-exact) needed no code fixes — round 1 caught 2 write-up-only issues,
both independently verified against committed data before fixing: a
semantic regression back to "natural threshold"/"biologically correct"
language that PR #30's own design review had already removed, and a real
factual sign error (a P-enrichment finding was claimed consistent with a
"weak-but-positive" revCSC↔P correlation that is actually negative per
PR #27's committed data). **PR #31 merged.**

**Real, honest result** (`docs/STEP6_TERTIARY_ANALYSIS_RESULTS.md`): the
same bimodal composition pattern from the secondary analysis holds across
the *entire* atlas, not just the revCSC-high subset — confirming it is
not an artifact of that subsetting. F-supported biology is, if anything,
marginally *less* concentrated in revCSC-high cells (52.9%) than in the
general population (56.7%); P-supported biology is only modestly more
concentrated there (8.7% vs. 8.0%). The SI-over-Colon regional skew gets
*stronger* population-wide (4.7× vs. 3× within the revCSC-high cohort).
The two revCSC cross-tabs (operational top-10%-rank cohort vs.
pre-specified percentile≥90 matched-null support) give genuinely
different, both-real pictures of the same 53,415 P-supported cells (89%
"outside" the fixed-size cohort vs. a near-even 49.5%/50.5% split at the
matched-null threshold) — the concrete real-data payoff of PR #30's
design fix. **This closes Step 6's full 3-analysis plan
(primary/secondary/tertiary) on `CRC_single_cell_atlas_2025`.**

### Step 6 Q1/Q2/Q4 cross-reference — CLOSED (PR #32)

`docs/Q1_Q2_Q4_CROSS_REFERENCE.md` maps this project's already-frozen,
already-reviewed work onto Q1/Q2/Q4 of the external 6-question framework
(`2026-GPT-TWEAKR-Oncofetal.md`, located in Zhixins-KB this session, read
in full for the first time 2026-08-16 — previously only referenced by
filename). Pure synthesis, no new compute. Went through 4 real review
rounds — every round caught a genuine wording/provenance issue (D/F/P
"orthogonal" overclaim vs. their actual mutually-exclusive-partition
construction; an untested inferential jump claiming the revCSC/M11
shared structure is "non-developmental"; Q4's Signature 4 mis-framed as
a decision made in response to the framework rather than a pre-existing
gap; a compressed mike_verzi summary that read as uniform success when
the real outcome is mixed; then 2 further population-scope/robustness-
count slips), all independently re-verified against the committed source
docs before fixing.

**Answers**: Q1 (D-shared/F-specific/P-specific from normal tissue) is
already delivered by Steps 4/4a at two resolutions (pan-organ, frozen;
gut-re-anchored, frozen, the primary CRC coordinate system) — with the
important nuance that the two resolutions use different evidence types
per component (threshold-based vs. DE-based) and gut F's external
validation is genuinely mixed (2 of 5 `mike_verzi` benchmark signatures
positively triangulate, 1 null, 2 discordant). Q2 (decompose
revCSC/M11 "Oncofetal") is already answered by Step 6 primary/secondary/
tertiary — `revCSC` does not decompose onto a single D/F/P axis; `M11`
concords with `revCSC` far more strongly than any D/F/P axis does, but
what explains that remaining concordance beyond the tested single-axis
comparisons is left explicitly open, not overclaimed. Q3 is cross-referenced
to the tertiary analysis's F×P quadrant table, not re-derived. Q4
(output signatures for scRNA/spatial) gets a table: 3 of 4 primary
signatures are frozen/reviewed matches; the 4th ("Consensus Oncofetal")
and the 2 secondary early/term-placenta modules are genuine gaps relative
to the framework — reported honestly as gaps, with `revCSC`+`M11`
noted as an existing reviewed alternative for the former, not a
substitute deliverable.

## What's next (not started)

- **2 additional CRC datasets** (`HTAN_CRC_progressive_plasticity`,
  `CRLM_NMP_ATLAS`) — extend the same primary/secondary/tertiary pipeline;
  not yet scoped. **This is the last substantive item for project
  completion.**
- **Q5/Q6 of the original 6-question framework** (macrophage/TWEAK-driven
  developmental program; independent functional consequences of the
  placental component) — **explicitly out of scope for project
  completion** (user-confirmed 2026-08-16): both require data this
  project does not have (spatial macrophage/TWEAK perturbation data;
  functional/invasion assay data). Documented here as a distinct future
  aim, not a blocker.

## Standing practices (apply to every future step, not just Step 4)

- Every real compute step runs on Argos via `qsub`, never locally — the
  Mac mini is for orchestration, git, and lightweight post-hoc
  verification of already-pulled results only.
- Every qsub job's output is pulled back and verified **byte-exact
  (md5)** before being trusted or written up — never assume a remote run
  succeeded just because the job exited.
- Every design decision or finding gets its own branch → commit → PR →
  submitted to the ChatGPT reviewer → iterated on `REQUEST_CHANGES` →
  merged only on `APPROVE`.
- Every completed/reviewed step updates `Worklog.md`'s progress tracker
  (weighted %, with the delta from the last reported number) and, for
  anything structural, this file.
- Never trust a script's own success message — recompute or
  independently re-derive key numbers (gene counts, sample counts,
  marker directions) before reporting them as fact.
