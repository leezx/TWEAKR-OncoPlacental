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
Steps 1–4 build. Applying it to real CRC Oncofetal cells (Step 5+) is
what actually answers the question — not done yet, see "What's next."

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

## What's next (not started)

- **The actual null-calibrated scoring/decomposition compute** across the
  4 inventoried CRC datasets, projecting revCSC + the locked gut D/F/P
  contract onto real CRC malignant cells — this is what actually answers
  the project's original Q1 with the anatomically-correct reference. Real,
  substantial engineering work, not yet started.
- **Q2–Q6** of the original 6-question framework — entirely unscoped,
  comes after the above.

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
