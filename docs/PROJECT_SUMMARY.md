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

## What's next (not started)

- **Tier-2 validation**: Tabula Sapiens (adult single-cell) independent
  check that P-specific/F-specific genes are genuinely absent across
  adult cell types — held out of signature definition specifically so
  this validation is not circular.
- **Apply the frozen D/F/P signature to real CRC Oncofetal cells**
  (single-cell/spatial) — this is the step that actually answers the
  project's original question.
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
