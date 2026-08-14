# Step 6a: Normal-Context D/F/P × Independent Fetal/Revival Signature Decomposition

## Purpose

**PR #18 review round 1 (REQUEST_CHANGES)** caught a real structural
confound in the original scoring plan: scoring HDMA in R
(`Seurat::AddModuleScore`) and Arutyunyan in Python
(`scanpy.tl.score_genes`) makes scoring-engine perfectly collinear with
dataset, which is perfectly collinear with the F-vs-P question being
tested — an apparent F/P difference could be a scoring-implementation
artifact, not biology. Fixed by eliminating cross-language randomness
structurally (freeze all stochastic gene-selection decisions once in
Python, both languages become pure deterministic arithmetic on the same
frozen assignments) plus an implementation-parity test — see "Scoring
implementation" below. Also adopted the reviewer's non-blocking wording
fix: the primary question no longer presupposes F-aligns/P-stays-distinct,
reworded to a genuinely open pre-registered question.

User redirected Step 6's execution order: before running the cancer-context
revCSC × CRC decomposition (Step 6b, `docs/STEP6_CRC_PROJECTION_DESIGN.md`,
design approved in PR #17), first check whether the frozen D/F/P signature
agrees with an **independent, normal-tissue-derived** fetal/revival/
regeneration reference — entirely within normal developmental data, no
cancer data touched yet.

User flagged the reference: `mike_verzi_fetal_signature.gmt`
(`/home/zz950/projects/ApcKO_multiomics/RA_output/`) — 5 published mouse
gene sets from normal-tissue intestinal-injury/regeneration and fetal
biology studies (not the cancer-context revCSC). **Scope boundary, per
direct user clarification**: use only this specific `.gmt` file (a static
gene-set list) from `ApcKO_multiomics` — do **not** pull any of that
project's own real single-cell/multiome data (confirmed to exist,
e.g. `D21_fetal_reprogram/scRNA`) without separate explicit authorization.
Instead, score these 5 signatures on **D/F/P's own construction data**
(HDMA + Arutyunyan placenta, both already in this project) — per direct
user instruction, avoiding any further reach into a sibling project's data
lake (the same category of concern flagged for the original M11 discovery).

## Real inventory: the 5 reference gene sets

Read directly from Argos (`wc -l` + `awk` field count, not assumed):

| Gene set | n genes (mouse) | Source |
|---|---|---|
| `REGENERATIVE_EPITHELIUM` | 203 | Wang et al., Cell 2019; Yui et al., Cell Stem Cell 2018 |
| `FETAL_INTESTINE_GENES` | 1,258 | Yui et al., Cell Stem Cell 2018 |
| `REVIVAL_STEM_CELL_GENES` | 236 | Ayyaz et al., Nature 2019 |
| `FETAL_SPHEROID_EPITHELIUM_GENES` | 270 | Mustata et al., Cell Reports 2013 |
| `YAP_SIGNALING_GENES` | 444 | Gregorieff et al., Nature 2015 |

**Note on a prior error, for the record**: earlier this session, revCSC
(the *cancer*-context CRC anchor in Step 6b) was incorrectly attributed to
Ayyaz et al. 2019 — that paper is real and relevant, but as the source of
`REVIVAL_STEM_CELL_GENES` *here*, a normal-tissue intestinal-injury
signature, not the cancer-context revCSC used in Step 6b (which the user
confirmed is from *An oncogenic phenoscape of colonic stem cell
polarization*). These are two distinct "revival stem cell" concepts from
different papers/contexts — flagging explicitly to avoid re-conflating them.

All 5 gene sets use mouse gene symbols (mixed-case convention, e.g.
`Pdgfra`, `Col1A1`) — same mouse→human orthology problem as revCSC, at
much larger scale (~2,400 raw gene entries across 5 sets, vs. revCSC's 32).

## Target data: D/F/P's own construction datasets (already in this project, no new data)

- **HDMA** (`Human Development Multiomic Atlas`), 7 organs already
  downloaded: Adrenal, Thyroid, Spleen, Thymus, Liver, Skin,
  StomachEsophagus. Real per-cell Seurat RNA objects:
  `/home/zz950/DATA/scRNAseq/HumanDevelopmentMultiomicAtlas/raw/
  per_organ_RNA_seurat/<Organ>_RNA_obj_clustered_final.rds`. This is the
  exact data F-developmental was built from (Step 4) — but note:
  F-developmental itself was built at **per-sample pseudobulk** resolution
  (`build_hdma_pseudobulk.R`), not per-cell. This step will score at
  **per-cell** resolution (the RDS objects have real per-cell counts) for
  a genuine continuous decomposition, then aggregate to per-sample summary
  stats for the correlation claim — matching the donor/sample-aware
  discipline used in Steps 4/5.
- **Arutyunyan placenta** (`Arutyunyan2023_MFI`), the same primary_tissue
  dataset P-developmental's trophoblast DE was built from:
  `/home/zz950/DATA/scRNAseq/Arutyunyan2023_MFI/raw/primary_tissue/
  adata_all_donors_all_cell_states_UPD_20230307.h5ad` — real per-cell
  h5ad, raw counts confirmed in Step 4.

**Same non-independence caveat as any signature scored on its own
construction data**: the mike_verzi gene sets themselves have zero
construction-time dependency on D/F/P (different species of origin,
different papers, different disease/model context — same category of
independence check already applied to M11/revCSC in Step 6b) — but scoring
them on the *exact cells* F-developmental was statistically carved from
means this is a **consistency check within the construction data**, not an
independent-cohort replication. That's an intentional, appropriate scope
for this step (per user's request) — a genuinely new external validation
(comparable to Step 5's Tabula Sapiens role for D/F/P itself) would need a
third-party normal-tissue dataset neither signature touched, which is not
what's being asked here.

## Plan

### 1. Mouse→human ortholog mapping (bulk, not one-by-one)

Revcsc's round-2 review this session (PR #17) already established the
right discipline — naive symbol-uppercasing is not verified orthology, and
even Ensembl's own single-gene `homology/symbol` REST endpoint was
unreliable/rate-limited for just 32 genes (many `HTTP 500`s, one gene
needed 2 minutes across retries). At ~2,400 raw gene entries, the same
one-by-one approach is impractical. Plan: use a **bulk** orthology
resource instead — Ensembl BioMart's bulk export (single query returning
mouse Ensembl ID ↔ human Ensembl ID ↔ orthology-type for all mouse genes
at once, or filtered to this gene list) rather than per-gene REST calls.
Same verification standard as revCSC: freeze a
`mike_verzi_human_FINAL.tsv` (mouse symbol, mouse Ensembl ID, human
ortholog symbol, human Ensembl ID, orthology type, inclusion decision),
primary = Compara-confirmed `one2one` only, `one2many`/ambiguous calls
reported as an extended/sensitivity variant only — same pre-registered
rule PR #17 round 3 established, applied from the start this time instead
of being caught by review.

### 2. Gene-overlap audit vs. D/F/P

Same discipline as the M11 and revCSC audits: intersect each of the 5
human-mapped gene sets against D-shared, F-specific (global + 7 lineage
modules), and P-specific *before* any scoring. Given `FETAL_INTESTINE_GENES`
alone is 1,258 genes, expect non-trivial overlap with F-specific
(2,504 genes) — report it plainly and apply the same overlap-exclusion
contract (drop shared genes from the mike_verzi score per pairwise
comparison) rather than assuming it will be small like revCSC's was.

### 3. Scoring implementation — one unified, language-independent contract (revised, PR #18 review round 1)

**Original plan (flawed, caught by review)**: score HDMA in R via
`Seurat::AddModuleScore` and Arutyunyan in Python via
`scanpy.tl.score_genes`, both funneled into "the same" null-calibrated
percentile. **Reviewer correctly identified a structural confound**: the
two functions' control-gene sampling, expression binning, and default
parameters are not guaranteed equivalent — and because HDMA is scored in
R and Arutyunyan in Python, *scoring engine is perfectly collinear with
dataset*, which is perfectly collinear with the F-vs-P biological
question. Any apparent F-vs-P difference could be a scoring-implementation
artifact, not biology. "Same signatures, same method, two datasets" was
not actually true — same *signatures*, different *methods*.

**Fix: eliminate cross-language randomness structurally, not by trying to
replicate one ecosystem's RNG bit-for-bit in the other.** All stochastic
decisions (detectability binning, control-gene sampling, the 500 null
gene panels) are computed **once, in Python**, per dataset, and frozen to
plain gene-list files — reusing the exact detectability-decile-matched
permutation-null code already built and reviewer-approved in Step 5. R and
Python then both become pure, deterministic arithmetic engines that
consume the *same precomputed assignments* — no independent sampling in
either language, so there is no RNG-parity problem to solve:

1. **Python (once per dataset — HDMA org-by-org, and Arutyunyan)**:
   compute each gene's detectability decile from that dataset's own
   expression distribution; for every scored gene (all mike_verzi/D/F/P
   genes), sample its control-gene set from the matching detectability
   decile; draw the 500 null gene panels. Freeze all of this as explicit
   gene-ID lists (`control_genes_<dataset>_<target_gene>.txt`-equivalent
   TSV, `null_panel_<dataset>_permXXX.txt`) — no randomness left for
   either language to redo.
2. **Observed/null module score, identically defined in both languages**:
   `score = mean(log1p-normalized expression, target genes) -
   mean(log1p-normalized expression, that gene's precomputed control set)`
   per cell — the R script and the Python script each just read the
   frozen gene-list files and compute this mean-difference; there is no
   remaining implementation choice for either side to diverge on.
3. **Empirical percentile**: each cell's observed score ranked against
   its own signature's 500 frozen null-panel scores (same `(n_ge+1)/
   (n_perm+1)` convention as Step 5) — computed identically in both
   languages from the same frozen null-panel gene lists.
4. **Implementation-parity test, run once before trusting any real
   result**: push one small shared toy expression matrix (same values,
   same gene IDs) through both the R and the Python arithmetic with the
   same frozen control/null gene-list assignment; confirm outputs match
   within a pre-stated numerical tolerance (`1e-9`, floating-point-only
   slack). This is the concrete, checkable artifact the reviewer asked
   for, not just an assertion of equivalence.
5. Both tracks score: the 5 mike_verzi signatures (overlap-excluded per
   comparison), D-shared, F-specific (global + 7 lineage modules where
   applicable), P-specific — genuinely the *same* method now, not just
   the same signatures.

### 4. Continuous decomposition + patient/sample validation

Per-cell null-calibrated percentile scores correlated (mike_verzi
signature × D/F/P signature), continuously, within each dataset (HDMA
per-organ, Arutyunyan placenta) — not binarized up front. Aggregated to
per-sample/per-donor summaries to confirm no single sample drives the
correlation (same discipline as Steps 4/5/6b).

**Primary question, reframed per review (non-blocking suggestion, adopted
to avoid confirmation framing)**: the original framing ("does F align in
HDMA while P stays distinct in placenta") presupposes the answer these
construction datasets would be primed to give. Reworded to a genuinely
open, pre-registered question: ***do independently-defined
fetal/revival/regenerative programs (the 5 mike_verzi signatures)
preferentially associate with F-developmental, P-developmental,
shared-D, or multiple developmental axes at once?*** If, say,
`REVIVAL_STEM_CELL_GENES` turns out to strongly associate with *both* F
and P, that is a real, reportable finding (a genuinely shared
regenerative axis) — not a "failure" of the expected F/P separation.

## Open items before compute

1. Confirm BioMart bulk-query mechanics work reliably for ~2,400 mouse
   gene IDs in one or a few batched calls (test before committing to it as
   the primary method; fall back to a smaller batched REST-based approach
   only if BioMart bulk export proves unworkable).
2. **Superseded by the unified-contract fix above (PR #18 round 1)**:
   no longer building separate `Seurat::AddModuleScore`/`scanpy.tl.score_genes`
   wrappers. Instead: (a) write the Python detectability-binning +
   control-gene-sampling + null-panel-freezing code (extends Step 5's
   existing permutation-null code to emit explicit gene-list files rather
   than only summary statistics), (b) write a minimal R script that only
   reads those frozen files and computes the mean-difference arithmetic —
   deliberately not using `AddModuleScore`'s own internal sampling at all,
   (c) run the implementation-parity test on a shared toy matrix before
   trusting any real HDMA/Arutyunyan result.
3. Per-organ HDMA cell-count/QC check (how many cells per organ actually
   pass whatever QC the pre-clustered RDS objects already encode) before
   assuming every organ supports a stable score.

## What this step is not

Not a re-opening of Step 4's frozen D/F/P cutoffs. Not a replacement for
Step 6b (revCSC × CRC) — this runs first, Step 6b (already-approved design,
PR #17) runs after. Not a pull of any `ApcKO_multiomics` single-cell/
multiome data beyond the one `.gmt` gene-set file explicitly flagged by the
user.

Submitting this design for review before running any real ortholog-mapping
or scoring compute.
