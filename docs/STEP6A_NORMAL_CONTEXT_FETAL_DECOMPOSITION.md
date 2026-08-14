# Step 6a: Normal-Context D/F/P × Independent Fetal/Revival Signature Decomposition

## Purpose

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

### 3. Scoring implementation (two tracks — data formats differ)

- **HDMA (7 organs, `.rds`/Seurat)**: score in R via
  `Seurat::AddModuleScore` (the RDS objects are pre-clustered Seurat
  objects; converting to h5ad first is an unnecessary detour). Build the
  same expression-matched-null-calibration wrapper in R that Step 5/6b use
  in Python (detectability-decile-matched random gene panels, N=500
  permutations, empirical percentile as primary common scale).
- **Arutyunyan placenta (`.h5ad`)**: score in Python via
  `scanpy.tl.score_genes`, same null-calibration code already built for
  Step 6b, reused as-is.
- Both tracks score: the 5 mike_verzi signatures (overlap-excluded per
  comparison), D-shared, F-specific (global + 7 lineage modules where
  applicable), P-specific — same signatures, same method, two datasets.

### 4. Continuous decomposition + patient/sample validation

Per-cell null-calibrated percentile scores correlated (mike_verzi
signature × D/F/P signature), continuously, within each dataset (HDMA
per-organ, Arutyunyan placenta) — not binarized up front. Aggregated to
per-sample/per-donor summaries to confirm no single sample drives the
correlation (same discipline as Steps 4/5/6b). Primary question: **does
D/F/P's own fetal-somatic (F) signal align with the independent
normal-tissue regenerative/fetal-like reference within HDMA, while
P-specific (placental) stays comparatively distinct within Arutyunyan** —
consistent with D/F/P's own internal F vs. P separation — **or does the
mike_verzi reference blur that separation**, which would itself be an
important finding about how "generically fetal-like" vs.
"developmentally-specific" these signatures really are.

## Open items before compute

1. Confirm BioMart bulk-query mechanics work reliably for ~2,400 mouse
   gene IDs in one or a few batched calls (test before committing to it as
   the primary method; fall back to a smaller batched REST-based approach
   only if BioMart bulk export proves unworkable).
2. R `Seurat::AddModuleScore` null-calibration wrapper does not exist yet
   in this project (Python version does, from Step 5/6b) — needs to be
   written and validated to give comparable percentile semantics to the
   Python version before cross-dataset comparison is trusted.
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
