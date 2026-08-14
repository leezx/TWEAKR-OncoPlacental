# Final D-shared / F-specific / P-specific gene set assembly

**Superseded in part — see `whole_body_disagreement_audit.md`.** The PR #13
reviewer correctly caught that this doc's original whole-body choice
(picking GTEx as an unjustified "primary" reference, 158 genes, with HPA
only reported alongside) needed real justification given GTEx-pass and
HPA-pass sets only overlapped by about half. The disagreement audit found
that most of that apparent overlap gap was itself a coverage-gap artifact
(only 815 of 1,007 genes are measured by both platforms at all), and of
the true disagreement population, real genuine disagreement is much
smaller (90/815 = 11%). **Final P-developmental uses the 84-gene
both-platforms-agree set as primary** (not the original 158-gene
GTEx-only set), with a 174-gene extended tier reported separately. The
"Final assembly" numbers directly below are from the superseded GTEx-only
version — see `whole_body_disagreement_audit.md`'s "Revised final D/F/P
assembly" section for the corrected numbers (D-shared=6, F-specific=2,504,
P-specific=78 using the primary-84/union combination).

Both halves of the D/F/P partition are now frozen at the per-dataset/
per-organ evidence level: P-developmental's `replicated_in_placenta`
(PR #10) and F-developmental's `elevated` + organ-matched `adult_excluded`
(PR #12). Per `docs/STEP4_DFP_DESIGN.md`, two things were still explicitly
left open before final assembly, both resolved here with real data:

1. **P-developmental's `adult_excluded(gene, whole_body)`** — part of the
   design doc's compound definition (`replicated_in_placenta AND
   adult_excluded(whole_body)`) but never frozen.
2. **F-developmental's organ-specific vs. P-developmental's global scope**
   — the design doc's reviewer note explicitly deferred union-vs-consensus
   to "once real per-organ results exist."

## 1. Whole-body `adult_excluded` calibration (new)

Computed with the same not-detected-floor + percentile-among-detected
design already used for organ-matched exclusion (PR #11/#12), across
GTEx's 68 tissues (primary reference) with HPA's 39 non-placenta tissues
reported alongside (fills GTEx's Thymus gap):

| Cutoff | Quorum | GTEx pass (whole-body) | → P-developmental via GTEx | HPA pass | → P-developmental via HPA | Both agree |
|---|---|---|---|---|---|---|
| 10 | all 68 | 39,104 | 57 | 1,972 | 50 | 25 |
| 10 | all but 1 | 47,668 | 116 | 3,452 | 91 | 55 |
| 10 | all but 2 | 49,706 | 139 | 4,242 | 123 | 73 |
| 25 | all 68 | 44,705 | 89 | 2,829 | 84 | 37 |
| 25 | all but 1 | 51,834 | **158** | 4,744 | 160 | 84 |
| 25 | all but 2 | 53,456 | 180 | 5,738 | 207 | 105 |

**Marker validation** (GTEx, pct=25, all-but-1 — the combination used
below): ERVFRD-1, CSH1, CSH2, PSG1, PSG3 fail whole-body exclusion in
**0/68** GTEx tissues (essentially universally excluded, as expected for
core placental hormone/fusion genes); CGA fails in only 1/68. By contrast
GATA3 (22/68), KRT7 (26/68), HLA-G (20/68) fail more broadly, consistent
with their known non-placenta-restricted biology (established already in
the P-developmental DE sanity check) — the whole-body criterion correctly
separates genuinely placenta-restricted genes from broader-function ones.

**Chosen for this assembly**: `pct_cut=25`, `quorum=all_but_1` (consistent
with the organ-matched `adult_excl_pct=25` already frozen for
F-developmental; "all but 1" tolerates a single outlier adult tissue
without over-pruning, same logic as the organ-matched design). **This
specific combination is proposed here, not yet separately frozen by
review** — flagged for the reviewer alongside the rest of this assembly.

## 2. F-developmental: organ-specific sets, union vs. consensus

Real per-organ F-developmental gene sets (frozen `elev75/adult_excl25/quorum1.0`):

| Organ | Genes |
|---|---|
| Adrenal | 681 |
| Thyroid | 684 |
| Spleen | 1,090 |
| Thymus | 1,191 |
| Liver | 800 |
| Skin | 1,128 |
| Stomach | 750 |

**Union** (any organ): 2,510 genes. **Consensus** (≥2 organs): 1,339
genes. Organ-count distribution of the union: 1,171 genes appear in
exactly 1 organ (organ-specific), 104 genes appear in **all 7 organs**
(broadly conserved fetal-somatic program) — a real, interpretable
structure, not a flat/uninformative distribution.

**Used union for the primary assembly below** (more permissive,
consistent with F-developmental's own per-organ evidence already being a
real statistical criterion — requiring ≥2 organs would double-penalize
organ-restricted developmental genes that are real in one organ but
biologically irrelevant to others). Consensus reported as an alternate,
stricter option.

## 3. Final D/F/P gene sets

Using `P-developmental` (1,007 `replicated_in_placenta` genes ∩ 158
passing whole-body `adult_excluded`) and `F-developmental` (union, 2,510
genes):

| Set | Union-based | Consensus-based (≥2 organs) |
|---|---|---|
| **D-shared** | 14 | 8 |
| **F-specific** | 2,496 | 1,331 |
| **P-specific** | 144 | 150 |

Full gene lists: `P_developmental.txt`, `F_developmental_union.txt`,
`F_developmental_consensus2plus.txt`, `D_shared_{union,consensus}.txt`,
`F_specific_{union,consensus}.txt`, `P_specific_{union,consensus}.txt`,
and per-organ `F_developmental_<Organ>.txt`.

## Validation: does the partition make biological sense?

**All 6 core placental hormone/fusion markers used throughout this
project's DE sanity checks land exactly where they should** — direct
verification, not inferred:

| Marker | P-developmental | F-developmental (any organ) |
|---|---|---|
| ERVFRD-1 | **True** | False |
| CGA | **True** | False |
| CSH1 | **True** | False |
| CSH2 | **True** | False |
| PSG1 | **True** | False |
| PSG3 | **True** | False |
| GATA3 | False | False |
| KRT7 | False | False |
| HLA-G | False | False |

The 6 genuinely placenta-restricted markers are P-developmental and never
leak into F-developmental — exactly the expected result, since they're
trophoblast-lineage genes with no reason to be broadly elevated across
fetal Adrenal/Thyroid/Spleen/Thymus/Liver/Skin/Stomach. GATA3/KRT7/HLA-G
(broader-function genes, already shown to have real adult expression
elsewhere) correctly fail P-developmental's whole-body exclusion and also
don't happen to qualify for F-developmental — consistent, not
contradictory.

**Gene-level spot check on D-shared** (14 genes, union basis): includes
`TRIM71` (a well-established stemness/developmental gene, LIN28 pathway
partner, biologically sensible as genuinely shared fetal/placental
biology) and `PLAC4` (a known placenta-associated developmental gene) —
plausible real D-shared candidates, not junk.

**Gene-level spot check on P-specific** (144 genes, union basis): includes
`CGA`, `CGB3`, `CGB5` (hCG subunits), `CSH1`/`CSH2`/`CSHL1` (placental
lactogens), `DLX4` (trophoblast transcription factor) — textbook placental
genes, correctly isolated from the fetal-somatic side.

## What remains open

- The `pct_cut=25, all_but_1` whole-body choice is proposed, not yet
  independently frozen by review (organ-matched 25 was frozen; whole-body
  wasn't explicitly revisited).
- Union vs. consensus for F-developmental's cross-organ scope — union used
  as primary here, consensus reported as the stricter alternative; not
  independently frozen by review either.
- This is the CRC Oncofetal signature's D/F/P **candidate gene sets** —
  not yet applied to actual CRC Oncofetal cells (Step 4's next phase,
  Worklog phase I, answering the project's original Q1).
