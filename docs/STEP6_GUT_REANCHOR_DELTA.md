# Step 6 re-anchor: substitute Gut-specific D/F/P for pan-organ D/F/P

Per `Worklog.md`/`docs/PROJECT_SUMMARY.md`'s "what's next": re-run the
already-designed, fully-approved Step 6 CRC decomposition
(`docs/STEP6_CRC_PROJECTION_DESIGN.md`, PR #16/#17, ~7 review rounds,
compute never run) against the new Gut-specific D/F/P instead of the
pan-organ one. **The method is unchanged and does not need re-review** —
Layer 1/2 structure, revCSC as the primary independent Oncofetal anchor
(27-gene primary/1-gene extended-only sensitivity addition), null-calibrated
percentile scoring, overlap-exclusion contract, donor/study-aware
aggregation, and the primary/secondary/tertiary analysis structure all
carry over exactly as approved. **What changes is only which D/F/P gene
sets are Layer 2** — this went through 3 real review rounds to get right
(see corrections below); the final, locked mapping is:

| Original (pan-organ, PR #16/#17 approved design) | Now (gut-specific, Step 4a) |
|---|---|
| `D_shared_FINAL.txt` (6 genes, **one global axis**) | `D_Gut-shared.txt` (8 genes, computed = `D_Colon-shared ∪ D_SI-shared`) — **the global D axis** |
| `F_specific_FINAL.txt` (2,504, **one global axis**, coarse/secondary summary) | `F_Gut-specific.txt` (2,192 genes, computed = `(F_Colon-developmental ∪ F_SI-developmental) \ P_developmental`) — **the global F axis**, coarse/secondary summary |
| 7 per-organ `F_developmental_<Organ>.txt` lineage modules (each ∩ `F_specific_FINAL`) — primary, region-resolved interpretation, **F's unique per-organ feature** | `F_Colon-specific.txt` (1,451) + `F_SI-specific.txt` (1,452) — primary, region-resolved interpretation |
| `P_specific_FINAL.txt` (78 genes, **one global axis**) | `P_Gut-specific.txt` (76 genes, computed = `P_Colon-specific ∩ P_SI-specific`) — **the global P axis** |

`D_Colon-shared`/`D_SI-shared`/`P_Colon-specific`/`P_SI-specific` remain
available as secondary, regional descriptive views — **not** independent
replication panels (see round-3 correction below).

This is a strictly closer match to CRC's actual tissue of origin (gut
epithelium) than the pan-organ set ever was — the entire reason for the
Step 4a redirect.

## Corrections (3 real review rounds before this contract was right)

**Round 1**: the first version wrongly used the raw
`F_Colon-developmental`/`F_SI-developmental` (still containing the
D-overlap genes) as the Layer 2 F input. Checked directly against
`docs/STEP4_DFP_DESIGN.md` ("F-specific = F-developmental AND NOT
P-developmental") and fixed to `F_Colon-specific`/`F_SI-specific` — the
real, already-built P-deduplicated sets (1,451/1,452 vs. 1,456/1,456).

**Round 2**: the original design scored BOTH a global `F_specific_FINAL`
(coarse/secondary summary) AND the 7 per-organ lineage modules (primary) —
round 1's fix mapped only the lineage-module role, silently dropping the
global-summary role. Fixed: `F_Gut-specific` computed (not hardcoded) as
`(F_Colon-developmental ∪ F_SI-developmental) \ P_developmental`.

**Round 3 (two issues)**: (a) the approved design has exactly **one**
global D axis and **one** global P axis — only F was ever
organ/lineage-resolved. Mapping `D_shared_FINAL`/`P_specific_FINAL`
directly to two *regional* axes each (`D_Colon-shared`+`D_SI-shared`,
`P_Colon-specific`+`P_SI-specific`) silently turned one global axis into
two regional ones, which the original design never had for D/P. Fixed:
`D_Gut-shared` = `D_Colon-shared ∪ D_SI-shared` (8 unique genes — the two
sets share only `TRIM71`) and `P_Gut-specific` = `P_Colon-specific ∩
P_SI-specific` (76 genes — equivalent to `P_developmental \
(F_Colon-developmental ∪ F_SI-developmental)`, i.e. genes P-developmental
but fetal-developmental in *neither* region) are now the correct global
D/P axes, both computed by the script, not hardcoded. This also
directly matches PR #24's own finding that `P_Colon-specific`/
`P_SI-specific` already share 76 of their 79/80 genes and must not be
treated as two independent panels — the regional P sets were never
meant to stand in as two replications of a global P axis.
(b) the overlap test was gated on BioMart `authoritative_symbol`
resolution, leaving 44–65 genes per set `NOT_TESTABLE` even though the
underlying GCA `gene_id` (Ensembl) is known and unique for every gene
(zero duplicate Ensembl `gene_id`s across all 33,538 GCA variables, per
PR #21's audit) and revCSC's own frozen provenance table
(`revCSC_human_FINAL.tsv`) already carries a verified human Ensembl ID
for every primary/extended member. Fixed: overlap is now tested via GCA
`gene_id` == revCSC `human_ensembl_id` (Ensembl-ID identity), with
`authoritative_symbol` retained only for human-readable reporting — same
lesson already applied to GTEx matching in PR #23/#24. Result: **100%
Ensembl-ID resolution** for every gut set tested (0 `n_no_gene_id`
anywhere), vs. the prior symbol-gated method's 44–65 unresolved genes
per set. The overlap conclusion is numerically unchanged (as expected —
BioMart-symbol failure was never evidence the underlying gene identity
was wrong).

## Re-run: revCSC × Gut-D/F/P gene-overlap audit (done, this PR)

Same discipline as the original `revcsc_dfp_overlap_audit.py` (PR
#16/#17), against the corrected gut sets. Uses the **final frozen**
ortholog-audited revCSC provenance table (`revCSC_human_FINAL.tsv`) —
not re-derived from the raw CSC table, matching PR #17 round 2's frozen
labeling.

**Result** (`results/06_crc_projection/revcsc_gut_overlap_audit/revcsc_gut_dfp_overlap_audit.tsv`):

| Gut set | Role | n genes | n Ensembl-resolved | Overlap with revCSC primary (27) |
|---|---|---|---|---|
| `F_Gut-specific` | **global F** | 2,192 | 2,192 (100%) | **2** (`CLU`, `ASS1`) |
| `F_Colon-specific` | regional F (primary) | 1,451 | 1,451 (100%) | 1 (`CLU`) |
| `F_SI-specific` | regional F (secondary) | 1,452 | 1,452 (100%) | 1 (`ASS1`) |
| `D_Gut-shared` | **global D** | 8 | 8 (100%) | 0 |
| `D_Colon-shared` | regional D (descriptive) | 5 | 5 (100%) | 0 |
| `D_SI-shared` | regional D (descriptive) | 4 | 4 (100%) | 0 |
| `P_Gut-specific` | **global P** | 76 | 76 (100%) | 0 |
| `P_Colon-specific` | regional P (descriptive) | 79 | 79 (100%) | 0 |
| `P_SI-specific` | regional P (descriptive) | 80 | 80 (100%) | 0 |
| `F_Gut-core` | tertiary concordance/core (not a Layer 2 input) | 712 | 712 (100%) | 0 |

**Overlap is negligible**, matching the original pan-organ audit's
finding (there, F-specific global shared 2/2,504 genes with revCSC; here,
the global `F_Gut-specific` shares exactly the same 2 — the union of both
regions' 1-gene overlaps, as expected). `CLU` (clusterin, a
stress-response/chaperone gene) and `ASS1` (urea-cycle enzyme, also a
known stress/metabolic-reprogramming marker) are both plausible generic
"activated" markers, not developmentally organ-specific — same
interpretation as the original audit's `ACTA1`/`ANKRD1` finding. All D
and P axes (global and regional) show **zero** overlap.

**Overlap-exclusion contract** (per the original design, extended to the
new global axes): revCSC↔`F_Gut-specific` (global) drops both `CLU` and
`ASS1` (25-gene overlap-excluded score); revCSC↔`F_Colon-specific` drops
`CLU` only (26 genes); revCSC↔`F_SI-specific` drops `ASS1` only (26
genes); all D/P comparisons (global and regional) need no exclusion
(zero overlap). Full (non-excluded) revCSC score retained as a
sensitivity check throughout, per the original design.

## What's next (not in this PR — a separate, larger compute step)

The actual null-calibrated scoring/decomposition compute (Layer 2:
projecting revCSC + gut D/F/P onto real CRC malignant cells across the 4
inventoried datasets, primary/secondary/tertiary analyses,
donor/study-aware aggregation) is real, substantial engineering work —
scoring up to 665,473 cells (`CRC_single_cell_atlas_2025`, primary
dataset) plus 3 secondary/tertiary datasets, building a fresh
cross-dataset gene-ID mapping (each CRC dataset uses its own Ensembl ID
space; the Ensembl-ID-primary contract locked in this PR extends
directly). This will be built and submitted as its own PR, on top of
this one, rather than folded in here — keeps this PR's real, reviewable
diff small (a gene-overlap audit) and the larger compute's PR focused
only on the scoring pipeline itself.
