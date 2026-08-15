# Step 6 re-anchor: substitute Gut-specific D/F/P for pan-organ D/F/P

Per `Worklog.md`/`docs/PROJECT_SUMMARY.md`'s "what's next": re-run the
already-designed, fully-approved Step 6 CRC decomposition
(`docs/STEP6_CRC_PROJECTION_DESIGN.md`, PR #16/#17, ~7 review rounds,
compute never run) against the new Gut-specific D/F/P instead of the
pan-organ one. **The method is unchanged and does not need re-review** —
Layer 1/2 structure, revCSC as the primary independent Oncofetal anchor
(27-gene primary/28-gene extended), null-calibrated percentile scoring,
overlap-exclusion contract, donor/study-aware aggregation, and the
primary/secondary/tertiary analysis structure all carry over exactly as
approved. **What changes is only which D/F/P gene sets are Layer 2**:

| Original (pan-organ, PR #16/#17 approved design) | Now (gut-specific, Step 4a) |
|---|---|
| `D_shared_FINAL.txt` (6 genes) | `D_Colon-shared.txt` (5) + `D_SI-shared.txt` (4) |
| `F_specific_FINAL.txt` (2,504, global, **already `F ∧ ¬P`**, coarse/secondary summary) | `F_Gut-specific.txt` (2,192, computed — see below) — the coarse/global summary role |
| 7 per-organ `F_developmental_<Organ>.txt` lineage modules (each further intersected with `F_specific_FINAL` to stay P-deduplicated) — primary, region-resolved interpretation | `F_Colon-specific.txt` (1,451) + `F_SI-specific.txt` (1,452) — primary, region-resolved interpretation |
| `P_specific_FINAL.txt` (78 genes) | `P_Colon-specific.txt` (79) + `P_SI-specific.txt` (80) |

**Round-1 correction**: the first version of this doc/script wrongly
listed `F_Colon-developmental`/`F_SI-developmental` as the Layer 2 F
input. Checked directly against `docs/STEP4_DFP_DESIGN.md`
("F-specific = F-developmental AND NOT P-developmental") and confirmed
`F_Colon-specific`/`F_SI-specific` are real, already-built, smaller
sets (1,451/1,452 vs. 1,456/1,456 — exactly `F_{region}-developmental`
minus the 5/4 `D_{region}-shared` genes) that correctly preserve the
D/F/P mutual-exclusivity the original Step 6 design explicitly required
("organ module ∩ F-specific, to stay within the frozen, P-deduplicated
signature"). Using the raw `-developmental` sets would have silently
broken that property for the gut re-anchor. Fixed in both this doc and
the overlap-audit script; re-run confirmed the numeric overlap
conclusion is unaffected (neither `CLU` nor `ASS1`, the two genes found
to overlap revCSC, is among the tiny `D_Colon-shared`/`D_SI-shared`
sets that got removed).

**Round-2 correction**: round-1's fix mapped only the lineage-module role
of the original F arm, silently dropping the coarse/global
`F_specific_FINAL` role the original design also explicitly scored
("score both the global F-specific set ... and each of the 7 per-organ F
modules separately ... with the global F score used only as a coarse,
secondary summary"). Fixed: `F_Gut-specific` is now **computed by the
script** (not hardcoded, not a pre-existing file) as `(F_Colon-developmental
∪ F_SI-developmental) \ P_developmental` — the exact gut-scoped analogue
of `F_specific_FINAL`'s own definition, just unioned across both regions
first. `F_Gut-core` (Colon∩SI, **not** P-deduplicated) remains a separate
tertiary concordance/core summary — it does not substitute for either the
global or region-resolved F role.

This is a strictly closer match to CRC's actual tissue of origin (gut
epithelium) than the pan-organ set ever was — the entire reason for the
Step 4a redirect.

## Re-run 1: revCSC × Gut-D/F/P gene-overlap audit (done, this PR)

Same discipline as the original `revcsc_dfp_overlap_audit.py` (PR
#16/#17), against the gut sets instead. Uses the **final frozen**
ortholog-audited revCSC sets (`revCSC_symbols.primary27.txt`,
`revCSC_symbols.extended28.txt`) — not re-derived from the raw CSC
table, matching PR #17 round 2's frozen labeling.

**Cross-dataset gene-ID contract**: the gut F/D/P gene lists are raw Gut
Cell Atlas `var_name`s (may include anndata duplicate-symbol artifacts,
e.g. `IGF2-1` not `IGF2` — the same class of issue this project has
already found and fixed twice, PR #21 and PR #23/#24). Overlap is
computed via each gut gene's `authoritative_symbol` (PR #21's own
`var_id_map.tsv`), not the raw `var_name` string — same contract as the
Step 4a adult-validation work.

**Result** (`results/06_crc_projection/revcsc_gut_overlap_audit/revcsc_gut_dfp_overlap_audit.tsv`):

| Gut set | n genes | n resolved | Overlap with revCSC primary (27) |
|---|---|---|---|
| `F_Colon-specific` | 1,451 | 1,407 | 1 (`CLU`) |
| `F_SI-specific` | 1,452 | 1,414 | 1 (`ASS1`) |
| `D_Colon-shared` | 5 | 5 | 0 |
| `D_SI-shared` | 4 | 4 | 0 |
| `P_Colon-specific` | 79 | 79 | 0 |
| `P_SI-specific` | 80 | 80 | 0 |
| `F_Gut-core` | 712 | 695 | 0 (reported for completeness only — not a Layer 2 input, see note above the mapping table) |
| `F_Gut-specific` (computed, coarse/global) | 2,192 | 2,127 | **2** (`CLU`, `ASS1`) |

**Overlap is negligible**, matching the original pan-organ audit's
finding (there, F-specific global shared 2/2,504 genes with revCSC;
here, each region shares exactly 1, and the global `F_Gut-specific`
union — as expected — shares exactly the union of both, `CLU` + `ASS1`).
`CLU` (clusterin, a stress-response/chaperone gene) and `ASS1`
(urea-cycle enzyme, also a known stress/metabolic-reprogramming marker)
are both plausible generic "activated" markers, not developmentally
organ-specific — same interpretation as the original audit's
`ACTA1`/`ANKRD1` finding. **Same overlap-exclusion contract as the
original design**: primary revCSC↔`F_Colon-specific` uses an
overlap-excluded revCSC score (drop `CLU`, 26 genes) for that one
comparison; revCSC↔`F_SI-specific` drops `ASS1` (26 genes);
revCSC↔`F_Gut-specific` (the coarse/global comparison) drops both (25
genes); all other revCSC↔gut-D/F/P comparisons need no exclusion (zero
overlap). Full (non-excluded) revCSC score retained as a sensitivity
check throughout, per the original design.

## What's next (not in this PR — a separate, larger compute step)

The actual null-calibrated scoring/decomposition compute (Layer 2:
projecting revCSC + gut D/F/P onto real CRC malignant cells across the 4
inventoried datasets, primary/secondary/tertiary analyses,
donor/study-aware aggregation) is real, substantial engineering work —
scoring up to 665,473 cells (`CRC_single_cell_atlas_2025`, primary
dataset) plus 3 secondary/tertiary datasets, building a fresh
cross-dataset gene-ID mapping (each CRC dataset uses its own Ensembl ID
space; gut D/F/P genes need the same `var_id_map.tsv`-based resolution
used above). This will be built and submitted as its own PR, on top of
this one, rather than folded in here — keeps this PR's real, reviewable
diff small (a gene-overlap audit) and the larger compute's PR focused
only on the scoring pipeline itself.
