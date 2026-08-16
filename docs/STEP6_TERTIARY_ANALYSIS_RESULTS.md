# Step 6 tertiary analysis: full-atlas, revCSC-independent D/F/P landscape — results

Real compute, per the approved design (`docs/STEP6_TERTIARY_ANALYSIS_DESIGN.md`,
APPROVE after 4 review rounds). Applies Step 6 secondary's (PR #28/#29)
already-reviewed axis-supported-status methodology to **every one of the
665,473 malignant cells in `CRC_single_cell_atlas_2025`, unconditionally
on revCSC status** — the last item in Step 6's three-analysis plan, and
directly answers the axis-composition component of Q3 of the project's
original 6-question framework. All numbers below are real qsub output,
pulled back and verified byte-exact (md5, 5 files) before being trusted
or written up — same discipline as every prior step.

## Job provenance

| Stage | Job ID | Scope | Outcome |
|---|---|---|---|
| **Composition compute** | **3621204** | 665,473 cells, no new scoring (reuses PR #27's already-verified full-atlas D/F/P + revCSC percentiles) | **Completed cleanly, <1 min; used for this doc** |

No per-cell files produced — same "small deliverable" size norm as Step 6
secondary. All 5 output tables committed in full
(`results/06_crc_projection/tertiary_analysis_composition/`).

## Section 2: axis-supported status, full atlas (n=665,473, 529 donors, 36 studies)

| Category | Pooled % | Donor-unweighted mean % | Study-unweighted mean % |
|---|---|---|---|
| F_only | 47.8% | 44.1% | 49.0% |
| none | 38.8% | 41.4% | 37.1% |
| F+P | 4.4% | 4.1% | 4.8% |
| D+F | 4.1% | 3.2% | 3.7% |
| P_only | 3.1% | 4.8% | 3.4% |
| D_only | 1.3% | 1.9% | 1.6% |
| D+F+P | 0.4% | 0.3% | 0.3% |
| D+P | 0.1% | 0.3% | 0.2% |

**Real, honest read**: across the *entire* malignant population — not
just revCSC-high cells — the same bimodal pattern from Step 6 secondary
holds: either no axis clears its matched-null bar (38.8%) or
`F_Gut-specific` alone does (47.8%). **A genuinely new, non-circular
observation this analysis adds**: F-supported biology (F_only + D+F + F+P
+ D+F+P = 56.7% of the full atlas) is *not* concentrated in revCSC-high
cells — it is, if anything, marginally *more* common in the general
population than within Step 6 secondary's revCSC-high cohort (52.9% of
that 66,548-cell cohort). Same for `F_only` alone (47.8% full atlas vs.
42.6% within revCSC-high). This is consistent with, and independently
corroborates, the primary analysis's finding (PR #27) that revCSC↔F
correlations are weak — revCSC-high status is not meaningfully
concentrating F-axis evidence. P-supported biology (P_only + F+P + D+P +
D+F+P = 8.0% of the full atlas) shows the opposite, small effect: modestly
*enriched* within the revCSC-high cohort (8.7% of that cohort vs. 8.0%
full-atlas baseline) — a small, real signal consistent with the primary
analysis's weak-but-positive revCSC↔P correlation direction.

**Argmax distortion, quantified again on a different population (same
real pattern PR #29 first showed)**: the descriptive-only coarse argmax
assigns F to 75.4% of all cells — 18.7 points higher than the true
any-F-supported total (56.7%) — confirming this is a structural property
of the unconditional-argmax construction itself, not an artifact specific
to the revCSC-high subset PR #29 originally caught it in.

**Step B regional refinement** (within the 56.7% F-assigned cells):
`SI_biased` (18.8% of the full atlas, 33.1% of F-assigned cells) is
**~4.7× more common than** `Colon_biased` (4.0% of the full atlas, 7.1%
of F-assigned cells) — an even *more* pronounced SI-over-Colon skew than
Step 6 secondary found within the revCSC-high cohort specifically (~3×
there). This strengthens, rather than weakens, the earlier flagged
observation: across the entire CRC malignant population, not just
revCSC-high cells, F-assigned cells skew toward the *secondary* regional
axis (`F_SI-specific`) more than the *primary* one (`F_Colon-specific`,
per Step 4a's locked hierarchy) — still not investigated further here,
flagged again for future scrutiny.

## Section 3: Q3 F×P quadrant table (full atlas)

| Quadrant | Pooled % |
|---|---|
| Fetal-high / Placenta-low | 51.9% |
| double-low | 40.1% |
| Fetal-low / Placenta-high | 3.2% |
| Fetal-high / Placenta-high | 4.8% |

Derived from the identical `F_Gut-specific`/`P_Gut-specific`
axis-supported flags used in §2 (same percentile≥90 threshold, no
separate pick). **Per the design's explicit, review-established
boundary, this table shows quadrant occupancy / axis-defined
composition — it is not, on its own, evidence that these four quadrants
are formally separable biological states** (a continuous F×P gradient
thresholded at the same bar produces an identical-looking table; a
genuine separability claim would need a pre-specified clustering/mixture
model, out of scope here).

## Section 4: revCSC cross-tabs — two different questions, two different answers

**This is the concrete, real-data payoff of round 1's review fix** (which
required splitting one conflated cross-tab into two correctly-labeled
ones) — the two tables give genuinely different, both-real pictures of
the same underlying data:

**Table 4a — axis-supported status × operational revCSC-high cohort**
(top 10% by cross-cell z-score rank, PR #28/#29's exact construction).
Of the 53,415 P-supported cells (any of P_only/D+P/F+P/D+F+P) in the
full atlas, only **5,767 (10.8%)** fall inside the top-10% operational
revCSC-high cohort — **47,648 (89.2%) are outside it**. Read naively,
this looks like strong support for "a P-high population exists outside
what revCSC captures."

**Table 4b — axis-supported status × revCSC matched-null support**
(`revCSC_primary27_minus_CLU_ASS1_percentile`≥90, the same
null-calibrated-evidence semantics used for D/F/P — not a cross-cell
rank). Of the same 53,415 P-supported cells, **26,424 (49.5%) are
revCSC-supported and 26,991 (50.5%) are revCSC-not-supported** — an
almost even split, not the dramatic 89%-outside picture Table 4a
suggests.

**The honest reconciliation, stated directly rather than left as an
apparent contradiction**: Table 4a's "89% outside" number is real, but it
is substantially a mechanical consequence of the operational cohort being
constructed to be exactly 10% of all cells by design — *any* biological
category with a roughly-uniform-or-worse-than-uniform relationship to
revCSC rank will show most of its members "outside" a fixed top-10% cut,
independent of whether real revCSC evidence for those cells is high or
low. Table 4b, which uses revCSC's own natural (and much larger, ~51%
of the atlas) evidence threshold, shows the biologically correct picture:
P-supported cells are **not** concentrated among revCSC-evidence cells,
but they are not concentrated *away* from them either — an essentially
even split. **This same near-even-split pattern holds for every
axis-supported category, including `none`** (133,416 revCSC-supported vs.
124,681 revCSC-not-supported) — reinforcing, on the full 665,473-cell
population this time, the primary analysis's core finding that
revCSC-evidence status is only weakly related to gut-developmental
axis-supported status. Table 4a's cohort-membership framing and Table
4b's evidence-based framing are both real, correctly-labeled, and
answer different questions — reporting only one would have been
misleading in opposite directions.

Full per-category, donor/study-aware detail for both tables is committed
in full (`results/06_crc_projection/tertiary_analysis_composition/`).

## What this shows, honestly

**The full-atlas, revCSC-independent landscape does not surface a large,
qualitatively new developmental state that the revCSC-high framework
entirely misses.** The same bimodal (none/F_only-dominated) axis-supported
pattern found within the revCSC-high subset (Step 6 secondary) holds
across the *entire* malignant population — if anything, F-supported
biology is *marginally less* concentrated in revCSC-high cells than in
the general population, and P-supported biology is only *modestly* more
concentrated there. **This is itself a real, informative negative
result**: it means Step 6 secondary's revCSC-high-cohort composition
findings were not an artifact of restricting to that cohort — they
generalize to the full atlas. The SI-over-Colon regional skew,
independently confirmed here at an even larger magnitude (4.7× vs. 3×),
is the one finding that gets *stronger*, not weaker, when checked
population-wide.

**What this does NOT show**: formal statistical separability of the F×P
quadrants (see §3's explicit boundary); a P-high population that revCSC
"does not capture" in any absolute sense (Table 4a's 89%-outside number
is a real fact about the fixed-size operational cohort, not evidence that
revCSC evidence is absent for those cells — Table 4b shows it usually
isn't); anything about the mechanism behind the F/P split (gene-set size,
biological program identity, or otherwise — not tested here).

## What this does NOT do

Same explicit scope boundary as the approved design: still
`CRC_single_cell_atlas_2025` only (the 2 additional CRC datasets,
`HTAN_CRC_progressive_plasticity` and `CRLM_NMP_ATLAS`, are separate,
already-scoped future work); does not touch M11 (only defined on its own
297,307-cell subset); does not re-derive or re-score any frozen gut
D/F/P or revCSC gene set.

## Review history

Submitting for review before merge, same discipline as every prior step.
