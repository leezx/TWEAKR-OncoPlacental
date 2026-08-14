# F-developmental multi-organ marker panel check — revised cutoff proposal

Per the PR #12 reviewer's explicit request: AFP alone can discriminate
`adult_excl_pct` (10 vs 25) but cannot discriminate `elevated_pct` (75 vs
90) — its within-organ percentile is 99.7, comfortably clearing both. A
real, predefined, cross-organ positive-control panel was needed before
freezing `elevated_pct`.

## Panel (predefined before looking at results)

**Cross-organ pan-fetal/oncofetal reactivation genes** (applied to all 7
organs): `DLK1`, `IGF2`, `H19`, `LIN28B`, `PEG10` — well-established
imprinted/fetal-reactivation genes from the literature, expected to be
reduced in most adult tissues. **Confidence caveat, stated up front**:
this panel is well-supported for broad fetal-vs-adult biology but is NOT
organ-specific developmental biology — genuinely organ-restricted fetal
markers are far better established for liver than for Adrenal/Thyroid/
Skin/Stomach/Spleen/Thymus.

**Liver-specific additions**: `AFP`, `GPC3` — established liver oncofetal
markers.

## Results: this is a real correction to the original proposal

| Combo | Retention (raw) | Retention rate | Mean candidate count |
|---|---|---|---|
| elev75 / adult10 | 2/37 | 5.4% | 394 |
| elev75 / adult25 | 4/37 | 10.8% | 666 |
| elev90 / adult10 | 0/37 | 0.0% | 129 |
| elev90 / adult25 | 1/37 | 2.7% | 231 |

**Raising `elevated_pct` from 75→90 loses 3 of the 4 real marker hits.**
At elev75/adult25, four organ×gene pairs pass: Thyroid-DLK1, Skin-PEG10,
Thymus-IGF2, Liver-AFP. At elev90/adult25, only Liver-AFP survives — DLK1
(Thyroid, percentile 79.8), PEG10 (Skin, percentile 79.1), and IGF2
(Thymus, percentile 75.7) all fall below the 90th-percentile bar. **This
directly contradicts the original proposal's `elevated_pct=90`** — the
evidence argues for **`elevated_pct=75`** instead, since 90 was chosen on
"comfortably above background" reasoning that the reviewer correctly
identified as circular (percentile is organ-internal by definition; 90
being stricter than 75 is definitional, not external evidence).

**`adult_excl_pct=25` remains necessary**: AFP fails at `adult_excl_pct=10`
regardless of `elevated_pct` (HPA's within-liver percentile for AFP is
20.5, clearing 25 but not 10) — same result as the original AFP check,
now confirmed as part of a real panel rather than a single anecdote.

## Two additional real findings from this check

1. **GPC3 fails at every combo** despite being "elevated" (99.96th
   percentile in fetal liver) — it never clears `adult_excluded` at either
   10 or 25. This is a genuine negative case, not a bug: GPC3 is an
   oncofetal marker because its expression is *reactivated above adult
   baseline* in HCC, not because it's completely silent in normal adult
   liver — it retains enough low-level adult expression to fail a strict
   "absent everywhere" bar. Useful calibration data point: not every
   textbook oncofetal gene is a good `adult_excluded` case, and that's
   expected, not a flaw in the criterion.
2. **H19 has a cross-platform gene-coverage gap**: present in GTEx but
   **absent from HPA's gene panel entirely** (not merely undetected — not
   in the `rna_tissue_hpa.tsv` gene list at all). For organs with
   `GTEx+HPA` provenance, the current `adult_excluded_mask()` logic
   silently drops genes missing from either platform's shared index
   rather than falling back to single-platform evidence — meaning H19
   (and presumably any other GTEx-only gene) can never pass F-developmental
   in a GTEx+HPA organ regardless of its real GTEx evidence. This is a
   conservative behavior (both sources must agree) but was not an
   explicit design decision — flagged for reviewer input: keep as
   "both-required," or fall back to single-platform evidence when only
   one platform measures the gene at all.

## Per-organ marker-hit summary (any panel gene passing, elev75/adult25)

| Organ | Any marker retained? |
|---|---|
| Adrenal | No (0/5 panel genes pass at any tested combo) |
| Thyroid | Yes — DLK1 |
| Spleen | No (0/5) |
| Thymus | Yes — IGF2 (HPA-only provenance) |
| Liver | Yes — AFP (GPC3 fails, see above) |
| Skin | Yes — PEG10 |
| Stomach | No (0/5) |

Adrenal, Spleen, and Stomach show zero panel-gene hits under any tested
combo. Most likely explanation: the cross-organ panel just isn't
well-suited to these organs' biology (as flagged in the confidence
caveat above) rather than a defect in the cutoff — but flagged honestly,
not glossed over, since it could also mean `adult_excl_pct=25` is still
too strict for some organs.

## Round 2 fix: coverage-aware adult-exclusion evidence (real bug, not just an open question)

The reviewer correctly identified the H19 finding above as a real bug, not
a minor caveat: the original `adult_excluded_mask()` intersected GTEx's
and HPA's gene panels *before* evaluating, so a gene present in GTEx but
simply outside HPA's gene panel (H19: in GTEx's 73,321 genes, absent from
HPA's 20,151) was silently dropped from the evaluable set entirely —
conflating "platform doesn't measure this gene" with "gene failed the
test." Two different failure modes, same (wrong) result.

**Fixed**: `adult_excluded_mask()` now decides each gene's evidence
provenance individually — `GTEx+HPA` (both platforms measure it, both
must pass), `GTEx-only` or `HPA-only` (only one platform measures it,
evaluated on that platform's evidence alone), or absent from the
evaluable index entirely if neither platform measures it (honestly
unresolved, never silently marked pass or fail). Re-ran both the full
calibration and the marker panel with the fix.

**Verified**: H19 in Thymus now correctly shows
`UNCOVERED_BY_EITHER_PLATFORM` (Thymus has no GTEx column and H19 isn't
in HPA either — genuine absence of evidence) instead of being silently
dropped. H19 everywhere else now shows `GTEx-only` provenance and still
fails `adult_excluded` — but now for a real, attributable reason (real
GTEx evidence says H19 isn't sufficiently excluded in that tissue), not a
coverage artifact. **The 4 passing marker×organ pairs at elev75/adult25
are unchanged** (Thyroid-DLK1, Skin-PEG10, Thymus-IGF2, Liver-AFP) —
confirms the fix only grew the previously-undercounted candidate pools
(e.g. Adrenal 466→681, Thyroid 461→684 at elev75/adult25, from genes
correctly restored via single-platform evidence) without changing any
marker's individual verdict. Retention rate unchanged (4/37 = 10.8% at
elev75/adult25).

## Revised proposal (for reviewer sign-off)

**`elevated_pct=75`, `adult_excl_pct=25`** (quorum still not
discriminating, defaults to 1.0) — revised down from the original
`elevated_pct=90` proposal based on this panel's direct evidence that 90
loses real marker hits without a compensating gain shown by any positive
evidence. Candidate counts under this combo, after the coverage-aware fix
(from `f_developmental_calibration.tsv`), with provenance breakdown:

| Organ | Candidates | Provenance breakdown |
|---|---|---|
| Adrenal | 681 | GTEx+HPA=466; GTEx-only=210; HPA-only=5 |
| Thyroid | 684 | GTEx+HPA=461; GTEx-only=219; HPA-only=4 |
| Spleen | 1,090 | GTEx+HPA=754; GTEx-only=323; HPA-only=13 |
| Thymus | 1,191 | HPA-only=1,191 (no GTEx column; unaffected by this fix) |
| Liver | 800 | GTEx+HPA=516; GTEx-only=270; HPA-only=14 |
| Skin | 1,128 | GTEx+HPA=786; GTEx-only=331; HPA-only=11 |
| Stomach | 750 | GTEx+HPA=545; GTEx-only=200; HPA-only=5 |

Larger than both the original 90/25 proposal and the pre-fix 75/25
numbers, reflecting genes now correctly restored via single-platform
(GTEx-only/HPA-only) evidence that the intersection bug had silently
excluded. Thymus is unchanged (it never had a GTEx column to intersect
against, so the bug never affected it).

**Pending final reviewer sign-off** on both `elevated_pct=75,
adult_excl_pct=25` and the coverage-aware provenance model.
