# Step 4a design: external adult-expression audit of the gut D/F/P gene sets

Per `docs/PROJECT_SUMMARY.md`'s "What's next": check the frozen
gut-specific gene sets (`results/04a_dfp_gut/dfp_gut_gene_sets/`) against
adult references beyond the Gut Cell Atlas's own internal fetal-vs-adult
contrast. **Revised after round-1 REQUEST_CHANGES** — the reviewer found
4 real design-level problems (independence framing was factually wrong
for D/P, the F-arm's own scientific question was mis-stated, no
cross-dataset gene-ID mapping contract existed, and a claimed "reused"
permutation contract for D/P whole-body doesn't actually exist in the
approved Step 5 code) plus several lower-severity but real fixes. All
addressed below before any compute is run — still no qsub job submitted.

## What this actually tests (reframed, not "adult-negative validation")

**Two different scientific questions, kept explicitly separate — conflating them was round-1's core problem:**

1. **For `F_Colon-developmental`/`F_SI-developmental`/`F_Gut-core`
   (built solely from the Gut Cell Atlas's own fetal-vs-adult contrast,
   `FDR<0.05 & logFC>1`, fetal significantly *higher* than adult)**: the
   frozen definition does **not** require adult expression near zero —
   only that fetal is significantly higher. A gene can be real, correct
   fetal-up biology and still show real, detectable adult expression.
   Tabula Sapiens and GTEx are **adult-only** data (no fetal arm), so
   they structurally cannot re-test "is fetal > adult" — that question
   is already answered, inside the Gut Cell Atlas, by the frozen edgeR
   result. What TS/GTEx *can* answer is a genuinely different, still
   useful question: **how much residual adult expression does each
   fetal-up gene carry, externally, in independent adult tissue/cell-type
   data** — an adult-expression / adult-specificity **audit**, not a
   validate/refute test of the F sets themselves. A significant hit here
   is a **red flag worth manual scrutiny** (evidence the gene may not be
   very adult-specific, useful context for later CRC interpretation);
   non-significant is **inconclusive**, not proof of adult-negativity —
   this asymmetry must be stated wherever results are reported, not just
   here.
2. **For `D_Colon-shared`/`D_SI-shared`/`P_Colon-specific`/`P_SI-specific`
   (all built by intersecting/subtracting `F_{Colon,SI}-developmental`
   against the frozen `P_developmental_primary84.txt`)**: `P_developmental`
   itself already required whole-body `adult_excluded` evidence (GTEx +
   HPA) as part of its own Step 4 construction. So for these four sets,
   GTEx is **not** an independent adult reference — checking them against
   GTEx again is a **construction-consistency / provenance audit**: does
   the frozen membership still look internally consistent against the
   same kind of evidence that helped define `P_developmental` in the
   first place, not a new independent test. Tabula Sapiens, by contrast,
   **is** genuinely independent for these four sets too (never touched
   any part of the pan-organ or gut D/F/P construction) — a real
   held-out check.

`F_Gut-core` (`F_Colon-developmental ∩ F_SI-developmental`, verified
directly against `build_dfp_gut_gene_sets.py`) inherits **neither**
whole-body constraint — it is excluded from every whole-body-style check
below (round-1 finding, unchanged from the first fix).

## Data sources (both already on disk, verified in prior steps — not re-downloaded)

- **Tabula Sapiens `Large_Intestine`** (`datasets/TabulaSapiens/`, also on
  Argos): 13,680 cells, 2 donors (`TSP14`, `TSP2`), raw counts in
  `layers['raw_counts']` (confirmed 100% integer-valued, Step 5 inventory).
  `compartment == "epithelial"` cell types confirmed present, 6,220 cells
  (`enterocyte of epithelium of large intestine`, `immature enterocyte`,
  `mature enterocyte`, `paneth cell of epithelium of large intestine`,
  `large intestine goblet cell`) — matches `F_Colon-developmental`'s own
  epithelial-lineage-only construction population. **Real, stated gap**:
  no `Small_Intestine` organ file — `F_SI-developmental` gets no
  single-cell check from this source.
- **GTEx v11 median TPM** (`datasets/GTEx_v11_median_tpm/`): confirmed
  columns include `Colon_Sigmoid`, `Colon_Transverse` and
  `Small_Intestine_Terminal_Ileum`. **`Small_Intestine_Terminal_Ileum` is
  partial anatomical coverage** (terminal ileum only, not the whole small
  intestine) — a low result there supports "adult terminal ileum shows
  low bulk expression," not a claim about the entire adult small
  intestine; stated this way in reporting, not "GTEx fills the SI gap."

## Compute contract: cross-dataset gene-ID mapping (frozen before any compute — round-1 blocker)

The frozen `F_Colon-developmental`/`F_SI-developmental` gene lists contain
Gut Cell Atlas `var_name`s directly, including anndata duplicate-symbol
artifacts (e.g. `IGF2-1`, not canonical `IGF2` — the exact bug class PR #21
already found and fixed once for the marker panel). GTEx uses canonical
symbol + versioned Ensembl ID; Tabula Sapiens uses what Step 5 confirmed
are native/canonical symbols. A naive string-match of 1,456 raw GCA
`var_name`s against either external reference would silently mis-handle
every renamed/duplicate-suffixed gene — undercounting real matches as
"not found" and, worse, letting some slip into being counted as
"adult-negative" by default. **Frozen mapping rule, reusing PR #21's own
audit output (`results/04a_dfp_gut/gene_id_audit/var_id_map.tsv`), not a
new lookup:**

1. Every GCA `var_name` → `gene_id` (Ensembl) → `authoritative_symbol`
   (BioMart), exactly as already computed and verified in PR #21.
2. **GTEx**: match primarily by Ensembl ID (GTEx's own `Name` column,
   version-stripped); fall back to `authoritative_symbol` string match
   only if Ensembl match fails.
3. **Tabula Sapiens**: match by `authoritative_symbol` (not the raw GCA
   `var_name`) against TS's own `var_names` — avoids exactly the
   `IGF2-1`-vs-`IGF2` mismatch class.
4. Any GCA gene whose `authoritative_symbol` is unresolved (BioMart
   `found_in_biomart=False`) or whose mapped ID/symbol is absent from the
   external reference is marked **`NOT_TESTABLE`** — reported as its own
   explicit category, **never** silently folded into "adult not detected."
5. Every check's output reports `n_input` / `n_mapped` / `n_present_in_reference`
   / `n_not_testable` alongside the results, so the denominator behind
   every summary statistic is auditable.

## Method

### Check 1: Tabula Sapiens, epithelial-restricted, organ-matched (`F_Colon-developmental` only)

Same donor-aware pseudobulk contract as `tier2_validation.py` (PR #14/#15):
primary unit `(donor, cell_ontology_class)` pseudobulk from raw counts,
CPM per pseudobulk library's own size, minimum 20-cell floor at the
donor×cell-type level, cross-donor-consistent detection (CPM≥1 in **all**
eligible donors, ≥2 eligible donors) as the primary flag — reported and
interpreted per the reframed question above (a flag = adult-expression
red flag, not a pass/fail on the gene's validity). Restricted to
`compartment == "epithelial"` cell types (confirmed present, 6,220 cells,
2 donors — both must be eligible for cross-donor-consistent status,
labeled `single-donor evidence` otherwise, same as Step 5's Thymus/Liver
precedent).

**Optional permutation layer, only if run, with its estimand corrected**:
the original `background_detection_rate.py`'s detectability-matched null
was computed from *all* eligible adult cell types (whole-organ), which
would silently pull the matching covariate back to a whole-organ
estimand the moment the primary hypothesis is epithelial-restricted. If a
null is run here, its detectability covariate must be computed **only
from the eligible epithelial donor×cell-type pseudobulks**, not the
original whole-organ pool. Its p-value supports "adult expression is
anomalously high in this cell type relative to matched genes" (evidence
against specificity) — a non-significant result is **inconclusive**, not
evidence of adult-negativity (same asymmetry as above, restated here
because this is where it's easiest to misread a null result as a pass).

### Check 2: GTEx bulk, adult-expression audit for F / construction-consistency audit for D+P

Same not-detected-floor + percentile-among-detected mechanics as
`adult_excluded_percentile_audit.py`, relabeled **`below_bulk_expression_floor`**
throughout (not "not detected") — `TPM<1` is this project's own reporting
convention, not an official platform detection rule, and is not
comparable in strength to Tabula Sapiens' `CPM≥1` donor-level floor; the
two are never combined into one "negative" verdict.

- `F_Colon-developmental` vs. `Colon_Sigmoid` + `Colon_Transverse`
  (adult-expression audit, per the reframed F question above).
- `F_SI-developmental` vs. `Small_Intestine_Terminal_Ileum` (same, with
  the partial-anatomical-coverage caveat stated above).
- `D_Colon-shared`, `D_SI-shared`, `P_Colon-specific`, `P_SI-specific` vs.
  all 68 GTEx tissues (**construction-consistency audit**, not
  independent validation — corrected from round-1's mislabeling; now
  covers **both** D and P consistently, fixing round-1's inconsistency
  where P was checked by Check 3 but not Check 2). `F_Gut-core` excluded
  (no whole-body constraint in its own definition).

### Check 3: Tabula Sapiens, all 5 organs — D/P construction-consistency + genuinely-independent audit

**Reuses Step 5's actual primary contract, not an invented permutation
layer** (round-1 correction: the original whole-body D-shared/P-specific
check in `tier2_validation.py` is direct gene×cell-type×donor pseudobulk
evidence — cross-donor-consistent CPM≥1 flags — with **no** permutation
test attached; the 500-permutation matched null lives in a *different*
script, `background_detection_rate.py`, and only ever covered the
organ-matched F checks. Claiming D/P whole-body "reuses the same
permutation null" was itself a design error, now fixed by reusing what
the code actually does rather than what round-1's doc assumed it did).

`D_Colon-shared`, `D_SI-shared`, `P_Colon-specific`, `P_SI-specific` —
direct gene×cell-type×donor flags across all 5 organs (Liver, Skin,
Spleen, Thymus, Large_Intestine), no lineage restriction (matches these
sets' own whole-body definition). **Unlike Check 2, this is genuinely
independent** for all four sets (Tabula Sapiens never touched any part of
either D/F/P construction). If a permutation layer is wanted for this
check later, it needs its own fresh design (null pool, matching scope,
organ stratification, hypothesis unit, FDR family) — explicitly out of
scope here, not silently borrowed from an unrelated script.

**`P_Colon-specific`/`P_SI-specific` are not two independent panels**:
verified directly (`comm -12`), 76 of 79/80 genes are identical between
them (`D_Colon-shared`/`D_SI-shared` share only `TRIM71`, so the two
`P_*-specific` sets differ by only a handful of genes each). Running two
separate panel-level tests and presenting concordant results as "two
independent supports" would be pseudo-replication. **Fixed design**: the
gene-level compute runs once on the union of unique genes across
`P_Colon-specific ∪ P_SI-specific`, with each gene separately annotated
for `P_Colon-specific`/`P_SI-specific` membership; summaries can be
sliced by membership afterward but are explicitly labeled as
overlapping, not independent, evidence.

### Hypothesis family and FDR scope (explicit, was undefined in round-1)

If any permutation-based inferential test is run (Check 1 only, and only
if the optional null layer above is included), its BH-FDR family is
scoped to **exactly the organ-matched epithelial `F_Colon-developmental`
tests**, not pooled with Check 2/3's direct-evidence audits (which report
raw flags/percentiles, not p-values, and so have no FDR family of their
own). 500 permutations floor the smallest reportable empirical p at
~0.002 — if the eventual epithelial cell-type × gene hypothesis count
grows large enough that this resolution becomes limiting, permutation
count will be increased before finalizing, not left under-resolved.

### F-specific summary (derived, no extra compute)

`F_Colon-specific`/`F_SI-specific` (`F_{region}-developmental \
P_developmental`) are not separately tested against external data — no
new compute needed, since Check 1/2's results on `F_{region}-developmental`
already cover every gene they contain. But the write-up will explicitly
derive and report an `F_{region}-specific` slice of those results (which
genes in the audited set are F-specific vs. D-shared), so the full D/F/P
partition story is visible, not just the F-developmental precursor —
round-1 correctly flagged that skipping this entirely would leave the
partition-level narrative incomplete even though no extra compute is
required.

## What this does NOT do

- Does not re-open or re-run any part of the gut D/F/P's own construction
  (PR #19/#20/#21) or the mike_verzi external validation (PR #22) — this
  is strictly additive.
- Does not propose a new cutoff or gate any gene's membership in the
  frozen sets.
- Does not claim a non-significant result "validates" adult-negativity
  for F-arm checks (see reframing above) — inconclusive is reported as
  inconclusive.

## What's next after this design is approved

Real compute: `scripts/04a_dfp_gut/gut_adult_validation_tabula_sapiens.py`
(Checks 1+3) and `scripts/04a_dfp_gut/gut_adult_validation_gtex.py`
(Check 2), run via qsub, results pulled back byte-exact, written up
honestly (including any null or contradictory findings), submitted for
review, merged only after user confirmation.
