# Step 4a design: re-anchoring F-developmental on real fetal gut/colon data

Builds on `docs/STEP4_DFP_DESIGN.md` / `docs/STEP4_STATISTICAL_DESIGN.md` (the frozen D/F/P logic and its statistical design) and directly addresses a scope mismatch the user flagged after reviewing the mike_verzi enrichment result: HDMA's 7-organ F-developmental reference (Adrenal/Liver/Skin/Spleen/StomachEsophagus/Thymus/Thyroid — `datasets/HumanDevelopmentMultiomicAtlas/dataset.md`) has **zero gut/colon representation**, which is scientifically mismatched for a CRC project. This doc defines `F_Colon-developmental` and `F_SI-developmental` — built from real human fetal gut single-cell data — as the new primary F-axis for CRC-facing work, before any qsub compute runs, matching the review-before-compute discipline that produced `STEP4_DFP_DESIGN.md`.

**Revised per round-1 REQUEST_CHANGES**: the reviewer's core correction — since the Gut Cell Atlas itself spans fetal through adult intestine in one processed object (`epi_raw_counts02_v2.h5ad`, epithelium-only, full lifespan), the primary construction method should be a **same-atlas fetal-vs-adult epithelial pseudobulk DE**, not a mechanical port of HDMA's single-population "elevated-in-fetal percentile + external GTEx/HPA adult-percentile-exclusion" design. GTEx/HPA/Tabula Sapiens are demoted further, to **external adult-negative validation only** (they were adult-exclusion *co-inputs* to the definition in the first draft below; they're not that anymore). See "Revised primary construction" below — the first-draft percentile-based definitions are kept struck through/superseded rather than deleted, so the review history stays legible.

## What changes and what doesn't

- **P-developmental is unchanged.** It was already computed as a whole-body (not organ-scoped) program — `replicated_in_placenta(gene)` AND `adult_excluded(gene, whole_body)` — so it needs no gut-specific rebuild. It's reused as-is from the frozen Step 4 output (`results/04_dfp_signature/dfp_gene_sets/`).
- **F-developmental is rebuilt**, replacing HDMA's 7 pan-fetal organs with two real fetal gut regions as the primary reference:
  - `F_Colon-developmental` (large intestine) — **primary**, matches CRC's tissue of origin directly.
  - `F_SI-developmental` (small intestine) — **secondary, parallel**, run with the identical method for cross-region comparison.
  - A pooled `F_Gut-developmental` ("Gut-core" = Colon ∪ SI, or Colon ∩ SI as the stricter alternative) is reported only as a **tertiary summary**, and only if Colon and SI turn out concordant enough for pooling to be meaningful — not assumed up front.
- **D/F/P set operations are redefined** using the new gut-specific F, per the user's explicit formulas:
  - `D_Gut-shared = F_Gut-developmental ∩ P_developmental`
  - `F_Gut-specific = F_Gut-developmental \ P_developmental`
  - `P_Gut-specific = P_developmental \ F_Gut-developmental`
  - (Colon and SI versions computed the same way against the same unchanged `P_developmental`, before any Colon/SI pooling decision.)
- **The existing pan-organ HDMA-based D/F/P is not deleted.** It's demoted to a secondary "pan-fetal / non-gut developmental validation framework" — useful for asking "is this gene generically fetal, anywhere in the body" but no longer the primary CRC coordinate system.

## Reference data: Gut Cell Atlas (Elmentaite et al., Nature 2021; first-trimester fetal object originally Elmentaite et al., Developmental Cell 2020)

**Provenance correction (round-1 review)**: the 62,849-cell, 6–11 PCW `fetal_RAWCOUNTS_cellxgene.h5ad` object is first-trimester data **previously published in Elmentaite et al., Developmental Cell 2020**, which the Nature 2021 Space-Time Gut Cell Atlas paper integrated alongside its own newly-generated second-trimester (12–17 PCW) and adult (29–69 yr) samples. Citing it as "Elmentaite et al., Nature 2021" alone (as the first draft of this doc did) is imprecise — both citations are now given.

Three files, two download rounds:

| File | Cells | Genes | Age range | Role |
|---|---|---|---|---|
| `fetal_RAWCOUNTS_cellxgene.h5ad` | 62,849 | 33,694 | 6–11 PCW (first-trimester only, Dev Cell 2020 cohort) | Raw counts (verified integer, `X_integer_fraction=1.0`) — early-fetal developmental-stage sensitivity/validation arm (see revised layering below), not primary construction |
| `final_fetal_object_cellxgene.h5ad` | 62,849 | 26,757 (HVG subset) | 6–11 PCW | Normalized, has `X_pca`/`X_umap` — metadata cross-check only |
| `epi_raw_counts02_v2.h5ad` | 142,113 | 33,538 | Full lifespan: Adult 77,341, First trim 31,689, Second trim 20,495, Pediatric 8,157, Pediatric_IBD 4,398, Adult_MLN 23, Second trim_MLN 10 | Raw counts, epithelium-only (`category`: Epithelial 142,104 / Mesenchymal 9, negligible) — **primary construction data**, downloaded this round (job 3620693) and inventoried (job 3620698), size/md5 verified byte-exact against source: 666,604,440 bytes / `2a149b8cf04567569707e9d1fab27209` |

Real structure confirmed by direct inventory (not assumed):

- **Raw counts confirmed** (`X_integer_fraction=1.0`); gene IDs: symbols in `var_names` + `gene_ids`/`feature_types` columns in `var`.
- **Region** (`obs['Region']`, 5 categories): `SmallInt` 73,023 / `LargeInt` 47,020 / `REC` (rectum) 17,348 / `APD` (appendix) 4,689 / `lymph node` 33. Finer `Region code` (18 categories, e.g. `DUO`/`JEJ`/`ILE`/`CAE`/`SCL`/`TCL`/`ACL`/`FPIL`/`FTIL`/`FMIL`/`FLI`) available for sub-region granularity if needed later.
- **Age × Region cross-tab (the fact that actually determines the primary contrast)**: fetal samples (First trim + Second trim) only have `LargeInt` and `SmallInt` — **no `REC`/`APD` fetal samples exist in this atlas**. So "colon" for the primary fetal-vs-adult contrast must be `LargeInt` specifically (not `REC`, which is adult-only and therefore can't be part of a fetal-vs-adult comparison):
  | Age_group | APD | LargeInt | REC | SmallInt |
  |---|---|---|---|---|
  | Adult | 4,689 | 38,612 | 17,348 | 16,692 |
  | First trim | 0 | 4,927 | 0 | 26,762 |
  | Second trim | 0 | 3,481 | 0 | 17,014 |
- **Chemistry claim, checked and partially revised**: the reviewer's claim (Second trim + Adult share 10x 5′ chemistry, First trim is a separate mostly-3′ batch) is **directionally correct but not a clean/exclusive split** — real cross-tab: Adult is 39,738/77,341 (51.4%) 5′ vs 37,603 (48.6%) 3′ (a genuine mix, not uniformly 5′); Second trim is 17,055/20,495 (83.2%) 5′; First trim is 27,578/31,689 (87.0%) 3′. So Second trim does skew heavily 5′ like Adult, and First trim does skew heavily 3′ unlike both — the qualitative "First trim is the more different batch" conclusion holds, but "Second trim and Adult share chemistry" is an approximation, not exact — the primary DE model should account for `10X` chemistry as a covariate/blocking factor rather than assume it away by cohort choice alone.
- **Donor overlap with the already-downloaded fetal-only object, confirmed**: First trim's 10 donors here (`BRC2121`, `BRC2043`, `BRC2133`, `BRC2258`, `BRC2134`, `BRC2119`, `BRC2046`, `BRC2049`, `BRC2026`, `BRC2029`) directly match (modulo the `BRC` prefix) 9 of the 9 `Donor_id`s already inventoried in `fetal_RAWCOUNTS_cellxgene.h5ad` (2026/2029/2043/2046/2049/2119/2121/2133/2134), plus one new donor (`BRC2258`) not in that file. **This confirms real donor/data overlap between the two downloaded files** — expected, since both trace back to the same Dev Cell 2020 first-trimester cohort integrated into the Nature 2021 atlas, but worth stating explicitly rather than treating the two files as independent evidence of anything.

Real structure confirmed by direct inventory, not assumed:

- **Gene IDs**: symbols in `var_names`, verified real Ensembl IDs in `var['gene_ids']` (33,694 unique `ENSG...` IDs, 1:1 with 33,694 unique symbols, no collisions) — same dual-ID discipline as HDMA/placental datasets, enables clean cross-dataset alignment.
- **Region** (`obs['Organ']`, 3 categories): `duojejunum` (duodenum+jejunum combined) 21,592 cells, `ileum` 20,110 cells, `colon` 21,147 cells. SI = duojejunum ∪ ileum; LI = colon.
- **Lineage** (`obs['cell_type_group']`, 4 categories): `epithelium` 16,937 / `mesenchymal` 40,671 / `vasculature` 3,274 / `immune` 1,967. Finer `cell_name_detailed` (28 categories) includes `Colonocytes`, `Small intestinal Epi`, `Enterocyte`, `Early/G2M/S enterocytes`, `LGR5 stem`, `BEST4+ enterocyte`, `Secretory Epi`.
- **Donor/sample structure** (real replicate structure, checked directly — same discipline as Step 4's HDMA per-organ-per-sample table):

  | Region | Epithelium-only cells | Donors contributing | Samples |
  |---|---|---|---|
  | duojejunum | 7,793 | 9/9 | 9 |
  | ileum | 6,302 | 9/9 | 10 |
  | colon | 2,842 | 9/9 | 9 |

  All 9 donors contribute epithelial cells to all 3 regions — no region is donor-confounded. Colon epithelium is the smallest (2,842 cells; thinnest donor contributes 64 cells) but still has full donor coverage, matching the Step 4 precedent of proceeding with organs that have 3–7 samples (Thyroid/Spleen/Thymus/Skin had only 3 each).
- **Age**: `PCW` (F6.1–F10.2) and `CRL` (crown-rump length, 9 categories) both present as parallel age proxies — same 9 donors, consistent with the previously-verified "6-11 PCW" description of this atlas.

## Primary population: epithelial lineage

Per the user's explicit direction, `cell_type_group == 'epithelium'` is the primary subset for the fetal-gut side of this analysis — it matches the CRC malignant-epithelial-cell-state focus the whole D/F/P signature exists to serve. Mesenchymal/vasculature/immune populations in this atlas are not used for `F_Colon-developmental`/`F_SI-developmental` construction (they may be useful later for stroma-context questions, out of scope here).

## Revised primary construction: same-atlas fetal-vs-adult epithelial DE (round-1 review fix)

**The first-draft definitions below this point (§"Superseded: HDMA-style percentile port") are struck through, not deleted** — the reviewer's fix isn't a small patch, it changes the primary evidence model entirely, and the review history should stay visible in the doc rather than silently vanish.

The correction: `epi_raw_counts02_v2.h5ad` already contains matched fetal *and* adult epithelium in one processed object, from the same atlas/pipeline. Using this internal contrast — the same "positive evidence from a real DE, not an isolated percentile against an external platform" logic already used for P-developmental's trophoblast-vs-rest arm in `STEP4_STATISTICAL_DESIGN.md` §1 — is strictly stronger than HDMA's design (which had no choice but to reach for an external adult reference, because HDMA itself is fetal-only). Reusing HDMA's percentile-based recipe here would have been "matching Step 4's *code*" at the cost of not matching Step 4's *actual statistical intent* (positive evidence should come from the best real contrast available, and here a same-cohort contrast is available).

- **Primary comparison**: fetal-`LargeInt`-epithelium vs. adult-`LargeInt`-epithelium, donor pseudobulk, within `epi_raw_counts02_v2.h5ad`. **Not donor-paired** — confirmed via inventory: First-trim/Second-trim donors (`BRC2xxx` codes) and Adult donors (`A2x`/`A3x`/`T0xx` codes) are structurally disjoint (different individuals, as expected for fetal-vs-adult), so the model fits an unpaired two-group edgeR/DESeq2 quasi-likelihood design.
  - **Round-2 correction — pooling First+Second trimester "because it has more cells" was rejected on review, correctly**: cell count is not the relevant unit (donor/pseudobulk is), and cell count was never a valid reason to prefer one grouping over another anyway. A real **pseudobulk-unit design audit** was run instead (`scripts/04a_dfp_gut/gut_epi_pseudobulk_design_audit.py`, job 3620700) per the reviewer's required gate, before locking anything:
    - **Donor structure, confirmed**: 0/41 donors span >1 `10X` chemistry value (donor-level pseudobulk aggregation is safe — no donor internally mixes chemistries) and 0/41 span >1 `Age_group` (no donor age-contamination). 22/41 donors do span >1 `batch` value and 6/41 span >1 `Fraction` (FACS sort) value — expected technical replicate structure within a donor, not a biological confound, and doesn't block donor-level aggregation.
    - **Real donor counts, `LargeInt`**: Adult 7 donors (2×3′, 5×5′), First trim 10 donors (9×3′, 1×5′ — chemistry is nearly collinear with cohort here), Second trim 5 donors (2×3′, 3×5′).
    - **The reviewer's factual concern about the Nature 2021 paper's chemistry claim is directly confirmed by this donor-level table**: First trim is overwhelmingly one chemistry (9 of 10 donors 3′) — chemistry and cohort are *nearly* confounded in First trim, exactly the risk the reviewer flagged, and pooling First+Second trim would have imported that near-collinearity into the primary contrast. Second trim, by contrast, has real donors on both chemistries (2×3′, 3×5′) and so does Adult (2×3′, 5×5′) — **a genuinely non-collinear, checkable design**.
    - **Decision, driven by this real table, not assumed**: **Second-trimester vs. Adult `LargeInt` (5 vs. 7 donors, `~ 10X + Age_group`) is the primary DE contrast** — matches the reviewer's own suggested resolution ("if Second trimester vs Adult in a matched source/protocol stratum has enough donors, prioritize it as primary"). First-trimester `LargeInt` (10 donors) moves fully into the early-fetal developmental-stage sensitivity tier (layering below) — not pooled into the primary model. `SmallInt` mirrors this: Second trim 5 donors (2×3′, 3×5′) vs. Adult 5 donors (1×3′, 4×5′), First trim 11 donors (9×3′, 2×5′) as sensitivity. `10X` chemistry stays in the primary model as an explicit covariate (not perfectly balanced even within Second-trim-vs-Adult, but not collinear with `Age_group` either) — real design-matrix identifiability, not assumed adequate.
  - **Region**: `LargeInt` primary (not `REC` — fetal has zero `REC` samples, so rectum can't be part of a fetal-vs-adult contrast; `REC` stays available as an adult-only descriptive/exploratory note, never a D/F/P input), `SmallInt` secondary parallel.
  - **Population**: `category == 'Epithelial'` (142,104/142,113 cells) — unchanged in spirit from the first draft's `cell_type_group == 'epithelium'` choice, now applied to this atlas's own column name.
- **Continuous statistic**: `T_g = fetal-vs-adult-colon-epithelium DE statistic` (signed logFC or an equivalent effect-size+significance combination from the edgeR/DESeq2 fit) — this is the exact continuous ranking statistic the later preranked-GSEA validation plan needs, and it now falls directly out of the primary construction itself rather than needing to be separately invented later.
- **`F_Colon-developmental(gene)` = {g : logFC_fetal/adult(g) > c, FDR(g) < q}** — real `c`/`q` cutoffs to be calibrated against this dataset's actual DE-statistic distribution (marker-panel-check discipline, e.g. LGR5/OLFM4 for crypt stem, MUC2 for goblet, plus the DLK1/IGF2/H19/LIN28B/PEG10 cross-tissue panel used for HDMA's calibration in PR #12), not assumed to transfer from HDMA's `elevated_pct=75`/`adult_excl_pct=25` (those percentile parameters don't even apply to a DE-statistic-based definition — a genuinely new calibration step, not a reused number).
- **`F_SI-developmental(gene)`**: identical method, fetal-vs-adult SI epithelium contrast within the same atlas object.

### Revised role of GTEx/HPA/Tabula Sapiens: external adult-negative validation only

These stay in the project, but move one tier down from what the first draft proposed: they are **not** part of `F_Colon-developmental`/`F_SI-developmental`'s definition anymore (the same-atlas adult epithelium supplies that role now). Their role is now purely a downstream *sanity check* — "does a gene called `F_Colon-developmental` also stay low in fully independent, bulk, adult-population-level colon data (GTEx `Colon_Sigmoid`/`Colon_Transverse`, HPA `colon`/`rectum`) and SI data (GTEx `Small_Intestine_Terminal_Ileum`, HPA `small intestine`)?" — reported alongside the primary result, never substituted for it. Tabula Sapiens `TS_Large_Intestine.h5ad` stays at the same Tier-2 validation role it already had.

## Revised layering (round-1/round-2 review)

1. **Primary discovery**: `epi_raw_counts02_v2.h5ad`, **Second-trimester-fetal vs. Adult** epithelial donor-pseudobulk DE (`LargeInt` primary, 5 vs. 7 donors; `SmallInt` secondary parallel, 5 vs. 5 donors) — the specific contrast identified by the round-2 design audit as the one with a real, non-collinear chemistry/cohort structure.
2. **Early-fetal developmental-stage sensitivity/validation**: **First-trimester `LargeInt`/`SmallInt` epithelium within `epi_raw_counts02_v2.h5ad`** (10/11 donors respectively — moved here from the primary model per the round-2 design audit, since chemistry is near-collinear with cohort in First trim) **and** the separately-downloaded 6–11 PCW, 62,849-cell first-trimester-only fetal object (`fetal_RAWCOUNTS_cellxgene.h5ad`, known to share donors with this First-trimester cohort — see donor-overlap note above) — checks whether `F_Colon-developmental`/`F_SI-developmental` genes are also elevated at this earlier developmental stage, not primary evidence.
3. **External adult-negative validation**: GTEx colon/terminal-ileum, HPA colon/rectum/small-intestine, Tabula Sapiens Large Intestine — sanity-check only, per the demotion above.
4. **External fetal replication (later, deferred)**: GSE158702 (Fawkner-Corbett et al., Cell 2021), GSE95630/GSE103239 (Nature Cell Biology 2018) — not downloaded yet; a donor/dataset-provenance overlap audit against Gut Cell Atlas is required before either can be called "independent replication" (Gut Cell Atlas's own first-trimester arm is itself previously-published Dev Cell 2020 data, so overlap is a real, not hypothetical, risk to check).

## §ARCHIVED — Superseded: HDMA-style percentile port (first draft, kept for review-history legibility only, not used)

~~- **`elevated_in_fetal_gut_epithelium(gene, region)`**: gene's expression in per-sample epithelial pseudobulk for that region clears a within-region percentile threshold, computed from the region's own sample distribution — same `elevated_pct` logic as HDMA (frozen at 75 in PR #12).~~
~~- **`adult_excluded(gene, matched_organ)`**: gene's expression in the matched adult reference (GTEx/HPA) falls below that dataset's own internal percentile.~~
~~- **`F_Colon-developmental(gene)` = `elevated_in_fetal_gut_epithelium(gene, colon)` AND `adult_excluded(gene, colon-matched)`**~~
~~- **`F_SI-developmental(gene)` = `elevated_in_fetal_gut_epithelium(gene, duojejunum ∪ ileum)` AND `adult_excluded(gene, SI-matched)`**~~

The GTEx v11 whole-body file's `Small_Intestine_Terminal_Ileum` column and HPA's `small intestine` row (discovered this session, no new download needed) remain a real, useful finding — they're just repurposed into the "external adult-negative validation" tier above instead of the (now-superseded) primary-definition role.

## Non-circularity boundary (explicit, per user directive)

The 5 `mike_verzi_fetal_signature.gmt` mouse gene sets (YAP_SIGNALING_GENES, REVIVAL_STEM_CELL_GENES, FETAL_SPHEROID_EPITHELIUM_GENES, REGENERATIVE_EPITHELIUM, FETAL_INTESTINE_GENES) **must not be used anywhere in `F_Colon-developmental`/`F_SI-developmental` construction** — no gene-selection, no threshold-tuning, no candidate filtering. Correct causal order: build `F_Colon-developmental`/`F_SI-developmental` independently from the Gut Cell Atlas data first; only afterward use the (fully independent, human-orthology-verified) mike_verzi signatures as external validation, per the three-layer statistical validation plan below. This mirrors the primary/extended-set-ordering lesson from PR #17 round 3 (Ly6a) — outcome-dependent evidence must never influence which genes make it into a primary definition.

## What this design does NOT do yet

No qsub DE compute, no gene lists, no thresholds chosen — same discipline as `STEP4_DFP_DESIGN.md`. Explicitly left open, to decide against real data:

- ~~`epi_raw_counts02_v2.h5ad`'s real structure~~ — **resolved this round** (job 3620693 download + job 3620698 inventory): 142,113 cells, `Diagnosis`/`Age_group`/`Region`/`10X`/`Sample name` columns confirmed, fetal-vs-adult unpaired (disjoint donors), chemistry claim checked and refined, region resolved to `LargeInt`/`SmallInt` as the usable fetal-vs-adult-matched categories.
- ~~Which fetal cohort(s) form the primary DE contrast~~ — **resolved this round** (job 3620700, pseudobulk-unit design audit): Second-trimester vs. Adult (real, non-collinear chemistry structure) is primary; First-trimester moves to the sensitivity tier (near-collinear chemistry/cohort in that cohort). See "Round-2 correction" above.
- Real `c`/`q` cutoffs for `F_Colon-developmental`/`F_SI-developmental`'s DE-statistic-based definition — a new calibration against this dataset's actual Second-trim-vs-Adult DE distribution, not a reused HDMA percentile (the percentile parameters don't even apply to a different statistic family). Requires the actual edgeR/DESeq2 fit to run first, not part of this design doc.
- Whether the `10X` chemistry covariate in `~ 10X + Age_group` adequately absorbs the residual chemistry imbalance even within the Second-trim-vs-Adult contrast (not perfectly balanced: 2×3′/3×5′ vs. 2×3′/5×5′ for `LargeInt`), or whether a chemistry-stratified sensitivity analysis is also needed — to be checked against the real fit's diagnostics, not assumed adequate in advance.
- Whether First-trimester and Second-trimester `LargeInt`/`SmallInt` DE results (primary vs. sensitivity) turn out concordant — a real check to report, not assumed either way.
- Whether duojejunum and ileum need to stay split with a union/consensus reconciliation step, or can be pooled directly into one SI sample set — decide once the two regions' real DE-gene overlap is visible.
- Whether pooled `F_Gut-developmental` (Colon ∪ SI or ∩) is reported at all — contingent on Colon/SI concordance once both are built.
- The three-layer statistical validation (hypergeometric enrichment with CI/OR, preranked GSEA against the continuous fetal-vs-adult-colon DE statistic, size-and-expression-matched permutation nulls) against the 5 mike_verzi signatures — this is the next design doc after `F_Colon-developmental`/`F_SI-developmental` gene sets exist, not part of this one.
- Donor/dataset-provenance overlap audit between Gut Cell Atlas and the two other candidate fetal-gut datasets (GSE158702, GSE95630/GSE103239) — needed before any future "independent replication" claim across those atlases; deferred per the reviewer's round-1 note that this isn't a current blocker.
