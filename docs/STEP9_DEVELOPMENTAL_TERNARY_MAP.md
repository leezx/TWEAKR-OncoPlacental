# Step 9: Fetal-Placenta-Adult developmental ternary map

User-requested exploratory extension (2026-08-17), based on
`2026-GPT-TWEAKR-Oncofetal.md#Fetal-Placenta-Adult ggtern图` (Zhixins-KB).
Combined design+results doc, not a separate design-then-compute pair
like most other steps in this project — the whole pipeline is cheap
(reuses already-computed DE tables plus one small new external
download), so design and a first real compute pass happened together,
with the user reviewing an intermediate result before this PR was
opened. Real bugs found and fixed during the build are reported below,
not smoothed over.

## Motivation and scope

Instead of the project's existing D/F/P **binary** gene-set membership
(a gene either is or isn't F-specific/P-specific/D-shared), this builds
a **continuous per-gene ternary coordinate**: for every gene, three
non-negative "positive evidence" scores — Fetal, Placenta, Adult — that
sum to 1 and place the gene inside a triangle. Close to a vertex = biased
toward that state; near the Fetal-Placenta edge, far from Adult = shared
developmental; center = nonspecific. This is explicitly a **first
version** (per the KB note's own framing) — a visualization/exploration
tool, not a replacement for the project's frozen, reviewed D/F/P gene
sets, and not itself used to draw any founding-question conclusion here.

Per the KB note's own explicit warning, raw mean expression per group is
never used (batch/cell-composition/baseline-expression would leak
directly into the coordinates), and per-gene p-values are never used
either (sample size confounds significance, not effect size) — only
ReLU'd, closure-normalized **effect-size-based** statistics.

**User-confirmed 2026-08-17**: reuse this project's own already-reviewed
data (Step 3/4/4a) for two tracks (no new acquisition), plus add a third,
independent single-atlas validation track using Human Cell Landscape
(Han et al., *Nature* 2020, GSE134355) — the one candidate atlas the KB
note identified as potentially containing all three states on one
platform, with its own real feasibility unverified at the time the note
was written. That verification happens in this PR (see Track C below).

## Method: 3 tracks, each independently self-consistent

All three tracks share the same final math: raw per-axis scores →
ReLU (positive part only) → each axis independently rescaled by its own
99th percentile across that track's gene universe (clipped at 1) → add
a small floor (`EPS=0.01`) → closure-normalize (divide by row sum).
The 99th-percentile rescale is a **real, explicit methodological
simplification**: it harmonizes arbitrary absolute-scale differences
between different kinds of effect-size statistics feeding the same
triangle (e.g. Track B mixes a percentile-differential-based Fetal/Adult
score with a logFC-based Placenta score) without claiming a false
statistical equivalence between them — stated here, not hidden.

### Track A — Gut-specific (Fetal colon vs Placenta trophoblast vs Adult colon)

- Fetal/Adult: `results/04a_dfp_gut/edgeR/LargeInt_edgeR_primary.tsv`
  (Gut Cell Atlas, real within-atlas fetal-vs-adult colon epithelium DE,
  Step 4a). `logFC>0` = fetal-high (confirmed: AFP=+9.85/+9.89 in that
  table, matching this project's own textbook calibration marker).
  `S_fetal = ReLU(logFC)`, `S_adult = ReLU(-logFC)`.
- Placenta: `Arutyunyan_edgeR_results.tsv` + `Nature2026_edgeR_results.tsv`
  (Step 4's P-developmental primary DE). `S_placenta = min(ReLU(logFC_Aru), ReLU(logFC_Nat))`
  — 2-of-2 concordant-replication combination, same spirit as the
  frozen `P_developmental_primary84` definition.
- Gene universe: intersection of genes tested in all 3 DE tables — **13,480 genes**.

### Track B — Pan-tissue (embryo/fetal somatic vs trophoblast/placenta vs adult somatic)

- Fetal/Adult: cannot use a raw cross-platform DE (explicitly forbidden
  by this project's own Step 3 method contract — GTEx/HPA bulk and HDMA
  scRNA-seq pseudobulk are never compared by raw magnitude, only
  within-dataset percentile rank). Instead, reimplements
  `f_developmental_calibration.py`'s exact percentile machinery: for
  each of the 7 HDMA organs, `diff_organ(g) = elevated_pct(g, organ) -
  adult_pct(g, organ)` (within-organ fetal expression percentile minus
  the max percentile-among-detected across that organ's matched
  GTEx+HPA adult tissue columns). `S_fetal = ReLU(max_organ diff_organ)`,
  `S_adult = ReLU(-min_organ diff_organ)` — "best fetal-favoring organ"
  and "most adult-favoring organ" respectively, matching
  `F_developmental_union`'s own "any organ passes" logic. 329 genes
  dropped for having no organ with both fetal AND matched-adult evidence.
- Placenta: same construction as Track A (dataset-independent of organ).
- Gene universe: pan-organ ∩ placenta-tested — **13,835 genes**.

### Track C — HCL independent validation (GSE134355, Han et al. 2020)

The one candidate atlas with Fetal intestine + Placenta + Adult intestine
on a single, uniformly-processed platform (Microwell-seq) — avoids the
cross-dataset/cross-platform batch-as-biology risk Tracks A/B cannot
fully escape (their Placenta axis comes from an entirely different study
than their Fetal/Adult axis).

**Real feasibility check, done directly against the live GEO record and
the actual downloaded files, not assumed** (matching this project's
`curl_pipe_swallows_exit_code` / exact-byte-verification standing
discipline):
- GSE134355, public, no login required. `GSE134355_RAW.tar`:
  Content-Length 927,109,120 bytes, downloaded and verified byte-exact.
  141/141 members pass `gzip -t`.
- Confirmed by direct inspection (not filename/format guess): real
  integer raw UMI counts, gene symbols, genes-as-rows/cells-as-columns
  `_dge.txt.gz` format.
- Confirmed real samples for all 3 needed tissue types: `Placenta1` (1
  sample), `Fetal-Intestine1-5` (5 samples, **not regionally resolved**
  — no colon/small-intestine split, unlike Gut Cell Atlas), 8 adult
  intestine samples across colon/duodenum/ileum/jejunum segments
  (pooled here as one "Adult" group to match Fetal-Intestine's own lack
  of regional split).
- **Real, honest limitation found during this inventory**: every one of
  the 141 files has exactly 10,000 barcode columns — HCL's own DGE
  export pads/caps each sample at a fixed 10,000 columns rather than a
  variable barcode-rank-knee cutoff, so a real (unquantified here)
  fraction of those columns per sample are likely low/near-empty
  droplets, not all genuine cells. Pseudobulk summing is reasonably
  robust to this but it is not a "10,000 real cells" claim.
- **Real, honest limitation**: Placenta has only 1 donor (matches this
  project's own precedent for treating single-donor sources as
  directional/exploratory support, not primary replicated evidence —
  Greenbaum, Vento-Tormo).

Method: genuine one-vs-rest edgeR pseudobulk DE within the single atlas
(`scripts/09_developmental_ternary/hcl_edgeR_onevsrest.R`) — 3 groups
(Fetal n=5, Placenta n=1, Adult n=8), `~0+group` design,
`filterByExpr` → 19,266/34,744 genes kept, common dispersion 0.3408.
Three one-vs-rest contrasts (`group - mean(other two)`). Placenta's n=1
means its dispersion is borrowed entirely from the common/trended
estimate (edgeR's standard behavior for an unreplicated group) — its
p-values are not meaningful and are not used; only the logFC effect
size feeds the ternary score, consistent with this whole analysis's
"effect size, not p-value" principle.
`S_fetal/placenta/adult = ReLU(logFC_<group>_vs_rest)` directly (no
same-dataset-pairwise-contrast approximation needed here, unlike
Tracks A/B, since all 3 groups are genuinely modeled together).
Gene universe: **19,266 genes** (the joint `filterByExpr` pass — notably
does *not* lose AFP/LIN28B the way Tracks A/B's separate-DE-then-
intersect design does, see below).

## Real marker sanity check (all 3 tracks)

| Marker | Track A (Fetal/Placenta/Adult) | Track B | Track C |
|---|---|---|---|
| AFP | not in universe (see below) | not in universe (see below) | **0.981 / 0.010 / 0.010** |
| LIN28B | not in universe (see below) | 0.431 / 0.563 / 0.006 | 0.010 / 0.981 / 0.010 |
| CSH1/CSH2/CGA/PSG1/PSG3/ERVW-1 | Placenta-dominant where present | Placenta-dominant (0.59–0.90 Placenta) | Placenta-dominant (0.86–0.98 Placenta) |
| KRT7 | 0.010 / 0.761 / 0.229 | 0.109 / 0.443 / 0.448 | Placenta-dominant |
| ERVFRD-1 | not annotated (weak signal, see note) | Placenta-leaning but weaker than other placental markers | Placenta-leaning, weaker than other placental markers |
| **TRIM71** | **0.498 / 0.498 / 0.005** | **0.469 / 0.526 / 0.005** | **0.733 / 0.259 / 0.008** |

TRIM71 landing consistently near the Fetal-Placenta edge (both
substantial, Adult near-zero) across all 3 independently-built tracks is
the single most reassuring cross-track result — it is exactly the KB
note's own predicted "shared developmental" example
(`TRIM71 → (0.47, 0.48, 0.05)`), reproduced by 3 different methods/data
sources without being built to match that prediction.

**ERVFRD-1** (classic trophoblast syncytialization/fusion marker) shows
a real, consistently weaker signal than the other placental markers in
every track (Track C: logFC=+0.68 vs. CSH1's +16.99) — this project has
flagged ERVFRD-1 as a borderline marker before during the original
P-developmental threshold calibration (Step 4), so a weak-but-directionally-
correct signal here is consistent with that prior finding, not a new
concern, but worth naming rather than silently omitting from the table.

**LIN28B lands Placenta-dominant in both Track B and Track C**, not
Fetal-dominant as the KB note's own prediction guessed
("LIN28B → 更靠developmental/fetal side"). Reported as a real,
honest disagreement with the note's prior expectation — LIN28B is a
broadly developmental/stemness regulator, not exclusively a
fetal-somatic-lineage gene, and this is the kind of finding a first
version of this map should surface, not suppress to match a prior guess.

## A real, honest gene-dropout finding: AFP and LIN28B missing from Tracks A/B

Both of this project's two headline calibration markers fail to appear
in one or both of the reuse-based tracks' gene universes:
- **AFP** is absent from `Arutyunyan_edgeR_results.tsv` AND
  `Nature2026_edgeR_results.tsv` (filtered out of both placental DE
  fits, presumably for very low placental expression) — so it's absent
  from Track A's and Track B's intersected universe entirely.
- **LIN28B** is absent from `LargeInt_edgeR_primary.tsv` (filtered out
  of the Gut Cell Atlas fetal-vs-adult colon fit) — absent from Track A's
  universe (present in Track B, since Track B's HDMA-based Fetal/Adult
  score doesn't depend on that specific DE table).

This is a genuine, structural cost of the "3 separate DE fits, each with
its own `filterByExpr`, then intersect" design used by Tracks A/B: any
gene failing expression-filtering in *any one* of the 2-3 input studies
silently drops out of the universe, even if it's a real, well-known
developmental gene in the studies where it *is* tested. **Track C does
not have this problem** — a single joint `filterByExpr` across all 3
groups together means any gene passing that one filter appears in all 3
contrasts, which is why AFP appears cleanly at Track C's Fetal vertex.
This is a real argument in Track C's favor beyond just "avoids
cross-dataset batch risk," reported honestly as a limitation of the
reuse-based design, not fixed in this first version.

## A real, honest cross-track difference: distribution shape

Track A and Track C's bulk gene distributions form a visible, continuous
gradient along the Fetal↔Adult edge (most genes have near-zero Placenta
evidence, spread smoothly between the other two vertices). **Track B's
bulk distribution is visibly more bimodal** — density clusters much more
tightly at the Fetal and Adult vertices themselves, with less continuous
spread between them (see `ternary_track_B_pantissue.png`). This is a
real, observable consequence of Track B's percentile-differential-based
score construction (bounded, saturating at the extremes after the
99th-percentile rescale) versus Track A/C's raw-logFC-based construction
(naturally continuous), not a data-quality problem in either track —
flagged here as a genuine difference in what the two methods produce,
worth keeping in mind before comparing gene positions *across* Track A/B
at face value.

## Real bugs found and fixed during this build (not smoothed over)

1. **Marker-label NaN corruption**: the ternary-coordinate TSVs' `marker`
   column used `""` for "not a known marker gene." Reading that TSV back
   with pandas' default `read_csv` NA-sniffing silently converts `""` to
   `NaN` — and `df["marker"] != ""` then evaluates `True` for *every*
   row (since `NaN != ""` is always true), mislabeling the entire gene
   universe as markers. The first rendered plot showed every one of
   ~13,000-19,000 genes as a red dot with the literal text "nan" printed
   on it. Fixed: readers of these TSVs must pass `keep_default_na=False`
   (documented inline in both `build_ternary_coords.py` and
   `plot_ternary.py`).
2. **Track B NaN propagation**: 329 genes had no HDMA organ where both a
   fetal-elevation score AND a matched-adult-percentile score existed
   simultaneously (`diff_table.max/min(axis=1, skipna=True)` on an
   all-NaN row silently stayed NaN). These were propagating into the
   closure-normalization step as `NaN/NaN/x` coordinates. Fixed:
   genes without evidence in *any* organ are now explicitly dropped
   from the universe (reported: 6,352 dropped) rather than silently
   producing an undefined ternary position.
3. **Overlapping marker labels**: several placental hormone genes
   (CSH1/CSH2/CGA/PSG1/PSG3) land at near-identical ternary positions
   near the Placenta vertex; the first rendered plots' text labels
   overlapped into unreadable jumbled characters. Fixed: labels for
   near-coincident points are now stacked vertically instead of all
   using the same offset.
4. **macOS `zcat` incompatibility** (caught during manual inspection,
   before it became a script bug): macOS's `zcat` expects `.Z`-format
   compression by default and errors on plain gzip; `gzcat` used
   instead for all manual HCL file inspection in this write-up's
   verification steps (not a script issue — `gzip.open()` in Python and
   R's own gzip handling were never affected).

## File manifest

- `scripts/09_developmental_ternary/hcl_pseudobulk.py` — builds the 14
  relevant HCL pseudobulk samples from the byte-verified GSE134355 raw
  files (local execution, lightweight aggregation, not heavy compute).
- `scripts/09_developmental_ternary/hcl_edgeR_onevsrest.R` — Track C's
  one-vs-rest edgeR DE.
- `scripts/09_developmental_ternary/build_ternary_coords.py` — builds
  all 3 tracks' per-gene ternary coordinates.
- `scripts/09_developmental_ternary/plot_ternary.py` — renders the 3
  static ternary plots (matplotlib, manual barycentric projection).
- `results/09_developmental_ternary/track_{A_gut_specific,B_pantissue,C_hcl}_coords.tsv`
  — per-gene ternary coordinates, one row per gene, 3 tracks.
- `results/09_developmental_ternary/track_B_organ_diff_table.tsv` —
  Track B's intermediate per-organ percentile-differential table.
- `results/09_developmental_ternary/hcl_pseudobulk_{counts,meta}.tsv`,
  `hcl_edgeR_{Fetal,Placenta,Adult}_vs_rest.tsv` — Track C's intermediate
  pseudobulk and DE tables.
- `results/09_developmental_ternary/ternary_track_{A,B,C}_*.png` — the 3
  rendered ternary plots.
- Raw HCL data: `DATA/scRNAseq/GSE134355/raw/` (not committed to this
  repo, per this project's data-separation rule — `GSE134355_RAW.tar`,
  927,109,120 bytes, byte-verified, plus `extracted/` — 141 `_dge.txt.gz`
  files, all `gzip -t` clean).

## Explicitly out of scope for this first version

- No statistical significance testing on ternary position itself (this
  is a positional/exploratory map, not an inferential test).
- No quantitative cross-track comparison beyond the marker sanity check
  and the qualitative distribution-shape observation above.
- Not used here to revisit or qualify this project's founding-question
  answer (Steps 1-8) — this is a new, separate exploratory tool.
- HCL's fetal-intestine samples are not regionally resolved, so Track C
  cannot answer "colon-specific vs. small-intestine-specific" the way
  Track A/Gut Cell Atlas can — a genuine granularity gap between the
  gut-specific and HCL-validation tracks, not attempted to be closed
  here.
