# Raw-counts availability audit: consequence for the P-developmental quorum

`docs/STEP4_STATISTICAL_DESIGN.md` (PR #8, approved) named 3 "adequately-powered"
datasets for the primary trophoblast-vs-rest paired DE — Arutyunyan (n=17),
Nature2026 (n=23), VentoTormo (n=12) — with `replicated_in_placenta(gene)`
defined as passing a quorum (e.g. ≥2 of 3) across them. That count assumed all
three had raw counts available for the count-based edgeR/DESeq2 model. This
was not checked directly at design time. It has now been checked directly,
against the actual files, while building the pseudobulk matrices.

## Finding

| Dataset | Raw counts available? | Evidence |
|---|---|---|
| Arutyunyan (`primary_tissue/adata_all_donors_all_cell_states_UPD_20230307.h5ad`) | **Yes** — `.raw.X` | `X` is normalized (non-integer); `.raw.X` holds true integer counts. |
| Nature2026 (`scPlacenta_host.h5ad`) | **Yes** — via sibling file | `X` is normalized, no `.raw` slot. But `snRNA_raw_counts.h5ad` (same dataset directory) holds integer raw counts, and its `obs_names` match `scPlacenta_host.h5ad`'s 100% (191,735 / 191,735) — direct index alignment, no join needed. This also resolves Step 1 Finding #3 ("`snRNA_raw_counts.h5ad` has no usable annotation") as a side effect: the annotation is borrowed from `scPlacenta_host.obs` via the same index. |
| VentoTormo (`decidua-v3.h5ad`) | **No** | `X` is `normalize_total` + `log1p`'d and not reversible: `np.expm1(X)` is close to an integer for only ~1.4% of entries (checked directly, not assumed from the transform name). No `.raw` slot. No `layers`. No separate raw-counts file on disk — only `link.md` / `source_notes.md` pointing to external ArrayExpress accessions (E-MTAB-6701 / 6678 / 7304) that were never locally downloaded as raw data. |

Confirmed empirically (`np.allclose(X, np.round(X))` for integer-count check;
`np.expm1(X)` closeness-to-integer test for suspected log1p'd data), not
assumed from filenames or prior inventory notes.

## Consequence for the approved quorum design

`replicated_in_placenta(gene)` cannot currently use VentoTormo as one of its
required votes — there is no valid raw-count input for the primary
count-based (edgeR/DESeq2) model on this dataset. Effective state:

- **2 of 2** datasets with real raw counts (Arutyunyan, Nature2026) are
  usable for the primary paired edgeR QLF model.
- **VentoTormo is demoted** to the same role already defined for Greenbaum
  in `STEP4_STATISTICAL_DESIGN.md` §1: an optional directional-concordance
  booster (sign/rank agreement on its own normalized data), never one of the
  required primary votes, never sufficient alone.

## Options considered

1. **Accept 2-of-2 as the interim primary quorum** (both required, or a
   relaxed "at least 1 of 2 plus VentoTormo/Greenbaum directional support"
   rule) and proceed. No new downloads needed; consistent with the discipline
   already used for Greenbaum's n=3 exclusion.
2. **Pursue fresh VentoTormo raw counts from ArrayExpress** (E-MTAB-6701 /
   6678 / 7304) to restore a 3rd primary vote. Adds a new download +
   harmonization cycle, out of scope for the current threshold-audit task.
3. Some hybrid (e.g., revisit after signature freeze if 2-of-2 proves too
   permissive/restrictive in practice).

This audit documents the finding and its immediate consequence (VentoTormo
demoted, matching Greenbaum's existing secondary role) but does not itself
pick between options 1–3 for the long run — that is a design decision for
`STEP4_STATISTICAL_DESIGN.md` to make explicitly, flagged here for review.
**Recommendation: option 1** (2-of-2 primary quorum now), since it requires
no new data acquisition and keeps the project moving; option 2 can be
revisited later as an independent validation add-on if the 2-of-2 result set
looks like it needs a tie-breaker.
