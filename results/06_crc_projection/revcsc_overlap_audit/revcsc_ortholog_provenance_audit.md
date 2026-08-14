# revCSC mouse→human orthology provenance audit

Per PR #17 review round 1 (`REQUEST_CHANGES`): the pre-existing
`CSC_subtype_signatures.ensembl_mapping.tsv` table only proves that a
naive uppercased mouse symbol *matches* an HGNC-approved human symbol —
that is a **symbol-mapping success**, not a verified **mouse→human
orthology relationship**. Before revCSC could be promoted to the primary
Layer 1 Oncofetal anchor, each of its 32 distinct mouse gene symbols
(cluster=revCSC in that table) was checked directly against **Ensembl
Compara's own computed orthology calls** (`homology/symbol` REST endpoint),
not assumed from symbol capitalization alone.

**Source correction**: an earlier draft of `STEP6_CRC_PROJECTION_DESIGN.md`
mischaracterized revCSC's source paper as "mouse-intestinal-regeneration-
derived" (implicitly Ayyaz et al., Nature, 2019 — a normal-tissue
regeneration paper). Per direct user correction, the actual source is
*An oncogenic phenoscape of colonic stem cell polarization* — a CRC/colonic
stem-cell-polarization study. This project has not independently
re-verified the original mouse gene list against that paper's own
supplementary table; the 32 symbols audited here are what was already
extracted into the pre-existing `CSC_subtype_signatures.ensembl_mapping.tsv`
on Argos, taken as the input to this orthology audit (a documented scope
limitation, not silently assumed complete).

## Method

For each of the 32 distinct mouse gene symbols, queried
`https://rest.ensembl.org/homology/symbol/mus_musculus/<symbol>` (target
species human, orthologues only) — a network call to a public API, run
locally, same category as Step 2's `query_biotype.py` precedent (not
"analysis of local biological data," so not an Argos-only-compute item).
Genes returning `HTTP 500`/timeouts were retried (Ensembl's public REST
endpoint was intermittently overloaded during this session); genes with a
stable, repeated empty-homology response across independent retries were
recorded as `no_ortholog_found`, not left ambiguous.

## Result: 3 real corrections + 1 wrong-gene fix, found by not trusting the naive mapping

| Outcome | n genes | Genes |
|---|---|---|
| `ortholog_one2one`, confirmed, unambiguous | 27 | (see `revCSC_human_FINAL.tsv`) |
| `ortholog_one2many`, ambiguous but resolvable to the pre-validated target | 1 | `Ly6a` → `LY6A`/`ENSG00000291309` |
| **No Compara-confirmed ortholog — excluded** (previously included via naive symbol match) | 3 | `Cldn4`, `Ctsl`, `Sprr1a` |
| No ortholog at all (already known/expected) | 1 | `Ctla2a` (also failed the original symbol-mapping step) |

**One outright wrong Ensembl ID found and corrected**: the pre-existing
table mapped `Ccn1` → `ENSG00000145386`, which is actually **`CCNA2`
(Cyclin A2)** — a completely unrelated cell-cycle gene, confirmed by direct
`lookup/id` query (`display_name: CCNA2`, `description: "cyclin A2"`).
Ensembl Compara's own `ortholog_one2one` call for `Ccn1` is
`ENSG00000142871` (confirmed by direct lookup: `display_name: CCN1`,
`description: "cellular communication network factor 1"`) — used instead.
Had this gone uncorrected, the revCSC score would have silently included
`CCNA2`'s expression as if it were part of the Oncofetal signature.

**Three genes dropped that symbol-matching had silently kept**: `Cldn4`,
`Ctsl`, `Sprr1a` all have a valid-looking HGNC-approved human symbol match
(`CLDN4`, `CTSL`, `SPRR1A`), but Ensembl Compara reports **zero** computed
mouse→human orthologs for any of the three (empty `homologies` list,
confirmed stable across 2–3 independent re-queries each, not a one-off
API hiccup). Plausible explanation for at least two of these (not
asserted as proven): `CTSL`/`SPRR1A` both belong to gene families with
known lineage-specific paralog expansion in mouse and human (cathepsins,
small-proline-rich proteins), which can defeat Compara's ortholog-calling
confidence even when a biological relationship plausibly exists — but
without a confirmed 1:1 call, none of the three are included in a
signature meant to stand as a rigorous Oncofetal anchor.

**Net effect on the frozen scoring set**: revCSC shrinks from the
previously-reported "31 mapped genes" to a corrected, ortholog-verified
**28-gene primary set** (27 confirmed `one2one` + `Ly6a` flagged
`one2many`), with a **27-gene strict `one2one`-only sensitivity set**
(dropping `Ly6a`) reported alongside every primary result per the
pre-registered inclusion rule below.

## Pre-registered inclusion rule (decided here, before any scoring compute — not post-hoc)

- **Primary revCSC score**: the 28-gene set (27 `one2one` + `Ly6a`).
- **Sensitivity variant, reported alongside every primary result, not
  optionally**: the 27-gene strict `one2one`-only set (`Ly6a` dropped).
  Any material difference between the two is itself worth reporting, not
  discarded.
- Genes with no Compara-confirmed ortholog (`Ctla2a`, `Cldn4`, `Ctsl`,
  `Sprr1a`) are excluded from both variants.
- All future scoring, overlap-audit, or write-up steps reference
  `revCSC_human_FINAL.tsv` as the single frozen artifact — not the
  historical `CSC_subtype_signatures.ensembl_mapping.tsv` directly.

## No re-run needed for the D/F/P gene-overlap audit

Checked directly: none of the 4 excluded genes (`Cldn4`, `Ctsl`, `Sprr1a`,
`Ctla2a`) or the corrected `Ccn1`/`Ly6a` appear in D-shared, F-specific
(global or any of the 7 lineage modules), or P-specific. The previously
reported overlap result (D=0, F-specific global=2 [`ACTA1`, `ANKRD1`],
P=0) is unaffected by this correction and does not need to be re-run.

## Files

- `revCSC_human_FINAL.tsv` — the frozen artifact: one row per mouse gene,
  with human ortholog symbol/Ensembl ID, orthology type, inclusion
  decision, and an explicit note for every corrected/excluded gene.
- `revcsc_mouse_human_ortholog_audit.tsv` — raw Ensembl Compara query
  results (pre-consolidation), kept for traceability.
- `revCSC_symbols.primary28.txt` / `revCSC_symbols.strict_one2one27.txt` —
  the two scoring-set gene lists per the inclusion rule above.
