# Replicate-structure audit: which placental datasets can run the trophoblast-vs-rest DE

Resolves the first open item in `docs/STEP4_DFP_DESIGN.md`: "whether all 5 placental datasets have enough donor/sample-level replicate structure for within-dataset DE, or whether some need to be pseudobulk-per-cluster instead — needs a real inventory check, not an assumption."

**Final numbers (see "Round 2" below for the real donor×trophoblast-status cross-tab — read this, not the round-1 table's donor counts, which undercount how thin Greenbaum's replication actually is)**: Arutyunyan primary_tissue 17/18 donors eligible, Nature2026 scPlacenta_host 23/23, VentoTormo decidua-v3 12/12, **Greenbaum 3/3 — not 8** as round 1's marginal-count table below implies (round 1 counted all 8 donors in the full `metadata.csv`; the annotated subset actually usable for the DE only covers 3 of them). All 4 datasets remain usable, but treat Greenbaum's evidence as much weaker than the other three (n=3 vs. 12–23).

**Method**: re-read the per-dataset `obs_columns` / `obs_value_counts` (and Greenbaum's `metadata_value_counts` / `cluster_value_counts`) already captured in `results/01_inventory/*.json` (Step 1, verified via qsub on Argos) — no new compute needed, this is a metadata lookup against already-produced, already-verified inventory output.

## Result: 4 of 7 placental datasets can run the DE; 3 structurally cannot (not a data-quality problem)

| Dataset | Donor/sample column | # donors/samples | Trophoblast cells present? | Can run trophoblast-vs-rest DE? |
|---|---|---|---|---|
| Arutyunyan primary_tissue | `donor` | 18 (1,249–66,038 cells/donor) | Yes — `coarse_annot == Trophoblast`, 75,042 cells vs. 179,654 non-trophoblast (dS_uSMC/Myeloid/NK/T/Epithelial/fF/Endothelial) | **Yes** |
| Nature2026 scPlacenta_host | `sample_id` | 23 (1,882–19,449 cells/sample) | Yes — `major_class` SCT+VCT+EVT = 95,872 cells vs. 88,229 non-trophoblast (Epi/DSC/dNK/M/EC/FB/HB/PV/T/cDC/B/Ery) | **Yes** |
| VentoTormo decidua-v3 | `Fetus` | 12 (634–17,591 cells/fetus; **note**: some values have leading whitespace, e.g. `" F15"` vs `"F19"` — needs trimming before use, real data-cleanliness issue, not a missing-data one) | Yes — VCT+EVT+SCT = 14,366 cells vs. 50,510 non-trophoblast (dS/dNK/dM/Tcells/HB/fFB/etc.) | **Yes** |
| Greenbaum NatMed 2024 | `donor_id` | ~~8~~ **superseded by Round 2 below — actually 3** (3,108–6,889 cells/donor is the full `metadata.csv`'s 8 donors, but the cluster-level annotation subset that actually carries `cell_type`, ~1,923 cells, only covers 3 of them — see Round 2) | Yes — vCTB+STB+EVT+EVT-progenitor+STB-progenitor = 1,038 cells vs. 885 non-trophoblast in the cluster file | **Yes, but on the smaller clustered subset (3 donors, not 8 — see Round 2)** |
| Arutyunyan organoid PTO | `donor` present | — | **No** — organoid culture, 100% trophoblast-lineage by construction (Step 1: "Pure trophoblast by construction... all values are VCT/EVT/SCT-lineage subtypes, no filtering needed") | **No — structurally excluded**, not a replicate-structure problem: there is no non-trophoblast population in the dataset to contrast against |
| Arutyunyan organoid TSC | `donor` present | — | No, same reason | **No — structurally excluded** |
| Arutyunyan organoid Fig3 | `donor` present | — | No, same reason | **No — structurally excluded** |

## What this settles

- **4 usable datasets** for the trophoblast-vs-rest within-dataset DE that feeds `P-developmental`: Arutyunyan primary_tissue, Nature2026 scPlacenta_host, VentoTormo decidua-v3, Greenbaum. All 4 have real donor/sample-level replicate structure (8–23 donors each) — pseudobulk-per-donor DE (e.g. Wilcoxon on donor-level pseudobulk trophoblast vs. non-trophoblast values) is valid for all four, not just cluster-level pseudo-replication.
- **3 organoid datasets structurally cannot contribute** to this specific DE, regardless of how much donor replication they have — there's no non-trophoblast population inside them to compare against. This is expected given what they are (trophoblast organoid cultures), not a data problem. They may still be useful elsewhere (e.g. as an independent trophoblast-positive expression reference, outside the DE test itself) but that's a separate question from this audit.
- Revises `STEP4_DFP_DESIGN.md`'s placeholder quorum language ("≥3 of 5") to **≥3 of 4**, since only 4 datasets are actually eligible to vote.

## Round 2 (PR #7 review): donor × trophoblast-status cross-tab, not just marginal totals

Reviewer correctly caught that round 1 only checked marginal totals ("dataset has N donors" and "dataset has trophoblast + non-trophoblast cells overall") — not whether any single donor actually has both. If trophoblast status is confounded with donor identity (e.g. some donors contribute only trophoblast cells, others only non-trophoblast), a donor-level pseudobulk contrast is invalid regardless of how many nominal "donors"/"replicates" exist.

**Method**: `scripts/04_dfp_signature/donor_troph_crosstab.py`, run on Argos (backed-mode `obs` reads for the 3 h5ad files, direct CSV read for Greenbaum) — this is genuinely new compute (the marginal counts in round 1 came from Step 1's existing inventory JSONs, but the joint donor × cell-type breakdown wasn't captured there and had to be computed fresh). Full log: `results/04_dfp_signature/donor_troph_crosstab.txt`.

| Dataset | Donors with both groups | Total donors | Notes |
|---|---|---|---|
| Arutyunyan primary_tissue | 17 | 18 | Only `R1` (2,309 cells) has zero trophoblast cells — excluded from the paired contrast, doesn't affect the other 17 |
| Nature2026 scPlacenta_host | 23 | 23 | Clean — every sample has both groups |
| VentoTormo decidua-v3 | 12 | 12 | Clean **after** trimming the `Fetus` whitespace bug found in round 1 — confirms that fix was necessary (untrimmed, `" F15"` and a hypothetical clean `"F15"` would have silently split into two fake donors) |
| Greenbaum | 3 | 3 (**not 8**) | Found and fixed a second real bug in the round-1 script: joining `cluster.csv`'s `NAME` to `metadata.csv`'s `NAME` fails 100% of the time — the two files use incompatible barcode schemes (`cluster.csv`: `W9_AAACCAACACCTGCCT`; `metadata.csv`: `JS34#ACGTCAAGTTGCAATG-1`, a different sample-ID system). Fixed: `cluster.csv`'s own `NAME` already embeds the donor ID as the prefix before the last `_` — parsed directly, no join needed. But this reveals the annotated ~1,923-cell subset only covers **3 of the full metadata's 8 donors** (W8-2, W9, W11), not all 8 — a real, smaller-than-hoped replicate count for this one dataset, not a bug to fix further |

**Revised conclusion**: all 4 datasets are usable, but with real (not assumed) donor counts for the actual paired contrast — 17, 23, 12, and 3 respectively, not 18/23/12/8. Greenbaum's n=3 is thin; it can still vote in a `replicated_in_placenta` quorum but is the weakest-powered of the four and its result should be weighted accordingly when Step 4 actually runs the DE, not treated as equally strong evidence to the other three.

## Not yet resolved by this audit

- Nature2026's `snRNA_raw_counts` file remains excluded (Step 1 Finding #3 — no usable annotation), so it was never a candidate for this DE regardless of replicate structure.
- The exact statistical test and threshold for "elevated in trophoblast" within each dataset (Wilcoxon on donor-pseudobulk vs. some other test, effect-size cutoff) — this audit only confirms *that* a valid test is possible on 4 datasets, not *which* test/threshold to use, or how to weight Greenbaum's much smaller n=3. Left for the next Step 4 sub-task.
