# Replicate-structure audit: which placental datasets can run the trophoblast-vs-rest DE

Resolves the first open item in `docs/STEP4_DFP_DESIGN.md`: "whether all 5 placental datasets have enough donor/sample-level replicate structure for within-dataset DE, or whether some need to be pseudobulk-per-cluster instead — needs a real inventory check, not an assumption."

**Method**: re-read the per-dataset `obs_columns` / `obs_value_counts` (and Greenbaum's `metadata_value_counts` / `cluster_value_counts`) already captured in `results/01_inventory/*.json` (Step 1, verified via qsub on Argos) — no new compute needed, this is a metadata lookup against already-produced, already-verified inventory output.

## Result: 4 of 7 placental datasets can run the DE; 3 structurally cannot (not a data-quality problem)

| Dataset | Donor/sample column | # donors/samples | Trophoblast cells present? | Can run trophoblast-vs-rest DE? |
|---|---|---|---|---|
| Arutyunyan primary_tissue | `donor` | 18 (1,249–66,038 cells/donor) | Yes — `coarse_annot == Trophoblast`, 75,042 cells vs. 179,654 non-trophoblast (dS_uSMC/Myeloid/NK/T/Epithelial/fF/Endothelial) | **Yes** |
| Nature2026 scPlacenta_host | `sample_id` | 23 (1,882–19,449 cells/sample) | Yes — `major_class` SCT+VCT+EVT = 95,872 cells vs. 88,229 non-trophoblast (Epi/DSC/dNK/M/EC/FB/HB/PV/T/cDC/B/Ery) | **Yes** |
| VentoTormo decidua-v3 | `Fetus` | 12 (634–17,591 cells/fetus; **note**: some values have leading whitespace, e.g. `" F15"` vs `"F19"` — needs trimming before use, real data-cleanliness issue, not a missing-data one) | Yes — VCT+EVT+SCT = 14,366 cells vs. 50,510 non-trophoblast (dS/dNK/dM/Tcells/HB/fFB/etc.) | **Yes** |
| Greenbaum NatMed 2024 | `donor_id` | 8 (3,108–6,889 cells/donor in the full `metadata.csv`; **note**: cluster-level annotation file is a much smaller subsampled/clustered subset, ~1,923 cells total, not the full dataset) | Yes — vCTB+STB+EVT+EVT-progenitor+STB-progenitor = 1,038 cells vs. 885 non-trophoblast in the cluster file | **Yes, but on the smaller clustered subset** — the DE test's effective N is bounded by whichever cell set actually carries the `cell_type` annotation, not the full ~36K-cell RNA matrix |
| Arutyunyan organoid PTO | `donor` present | — | **No** — organoid culture, 100% trophoblast-lineage by construction (Step 1: "Pure trophoblast by construction... all values are VCT/EVT/SCT-lineage subtypes, no filtering needed") | **No — structurally excluded**, not a replicate-structure problem: there is no non-trophoblast population in the dataset to contrast against |
| Arutyunyan organoid TSC | `donor` present | — | No, same reason | **No — structurally excluded** |
| Arutyunyan organoid Fig3 | `donor` present | — | No, same reason | **No — structurally excluded** |

## What this settles

- **4 usable datasets** for the trophoblast-vs-rest within-dataset DE that feeds `P-developmental`: Arutyunyan primary_tissue, Nature2026 scPlacenta_host, VentoTormo decidua-v3, Greenbaum. All 4 have real donor/sample-level replicate structure (8–23 donors each) — pseudobulk-per-donor DE (e.g. Wilcoxon on donor-level pseudobulk trophoblast vs. non-trophoblast values) is valid for all four, not just cluster-level pseudo-replication.
- **3 organoid datasets structurally cannot contribute** to this specific DE, regardless of how much donor replication they have — there's no non-trophoblast population inside them to compare against. This is expected given what they are (trophoblast organoid cultures), not a data problem. They may still be useful elsewhere (e.g. as an independent trophoblast-positive expression reference, outside the DE test itself) but that's a separate question from this audit.
- Revises `STEP4_DFP_DESIGN.md`'s placeholder quorum language ("≥3 of 5") to **≥3 of 4**, since only 4 datasets are actually eligible to vote.

## Not yet resolved by this audit

- Nature2026's `snRNA_raw_counts` file remains excluded (Step 1 Finding #3 — no usable annotation), so it was never a candidate for this DE regardless of replicate structure.
- The exact statistical test and threshold for "elevated in trophoblast" within each dataset (Wilcoxon on donor-pseudobulk vs. some other test, effect-size cutoff) — this audit only confirms *that* a valid test is possible on 4 datasets, not *which* test/threshold to use. Left for the next Step 4 sub-task.
- Greenbaum's cluster-file-vs-full-matrix cell-count mismatch (1,923 vs. ~36,456 per Step 1's shape record) needs a closer look before running the actual DE — is the clustered subset representative, or does it need re-deriving from the full matrix with proper cell-type calls first?
