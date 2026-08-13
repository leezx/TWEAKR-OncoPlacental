# Step 3 method contract: cross-platform comparison rules

Locked down 2026-08-13, before any D-shared/F-specific/P-specific (D/F/P) signature construction starts, per ChatGPT reviewer feedback on PR #5 (round 1): the acquired reference datasets sit on genuinely different platforms and normalizations, and the comparison method has to be fixed explicitly rather than left implicit — otherwise "F-specific"/"P-specific" would silently absorb platform, sequencing-depth, cell-composition, and normalization differences instead of real biology.

## The platform problem

| Dataset | Platform | Unit |
|---|---|---|
| HDMA (fetal-somatic), placental scRNA-seq (trophoblast) | Single-cell / single-nucleus RNA-seq, Seurat counts | Raw/normalized UMI counts per cell |
| `GTEx_v11_median_tpm` | Bulk RNA-seq | Median TPM per tissue |
| `HPA_RNA_tissue_consensus` | Bulk RNA-seq (HPA-generated + consensus with other sources) | nTPM per tissue |

These three are **not comparable as raw magnitudes** — a gene's TPM in a GTEx tissue and its mean UMI count in an HDMA scRNA cluster are on different scales for reasons that have nothing to do with biology (bulk vs. single-cell capture efficiency, sequencing depth, cell-type composition of the bulk sample, normalization method). Directly computing fold-change or a merged DE model across them would contaminate the signature with these platform effects.

## The rule

**Developmental evidence and adult-exclusion evidence are computed as two separate, independent evidence axes — never merged into one cross-platform DE model.**

1. **Developmental-side evidence (F-specific / P-specific candidate identification)**: computed entirely *within* the scRNA-seq data — e.g. trophoblast vs. other placental cell types (all within the same placental scRNA-seq dataset(s)), or fetal organ vs. other fetal organs (all within HDMA). No GTEx/HPA values enter this step.

2. **Adult-exclusion evidence (the actual role of GTEx/HPA in this project)**: GTEx/HPA are used only to answer a binary/ranked question **within their own dataset** — "is this candidate gene still meaningfully expressed in the matched adult organ / across the adult body?" — using each dataset's own internal rank, percentile, or a pre-defined expression threshold computed from that dataset alone (e.g. "top X% of genes by TPM in this tissue" or "TPM > threshold Y"). GTEx values are never compared in raw magnitude against HPA values, and neither is ever compared in raw magnitude against HDMA/trophoblast expression.

3. **Combining the two**: a gene qualifies as F-specific or P-specific only if it satisfies *both* axes — real developmental evidence (from the scRNA-seq-internal comparison) **and** adult-depletion evidence (from the GTEx/HPA-internal rank/percentile/threshold check) — combined as a logical AND over two independently-computed evidence types, not as a single merged statistical model.

## Placenta in HPA — explicit exclusion from the adult negative reference

`HPA_RNA_tissue_consensus` includes a `placenta` row (documented in `link.md` as a "bonus... cross-check" row) alongside the 39 other genuinely-adult tissue rows. This was previously described inconsistently — the dataset doc called it "40 adult tissues" without carving placenta out, which would have silently let placenta count as part of the adult-negative background for P-specific gene calls (self-contradictory: using placenta expression as evidence *against* a gene being placenta-specific).

**Explicit rule**: `placenta` is excluded from every adult-exclusion-reference computation (GTEx has no placenta row at all, so this only affects HPA). It may only be used as a **positive cross-check** — e.g. confirming a P-specific candidate is indeed high in HPA's placenta row, as an independent sanity check against the project's own placental scRNA-seq calls, never as a negative/background tissue.

## What this contract does NOT do

It does not start Step 3's actual signature construction — no candidate gene lists, no thresholds chosen, no code implementing the above written yet. It only fixes the comparison boundary so that when Step 3 code is written, it can't accidentally reach for a cross-platform fold-change as a shortcut.
