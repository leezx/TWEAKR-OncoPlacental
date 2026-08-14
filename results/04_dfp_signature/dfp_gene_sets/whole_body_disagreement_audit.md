# GTEx/HPA whole-body disagreement audit

Per the PR #13 reviewer's request: the initial assembly used GTEx as an
unjustified "primary" whole-body reference (158 genes) with HPA (160
genes) only reported alongside — real platform sensitivity that needed
explaining before freezing, not glossing over.

## First finding: most of the apparent "50% disagreement" was a coverage-gap artifact, not real disagreement

A first pass at bucketing `replicated_in_placenta`'s 1,007 genes by
GTEx/HPA whole-body pass/fail repeated the **exact same bug already fixed
for organ-matched exclusion in PR #12** (the H19 case): genes not in a
platform's gene panel at all were being counted as that platform
"failing" them, rather than "no evidence." Caught this myself before
reporting it — same coverage-aware fix applied here.

**Real coverage breakdown**: of the 1,007 genes, only **815** are actually
measured by *both* GTEx and HPA. 77 are GTEx-only-covered (not in HPA's
panel at all), 1 is HPA-only-covered, and 114 are covered by **neither**
platform (present in the scRNA-seq datasets but absent from both bulk
references — a separate, real gap, not investigated further here).

**Of the 815 genes both platforms actually measure — the only valid
disagreement population** — the real numbers are much smaller than the
raw pass-count comparison suggested:

| | Count |
|---|---|
| Both pass | 84 |
| GTEx-pass / HPA-fail (real disagreement) | 15 |
| HPA-pass / GTEx-fail (real disagreement) | 75 |
| Both fail | 641 |

Genuine disagreement is 90/815 = **11%**, not the ~50% the raw 158-vs-160
numbers implied.

## What's driving the two disagreement directions

**GTEx-pass/HPA-fail (15 genes)**: `CGA`, `CSH1`, `CSH2`, `ADCY10`,
`ATP1A4`, `DPPA4`, `FAM72D`, `LRRC8E`, `LRTM1`, `NOXO1`, `P2RY6`, `RFESD`,
`RGPD2`, `SLC6A2`, `SOHLH2`. Checking which HPA tissues they fail in:
**testis appears repeatedly** (ATP1A4, DPPA4, RGPD2, SLC6A2, SOHLH2,
RFESD all fail partly via testis) — testis is a well-documented outlier
in human transcriptomics for genome-wide "leaky" low-level transcription
(open chromatin during spermatogenesis), a known confound independent of
this project. CGA fails specifically in parathyroid gland and stomach —
CGA is the common glycoprotein-hormone alpha subunit shared with LH/FSH/TSH,
so real (if minor) non-placental expression is expected, genuine biology,
not noise. CSH1/CSH2 fail in scattered, biologically unrelated tissues
(endometrium, lung, smooth muscle, spleen, stomach) — no coherent story,
more consistent with detection noise than real cross-tissue expression.

**HPA-pass/GTEx-fail (75 genes)**: includes `HLA-G`, `FBN2`, `HSD3B1`,
`SIGLEC6` (in HPA's own official 65-gene set — real placental biology that
HPA's smaller panel happens to correctly retain) but the majority are
genes with clearly broad, non-placental adult biology that GTEx's larger
68-tissue panel catches and HPA's 39-tissue panel simply has less power to
detect: DNA-repair/meiosis genes (`BRCA2`, `POLQ`, `XRCC2`, `DMC1`,
`RDM1`, `TICRR`, `SYCP2L`), and genes failing in dozens of GTEx tissues
(`ANKRD24` fails in 27/68, `CDH24` in 57/68, `C2orf72` in 33/68) —
essentially ubiquitous expression, clearly not placenta-restricted.

**Interpretation**: neither direction of disagreement is obviously "the
other platform is wrong" across the board — some GTEx-only passes look
like real placental biology missed by HPA's narrower panel (CGA via
minor cross-reactivity) mixed with likely noise (CSH1/CSH2, testis-driven
cases); some HPA-only passes look like real broadly-expressed genes GTEx's
larger panel correctly caught (DNA-repair genes) mixed with real placental
biology HPA's official classification independently confirms (HLA-G,
FBN2, HSD3B1, SIGLEC6). **Not confident enough to write a bespoke
combination rule** (e.g. "ignore testis-driven HPA failures") without
further validation — following the reviewer's own stated fallback for
exactly this situation.

## Resolution: both-pass primary + explicit extended tier (not frozen blindly, not silently merged)

- **`P_developmental_primary84.txt`** (84 genes): the both-pass set — both
  platforms measure the gene and both agree it's whole-body adult-excluded.
  Highest confidence, 31% overlap with HPA's own official 65-gene placenta
  classification (vs. ~0.3% background rate — massive enrichment).
- **`P_developmental_extended174.txt`** (174 genes = 84 + 15 + 75): adds
  the genuine single-platform-disagreement genes as a lower-confidence
  extension, kept explicitly separate, not merged into the primary tier.

**Marker re-check against the primary 84-gene tier**: ERVFRD-1, PSG1, PSG3
remain (core placental genes, both platforms agree). CGA, CSH1, CSH2,
HLA-G move to the extended-only tier (real biology, but not both-platform
agreement) — GATA3/KRT7 remain excluded from both tiers (already
established as broadly-expressed, non-placenta-restricted).

## Revised final D/F/P assembly

| P-developmental | F-developmental | D-shared | F-specific | P-specific |
|---|---|---|---|---|
| Primary (84) | Union (2,510) | **6** | 2,504 | 78 |
| Primary (84) | Consensus (1,339) | 3 | 1,336 | 81 |
| Extended (174) | Union (2,510) | 24 | 2,486 | 150 |
| Extended (174) | Consensus (1,339) | 14 | 1,325 | 160 |

**Primary recommendation: P=primary(84), F=union(2,510)** → D-shared=6
(`IL1RAPL2`, `KCNH5`, `PCDH11X`, `TMC1`, `TRIM71`, `ZNF730` — TRIM71
survives even this stricter cut, still a plausible core shared gene),
F-specific=2,504, P-specific=78.

Full gene lists: `P_developmental_primary84.txt`,
`P_developmental_extended174.txt`, `D_shared_FINAL.txt`,
`F_specific_FINAL.txt`, `P_specific_FINAL.txt` (all P=primary84,
F=union combination); disagreement bucket files
`whole_body_{both_pass_84,gtex_only,hpa_only,both_fail}.txt`.
