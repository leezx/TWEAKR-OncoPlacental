# Dataset: Tabula Sapiens adult single-cell atlas (organ-matched subset)

- Source: Tabula Sapiens Consortium, Figshare (DOI 10.6084/m9.figshare.14267219)
- Role in project: Tier-2 adult reference — deliberately held out of Step 3's D/F/P signature *definition*, used only as independent post-freeze validation ("do P-specific genes stay absent across real adult cell types, not just bulk tissue averages").

## Scope decision

Downloaded only 5 of the ~24 per-organ files: `TS_{Liver,Skin,Spleen,Thymus,Large_Intestine}.h5ad.zip` (~2.85GB) — the 4 that match HDMA organs plus Large_Intestine for the CRC-relevant adult colon check. Skipped the full unified atlas (15.6GB) and the compartment-level files (epithelial/stromal/immune/endothelial, ~15GB) as not needed for this project's scope.

**Confirmed gap**: Tabula Sapiens has no Adrenal, Thyroid, or Stomach/Esophagus organ file — checked against the full Figshare file listing, not assumed. GTEx/HPA (Tier 1) remain the only adult reference for those 3 HDMA organs; no cell-type-level validation available for them from this source.

## What's on disk

`DATA/1.Databases/TabulaSapiens/raw/TS_{Liver,Skin,Spleen,Thymus,Large_Intestine}.h5ad.zip`

## A download bug worth flagging for anyone reusing this pattern

Figshare's `download_url` 302-redirects to a **pre-signed S3 URL that expires in 10 seconds**. `curl -o file url` without `-L` "succeeds" (exit 0) while writing a 0-byte file — it just saves the empty redirect response instead of following it. Not a loud failure; caught only because every download here is checked against Figshare's own `size`/`supplied_md5`, not by trusting curl's exit code. Fix: always use `curl -sSL`.

## Verification

Every file's size and MD5 checked against Figshare's own API-supplied values before being marked OK.

## Full record

`DATA/1.Databases/TabulaSapiens/link.md`
