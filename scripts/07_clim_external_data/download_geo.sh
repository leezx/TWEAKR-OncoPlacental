#!/bin/bash
# Step 7 inventory: download the public scRNA-seq + bulk RNA-seq/microarray
# cohorts locked in docs/STEP7_CLIM_DATA_ACQUISITION_DESIGN.md (PR #35,
# APPROVE after 3 review rounds).
#
# Runs on Argos, argos-qsub1, plain network I/O (not qsub-wrapped --
# I/O-bound, same precedent as this project's original HDMA/Arutyunyan
# acquisition phase). Never piped through tail/head (standing
# curl_pipe_swallows_exit_code lesson) -- every download's exit code is
# checked directly, then independently re-verified against the exact
# integer byte count from the HTTP Content-Length header (never GEO's
# rounded MB/GB webpage display).
#
# Usage: bash download_geo.sh <DATA_ROOT>
# e.g.:  bash download_geo.sh /home/zz950/DATA

set -uo pipefail  # NOT -e: a single failed download must not abort the
                   # whole manifest; each is checked and reported individually

DATA_ROOT="${1:-/home/zz950/DATA}"
LOG="download_geo.log"
echo "=== Step 7 GEO download, started $(date -u +%FT%TZ) ===" | tee -a "$LOG"

FAILED=0

# download_and_verify <url> <dest_path>
# Downloads via curl directly to disk (no pipe), checks curl's own exit
# code, then verifies the on-disk byte count against the Content-Length
# header fetched independently via -I.
download_and_verify() {
    local url="$1" dest="$2"
    local expect
    expect=$(curl -sI "$url" | grep -i '^content-length:' | tail -1 | tr -d '\r' | awk '{print $2}')
    if [ -z "$expect" ]; then
        echo "[FAIL] $dest -- could not fetch Content-Length for $url" | tee -a "$LOG"
        FAILED=$((FAILED+1))
        return 1
    fi
    if [ -f "$dest" ]; then
        local have
        have=$(stat -c%s "$dest" 2>/dev/null || stat -f%z "$dest" 2>/dev/null)
        if [ "$have" = "$expect" ]; then
            echo "[SKIP] $dest already present, $have bytes matches Content-Length -- not re-downloading" | tee -a "$LOG"
            return 0
        fi
    fi
    curl -sS --fail -o "$dest" "$url"
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "[FAIL] $dest -- curl exit code $rc for $url" | tee -a "$LOG"
        FAILED=$((FAILED+1))
        return 1
    fi
    local have
    have=$(stat -c%s "$dest" 2>/dev/null || stat -f%z "$dest" 2>/dev/null)
    if [ "$have" != "$expect" ]; then
        echo "[FAIL] $dest -- byte mismatch: expected $expect, got $have" | tee -a "$LOG"
        FAILED=$((FAILED+1))
        return 1
    fi
    echo "[OK] $dest -- $have bytes, exact match to Content-Length" | tee -a "$LOG"
    return 0
}

GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo"

# ---------------------------------------------------------------------
# 1. GSE231559 -- scRNA-seq, 26 samples (paper cites 9 CLiM + 6 primary
#    subset of these 26) -- full RAW.tar, clean 10x MTX loader pattern.
# ---------------------------------------------------------------------
d="$DATA_ROOT/scRNAseq/GSE231559/raw"
mkdir -p "$d"
download_and_verify \
    "$GEO_BASE/series/GSE231nnn/GSE231559/suppl/GSE231559_RAW.tar" \
    "$d/GSE231559_RAW.tar"

# ---------------------------------------------------------------------
# 2. GSE225857 -- scRNA-seq. ONLY GSM7058754/GSM7058755 are genuine
#    scRNA-seq (pooled immune/non-immune across all patients); the other
#    6 GSMs in this series are spatial (out of scope) -- selected by
#    exact GSM accession per the locked design, not the series RAW.tar
#    or file-shape heuristics.
# ---------------------------------------------------------------------
d="$DATA_ROOT/scRNAseq/GSE225857/raw"
mkdir -p "$d"
for gsm_file in \
    "GSM7058754:GSM7058754_immune_counts.txt.gz" \
    "GSM7058754:GSM7058754_immune_meta.txt.gz" \
    "GSM7058755:GSM7058755_non_immune_counts.txt.gz" \
    "GSM7058755:GSM7058755_non_immune_meta.txt.gz"; do
    gsm="${gsm_file%%:*}"; fname="${gsm_file##*:}"
    gsm_prefix="${gsm:0:7}nnn"
    download_and_verify \
        "$GEO_BASE/samples/${gsm_prefix}/${gsm}/suppl/${fname}" \
        "$d/${fname}"
done

# ---------------------------------------------------------------------
# 3. GSE285990 -- scRNA-seq. ONLY GSM8714595-GSM8714604 (P01_LM-P10_LM)
#    are the paper's cited human liver-metastasis cohort -- confirmed
#    human (Homo sapiens) during this round's inventory. GSM8714605-9
#    (NT_LM/RS_LM/RL_LM/WT_LM/WT_DTR_LM) and the series-level
#    KCs_gene_exp.txt.gz/WT_DTR_LM_* files are a SEPARATE Mus musculus
#    Kupffer-cell/FOLFOX mechanistic sub-study -- confirmed mouse
#    directly via each GSM's own record, out of scope (wrong species,
#    not the cited cohort). GSM8714610-1 (RS_LM_Stereo-seq/
#    RL_LM_Stereo-seq) are mouse Stereo-seq spatial data -- out of scope
#    on both species and modality grounds. None of the mouse/spatial
#    files are downloaded.
# ---------------------------------------------------------------------
d="$DATA_ROOT/scRNAseq/GSE285990/raw"
mkdir -p "$d"
GSE285990_GSMS="GSM8714595:P01_LM GSM8714596:P02_LM GSM8714597:P03_LM GSM8714598:P04_LM GSM8714599:P05_LM GSM8714600:P06_LM GSM8714601:P07_LM GSM8714602:P08_LM GSM8714603:P09_LM GSM8714604:P10_LM"
for gsm_pat in $GSE285990_GSMS; do
    gsm="${gsm_pat%%:*}"; pat="${gsm_pat##*:}"
    gsm_prefix="${gsm:0:7}nnn"
    for suffix in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do
        fname="${gsm}_${pat}_${suffix}"
        download_and_verify \
            "$GEO_BASE/samples/${gsm_prefix}/${gsm}/suppl/${fname}" \
            "$d/${fname}"
    done
done

# ---------------------------------------------------------------------
# 4. GSE131418 -- bulk microarray, 1,135-sample SuperSeries-scale series
#    (paper cites a 170-sample liver-met subset) -- full RAW.tar;
#    reconstructing the exact 170-sample cohort from clinical/treatment
#    metadata is the inventory script's job, not the download's.
# ---------------------------------------------------------------------
d="$DATA_ROOT/bulkRNAseq/GSE131418/raw"
mkdir -p "$d"
download_and_verify \
    "$GEO_BASE/series/GSE131nnn/GSE131418/suppl/GSE131418_RAW.tar" \
    "$d/GSE131418_RAW.tar"
# Also grab the clinical/annotation files needed for the 170-sample
# cohort reconstruction (not per-sample CEL data, small).
for fname in \
    GSE131418_GEO_submission_Recurrence_meta_data_V2_Updated.xls.gz \
    GSE131418_GEO_submission_stage4_meta_data_V2_Updated.xls.gz \
    GSE131418_Consortium_prim_met_GE_probe_level.txt.gz \
    GSE131418_MCC_prim_met_GE_probe_level.txt.gz; do
    download_and_verify \
        "$GEO_BASE/series/GSE131nnn/GSE131418/suppl/${fname}" \
        "$d/${fname}"
done

# ---------------------------------------------------------------------
# 5/6. GSE17536 (177 human CRC) + GSE17537 (55 human CRC) -- together
#      reconstruct the paper's cited n=232 primary-CRC cohort exactly.
# ---------------------------------------------------------------------
for acc_prefix in "GSE17536:GSE17nnn" "GSE17537:GSE17nnn"; do
    acc="${acc_prefix%%:*}"; prefix="${acc_prefix##*:}"
    d="$DATA_ROOT/bulkRNAseq/${acc}/raw"
    mkdir -p "$d"
    download_and_verify \
        "$GEO_BASE/series/${prefix}/${acc}/suppl/${acc}_RAW.tar" \
        "$d/${acc}_RAW.tar"
done

# ---------------------------------------------------------------------
# 7. GSE21510 -- bulk microarray, 148 samples total; paper cites n=146
#    (unresolved 148->146 discrepancy, per design -- identify during
#    inventory, not assumed). Full RAW.tar downloaded regardless.
# ---------------------------------------------------------------------
d="$DATA_ROOT/bulkRNAseq/GSE21510/raw"
mkdir -p "$d"
download_and_verify \
    "$GEO_BASE/series/GSE21nnn/GSE21510/suppl/GSE21510_RAW.tar" \
    "$d/GSE21510_RAW.tar"

echo "=== Step 7 GEO download finished $(date -u +%FT%TZ), $FAILED failure(s) ===" | tee -a "$LOG"
exit $FAILED
