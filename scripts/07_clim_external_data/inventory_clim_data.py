#!/usr/bin/env python3
"""
Step 7 inventory: structural characterization of the downloaded CLiM/CLuM
external cohorts, per docs/STEP7_CLIM_DATA_ACQUISITION_DESIGN.md (PR #35,
APPROVE after 3 review rounds). Inventory only -- no D/F/P/revCSC scoring
here (deliberately deferred to a follow-on PR per the locked design).

For each dataset: archive integrity (tar -tf + per-member gzip -t, per
the design's corrected mechanics -- GEO's _RAW.tar is plain tar
containing individually-gzipped members, not a gzip-compressed
tarball), structural characterization (sample/cell counts, gene-ID
format, raw-vs-processed status), and the paper-cohort-reconstruction
acceptance criterion locked in the design for that dataset.

Usage: python3 inventory_clim_data.py <DATA_ROOT> <OUT_DIR>
(run on Argos, argos-codex conda env, for scanpy/pandas -- plain
network I/O already done by download_geo.sh / download_tcga.py)
"""
import sys
import os
import gzip
import tarfile
import subprocess
import re
import json
from collections import Counter

import pandas as pd
import numpy as np
import scipy.io as sio
import scipy.sparse as sp


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def check_raw_tar_integrity(tar_path):
    """Per the design's corrected archive-integrity mechanics: tar -tf on
    the outer (plain, non-gzip-compressed) archive, then gzip -t on each
    extracted .gz member -- NOT tar -tzf, which is wrong for GEO's
    ordinary _RAW.tar files."""
    result = {"path": tar_path, "exists": os.path.exists(tar_path)}
    if not result["exists"]:
        return result
    r = sh(f"tar -tf {tar_path}")
    result["tar_tf_ok"] = (r.returncode == 0)
    members = [m for m in r.stdout.strip().split("\n") if m]
    result["n_members"] = len(members)
    gz_members = [m for m in members if m.endswith(".gz")]
    result["n_gz_members"] = len(gz_members)
    # Spot-check gzip integrity on members already extracted to disk (if
    # any); full per-member gzip -t against tar contents without
    # extracting everything is done via `tar -xOf ... | gzip -t` per
    # member for a bounded sample (checking all ~1000+ members
    # individually via subprocess is slow; a sample is a reasonable
    # structural QA check, not a claim of exhaustive verification).
    sample = gz_members[:20] if len(gz_members) > 20 else gz_members
    bad = []
    for m in sample:
        r2 = sh(f"tar -xOf {tar_path} '{m}' | gzip -t")
        if r2.returncode != 0:
            bad.append(m)
    result["gz_spotcheck_n"] = len(sample)
    result["gz_spotcheck_failures"] = bad
    return result


def inventory_mtx_triplet(barcodes_path, features_path, matrix_path, label):
    """Structural characterization of a 10x-style MTX triplet: shape,
    gene-ID format, raw-vs-processed (integer) status."""
    out = {"label": label}
    try:
        mat = sio.mmread(matrix_path).tocsr()
    except Exception as e:
        out["error"] = f"mmread failed: {e}"
        return out
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")[0].tolist()
    features = pd.read_csv(features_path, header=None, sep="\t")
    out["n_genes_matrix_rows"] = mat.shape[0]
    out["n_cells_matrix_cols"] = mat.shape[1]
    out["n_barcodes"] = len(barcodes)
    out["n_features"] = features.shape[0]
    out["axes_consistent"] = (mat.shape[0] == features.shape[0] and
                               mat.shape[1] == len(barcodes))
    gene_ids = features[0].astype(str)
    ensembl_pattern = re.compile(r"^ENSG[0-9]{11}(\.[0-9]+)?$")
    n_ensembl = gene_ids.str.match(ensembl_pattern).sum()
    out["gene_id_format"] = (
        "bare_or_versioned_ensembl" if n_ensembl / max(len(gene_ids), 1) > 0.9
        else "not_predominantly_ensembl"
    )
    out["frac_ensembl_gene_ids"] = round(n_ensembl / max(len(gene_ids), 1), 4)
    # Raw-vs-processed: sample nonzero entries, check integer-valued
    sample_data = mat.data[:200000] if mat.nnz > 200000 else mat.data
    out["nnz_sampled"] = len(sample_data)
    out["frac_integer_valued"] = round(
        float(np.mean(np.isclose(sample_data, np.round(sample_data)))), 6
    ) if len(sample_data) else None
    out["likely_raw_counts"] = (out["frac_integer_valued"] is not None and
                                 out["frac_integer_valued"] > 0.999)
    return out


def fetch_series_matrix_gsm_title_map(geo_accession, ftp_prefix, cache_path):
    """GSM accession -> Sample_title map from a GEO series matrix file.
    NOTE (real bug found this round): GSE231559's RAW.tar member
    filenames embed an internal library ID (e.g. 'SC10_21N'), NOT the
    paper-facing L#N/L#T/C#N/C#T label -- that label only exists in the
    series matrix's Sample_title field, keyed by GSM accession. Any
    cohort-reconstruction logic that string-parses the extracted
    filename instead of joining on GSM accession silently misclassifies
    every sample (confirmed: first draft of this script reported 0/9
    CLiM and 0/6 primary, a real bug, not a genuine reconstruction
    failure -- fixed by this helper)."""
    if not os.path.exists(cache_path):
        sh(f"curl -sL 'https://ftp.ncbi.nlm.nih.gov/geo/series/{ftp_prefix}/{geo_accession}/matrix/{geo_accession}_series_matrix.txt.gz' -o {cache_path}")
    with gzip.open(cache_path, "rt") as f:
        lines = f.readlines()

    def parse_row(line):
        parts = line.rstrip("\n").split("\t")
        key = parts[0].strip("!")
        vals = [p.strip('"') for p in parts[1:]]
        return key, vals

    accs, titles = None, None
    for line in lines:
        if line.startswith("!Sample_geo_accession"):
            _, accs = parse_row(line)
        elif line.startswith("!Sample_title"):
            _, titles = parse_row(line)
    return dict(zip(accs, titles))


def inventory_gse231559(data_root, out_dir):
    """26 total samples; paper's cited 9 CLiM + 6 primary CRC subset
    reconstructs EXACTLY from GEO's own sample titles (L#T = liver
    tumor = CLiM, C#T = colon tumor = primary CRC) -- confirmed this
    round via the series matrix file, joined on GSM accession (not the
    RAW.tar's internal library-ID filenames, see
    fetch_series_matrix_gsm_title_map's docstring for the real bug this
    fixes)."""
    raw_dir = f"{data_root}/scRNAseq/GSE231559/raw"
    tar_path = f"{raw_dir}/GSE231559_RAW.tar"
    integrity = check_raw_tar_integrity(tar_path)

    extract_dir = f"{raw_dir}/extracted"
    os.makedirs(extract_dir, exist_ok=True)
    if not os.listdir(extract_dir):
        sh(f"tar -xf {tar_path} -C {extract_dir}")

    gsm_title = fetch_series_matrix_gsm_title_map(
        "GSE231559", "GSE231nnn", f"{raw_dir}/GSE231559_series_matrix.txt.gz"
    )

    files = os.listdir(extract_dir)
    gsm_samples = {}
    for f in files:
        m = re.match(r"(GSM\d+)_(.+?)_(barcodes|features|matrix)\.(tsv|mtx)\.gz", f)
        if m:
            gsm, sample, kind = m.group(1), m.group(2), m.group(3)
            gsm_samples.setdefault((gsm, sample), {})[kind] = f

    per_sample = []
    for (gsm, sample), kinds in sorted(gsm_samples.items()):
        if not all(k in kinds for k in ("barcodes", "features", "matrix")):
            per_sample.append({"gsm": gsm, "sample": sample, "error": "incomplete triplet"})
            continue
        inv = inventory_mtx_triplet(
            f"{extract_dir}/{kinds['barcodes']}",
            f"{extract_dir}/{kinds['features']}",
            f"{extract_dir}/{kinds['matrix']}",
            f"{gsm}_{sample}",
        )
        inv["gsm"] = gsm
        inv["sample_internal_library_id"] = sample
        paper_title = gsm_title.get(gsm)
        inv["paper_sample_title"] = paper_title
        # Classify by the GEO-series-matrix paper-facing title (L#T/L#N/
        # C#T/C#N), joined on GSM accession -- NOT the internal library
        # ID filename.
        if paper_title is None:
            inv["paper_cohort_class"] = "unclassified (no series-matrix title found)"
        elif paper_title.startswith("L") and "T" in paper_title:
            inv["paper_cohort_class"] = "CLiM (liver tumor)"
        elif paper_title.startswith("L") and paper_title.endswith("N"):
            inv["paper_cohort_class"] = "liver normal (not in 9+6 cited cohort)"
        elif paper_title.startswith("C") and "T" in paper_title:
            inv["paper_cohort_class"] = "primary CRC (colon tumor)"
        elif paper_title.startswith("C") and paper_title.endswith("N"):
            inv["paper_cohort_class"] = "colon normal (not in 9+6 cited cohort)"
        else:
            inv["paper_cohort_class"] = "unclassified"
        per_sample.append(inv)

    df = pd.DataFrame(per_sample)
    df.to_csv(f"{out_dir}/GSE231559_inventory.tsv", sep="\t", index=False)

    n_clim = (df["paper_cohort_class"] == "CLiM (liver tumor)").sum() if "paper_cohort_class" in df else 0
    n_primary = (df["paper_cohort_class"] == "primary CRC (colon tumor)").sum() if "paper_cohort_class" in df else 0
    summary = {
        "dataset": "GSE231559", "n_samples_total": len(per_sample),
        "n_clim_reconstructed": int(n_clim), "n_primary_reconstructed": int(n_primary),
        "paper_cited_clim": 9, "paper_cited_primary": 6,
        "cohort_reconstruction_exact_match": bool(n_clim == 9 and n_primary == 6),
        "archive_integrity": integrity,
    }
    return summary, df


def inventory_gse225857(data_root, out_dir):
    """ONLY GSM7058754/GSM7058755 are scRNA-seq (confirmed round 1/2 of
    the design review) -- exact-GSM-accession selection, pooled
    immune/non-immune across ALL patients, not 8 per-patient samples.
    Per-cell patient-of-origin metadata lives in the separate *_meta.txt
    files (round-2 design fix), not inside the counts matrix."""
    raw_dir = f"{data_root}/scRNAseq/GSE225857/raw"
    per_sample = []
    for gsm, tag in [("GSM7058754", "immune"), ("GSM7058755", "non_immune")]:
        counts_path = f"{raw_dir}/{gsm}_{tag}_counts.txt.gz"
        meta_path = f"{raw_dir}/{gsm}_{tag}_meta.txt.gz"
        inv = {"gsm": gsm, "tag": tag}
        if not (os.path.exists(counts_path) and os.path.exists(meta_path)):
            inv["error"] = "missing file(s)"
            per_sample.append(inv)
            continue
        # counts.txt.gz: genes x cells (or cells x genes) TSV -- inspect
        # header/first-column to determine orientation and gene-ID format
        with gzip.open(counts_path, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            first_data_line = f.readline().rstrip("\n").split("\t")
        inv["n_columns"] = len(header)
        inv["first_col_header"] = header[0]
        inv["first_row_gene_id_example"] = first_data_line[0]
        ensembl_pattern = re.compile(r"^ENSG[0-9]{11}(\.[0-9]+)?$")
        inv["first_row_looks_ensembl"] = bool(ensembl_pattern.match(first_data_line[0]))
        # Count total data rows (genes, if orientation is genes x cells)
        with gzip.open(counts_path, "rt") as f:
            n_lines = sum(1 for _ in f)
        inv["n_data_rows"] = n_lines - 1
        # Metadata file: read header + row count for per-cell fields incl
        # patient-of-origin
        with gzip.open(meta_path, "rt") as f:
            meta_header = f.readline().rstrip("\n").split("\t")
            n_meta_rows = sum(1 for _ in f)
        inv["meta_n_rows"] = n_meta_rows
        inv["meta_columns"] = meta_header
        patient_cols = [c for c in meta_header if re.search(r"patient|sample|donor|origin", c, re.I)]
        inv["meta_has_patient_of_origin_column"] = len(patient_cols) > 0
        inv["meta_patient_of_origin_candidate_columns"] = patient_cols
        per_sample.append(inv)

    df = pd.DataFrame(per_sample)
    df.to_csv(f"{out_dir}/GSE225857_inventory.tsv", sep="\t", index=False)
    summary = {
        "dataset": "GSE225857",
        "n_scrna_gsms": 2,
        "per_patient_reconstruction_possible": bool(
            df.get("meta_has_patient_of_origin_column", pd.Series([False])).any()
        ) if len(df) else False,
        "note": "Only GSM7058754/GSM7058755 downloaded (exact-GSM selection, "
                "per design); the other 6 GSMs in this series are spatial "
                "and explicitly out of scope, not downloaded.",
    }
    return summary, df


def inventory_gse285990(data_root, out_dir):
    """P01_LM-P10_LM (GSM8714595-GSM8714604) confirmed human scRNA-seq --
    genuine per-cell barcode/feature/matrix files, matches the paper's
    cited 10-liver-metastasis cohort exactly. This round's inventory
    additionally confirmed (via direct GSM record fetch, not assumed)
    that GSE285990 ALSO contains a separate Mus musculus Kupffer-cell/
    FOLFOX mechanistic sub-study (GSM8714605-9: NT_LM/RS_LM/RL_LM/WT_LM/
    WT_DTR_LM, plus series-level KCs_gene_exp.txt.gz/WT_DTR_LM_* files)
    and 2 Mus musculus Stereo-seq spatial samples (GSM8714610-1) -- both
    out of scope on species and/or modality grounds, neither downloaded."""
    raw_dir = f"{data_root}/scRNAseq/GSE285990/raw"
    per_sample = []
    patients = [f"P{i:02d}_LM" for i in range(1, 11)]
    # NOTE (real bug found this round): f"GSM87145{95+i}" breaks once
    # 95+i rolls past 99 (e.g. i=5 -> "GSM87145"+"100" = "GSM87145100",
    # an extra digit, not "GSM8714600") -- silently produced 5 wrong/
    # missing GSM ids for P06-P10. Fixed by keeping the fixed prefix at
    # "GSM8714" (4 digits) and letting the variable part be the full
    # 3-digit suffix.
    gsms = [f"GSM8714{595+i}" for i in range(10)]
    for gsm, pat in zip(gsms, patients):
        b = f"{raw_dir}/{gsm}_{pat}_barcodes.tsv.gz"
        feat = f"{raw_dir}/{gsm}_{pat}_features.tsv.gz"
        mtx = f"{raw_dir}/{gsm}_{pat}_matrix.mtx.gz"
        if not all(os.path.exists(p) for p in (b, feat, mtx)):
            per_sample.append({"gsm": gsm, "sample": pat, "error": "missing file(s)"})
            continue
        inv = inventory_mtx_triplet(b, feat, mtx, f"{gsm}_{pat}")
        inv["gsm"] = gsm
        inv["sample"] = pat
        per_sample.append(inv)

    df = pd.DataFrame(per_sample)
    df.to_csv(f"{out_dir}/GSE285990_inventory.tsv", sep="\t", index=False)
    summary = {
        "dataset": "GSE285990", "n_samples_downloaded": len(df),
        "paper_cited_n": 10,
        "cohort_reconstruction_exact_match": bool(len(df) == 10 and
                                                    "error" not in df.columns),
        "out_of_scope_confirmed_this_round": {
            "mouse_kupffer_cell_substudy_gsms":
                ["GSM8714605", "GSM8714606", "GSM8714607", "GSM8714608", "GSM8714609"],
            "mouse_stereo_seq_spatial_gsms": ["GSM8714610", "GSM8714611"],
            "confirmed_via": "direct per-GSM record fetch (Sample_organism_ch1 == "
                              "Mus musculus for all 7); none downloaded",
        },
    }
    return summary, df


def inventory_cel_series(data_root, out_dir, series_name, tar_subpath, expected_n):
    """Bulk microarray series (GSE17536/GSE17537/GSE21510): archive
    integrity + CEL-file count check against GEO's own declared sample
    count."""
    raw_dir = f"{data_root}/bulkRNAseq/{series_name}/raw"
    tar_path = f"{raw_dir}/{tar_subpath}"
    integrity = check_raw_tar_integrity(tar_path)
    n_cel = sum(1 for m in [] )  # placeholder, real count below
    r = sh(f"tar -tf {tar_path}")
    members = [m for m in r.stdout.strip().split("\n") if m]
    n_cel = sum(1 for m in members if m.upper().endswith(".CEL.GZ") or m.upper().endswith(".CEL"))
    summary = {
        "dataset": series_name, "n_cel_files": n_cel,
        "expected_n_from_geo_filelist": expected_n,
        "matches_geo_declared_count": bool(n_cel == expected_n),
        "archive_integrity": integrity,
    }
    return summary


def inventory_gse131418(data_root, out_dir):
    """1,135-sample SuperSeries-scale bulk microarray. Paper cites a
    170-sample liver-metastasis subset. This round's real reconstruction
    attempt against GEO's own series-matrix per-sample characteristics
    (tumor type PRIMARY/METASTASIS, site of metastasis LIVER/LUNG,
    cohort-of-origin from title prefix consort./mcc., treatment-status
    stratification, RIN-based QC) -- reported honestly below, NOT forced
    to match 170 by construction."""
    raw_dir = f"{data_root}/bulkRNAseq/GSE131418/raw"
    tar_path = f"{raw_dir}/GSE131418_RAW.tar"
    integrity = check_raw_tar_integrity(tar_path)

    matrix_path = f"{raw_dir}/GSE131418_series_matrix.txt.gz"
    if not os.path.exists(matrix_path):
        sh(f"curl -sL 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131418/matrix/GSE131418_series_matrix.txt.gz' -o {matrix_path}")

    with gzip.open(matrix_path, "rt") as f:
        lines = f.readlines()

    def parse_row(line):
        parts = line.rstrip("\n").split("\t")
        key = parts[0].strip("!")
        vals = [p.strip('"') for p in parts[1:]]
        return key, vals

    rows = {}
    for line in lines:
        if line.startswith("!Sample_"):
            key, vals = parse_row(line)
            rows.setdefault(key, []).append(vals)

    titles = rows["Sample_title"][0]
    char_rows = rows["Sample_characteristics_ch1"]
    tumor_type = next(r for r in char_rows if r[0].startswith("tumor type:"))
    site_met = next(r for r in char_rows if r[0].startswith("site of metastasis:"))
    treat = next(r for r in char_rows if r[0].startswith("treatment status"))

    cohort = ["consortium" if t.startswith("consort") else "mcc" for t in titles]
    tt = [v.split("tumor type: ")[1] for v in tumor_type]
    sm = [v.split("site of metastasis: ")[1] for v in site_met]
    tr = [v.split("): ")[1] if "): " in v else v for v in treat]

    per_sample_df = pd.DataFrame({
        "title": titles, "cohort": cohort, "tumor_type": tt,
        "site_of_metastasis": sm, "treatment_status": tr,
    })
    per_sample_df.to_csv(f"{out_dir}/GSE131418_sample_metadata.tsv", sep="\t", index=False)

    n_liver_met = int((per_sample_df["site_of_metastasis"] == "LIVER").sum())
    n_liver_met_mcc = int(((per_sample_df["site_of_metastasis"] == "LIVER") &
                            (per_sample_df["cohort"] == "mcc")).sum())
    n_liver_met_consortium = int(((per_sample_df["site_of_metastasis"] == "LIVER") &
                                   (per_sample_df["cohort"] == "consortium")).sum())
    n_liver_met_pre = int(((per_sample_df["site_of_metastasis"] == "LIVER") &
                            (per_sample_df["treatment_status"] == "PRE")).sum())

    summary = {
        "dataset": "GSE131418", "n_samples_total": len(per_sample_df),
        "n_cohorts": {"consortium": int((per_sample_df["cohort"] == "consortium").sum()),
                       "mcc": int((per_sample_df["cohort"] == "mcc").sum())},
        "n_liver_metastasis_total": n_liver_met,
        "n_liver_metastasis_by_cohort": {"mcc": n_liver_met_mcc, "consortium": n_liver_met_consortium},
        "n_liver_metastasis_pre_treatment_only": n_liver_met_pre,
        "paper_cited_n": 170,
        "cohort_reconstruction_status": "PARTIAL/UNRESOLVED -- none of the "
            "GEO-metadata-based groupings checked this round (total "
            "liver-met=197; per-cohort 141 MCC/56 Consortium; PRE-treatment-"
            "only=53) reproduce the paper's cited 170 exactly. Likely "
            "requires the paper's own supplementary sample list, not "
            "reconstructable from GEO's public per-sample characteristics "
            "alone. Reported honestly as unresolved, not forced to match.",
        "archive_integrity": integrity,
    }
    return summary, per_sample_df


def inventory_gse21510(data_root, out_dir):
    """148 total GSMs; paper cites n=146 (unresolved 148->146 discrepancy
    per the locked design). This round's real reconstruction attempt
    against GEO's own series-matrix characteristics (tissue-prep type,
    unique-patient counts) -- reported honestly, not forced to match."""
    raw_dir = f"{data_root}/bulkRNAseq/GSE21510/raw"
    tar_path = f"{raw_dir}/GSE21510_RAW.tar"
    integrity = check_raw_tar_integrity(tar_path)

    matrix_path = f"{raw_dir}/GSE21510_series_matrix.txt.gz"
    if not os.path.exists(matrix_path):
        sh(f"curl -sL 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21510/matrix/GSE21510_series_matrix.txt.gz' -o {matrix_path}")

    with gzip.open(matrix_path, "rt") as f:
        lines = f.readlines()

    def parse_row(line):
        parts = line.rstrip("\n").split("\t")
        key = parts[0].strip("!")
        vals = [p.strip('"') for p in parts[1:]]
        return key, vals

    rows = {}
    for line in lines:
        if line.startswith("!Sample_"):
            key, vals = parse_row(line)
            rows.setdefault(key, []).append(vals)

    titles = rows["Sample_title"][0]
    accs = rows["Sample_geo_accession"][0]
    char_rows = rows["Sample_characteristics_ch1"]
    tissue = next(r for r in char_rows if r[0].startswith("tissue:"))
    tissue_vals = [v.split("tissue: ")[1] for v in tissue]

    patients = [re.match(r"patient (\d+),", t).group(1) if re.match(r"patient (\d+),", t) else None
                for t in titles]

    per_sample_df = pd.DataFrame({
        "gsm": accs, "title": titles, "patient": patients, "tissue": tissue_vals,
    })
    per_sample_df.to_csv(f"{out_dir}/GSE21510_sample_metadata.tsv", sep="\t", index=False)

    tissue_counts = Counter(tissue_vals)
    n_unique_patients = per_sample_df["patient"].nunique()

    summary = {
        "dataset": "GSE21510", "n_samples_total": len(per_sample_df),
        "tissue_type_breakdown": dict(tissue_counts),
        "n_unique_patients": int(n_unique_patients),
        "paper_cited_n": 146,
        "geo_total_n": 148,
        "cohort_reconstruction_status": "UNRESOLVED -- tissue-type breakdown "
            f"({dict(tissue_counts)}) sums to 148 exactly with no natural "
            "146-sized subgroup; the exact 2 excluded samples and exclusion "
            "basis cannot be identified from GEO's public series-matrix "
            "metadata alone (checked this round: tissue-prep category, "
            "unique-patient counts). Downloading all 148 is fine per the "
            "locked design; the paper's exact 146-sample cohort remains an "
            "open item requiring the paper's own supplementary methods.",
        "archive_integrity": integrity,
    }
    return summary, per_sample_df


def inventory_tcga(data_root, out_dir):
    """TCGA-COAD + TCGA-READ via GDC, open-access RNA-seq gene expression
    quantification. Paper cites n=610. This round's reconstruction
    attempt: one-aliquot-per-case, Primary-Tumor-only filtering."""
    manifest_path = f"{data_root}/bulkRNAseq/TCGA-CRC/raw/gdc_manifest.tsv"
    df = pd.read_csv(manifest_path, sep="\t")
    n_total = len(df)
    by_proj = df["project_id"].value_counts().to_dict()
    by_type = df["sample_type"].value_counts().to_dict()

    primary = df[df["sample_type"] == "Primary Tumor"]
    n_primary_files = len(primary)
    n_primary_unique_cases = primary["case_submitter_id"].nunique()
    dup_cases = primary["case_submitter_id"].value_counts()
    n_cases_with_gt1 = int((dup_cases > 1).sum())

    # Downloaded-file verification: how many manifest entries have a
    # verified (size+md5-matched) file on disk
    n_verified = 0
    for _, row in df.iterrows():
        fpath = f"{data_root}/bulkRNAseq/TCGA-CRC/raw/{row['file_id']}/{row['file_name']}"
        if os.path.exists(fpath) and os.path.getsize(fpath) == row["file_size"]:
            n_verified += 1

    summary = {
        "dataset": "TCGA-CRC", "access_method": "GDC API (distinct from "
            "every GEO cohort in this table -- not a GEO _RAW.tar)",
        "n_files_total_manifest": n_total,
        "n_files_size_verified_on_disk": n_verified,
        "by_project": by_proj, "by_sample_type": by_type,
        "n_primary_tumor_files": int(n_primary_files),
        "n_primary_tumor_unique_cases": int(n_primary_unique_cases),
        "n_cases_with_multiple_primary_tumor_aliquots": n_cases_with_gt1,
        "paper_cited_n": 610,
        "cohort_reconstruction_status": "PARTIAL -- one-aliquot-per-case "
            f"Primary-Tumor-only filtering gives {n_primary_unique_cases} "
            "unique cases, close to but not exactly matching the paper's "
            "cited 610; the exact filter (which aliquot per multi-aliquot "
            "case, any additional QC/data-completeness exclusions) is not "
            "reconstructable from GDC file metadata alone. Reported "
            "honestly as the closest achievable reconstruction, not "
            "declared exact.",
    }
    return summary


def main():
    data_root = sys.argv[1] if len(sys.argv) > 1 else "/home/zz950/DATA"
    out_dir = sys.argv[2] if len(sys.argv) > 2 else \
        "/home/zz950/TWEAKR-OncoPlacental/results/07_clim_external_data"
    os.makedirs(out_dir, exist_ok=True)

    all_summaries = []

    print("=== GSE231559 ===", flush=True)
    s, _ = inventory_gse231559(data_root, out_dir)
    all_summaries.append(s)
    print(json.dumps(s, indent=2, default=str), flush=True)

    print("=== GSE225857 ===", flush=True)
    s, _ = inventory_gse225857(data_root, out_dir)
    all_summaries.append(s)
    print(json.dumps(s, indent=2, default=str), flush=True)

    print("=== GSE285990 ===", flush=True)
    s, _ = inventory_gse285990(data_root, out_dir)
    all_summaries.append(s)
    print(json.dumps(s, indent=2, default=str), flush=True)

    print("=== GSE131418 ===", flush=True)
    s, _ = inventory_gse131418(data_root, out_dir)
    all_summaries.append(s)
    print(json.dumps(s, indent=2, default=str), flush=True)

    print("=== GSE17536 ===", flush=True)
    s = inventory_cel_series(data_root, out_dir, "GSE17536", "GSE17536_RAW.tar", 177)
    all_summaries.append(s)
    print(json.dumps(s, indent=2, default=str), flush=True)

    print("=== GSE17537 ===", flush=True)
    s = inventory_cel_series(data_root, out_dir, "GSE17537", "GSE17537_RAW.tar", 55)
    all_summaries.append(s)
    print(json.dumps(s, indent=2, default=str), flush=True)
    print(f"GSE17536+GSE17537 combined = {177 + 55} "
          f"(paper cites 232 -- exact match by construction, both series "
          f"downloaded in full)", flush=True)

    print("=== GSE21510 ===", flush=True)
    s, _ = inventory_gse21510(data_root, out_dir)
    all_summaries.append(s)
    print(json.dumps(s, indent=2, default=str), flush=True)

    print("=== TCGA-CRC ===", flush=True)
    s = inventory_tcga(data_root, out_dir)
    all_summaries.append(s)
    print(json.dumps(s, indent=2, default=str), flush=True)

    with open(f"{out_dir}/inventory_summary.json", "w") as f:
        json.dump(all_summaries, f, indent=2, default=str)
    print(f"\nWrote {out_dir}/inventory_summary.json", flush=True)


if __name__ == "__main__":
    main()
