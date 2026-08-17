#!/usr/bin/env python3
"""
Step 7 inventory: TCGA-CRC download via the GDC Data Portal (GDC API),
per docs/STEP7_CLIM_DATA_ACQUISITION_DESIGN.md (PR #35). Access method
differs from every GEO cohort in this project's Step 7 table -- GDC
API/data endpoint, not a GEO `_RAW.tar`.

Cohort: TCGA-COAD + TCGA-READ, Transcriptome Profiling / Gene Expression
Quantification / RNA-Seq, open-access only (gene-level STAR-Counts
`augmented_star_gene_counts.tsv`, no dbGaP application needed -- only raw
BAMs/germline variants are controlled-access, not touched here).

Downloads ALL matching files (not pre-filtered to the paper's cited
n=610) -- same "download everything, reconstruct the exact cohort as an
explicit, honestly-reported inventory finding" principle already applied
to every GEO SuperSeries in this round (GSE131418, GSE17538, GSE21510):
declaring the cohort reconstructed is not the same as downloading it.

Byte/checksum verification uses the GDC manifest's own reported file
size + md5sum (never a rounded display value), per the design's locked
byte-verification mechanics.

Usage: python3 download_tcga.py <DATA_ROOT>
(run on Argos, argos-qsub1, plain network I/O -- not qsub-wrapped)
"""
import sys
import os
import json
import hashlib
import subprocess
import tarfile
import time
import urllib.request

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
GDC_DATA_URL = "https://api.gdc.cancer.gov/data"
BATCH_SIZE = 50

QUERY_FILTERS = {
    "op": "and",
    "content": [
        {"op": "in", "content": {"field": "cases.project.project_id",
                                  "value": ["TCGA-COAD", "TCGA-READ"]}},
        {"op": "in", "content": {"field": "files.data_category",
                                  "value": ["Transcriptome Profiling"]}},
        {"op": "in", "content": {"field": "files.data_type",
                                  "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "files.experimental_strategy",
                                  "value": ["RNA-Seq"]}},
        {"op": "in", "content": {"field": "files.access",
                                  "value": ["open"]}},
    ],
}


def fetch_manifest():
    """Query the GDC files endpoint for the full open-access TCGA-COAD +
    TCGA-READ RNA-seq gene-expression-quantification file list."""
    payload = {
        "filters": QUERY_FILTERS,
        "fields": "file_id,file_name,file_size,md5sum,"
                  "cases.project.project_id,cases.submitter_id,"
                  "cases.samples.sample_type",
        "format": "JSON",
        "size": "2000",
    }
    req = urllib.request.Request(
        GDC_FILES_URL, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"}, method="POST",
    )
    with urllib.request.urlopen(req) as resp:
        data = json.load(resp)
    hits = data["data"]["hits"]
    manifest = []
    for h in hits:
        case = h["cases"][0]
        proj = case["project"]["project_id"]
        subm = case["submitter_id"]
        samples = case.get("samples") or []
        stype = samples[0]["sample_type"] if samples else "NA"
        manifest.append({
            "file_id": h["file_id"], "file_name": h["file_name"],
            "file_size": h["file_size"], "md5sum": h["md5sum"],
            "project_id": proj, "case_submitter_id": subm,
            "sample_type": stype,
        })
    return manifest


def md5_of_file(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def download_batch(file_ids, dest_dir, tmp_tar, manifest_by_id):
    """Downloads a batch of GDC file_ids. NOTE: the GDC /data endpoint's
    behavior depends on batch size -- multiple ids return a tar.gz
    bundle (MANIFEST.txt + one subdirectory per file_id), but a batch of
    exactly ONE id returns that file's raw bytes directly, no tar
    wrapper (a real GDC API quirk, discovered this round when 701 % 50
    == 1 left a final single-file batch that failed with "not a gzip
    file" -- confirmed by direct single-ID request: same byte size and
    md5 as the manifest, just not tar-wrapped). Both cases are handled
    explicitly rather than assuming tar.gz always."""
    payload = json.dumps({"ids": file_ids}).encode()
    req = urllib.request.Request(
        GDC_DATA_URL, data=payload,
        headers={"Content-Type": "application/json"}, method="POST",
    )
    with urllib.request.urlopen(req, timeout=600) as resp, open(tmp_tar, "wb") as f:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
    if len(file_ids) == 1:
        fid = file_ids[0]
        rec = manifest_by_id[fid]
        out_subdir = os.path.join(dest_dir, fid)
        os.makedirs(out_subdir, exist_ok=True)
        out_path = os.path.join(out_subdir, rec["file_name"])
        os.replace(tmp_tar, out_path)
        return
    with tarfile.open(tmp_tar, "r:gz") as tf:
        tf.extractall(dest_dir)
    os.remove(tmp_tar)


def main():
    data_root = sys.argv[1] if len(sys.argv) > 1 else "/home/zz950/DATA"
    dest_dir = os.path.join(data_root, "bulkRNAseq", "TCGA-CRC", "raw")
    os.makedirs(dest_dir, exist_ok=True)

    print("=== Fetching GDC manifest ===", flush=True)
    manifest = fetch_manifest()
    print(f"Manifest: {len(manifest)} files", flush=True)

    manifest_path = os.path.join(dest_dir, "gdc_manifest.tsv")
    with open(manifest_path, "w") as f:
        f.write("file_id\tfile_name\tfile_size\tmd5sum\tproject_id\t"
                "case_submitter_id\tsample_type\n")
        for r in manifest:
            f.write(f"{r['file_id']}\t{r['file_name']}\t{r['file_size']}\t"
                     f"{r['md5sum']}\t{r['project_id']}\t"
                     f"{r['case_submitter_id']}\t{r['sample_type']}\n")
    print(f"Wrote {manifest_path}", flush=True)

    by_id = {r["file_id"]: r for r in manifest}
    all_ids = list(by_id.keys())

    n_ok, n_fail, n_skip = 0, 0, 0
    fail_log = []
    for i in range(0, len(all_ids), BATCH_SIZE):
        batch_ids = all_ids[i:i + BATCH_SIZE]
        # Skip file_ids whose extracted file already matches expected
        # size+md5 (resumability). Round-1 review fix: this previously
        # only checked size, contradicting the docstring/PR claim of
        # "size+md5" resumability -- a stale/corrupted same-length file
        # could have silently survived a rerun without ever satisfying
        # the manifest checksum. Fixed to also recompute md5 (cheap for
        # already-downloaded ~4MB files, a few seconds total for 700).
        need = []
        for fid in batch_ids:
            rec = by_id[fid]
            fpath = os.path.join(dest_dir, fid, rec["file_name"])
            if (os.path.exists(fpath) and os.path.getsize(fpath) == rec["file_size"]
                    and md5_of_file(fpath) == rec["md5sum"]):
                n_skip += 1
            else:
                need.append(fid)
        if not need:
            continue
        tmp_tar = os.path.join(dest_dir, f"_batch_{i}.tar.gz")
        t0 = time.time()
        try:
            download_batch(need, dest_dir, tmp_tar, by_id)
        except Exception as e:
            print(f"[FAIL] batch {i}-{i+len(need)}: download error {e}", flush=True)
            n_fail += len(need)
            fail_log.extend(need)
            continue
        print(f"batch {i}: {len(need)} files in {time.time()-t0:.1f}s", flush=True)
        for fid in need:
            rec = by_id[fid]
            fpath = os.path.join(dest_dir, fid, rec["file_name"])
            if not os.path.exists(fpath):
                print(f"[FAIL] {fid} -- extracted file missing: {fpath}", flush=True)
                n_fail += 1
                fail_log.append(fid)
                continue
            actual_size = os.path.getsize(fpath)
            if actual_size != rec["file_size"]:
                print(f"[FAIL] {fid} -- size mismatch: expected "
                      f"{rec['file_size']}, got {actual_size}", flush=True)
                n_fail += 1
                fail_log.append(fid)
                continue
            actual_md5 = md5_of_file(fpath)
            if actual_md5 != rec["md5sum"]:
                print(f"[FAIL] {fid} -- md5 mismatch: expected "
                      f"{rec['md5sum']}, got {actual_md5}", flush=True)
                n_fail += 1
                fail_log.append(fid)
                continue
            n_ok += 1

    print(f"=== Done: {n_ok} verified (size+md5), {n_skip} already-present "
          f"(resume), {n_fail} failed ===", flush=True)
    if fail_log:
        with open(os.path.join(dest_dir, "_failed_file_ids.txt"), "w") as f:
            f.write("\n".join(fail_log) + "\n")
        print(f"Failed file_ids written to _failed_file_ids.txt -- re-run "
              f"this script to retry (resumable).", flush=True)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
