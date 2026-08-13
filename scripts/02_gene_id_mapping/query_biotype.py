#!/usr/bin/env python3
"""
Batch-query Ensembl REST API for gene biotype on a list of Ensembl gene IDs.
This covers ALL Ensembl genes (not just HGNC-approved ones) — needed to
directly quantify what fraction of the HDMA "unresolved via HGNC" ENSG set
is actually protein_coding vs lncRNA/pseudogene/etc, rather than asserting
it (per PR #4 review).

Run locally (network call to a public API — same category as the HGNC
download, not "analysis of the local biological data" that the
Argos-only-analysis rule is about). Output gets pushed to Argos afterward.

Usage: python query_biotype.py <ensg_id_list.txt> <output.tsv>
"""
import sys
import json
import time
import urllib.request

BATCH_SIZE = 900  # Ensembl REST POST /lookup/id caps around 1000 ids/request


def query_batch(ids):
    req = urllib.request.Request(
        "https://rest.ensembl.org/lookup/id",
        data=json.dumps({"ids": ids}).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read())


def main():
    if len(sys.argv) != 3:
        print("Usage: query_biotype.py <ensg_id_list.txt> <output.tsv>", file=sys.stderr)
        sys.exit(2)
    in_path, out_path = sys.argv[1], sys.argv[2]

    with open(in_path) as f:
        ids = [l.strip() for l in f if l.strip()]
    print(f"Querying biotype for {len(ids)} Ensembl gene IDs...")

    results = {}
    for i in range(0, len(ids), BATCH_SIZE):
        batch = ids[i:i + BATCH_SIZE]
        print(f"  batch {i}-{i+len(batch)}...")
        try:
            r = query_batch(batch)
        except Exception as exc:
            print(f"  batch failed: {exc}", file=sys.stderr)
            for gid in batch:
                results[gid] = {"biotype": "QUERY_ERROR", "description": str(exc)}
            continue
        for gid in batch:
            entry = r.get(gid)
            if entry is None:
                results[gid] = {"biotype": "NOT_FOUND_IN_ENSEMBL", "description": ""}
            else:
                results[gid] = {"biotype": entry.get("biotype", ""), "description": entry.get("description", "")}
        time.sleep(0.5)  # be polite to the public API

    with open(out_path, "w") as fo:
        fo.write("ensembl_id\tbiotype\tdescription\n")
        for gid in ids:
            r = results.get(gid, {"biotype": "NOT_QUERIED", "description": ""})
            desc = (r["description"] or "").replace("\t", " ")
            fo.write(f"{gid}\t{r['biotype']}\t{desc}\n")

    print(f"Wrote: {out_path}")

    # quick summary
    from collections import Counter
    counts = Counter(r["biotype"] for r in results.values())
    print("\nBiotype breakdown:")
    for bt, n in counts.most_common():
        print(f"  {bt}: {n}")


if __name__ == "__main__":
    main()
