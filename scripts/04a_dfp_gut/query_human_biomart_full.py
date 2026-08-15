#!/usr/bin/env python3
"""
Bulk-fetch BioMart's authoritative external_gene_name for every human
Ensembl gene ID, per PR #21 round-3 review: the previous duplicate-symbol
detection heuristic (strip a trailing "-\\d+" from var_name and treat
matches as anndata-uniquification artifacts) is unreliable, because many
real HGNC symbols legitimately end in a dash-number as part of their own
canonical name (e.g. immunoglobulin/TCR variable-region genes like
"IGHVIII-5" -- exactly the class of gene the previous heuristic flagged
as "mostly Ig/TCR families", which the reviewer correctly pointed out is
the heuristic's own likely failure mode, not confirmation of anything).

This is the same bulk-BioMart pattern already used successfully in this
project (mouse_biomart_full.tsv, ~20s for tens of thousands of genes) --
one query for the whole species dataset, not per-gene lookups.

Usage: python3 query_human_biomart_full.py <out_tsv>
(pure public-API gene-ID query, run locally per project's narrow
exception -- no qsub needed)
"""
import sys
import time
import urllib.request
import urllib.parse

OUT = sys.argv[1] if len(sys.argv) > 1 else "human_biomart_full.tsv"

BIOMART_URL = "https://jun2026.archive.ensembl.org/biomart/martservice"

QUERY_XML = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id" />
    <Attribute name="external_gene_name" />
  </Dataset>
</Query>"""


def fetch(max_retries=6):
    data = urllib.parse.urlencode({"query": QUERY_XML}).encode()
    for attempt in range(1, max_retries + 1):
        try:
            req = urllib.request.Request(BIOMART_URL, data=data, method="POST")
            with urllib.request.urlopen(req, timeout=120) as resp:
                body = resp.read().decode("utf-8")
            if body.strip().startswith("<html") or "Query ERROR" in body[:500]:
                raise RuntimeError(f"BioMart returned an error page: {body[:300]}")
            return body
        except Exception as exc:
            print(f"attempt {attempt}/{max_retries} failed: {exc}", flush=True)
            if attempt == max_retries:
                raise
            time.sleep(5)


if __name__ == "__main__":
    print("Querying full hsapiens_gene_ensembl table (ensembl_gene_id + external_gene_name)...", flush=True)
    body = fetch()
    lines = body.strip().split("\n")
    print(f"Got {len(lines)-1} human gene rows.", flush=True)
    with open(OUT, "w") as f:
        f.write(body)
    print(f"Wrote {OUT}", flush=True)
