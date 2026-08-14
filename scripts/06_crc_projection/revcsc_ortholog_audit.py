#!/usr/bin/env python3
"""
revCSC mouse->human orthology provenance audit, per PR #17 review round 1
(REQUEST_CHANGES): "31 successfully mapped to a human Ensembl ID" only
proves symbol-matching succeeded against the HGNC table, not that the
mouse->human ORTHOLOGY relationship is biologically correct (one-to-one,
not a paralog/false-positive symbol-capitalization match).

Locks the 32 distinct mouse gene symbols already extracted from
CSC_subtype_signatures.ensembl_mapping.tsv (cluster=revCSC) as the input
mouse gene list (this project did not re-derive this list from the Ayyaz
et al. 2019 Nature paper's own supplementary table -- that re-derivation
is out of scope here and flagged as a documented limitation, not silently
assumed done), then queries Ensembl's own computed orthology calls
(Compara, via the public homology REST API -- same "network call to a
public API, run locally" exception as Step 2's query_biotype.py) for each
mouse gene symbol, recording:
  - human ortholog gene symbol(s) + Ensembl ID(s) per Ensembl Compara
  - orthology type: ortholog_one2one / ortholog_one2many /
    ortholog_many2many / NO_ORTHOLOG_FOUND / MOUSE_GENE_NOT_FOUND
  - whether the existing table's naive-uppercase-mapped human symbol
    actually matches Ensembl Compara's own orthology call

Usage: python3 revcsc_ortholog_audit.py <out_dir>
(run locally -- network call, not local-data analysis, same category as
Step 2's Ensembl biotype lookup)
"""
import sys
import os
import json
import time
import urllib.request
import urllib.error

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(OUT_DIR, exist_ok=True)

# The 32 distinct mouse gene symbols extracted from
# CSC_subtype_signatures.ensembl_mapping.tsv, cluster=revCSC (as pulled
# directly from Argos this session -- not re-derived from the original
# paper's supplementary table, see docstring).
MOUSE_SYMBOLS = [
    "Pyy", "Ctla2a", "Clu", "Acta1", "F3", "Cldn4", "Basp1", "Marcksl1",
    "Sfn", "Tmsb4x", "Ankrd1", "Sox4", "Pmepa1", "Anxa1", "Prdx2",
    "Tuba1a", "Sprr1a", "Krt18", "Fn1", "Ctse", "Ctsl", "Tpm1",
    "Tnfrsf12a", "Ass1", "Ly6a", "Areg", "Ccn1", "Ccn2", "Ecm1", "Sox9",
    "Cd44", "Itga2",
]

# The existing table's naive-uppercase-mapped human symbol/Ensembl ID, for
# cross-checking against Ensembl Compara's own orthology call.
EXISTING_TABLE_MAP = {
    "Pyy": ("PYY", "ENSG00000131096"),
    "Ctla2a": ("CTLA2A", None),
    "Clu": ("CLU", "ENSG00000120885"),
    "Acta1": ("ACTA1", "ENSG00000143632"),
    "F3": ("F3", "ENSG00000117525"),
    "Cldn4": ("CLDN4", "ENSG00000189143"),
    "Basp1": ("BASP1", "ENSG00000176788"),
    "Marcksl1": ("MARCKSL1", "ENSG00000175130"),
    "Sfn": ("SFN", "ENSG00000175793"),
    "Tmsb4x": ("TMSB4X", "ENSG00000205542"),
    "Ankrd1": ("ANKRD1", "ENSG00000148677"),
    "Sox4": ("SOX4", "ENSG00000124766"),
    "Pmepa1": ("PMEPA1", "ENSG00000124225"),
    "Anxa1": ("ANXA1", "ENSG00000135046"),
    "Prdx2": ("PRDX2", "ENSG00000167815"),
    "Tuba1a": ("TUBA1A", "ENSG00000167552"),
    "Sprr1a": ("SPRR1A", "ENSG00000169474"),
    "Krt18": ("KRT18", "ENSG00000111057"),
    "Fn1": ("FN1", "ENSG00000115414"),
    "Ctse": ("CTSE", "ENSG00000196188"),
    "Ctsl": ("CTSL", "ENSG00000135047"),
    "Tpm1": ("TPM1", "ENSG00000140416"),
    "Tnfrsf12a": ("TNFRSF12A", "ENSG00000006327"),
    "Ass1": ("ASS1", "ENSG00000130707"),
    "Ly6a": ("LY6A", "ENSG00000291309"),
    "Areg": ("AREG", "ENSG00000109321"),
    "Ccn1": ("CCN1", "ENSG00000145386"),
    "Ccn2": ("CCN2", "ENSG00000118523"),
    "Ecm1": ("ECM1", "ENSG00000143369"),
    "Sox9": ("SOX9", "ENSG00000125398"),
    "Cd44": ("CD44", "ENSG00000026508"),
    "Itga2": ("ITGA2", "ENSG00000164171"),
}


def query_homology(mouse_symbol, max_retries=6):
    url = (
        f"https://rest.ensembl.org/homology/symbol/mus_musculus/{mouse_symbol}"
        f"?type=orthologues;target_species=homo_sapiens"
    )
    last_err = "unknown"
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(
                url, headers={"Content-Type": "application/json"}
            )
            with urllib.request.urlopen(req, timeout=20) as resp:
                body = resp.read()
            data = json.loads(body)
            return data
        except urllib.error.HTTPError as e:
            if e.code == 400:
                return {"error": "MOUSE_GENE_NOT_FOUND", "detail": str(e)}
            # 500/503 observed to be transient server-side load, not a real
            # data-not-found signal -- retry with a flat delay.
            last_err = f"HTTP {e.code}"
            print(f"    ({mouse_symbol}: HTTP {e.code}, retry {attempt+1}/{max_retries})", flush=True)
            time.sleep(3)
        except Exception as e:
            last_err = str(e)
            print(f"    ({mouse_symbol}: {e}, retry {attempt+1}/{max_retries})", flush=True)
            time.sleep(3)
    return {"error": "QUERY_FAILED", "detail": last_err}


def main():
    rows = []
    for sym in MOUSE_SYMBOLS:
        print(f"Querying {sym}...", flush=True)
        data = query_homology(sym)
        existing_sym, existing_ens = EXISTING_TABLE_MAP[sym]

        if "error" in data:
            rows.append({
                "mouse_symbol": sym,
                "orthology_type": data["error"],
                "human_ortholog_symbol": "",
                "human_ensembl_id": "",
                "matches_existing_table": "N/A",
                "existing_table_symbol": existing_sym,
                "existing_table_ensembl": existing_ens or "",
                "source": "Ensembl Compara (homology REST API)",
            })
            time.sleep(1)
            continue

        homology_entries = []
        try:
            homology_entries = data["data"][0]["homologies"]
        except (KeyError, IndexError):
            pass

        if not homology_entries:
            rows.append({
                "mouse_symbol": sym,
                "orthology_type": "NO_ORTHOLOG_FOUND",
                "human_ortholog_symbol": "",
                "human_ensembl_id": "",
                "matches_existing_table": "N/A" if existing_ens is None else "NO (Compara found none)",
                "existing_table_symbol": existing_sym,
                "existing_table_ensembl": existing_ens or "",
                "source": "Ensembl Compara (homology REST API)",
            })
            time.sleep(1)
            continue

        for h in homology_entries:
            target_id = h["target"]["id"]
            otype = h["type"]
            matches = "N/A"
            if existing_ens is not None:
                matches = "YES" if target_id == existing_ens else "DIFFERENT_TARGET"
            rows.append({
                "mouse_symbol": sym,
                "orthology_type": otype,
                "human_ortholog_symbol": "",  # Compara returns Ensembl ID, not symbol, via this endpoint
                "human_ensembl_id": target_id,
                "matches_existing_table": matches,
                "existing_table_symbol": existing_sym,
                "existing_table_ensembl": existing_ens or "",
                "source": "Ensembl Compara (homology REST API)",
            })
        time.sleep(1)

    out_path = f"{OUT_DIR}/revcsc_mouse_human_ortholog_audit.tsv"
    with open(out_path, "w") as f:
        cols = ["mouse_symbol", "orthology_type", "human_ensembl_id",
                "matches_existing_table", "existing_table_symbol",
                "existing_table_ensembl", "source"]
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"\nWrote {out_path}")

    # Summary
    from collections import Counter
    per_gene_types = {}
    for r in rows:
        per_gene_types.setdefault(r["mouse_symbol"], []).append(r["orthology_type"])
    type_summary = Counter()
    for sym, types in per_gene_types.items():
        if len(types) == 1:
            type_summary[types[0]] += 1
        else:
            type_summary["multiple_calls_" + "_or_".join(sorted(set(types)))] += 1
    print("\nPer-gene orthology-call summary:")
    for k, v in type_summary.most_common():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
