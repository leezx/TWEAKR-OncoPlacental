#!/usr/bin/env python3
"""
Freeze revCSC_human_FINAL from the Ensembl Compara orthology audit
(revcsc_mouse_human_ortholog_audit.tsv, re-verified/corrected below for the
2 genes that failed all retries in that run and the 1 gene with a wrong
Ensembl ID in the pre-existing table). Per PR #17 review round 1
(REQUEST_CHANGES): freeze this as the one artifact all future scoring and
overlap-audit steps reference -- not the historical
CSC_subtype_signatures.ensembl_mapping.tsv directly.

Pre-registered inclusion rule (revised per PR #17 review round 2 --
REQUEST_CHANGES): the round-1 version of this script put the 28-gene set
(27 one2one + Ly6a) in "primary" because Ly6a's specific human target
(ENSG00000291309) had already tested out as reasonable against real
CRC/NMF/Jaccard data -- reviewer correctly flagged this as
outcome-dependent evidence leaking into the definition of what is supposed
to be an orthology-only-driven primary anchor. Fixed:

  - PRIMARY revCSC human signature = the 27 Compara-confirmed
    ortholog_one2one genes only. Selection is driven purely by orthology
    provenance, with zero input from any downstream CRC/NMF result.
  - EXTENDED/sensitivity variant = the 28-gene set adding Ly6a
    (ortholog_one2many, ambiguous target choice). Reported alongside every
    primary result, not the other way around -- the ambiguous, CRC-data-
    informed gene is the one that needs the sensitivity label, not the
    unambiguous 27.
  - Genes with no Compara-confirmed ortholog are excluded from both.

Documented limitation, unchanged from round 1: the 32 mouse gene symbols
audited here are what was already extracted into the pre-existing
CSC_subtype_signatures.ensembl_mapping.tsv on Argos -- this project has
not independently re-derived the revCSC mouse gene list from the original
paper's own supplementary material gene-by-gene.
"""
import csv

rows = [
    # mouse_symbol, human_symbol, human_ensembl_id, orthology_type, decision, note
    ("Pyy", "PYY", "ENSG00000131096", "ortholog_one2one", "include_primary", ""),
    ("Clu", "CLU", "ENSG00000120885", "ortholog_one2one", "include_primary", ""),
    ("Acta1", "ACTA1", "ENSG00000143632", "ortholog_one2one", "include_primary", ""),
    ("F3", "F3", "ENSG00000117525", "ortholog_one2one", "include_primary", ""),
    ("Basp1", "BASP1", "ENSG00000176788", "ortholog_one2one", "include_primary", ""),
    ("Marcksl1", "MARCKSL1", "ENSG00000175130", "ortholog_one2one", "include_primary", ""),
    ("Sfn", "SFN", "ENSG00000175793", "ortholog_one2one", "include_primary", ""),
    ("Tmsb4x", "TMSB4X", "ENSG00000205542", "ortholog_one2one", "include_primary", ""),
    ("Ankrd1", "ANKRD1", "ENSG00000148677", "ortholog_one2one", "include_primary", ""),
    ("Sox4", "SOX4", "ENSG00000124766", "ortholog_one2one", "include_primary", ""),
    ("Pmepa1", "PMEPA1", "ENSG00000124225", "ortholog_one2one", "include_primary", ""),
    ("Anxa1", "ANXA1", "ENSG00000135046", "ortholog_one2one", "include_primary", ""),
    ("Prdx2", "PRDX2", "ENSG00000167815", "ortholog_one2one", "include_primary", ""),
    ("Tuba1a", "TUBA1A", "ENSG00000167552", "ortholog_one2one", "include_primary", ""),
    ("Krt18", "KRT18", "ENSG00000111057", "ortholog_one2one", "include_primary", ""),
    ("Fn1", "FN1", "ENSG00000115414", "ortholog_one2one", "include_primary", ""),
    ("Ctse", "CTSE", "ENSG00000196188", "ortholog_one2one", "include_primary", ""),
    ("Tpm1", "TPM1", "ENSG00000140416", "ortholog_one2one", "include_primary",
     "re-queried after initial run timed out 6/6 retries; confirmed stable one2one"),
    ("Tnfrsf12a", "TNFRSF12A", "ENSG00000006327", "ortholog_one2one", "include_primary", ""),
    ("Ass1", "ASS1", "ENSG00000130707", "ortholog_one2one", "include_primary", ""),
    ("Areg", "AREG", "ENSG00000109321", "ortholog_one2one", "include_primary", ""),
    ("Ccn2", "CCN2", "ENSG00000118523", "ortholog_one2one", "include_primary", ""),
    ("Ecm1", "ECM1", "ENSG00000143369", "ortholog_one2one", "include_primary", ""),
    ("Sox9", "SOX9", "ENSG00000125398", "ortholog_one2one", "include_primary", ""),
    ("Cd44", "CD44", "ENSG00000026508", "ortholog_one2one", "include_primary", ""),
    ("Itga2", "ITGA2", "ENSG00000164171", "ortholog_one2one", "include_primary", ""),
    ("Ccn1", "CCN1", "ENSG00000142871", "ortholog_one2one", "include_primary",
     "CORRECTED: the pre-existing CSC_subtype_signatures.ensembl_mapping.tsv "
     "mapped Ccn1 to ENSG00000145386, which Ensembl lookup/id confirms is "
     "actually CCNA2 (Cyclin A2) -- an unrelated gene, not CCN1/CYR61. "
     "Ensembl Compara's own one2one call (ENSG00000142871) is used instead, "
     "confirmed by direct lookup: display_name=CCN1, "
     "description='cellular communication network factor 1'."),
    ("Ly6a", "LY6A", "ENSG00000291309", "ortholog_one2many", "include_extended_sensitivity_only",
     "Compara reports ortholog_one2many for mouse Ly6a (part of the "
     "species-specific-expanded Ly6 family); ENSG00000291309 is the specific "
     "human target that also matches the pre-existing table's naive-symbol "
     "mapping AND already tested out as reasonable against real CRC/NMF/"
     "Jaccard data -- per PR #17 review round 2, that downstream-data "
     "agreement is exactly why this gene must NOT be in the primary set: "
     "letting outcome-dependent evidence pick among ambiguous ortholog "
     "targets would compromise the primary anchor's independence. Kept "
     "only in the 28-gene EXTENDED/sensitivity variant, reported alongside "
     "every primary (27-gene) result, not the reverse."),
    ("Ctla2a", "", "", "no_ortholog_found", "exclude",
     "No human Ensembl ID resolves via HGNC-based symbol mapping (the "
     "pre-existing table's own failure) and Ensembl Compara returns no "
     "computed mouse->human ortholog for Ctla2a. Excluded from all scoring."),
    ("Cldn4", "", "", "no_ortholog_found", "exclude",
     "CORRECTED: the pre-existing table mapped Cldn4 to CLDN4/ENSG00000189143 "
     "via naive symbol-case matching, but Ensembl Compara reports ZERO "
     "computed orthologs for mouse Cldn4 (empty homologies list, confirmed "
     "stable across 3 independent re-queries). Symbol-matching succeeded; "
     "orthology is not Compara-confirmed. Excluded from primary scoring."),
    ("Ctsl", "", "", "no_ortholog_found", "exclude",
     "CORRECTED: the pre-existing table mapped Ctsl to CTSL/ENSG00000135047 "
     "via naive symbol-case matching, but Ensembl Compara reports ZERO "
     "computed orthologs for mouse Ctsl (empty homologies list, confirmed "
     "stable across 2 independent re-queries) -- plausible given the "
     "cathepsin-L paralog family's complexity in both species. Excluded "
     "from primary scoring."),
    ("Sprr1a", "", "", "no_ortholog_found", "exclude",
     "CORRECTED: the pre-existing table mapped Sprr1a to SPRR1A/ENSG00000169474 "
     "via naive symbol-case matching, but Ensembl Compara reports ZERO "
     "computed orthologs for mouse Sprr1a (empty homologies list, confirmed "
     "on re-query) -- plausible given the small-proline-rich-protein gene "
     "family's known lineage-specific expansion in both species. Excluded "
     "from primary scoring."),
]

out_path = "revCSC_human_FINAL.tsv"
with open(out_path, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["mouse_symbol", "human_ortholog_symbol", "human_ensembl_id",
                "orthology_type", "inclusion_decision", "notes", "orthology_source"])
    for r in rows:
        w.writerow(list(r) + ["Ensembl Compara (homology REST API), re-verified "
                               "2026-08-14; cross-checked against the pre-existing "
                               "CSC_subtype_signatures.ensembl_mapping.tsv"])

primary = [r for r in rows if r[4] == "include_primary"]
extended = [r for r in rows if r[4] in ("include_primary", "include_extended_sensitivity_only")]
excluded = [r for r in rows if r[4] == "exclude"]

print(f"Total mouse genes audited: {len(rows)}")
print(f"PRIMARY scoring set (Compara-confirmed one2one only): {len(primary)} genes")
print(f"  -> {sorted(r[1] for r in primary)}")
print(f"EXTENDED/sensitivity set (primary + ambiguous one2many Ly6a): {len(extended)} genes")
print(f"Excluded (no Compara-confirmed ortholog): {len(excluded)} genes -> {[r[0] for r in excluded]}")
print(f"\nWrote {out_path}")
