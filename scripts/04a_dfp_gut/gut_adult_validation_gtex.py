#!/usr/bin/env python3
"""
Step 4a external adult-expression audit, Check 2 (GTEx bulk).
Per docs/STEP4A_GUT_ADULT_VALIDATION.md (PR #23, APPROVE after 3 review
rounds -- design doc has the full rationale for every choice below).

Two different scientific questions, kept separate throughout (round-1
correction -- conflating them was the design's original core problem):

  F_Colon-developmental / F_SI-developmental: adult-expression /
  adult-specificity AUDIT. These sets require fetal significantly higher
  than adult (FDR<0.05 & logFC>1), not near-zero adult expression, so a
  detected hit here is a red flag worth manual scrutiny, NOT evidence the
  gene is a false positive; below-floor is reported as-is, not claimed as
  "validated adult-negative."

  D_Colon-shared / D_SI-shared / (P_Colon-specific union P_SI-specific):
  CONSTRUCTION-CONSISTENCY audit. These sets are built from the frozen
  P_developmental_primary84.txt, which itself already required whole-body
  GTEx+HPA adult_excluded evidence (pct_cut=25, quorum=all_but_1, the
  actual frozen Step 4 criterion -- reused here unchanged, not reinvented)
  as part of its own construction. This re-checks whether the current
  frozen gut D/P membership still looks consistent with that same
  evidence, not an independent validation.

Gene-ID mapping contract (round-1/round-2 blockers, frozen before this
script was written -- see design doc "Compute contract" section):
  - Every GCA var_name -> gene_id (Ensembl) -> authoritative_symbol
    (BioMart), from PR #21's own results/04a_dfp_gut/gene_id_audit/var_id_map.tsv.
  - GTEx: TESTABLE if gene_id is found in GTEx's own (already
    version-stripped) ensembl_id column -- no BioMart symbol required for
    this path. Only if that fails does it fall back to authoritative_symbol
    string match against GTEx's symbol-collapsed matrix. NOT_TESTABLE only
    if both fail.
  - GTEx's approved percentile reference distribution
    (gtex.groupby("symbol")[tissue_cols].max(), the exact object already
    calibrated in adult_excluded_percentile_audit.py) is reused completely
    unchanged. Ensembl ID is used ONLY to resolve GCA gene -> correct GTEx
    symbol identity (fixing IGF2-1-vs-IGF2-class mismatches); it never
    builds a new Ensembl-level percentile universe.

Usage: python3 gut_adult_validation_gtex.py <out_dir>
(run on Argos, argos-codex env -- pure pandas, no qsub-specific deps, but
kept on the same qsub-only compute discipline as everything else)
"""
import os
import sys
import pandas as pd
import numpy as np

ROOT = "/home/zz950/TWEAKR-OncoPlacental"
GTEX_TSV = "/home/zz950/DATA/1.Databases/GTEx_v11_median_tpm/processed/v0.1/gtex_v11_median_tpm_clean.tsv"
DFP_GUT_DIR = f"{ROOT}/results/04a_dfp_gut/dfp_gut_gene_sets"
VAR_ID_MAP = f"{ROOT}/results/04a_dfp_gut/gene_id_audit/var_id_map.tsv"
OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/results/04a_dfp_gut/adult_validation"
os.makedirs(OUT_DIR, exist_ok=True)

NOT_DETECTED_FLOOR = 1.0  # TPM, below_bulk_expression_floor -- project convention, see design doc
WHOLE_BODY_PCT_CUT = 25   # the actual frozen P_developmental whole-body criterion, reused unchanged
WHOLE_BODY_ALLOWED_FAIL = 1  # quorum=all_but_1


def load_gene_list(path):
    with open(path) as f:
        return [l.strip() for l in f if l.strip()]


def load_var_id_map():
    df = pd.read_csv(VAR_ID_MAP, sep="\t")
    return df.set_index("var_name")


def detected_percentile(wide_df):
    """Per-tissue: not_detected mask (TPM<floor) and percentile rank
    (0-100, min-rank tie-break) computed ONLY among detected values --
    exact method from adult_excluded_percentile_audit.py, unchanged."""
    not_detected = wide_df < NOT_DETECTED_FLOOR
    pct = pd.DataFrame(index=wide_df.index, columns=wide_df.columns, dtype=float)
    for col in wide_df.columns:
        detected_vals = wide_df.loc[~not_detected[col], col]
        if len(detected_vals) > 0:
            pct.loc[detected_vals.index, col] = detected_vals.rank(pct=True, method="min") * 100
    return not_detected, pct


def resolve_genes(var_names, var_map, ensembl_to_symbol, gtex_symbol_index):
    """Cross-dataset gene-ID mapping contract, GTEx-specific NOT_TESTABLE
    rule (round-2 fix): TESTABLE via Ensembl gene_id match (no BioMart
    symbol required), falling back to authoritative_symbol only if the
    Ensembl lookup fails. Returns a DataFrame indexed by var_name with
    columns: gene_id, resolved_symbol, resolution_method, testable."""
    rows = []
    for vn in var_names:
        if vn not in var_map.index:
            rows.append({"var_name": vn, "gene_id": None, "resolved_symbol": None,
                         "resolution_method": "NOT_IN_VAR_ID_MAP", "testable": False})
            continue
        rec = var_map.loc[vn]
        gene_id = rec["gene_id"]
        auth_symbol = rec["authoritative_symbol"] if rec["found_in_biomart"] else None

        if gene_id in ensembl_to_symbol:
            sym = ensembl_to_symbol[gene_id]
            rows.append({"var_name": vn, "gene_id": gene_id, "resolved_symbol": sym,
                         "resolution_method": "ensembl", "testable": sym in gtex_symbol_index})
        elif auth_symbol is not None and auth_symbol in gtex_symbol_index:
            rows.append({"var_name": vn, "gene_id": gene_id, "resolved_symbol": auth_symbol,
                         "resolution_method": "symbol_fallback", "testable": True})
        else:
            rows.append({"var_name": vn, "gene_id": gene_id, "resolved_symbol": None,
                         "resolution_method": "NOT_TESTABLE", "testable": False})
    return pd.DataFrame(rows).set_index("var_name")


def main():
    print("=== Loading GTEx v11 median TPM ===", flush=True)
    gtex = pd.read_csv(GTEX_TSV, sep="\t")
    tissue_cols = [c for c in gtex.columns if c not in ("ensembl_id", "ensembl_id_versioned", "symbol")]
    print(f"GTEx: {gtex.shape[0]} rows x {len(tissue_cols)} tissues", flush=True)

    # Ensembl ID -> symbol map (pre-dedup, for target identity resolution ONLY)
    ensembl_to_symbol = gtex.drop_duplicates("ensembl_id").set_index("ensembl_id")["symbol"].to_dict()

    # The approved, unchanged, symbol-collapsed reference distribution
    gtex_wide = gtex.groupby("symbol")[tissue_cols].max()
    print(f"GTEx after symbol dedup (max per symbol): {gtex_wide.shape[0]} genes -- "
          f"this exact matrix is the approved percentile universe, unmodified", flush=True)
    not_detected, pct = detected_percentile(gtex_wide)

    var_map = load_var_id_map()

    # ============ F-arm: adult-expression / adult-specificity audit ============
    f_checks = [
        ("F_Colon-developmental", ["Colon_Sigmoid", "Colon_Transverse"], "organ_matched"),
        ("F_SI-developmental", ["Small_Intestine_Terminal_Ileum"], "organ_matched_partial_coverage"),
    ]
    f_rows = []
    mapping_summary_rows = []
    for set_name, cols, coverage_note in f_checks:
        genes = load_gene_list(f"{DFP_GUT_DIR}/{set_name}.txt")
        resolved = resolve_genes(genes, var_map, ensembl_to_symbol, set(gtex_wide.index))
        n_input, n_testable = len(genes), int(resolved["testable"].sum())
        n_ensembl = int((resolved["resolution_method"] == "ensembl").sum())
        n_symbol_fb = int((resolved["resolution_method"] == "symbol_fallback").sum())
        n_not_testable = n_input - n_testable
        print(f"{set_name}: n_input={n_input}, n_mapped_ensembl={n_ensembl}, "
              f"n_mapped_symbol_fallback={n_symbol_fb}, n_present_in_reference={n_testable}, "
              f"n_not_testable={n_not_testable}", flush=True)
        mapping_summary_rows.append({"gene_set": set_name, "n_input": n_input,
                                      "n_mapped_ensembl": n_ensembl, "n_mapped_symbol_fallback": n_symbol_fb,
                                      "n_present_in_reference": n_testable, "n_not_testable": n_not_testable})

        for vn in genes:
            rec = resolved.loc[vn]
            if not rec["testable"]:
                f_rows.append({"gene_set": set_name, "var_name": vn, "resolved_symbol": None,
                               "coverage": coverage_note, "tissue": None, "status": "NOT_TESTABLE"})
                continue
            sym = rec["resolved_symbol"]
            for c in cols:
                if not_detected.loc[sym, c]:
                    status = "below_bulk_expression_floor"
                    pct_val = None
                else:
                    status = "detected"
                    pct_val = round(float(pct.loc[sym, c]), 2)
                f_rows.append({"gene_set": set_name, "var_name": vn, "resolved_symbol": sym,
                               "coverage": coverage_note, "tissue": c, "status": status,
                               "percentile_among_detected": pct_val,
                               "resolution_method": rec["resolution_method"]})

    f_df = pd.DataFrame(f_rows)
    f_path = f"{OUT_DIR}/gtex_F_adult_expression_audit.tsv"
    f_df.to_csv(f_path, sep="\t", index=False)
    print(f"\nWrote {f_path} ({len(f_df)} rows)", flush=True)

    # ============ D/P-arm: construction-consistency audit ============
    d_colon = load_gene_list(f"{DFP_GUT_DIR}/D_Colon-shared.txt")
    d_si = load_gene_list(f"{DFP_GUT_DIR}/D_SI-shared.txt")
    p_colon = set(load_gene_list(f"{DFP_GUT_DIR}/P_Colon-specific.txt"))
    p_si = set(load_gene_list(f"{DFP_GUT_DIR}/P_SI-specific.txt"))
    p_union = sorted(p_colon | p_si)
    print(f"\nD_Colon-shared={len(d_colon)}, D_SI-shared={len(d_si)}, "
          f"P_Colon-specific={len(p_colon)}, P_SI-specific={len(p_si)}, "
          f"P union (unique, deduplicated -- not two independent panels, round-1 fix)={len(p_union)}", flush=True)

    dp_sets = [("D_Colon-shared", d_colon), ("D_SI-shared", d_si), ("P_union", p_union)]
    dp_rows = []
    for set_name, genes in dp_sets:
        resolved = resolve_genes(genes, var_map, ensembl_to_symbol, set(gtex_wide.index))
        n_input, n_testable = len(genes), int(resolved["testable"].sum())
        n_not_testable = n_input - n_testable
        mapping_summary_rows.append({"gene_set": set_name, "n_input": n_input,
                                      "n_mapped_ensembl": int((resolved["resolution_method"] == "ensembl").sum()),
                                      "n_mapped_symbol_fallback": int((resolved["resolution_method"] == "symbol_fallback").sum()),
                                      "n_present_in_reference": n_testable, "n_not_testable": n_not_testable})
        print(f"{set_name}: n_input={n_input}, n_present_in_reference={n_testable}, "
              f"n_not_testable={n_not_testable}", flush=True)

        for vn in genes:
            rec = resolved.loc[vn]
            in_p_colon = vn in p_colon
            in_p_si = vn in p_si
            if not rec["testable"]:
                dp_rows.append({"gene_set": set_name, "var_name": vn, "resolved_symbol": None,
                                "in_P_Colon_specific": in_p_colon, "in_P_SI_specific": in_p_si,
                                "status": "NOT_TESTABLE"})
                continue
            sym = rec["resolved_symbol"]
            # reuse the exact frozen P_developmental whole-body criterion (pct_cut=25, quorum=all_but_1)
            fails = ((~not_detected.loc[sym, tissue_cols]) & (pct.loc[sym, tissue_cols] > WHOLE_BODY_PCT_CUT))
            n_fail = int(fails.sum())
            still_passes = n_fail <= WHOLE_BODY_ALLOWED_FAIL
            detected_pcts = pct.loc[sym, tissue_cols][~not_detected.loc[sym, tissue_cols]].dropna()
            dp_rows.append({
                "gene_set": set_name, "var_name": vn, "resolved_symbol": sym,
                "in_P_Colon_specific": in_p_colon, "in_P_SI_specific": in_p_si,
                "status": "tested", "n_tissues_failing_pct25": n_fail,
                "n_tissues_total": len(tissue_cols),
                "still_passes_whole_body_adult_excluded_pct25_allbut1": still_passes,
                "median_percentile_where_detected": round(float(detected_pcts.median()), 2) if len(detected_pcts) else None,
                "max_percentile_where_detected": round(float(detected_pcts.max()), 2) if len(detected_pcts) else None,
            })

    dp_df = pd.DataFrame(dp_rows)
    dp_path = f"{OUT_DIR}/gtex_DP_construction_consistency_audit.tsv"
    dp_df.to_csv(dp_path, sep="\t", index=False)
    print(f"\nWrote {dp_path} ({len(dp_df)} rows)", flush=True)

    n_still_pass = int(dp_df["still_passes_whole_body_adult_excluded_pct25_allbut1"].fillna(False).sum())
    n_tested = int((dp_df["status"] == "tested").sum())
    print(f"\nD/P construction-consistency: {n_still_pass}/{n_tested} tested genes still pass the frozen "
          f"whole-body adult_excluded(pct_cut=25, quorum=all_but_1) criterion", flush=True)

    mapping_df = pd.DataFrame(mapping_summary_rows)
    mapping_path = f"{OUT_DIR}/gtex_gene_id_mapping_summary.tsv"
    mapping_df.to_csv(mapping_path, sep="\t", index=False)
    print(f"Wrote {mapping_path}", flush=True)


if __name__ == "__main__":
    main()
