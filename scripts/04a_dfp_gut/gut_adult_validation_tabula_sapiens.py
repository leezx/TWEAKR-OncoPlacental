#!/usr/bin/env python3
"""
Step 4a external adult-expression audit, Checks 1+3 (Tabula Sapiens).
Per docs/STEP4A_GUT_ADULT_VALIDATION.md (PR #23, APPROVE after 3 review
rounds -- design doc has the full rationale for every choice below).

Check 1: Large_Intestine, epithelial-restricted, organ-matched against
F_Colon-developmental only. Primary evidence = direct donor-aware
cross-donor-consistent CPM>=1 flags (reused unchanged from Step 5's
tier2_validation.py contract, PR #14/#15). Secondary evidence = a
pre-registered (not optional, round-2 fix) permutation layer, reusing
background_detection_rate.py's matched-null machinery but with BOTH the
detectability-matching covariate AND the hit matrix restricted to
epithelial donor x cell-type pseudobulks only (round-2 fix -- using the
whole-organ pool here would silently pull the matching estimand back to
a whole-organ question the moment the primary hypothesis is
epithelial-restricted). Hypothesis unit = one empirical p per eligible
epithelial cell type for the WHOLE F_Colon-developmental panel's
aggregate detection rate (verified directly against
background_detection_rate.py -- NOT one per gene x cell-type).

Check 3: all 5 Tabula Sapiens organs, no lineage restriction, D_Colon-
shared / D_SI-shared / (P_Colon-specific union P_SI-specific) --
construction-consistency audit for D/P (they inherit whole-body
adult_excluded evidence from P_developmental's own construction) but
genuinely independent for Tabula Sapiens specifically (never touched any
part of either D/F/P construction). Direct gene x cell-type x donor
flags ONLY -- reuses Step 5's actual primary contract, no permutation
(round-1 correction: the original whole-body check never had a
permutation layer; inventing one here would be a new, undeclared
statistical design).

Both checks always report: significant/flagged = adult-expression red
flag worth manual scrutiny for F-arm (NOT proof of a false positive);
non-flagged = inconclusive, NOT proof of adult-negativity. D/P flags
here are additionally a genuinely-independent held-out check (unlike
GTEx's construction-consistency framing for D/P).

Gene-ID mapping: GCA var_name -> authoritative_symbol (BioMart, from PR
#21's var_id_map.tsv), matched against Tabula Sapiens' own canonical
var_names (symbol-only alignment, no Ensembl IDs in this dataset).
Unresolved BioMart symbol or symbol absent from TS = NOT_TESTABLE,
reported explicitly, never silently folded into "adult not detected."

Usage: python3 gut_adult_validation_tabula_sapiens.py <out_dir>
(run on Argos, argos-codex env)
"""
import os
import sys
import zipfile
import tempfile
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

ROOT = "/home/zz950/TWEAKR-OncoPlacental"
TS_DIR = "/home/zz950/DATA/1.Databases/TabulaSapiens/raw"
DFP_GUT_DIR = f"{ROOT}/results/04a_dfp_gut/dfp_gut_gene_sets"
VAR_ID_MAP = f"{ROOT}/results/04a_dfp_gut/gene_id_audit/var_id_map.tsv"
OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/results/04a_dfp_gut/adult_validation"
os.makedirs(OUT_DIR, exist_ok=True)

ORGANS = ["Liver", "Skin", "Spleen", "Thymus", "Large_Intestine"]
MIN_CELLS = 20
NOT_DETECTED_FLOOR = 1.0
N_PERM = 500
N_BINS = 10
SEED = 20260815  # fixed, stated for reproducibility, same convention as PR #22


def load_gene_list(path):
    with open(path) as f:
        return [l.strip() for l in f if l.strip()]


def load_var_id_map():
    return pd.read_csv(VAR_ID_MAP, sep="\t").set_index("var_name")


def resolve_genes_ts(var_names, var_map, ts_var_names):
    """TS-specific NOT_TESTABLE rule (unchanged from round-1): symbol-only
    alignment, authoritative_symbol unresolved or absent from TS = NOT_TESTABLE."""
    ts_set = set(ts_var_names)
    rows = []
    for vn in var_names:
        if vn not in var_map.index:
            rows.append({"var_name": vn, "resolved_symbol": None, "testable": False})
            continue
        rec = var_map.loc[vn]
        auth_symbol = rec["authoritative_symbol"] if rec["found_in_biomart"] else None
        if auth_symbol is not None and auth_symbol in ts_set:
            rows.append({"var_name": vn, "resolved_symbol": auth_symbol, "testable": True})
        else:
            rows.append({"var_name": vn, "resolved_symbol": None, "testable": False})
    return pd.DataFrame(rows).set_index("var_name")


def load_organ(organ):
    path = f"{TS_DIR}/TS_{organ}.h5ad.zip"
    with zipfile.ZipFile(path) as z:
        h5ad_name = [n for n in z.namelist() if n.endswith(".h5ad")][0]
        with tempfile.TemporaryDirectory() as tmp:
            z.extract(h5ad_name, tmp)
            adata = ad.read_h5ad(os.path.join(tmp, h5ad_name))
    return adata


def build_donor_celltype_pseudobulk(adata, cell_mask=None):
    """Same contract as Step 5's tier2_validation.py. cell_mask, if given,
    restricts to a boolean subset of obs (e.g. compartment=='epithelial')
    before building groups -- this is how Check 1's epithelial restriction
    is applied, not a post-hoc filter of the resulting pseudobulk."""
    if cell_mask is not None:
        adata = adata[cell_mask].copy()
    counts = adata.layers["raw_counts"]
    donor = adata.obs["donor"].astype(str).values
    celltype = adata.obs["cell_ontology_class"].astype(str).values
    group = pd.Series(donor).str.cat(pd.Series(celltype), sep="||").values
    uniq_groups = sorted(set(group))
    group_idx = {g: i for i, g in enumerate(uniq_groups)}
    row_idx = np.arange(len(group))
    col_idx = np.array([group_idx[g] for g in group])
    indicator = sparse.csr_matrix(
        (np.ones(len(group)), (row_idx, col_idx)),
        shape=(len(group), len(uniq_groups)),
    )
    pb = indicator.T @ counts
    pb = np.asarray(pb.todense()) if sparse.issparse(pb) else np.asarray(pb)
    pb_df = pd.DataFrame(pb.T, index=adata.var_names, columns=uniq_groups)
    n_before = pb_df.shape[0]
    pb_df = pb_df.groupby(pb_df.index).sum()
    n_dupes = n_before - pb_df.shape[0]
    if n_dupes > 0:
        print(f"  collapsed {n_dupes} duplicate gene-symbol rows by summing", flush=True)
    n_cells = pd.Series(group).value_counts().reindex(uniq_groups).values
    meta = pd.DataFrame({
        "sample": uniq_groups,
        "donor": [g.split("||")[0] for g in uniq_groups],
        "cell_type": [g.split("||")[1] for g in uniq_groups],
        "n_cells": n_cells,
    })
    return pb_df, meta


def gene_celltype_report(pb_cpm, meta, genes, min_cells=MIN_CELLS, floor=NOT_DETECTED_FLOOR):
    """Direct gene x cell-type x donor flags -- exact Step 5 primary contract."""
    eligible = meta[meta.n_cells >= min_cells]
    rows = []
    for ct, sub in eligible.groupby("cell_type"):
        samples = sub["sample"].values
        n_donors = len(samples)
        if n_donors == 0:
            continue
        sub_cpm = pb_cpm.loc[[g for g in genes if g in pb_cpm.index], samples]
        for gene in sub_cpm.index:
            vals = sub_cpm.loc[gene]
            detected = vals >= floor
            n_detected = int(detected.sum())
            frac = n_detected / n_donors
            cross_donor_consistent = bool(n_donors >= 2 and n_detected == n_donors)
            rows.append({
                "cell_type": ct, "gene": gene,
                "n_eligible_donors": n_donors, "n_donors_detected": n_detected,
                "donor_detection_fraction": round(frac, 3),
                "median_cpm": round(float(vals.median()), 2),
                "max_cpm": round(float(vals.max()), 2),
                "single_donor_evidence": n_donors == 1,
                "cross_donor_consistent_hit": cross_donor_consistent and n_detected > 0,
            })
    return pd.DataFrame(rows)


def compute_hit_matrix_and_detectability(pb_counts, meta, min_cells=MIN_CELLS, floor=NOT_DETECTED_FLOOR):
    """Exact background_detection_rate.py contract, applied to whatever
    pb_counts/meta is passed in -- for Check 1 this is ALREADY
    epithelial-restricted at the pseudobulk-construction stage, so the
    detectability covariate and hit matrix are epithelial-only by
    construction (round-2 fix), not by a post-hoc filter."""
    lib_size = pb_counts.sum(axis=0)
    pb_cpm = pb_counts.div(lib_size, axis=1) * 1e6
    eligible_meta = meta[meta.n_cells >= min_cells]
    eligible_cols = eligible_meta["sample"].values
    detected = pb_cpm[eligible_cols] >= floor
    detectability = detected.mean(axis=1)

    hit_cols = {}
    n_donors_per_ct = {}
    for ct, sub in eligible_meta.groupby("cell_type"):
        cols = sub["sample"].values
        n_donors = len(cols)
        if n_donors == 0:
            continue
        sub_detected = detected[cols]
        n_detected = sub_detected.sum(axis=1)
        hit = (n_detected == n_donors) & (n_donors >= 2) & (n_detected > 0)
        hit_cols[ct] = hit
        n_donors_per_ct[ct] = n_donors
    hit_matrix = pd.DataFrame(hit_cols, index=pb_cpm.index)
    return hit_matrix, detectability, n_donors_per_ct


def matched_permutation_null(hit_matrix, detectability, f_genes, n_perm=N_PERM, n_bins=N_BINS, seed=SEED):
    """Exact background_detection_rate.py contract: ONE empirical p per
    (cell_type) for the whole F panel's aggregate detection rate (round-2
    corrected hypothesis unit -- NOT one per gene x cell-type)."""
    rng = np.random.default_rng(seed)
    f_genes_present = [g for g in f_genes if g in hit_matrix.index]

    try:
        bin_labels, _ = pd.qcut(detectability, q=n_bins, labels=False, retbins=True, duplicates="drop")
    except ValueError:
        bin_labels, _ = pd.qcut(detectability, q=min(n_bins, detectability.nunique()), labels=False, retbins=True, duplicates="drop")
    bin_labels = pd.Series(bin_labels, index=detectability.index)

    f_bins = bin_labels.loc[f_genes_present]
    bin_pools = {b: bin_labels.index[bin_labels == b].tolist() for b in bin_labels.dropna().unique()}

    observed = hit_matrix.loc[f_genes_present].mean(axis=0)

    null_rates = np.zeros((n_perm, hit_matrix.shape[1]))
    for i in range(n_perm):
        panel = []
        for b, genes_in_bin in f_bins.groupby(f_bins).groups.items():
            n_needed = len(genes_in_bin)
            pool = bin_pools.get(b, [])
            if len(pool) >= n_needed:
                draw = rng.choice(pool, size=n_needed, replace=False)
            else:
                draw = rng.choice(pool if len(pool) > 0 else hit_matrix.index.values, size=n_needed, replace=True)
            panel.extend(draw)
        null_rates[i, :] = hit_matrix.loc[panel].mean(axis=0).values

    null_median = np.median(null_rates, axis=0)
    n_ge = (null_rates >= observed.values[None, :]).sum(axis=0)
    empirical_p = (n_ge + 1) / (n_perm + 1)  # North et al. add-one correction

    result = pd.DataFrame({
        "cell_type": hit_matrix.columns,
        "observed_F_flag_rate": observed.values,
        "null_median": null_median,
        "empirical_p_value": empirical_p,
        "n_permutations": n_perm,
        "n_F_genes_matched": len(f_genes_present),
    })
    result["enrichment_vs_null_median"] = (result["observed_F_flag_rate"] / result["null_median"].replace(0, np.nan)).round(2)
    return result.round(4)


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    raw_q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(raw_q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


def run_check1(var_map):
    print("=== Check 1: Large_Intestine, epithelial-restricted, F_Colon-developmental ===", flush=True)
    adata = load_organ("Large_Intestine")
    epi_mask = (adata.obs["compartment"] == "epithelial").values
    n_epi = int(epi_mask.sum())
    print(f"Large_Intestine: {adata.n_obs} total cells, {n_epi} epithelial cells", flush=True)

    ts_var_names = adata.var_names.tolist()
    f_colon = load_gene_list(f"{DFP_GUT_DIR}/F_Colon-developmental.txt")
    resolved = resolve_genes_ts(f_colon, var_map, ts_var_names)
    n_input, n_testable = len(f_colon), int(resolved["testable"].sum())
    print(f"F_Colon-developmental: n_input={n_input}, n_present_in_reference={n_testable}, "
          f"n_not_testable={n_input - n_testable}", flush=True)
    testable_symbols = resolved.loc[resolved["testable"], "resolved_symbol"].tolist()

    pb_counts, meta = build_donor_celltype_pseudobulk(adata, cell_mask=epi_mask)
    del adata
    lib_size = pb_counts.sum(axis=0)
    pb_cpm = pb_counts.div(lib_size, axis=1) * 1e6
    n_eligible = int((meta.n_cells >= MIN_CELLS).sum())
    print(f"{n_eligible} donor x epithelial-cell-type samples pass >={MIN_CELLS} cells "
          f"(of {len(meta)} total)", flush=True)

    # primary: direct flags
    varname_of_symbol = {row["resolved_symbol"]: vn for vn, row in resolved.iterrows() if row["testable"]}
    direct = gene_celltype_report(pb_cpm, meta, testable_symbols)
    direct["var_name"] = direct["gene"].map(varname_of_symbol)
    direct_path = f"{OUT_DIR}/check1_tabula_sapiens_direct_flags.tsv"
    direct.to_csv(direct_path, sep="\t", index=False)
    n_flagged = int(direct.cross_donor_consistent_hit.sum())
    print(f"Check 1 direct evidence: {len(direct)} (gene, epithelial cell_type) pairs tested, "
          f"{n_flagged} cross-donor-consistent adult-expression red flags", flush=True)

    # secondary: pre-registered permutation, epithelial-only covariate + hit matrix
    hit_matrix, detectability, n_donors_per_ct = compute_hit_matrix_and_detectability(pb_counts, meta)

    # BUG FIX (found in PR #24 review): a cell type with only 1 eligible donor
    # can structurally NEVER produce a cross-donor-consistent hit (the flag
    # criterion requires n_donors>=2 by definition) -- for such a column,
    # hit_matrix is identically False for every gene, in BOTH the observed
    # data and every permutation draw, so the resulting "test" is a
    # degenerate p=1.0 with zero information, not a real non-significant
    # result. Including it in the BH-FDR family inflates the hypothesis
    # count with fake tests. Restrict the permutation/FDR family to cell
    # types with >=2 eligible donors -- the only ones where a
    # cross-donor-consistent hit is even structurally possible.
    single_donor_cts = [ct for ct, n in n_donors_per_ct.items() if n < 2]
    testable_cts = [ct for ct, n in n_donors_per_ct.items() if n >= 2]
    print(f"Cell types with only 1 eligible donor (structurally cannot produce a "
          f"cross-donor-consistent hit, EXCLUDED from the permutation/FDR family): "
          f"{single_donor_cts}", flush=True)
    hit_matrix_testable = hit_matrix[testable_cts]

    perm_result = matched_permutation_null(hit_matrix_testable, detectability, testable_symbols)
    perm_result["fdr_bh"] = bh_fdr(perm_result["empirical_p_value"].values).round(4)
    perm_result["nominal_p_lt_0.05"] = perm_result["empirical_p_value"] < 0.05
    perm_result["significant_fdr_0.05"] = perm_result["fdr_bh"] < 0.05
    perm_path = f"{OUT_DIR}/check1_tabula_sapiens_permutation_null.tsv"
    perm_result.to_csv(perm_path, sep="\t", index=False)
    print(f"Check 1 permutation layer ({N_PERM} perms, epithelial-only covariate): "
          f"{len(perm_result)} genuinely-testable (>=2-eligible-donor) epithelial cell types "
          f"tested (BH-FDR family) -- {len(single_donor_cts)} single-donor cell types excluded", flush=True)
    print(perm_result.to_string(index=False), flush=True)

    mapping_row = {"gene_set": "F_Colon-developmental (Check1, epithelial)", "n_input": n_input,
                   "n_present_in_reference": n_testable, "n_not_testable": n_input - n_testable}
    return mapping_row


def run_check3(var_map):
    print("\n=== Check 3: all 5 organs, D/P construction-consistency + independent audit ===", flush=True)
    d_colon = load_gene_list(f"{DFP_GUT_DIR}/D_Colon-shared.txt")
    d_si = load_gene_list(f"{DFP_GUT_DIR}/D_SI-shared.txt")
    p_colon = set(load_gene_list(f"{DFP_GUT_DIR}/P_Colon-specific.txt"))
    p_si = set(load_gene_list(f"{DFP_GUT_DIR}/P_SI-specific.txt"))
    p_union = sorted(p_colon | p_si)
    all_genes_union = sorted(set(d_colon) | set(d_si) | set(p_union))
    print(f"D_Colon-shared={len(d_colon)}, D_SI-shared={len(d_si)}, P union (dedup)={len(p_union)}, "
          f"total unique genes tested={len(all_genes_union)}", flush=True)

    # var_name -> authoritative_symbol resolution is organ-independent (comes from
    # PR #21's var_id_map.tsv / BioMart, not from any TS file) -- resolved ONCE here,
    # and used to map results (keyed by symbol) back to the original var_name for
    # correct gene-set/membership annotation (symbol can differ from var_name, e.g.
    # IGF2 vs IGF2-1 -- annotating by raw symbol string against the var_name-keyed
    # gene lists would silently mis-annotate exactly that class of gene).
    symbol_of = {}
    for vn in all_genes_union:
        if vn in var_map.index:
            rec = var_map.loc[vn]
            if rec["found_in_biomart"]:
                symbol_of[vn] = rec["authoritative_symbol"]
    varname_of_symbol = {s: vn for vn, s in symbol_of.items()}  # 1:1 within this gene union, verified below
    if len(varname_of_symbol) != len(symbol_of):
        print("WARNING: non-unique symbol->var_name mapping within the tested gene union -- "
              "check for a real many-to-one collision before trusting annotation", flush=True)

    reports = []
    testability_rows = []
    for organ in ORGANS:
        print(f"--- Loading {organ} ---", flush=True)
        adata = load_organ(organ)
        ts_var_names = adata.var_names.tolist()
        resolved = resolve_genes_ts(all_genes_union, var_map, ts_var_names)
        testable_symbols = resolved.loc[resolved["testable"], "resolved_symbol"].tolist()
        n_testable = len(testable_symbols)
        testability_rows.append({"organ": organ, "n_input": len(all_genes_union),
                                  "n_present_in_reference": n_testable,
                                  "n_not_testable": len(all_genes_union) - n_testable})

        pb_counts, meta = build_donor_celltype_pseudobulk(adata)
        del adata
        lib_size = pb_counts.sum(axis=0)
        pb_cpm = pb_counts.div(lib_size, axis=1) * 1e6
        report = gene_celltype_report(pb_cpm, meta, testable_symbols)
        report.insert(0, "organ", organ)
        reports.append(report)
        print(f"{organ}: n_testable={n_testable}/{len(all_genes_union)}", flush=True)

    combined = pd.concat(reports, ignore_index=True)
    combined["var_name"] = combined["gene"].map(varname_of_symbol)
    # BUG FIX (PR #24 round-2 review): a mutually-exclusive "D_Colon-shared
    # else D_SI-shared else P_union" gene_set label silently drops
    # membership for any gene in more than one set -- TRIM71 is in BOTH
    # D_Colon-shared and D_SI-shared (the only gene the two share) and
    # would only ever get the "D_Colon-shared" label. Doesn't change the
    # 54-flag headline (TRIM71 isn't among them), but the label alone
    # would silently misrepresent membership for any future
    # region-specific summary. Fixed: explicit non-exclusive membership
    # columns for all four original sets, same pattern already used for
    # P_Colon_specific/P_SI_specific; gene_set kept only as a coarse,
    # documented-as-non-exclusive label for quick grouping.
    combined["in_D_Colon_shared"] = combined["var_name"].isin(d_colon)
    combined["in_D_SI_shared"] = combined["var_name"].isin(d_si)
    combined["in_P_Colon_specific"] = combined["var_name"].isin(p_colon)
    combined["in_P_SI_specific"] = combined["var_name"].isin(p_si)
    combined["gene_set"] = combined["var_name"].apply(
        lambda vn: ("D_Colon-shared" if vn in d_colon else
                    "D_SI-shared" if vn in d_si else "P_union"))
    combined_path = f"{OUT_DIR}/check3_tabula_sapiens_direct_flags.tsv"
    combined.to_csv(combined_path, sep="\t", index=False)
    n_flagged = int(combined.cross_donor_consistent_hit.sum())
    print(f"\nCheck 3: {len(combined)} (organ, cell_type, gene) triples tested, "
          f"{n_flagged} cross-donor-consistent flags", flush=True)

    testability_df = pd.DataFrame(testability_rows)
    testability_path = f"{OUT_DIR}/check3_gene_id_mapping_summary.tsv"
    testability_df.to_csv(testability_path, sep="\t", index=False)
    print(f"Wrote {testability_path}", flush=True)
    return combined


def main():
    var_map = load_var_id_map()
    check1_mapping = run_check1(var_map)
    check3_df = run_check3(var_map)

    mapping_summary_path = f"{OUT_DIR}/tabula_sapiens_gene_id_mapping_summary.tsv"
    pd.DataFrame([check1_mapping]).to_csv(mapping_summary_path, sep="\t", index=False)
    print(f"\nWrote {mapping_summary_path}", flush=True)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
