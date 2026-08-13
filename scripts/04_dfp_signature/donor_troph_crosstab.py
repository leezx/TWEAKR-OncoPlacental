#!/usr/bin/env python3
"""
Cross-tabulate donor/sample x trophoblast-status for the 4 placental
datasets found eligible in results/04_dfp_signature/replicate_structure_audit.md,
per PR #7 review: having donors overall and having trophoblast+non-trophoblast
cells overall does NOT prove any single donor has both -- if trophoblast
status is confounded with donor identity, a donor-level pseudobulk DE is
not valid, regardless of how many "donors" or "replicates" nominally exist.

Run on Argos (obs-metadata-only for the h5ad files via backed='r'; direct
CSV reads for Greenbaum, which needs an explicit metadata.csv <-> cluster.csv
join since donor_id and cell_type live in different files).

Usage: python donor_troph_crosstab.py
"""
import anndata as ad
import pandas as pd

pd.set_option("display.max_rows", 30)


def report(name, xtab):
    xtab["n_non_troph"] = xtab["n_total"] - xtab["n_troph"]
    xtab["eligible"] = (xtab["n_troph"] > 0) & (xtab["n_non_troph"] > 0)
    print(f"=== {name} ===")
    print(xtab.sort_values("n_total", ascending=False))
    n_eligible = int(xtab["eligible"].sum())
    print(f"n_eligible_donors (both groups present): {n_eligible} / {len(xtab)}")
    print()
    return name, n_eligible, len(xtab)


def h5ad_crosstab(name, path, donor_col, ct_col, troph_values, trim_donor=False):
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs[[donor_col, ct_col]].copy()
    if trim_donor:
        obs[donor_col] = obs[donor_col].astype(str).str.strip()
    obs["is_troph"] = obs[ct_col].isin(troph_values)
    xtab = obs.groupby(donor_col)["is_troph"].agg(["sum", "count"])
    xtab.columns = ["n_troph", "n_total"]
    return report(name, xtab)


def greenbaum_crosstab():
    # NOTE: an earlier version of this script tried to join cluster.csv's
    # NAME to metadata.csv's NAME -- that fails 100% of the time because the
    # two files use incompatible barcode schemes (cluster.csv:
    # "W9_AAACCAACACCTGCCT" vs metadata.csv: "JS34#ACGTCAAGTTGCAATG-1", a
    # different sample-ID system entirely). Fixed: cluster.csv's own NAME
    # already embeds the donor ID as the prefix before the last "_" -- no
    # join needed, parse it directly.
    cluster = pd.read_csv(
        "/home/zz950/DATA/scRNAseq/Greenbaum_NatMed_2024/raw/SCP2601/other/humanplacenta_cluster.csv",
        skiprows=[1],
    )
    cluster["donor"] = cluster["NAME"].str.rsplit("_", n=1).str[0]
    troph = {"vCTB", "STB", "EVT", "EVT-progenitor", "STB-progenitor"}
    cluster["is_troph"] = cluster["cell_type"].isin(troph)
    xtab = cluster.groupby("donor")["is_troph"].agg(["sum", "count"])
    xtab.columns = ["n_troph", "n_total"]
    print(f"Greenbaum: {len(cluster)} annotated cells cover only "
          f"{cluster['donor'].nunique()} of the full metadata.csv's 8 donors "
          f"(the cluster-annotated subset is much smaller than the full "
          f"~36,456-cell RNA matrix, as already flagged in the replicate-"
          f"structure audit).")
    return report("Greenbaum (annotated ~1.9K-cell subset, donor parsed from cluster.csv NAME prefix)", xtab)


def main():
    results = []
    results.append(h5ad_crosstab(
        "Arutyunyan_primary_tissue",
        "/home/zz950/DATA/scRNAseq/Arutyunyan2023_MFI/raw/primary_tissue/adata_all_donors_all_cell_states_UPD_20230307.h5ad",
        "donor", "coarse_annot", {"Trophoblast"},
    ))
    results.append(h5ad_crosstab(
        "Nature2026_scPlacenta_host",
        "/home/zz950/DATA/scRNAseq/2026_human_maternal_fetal_Nature/raw/scPlacenta_host.h5ad",
        "sample_id", "major_class", {"SCT", "VCT", "EVT"},
    ))
    results.append(h5ad_crosstab(
        "VentoTormo_decidua_v3",
        "/home/zz950/DATA/scRNAseq/VentoTormo_Nature_2018/raw/decidua-v3.h5ad",
        "Fetus", "CellType", {"VCT", "EVT", "SCT"},
        trim_donor=True,
    ))
    results.append(greenbaum_crosstab())

    print("=== Summary ===")
    for name, n_eligible, n_total_donors in results:
        print(f"{name}: {n_eligible}/{n_total_donors} donors have both trophoblast and non-trophoblast cells")
    n_datasets_with_ge2_eligible = sum(1 for _, n_e, _ in results if n_e >= 2)
    print(f"\nDatasets with >=2 eligible donors (minimum for any within-dataset "
          f"donor-level contrast at all): {n_datasets_with_ge2_eligible} / {len(results)}")


if __name__ == "__main__":
    main()
