#!/usr/bin/env python3
"""
Step 8: D/F/P/revCSC scoring for the Step 7 scRNA-seq cohorts -- shared
loader/coverage-check machinery, per docs/STEP8_CLIM_SCRNA_SCORING_
DESIGN.md (PR #37, APPROVE after 3 review rounds). Imported by
clim_scrna_scoring_driver.py. Reuses crc_gut_scoring_core.py's scoring
machinery (score_all_panels etc.) unchanged -- this module only builds
correctly-axis-aligned `adata` objects matching that module's internal
contract (layers["counts"] = raw counts on the canonical gene axis, X =
normalize_total(1e4)+log1p of that same layer).

Locked scoring populations (design doc "Cell/sample-scoring-population
contract"): GSE231559 primary (9 CLiM + 6 primary-tumor samples, 15
total) + GSE231559 normal (11 samples, separate calibration pass);
GSE285990 (all 10 P01_LM-P10_LM samples); GSE225857 non-immune fraction
only (GSM7058755, 41,892 cells) -- the 196,473-cell immune fraction
(GSM7058754) is explicitly OUT OF SCOPE, not loaded here at all.

Gene-axis canonicalization contract (design doc, rounds 1-2 fixes,
locked in the correct order): canonicalize each reference to bare
Ensembl IDs (or map symbols->Ensembl for GSE225857) FIRST, assert no
within-reference duplicate collisions, THEN intersect across references
-- one shared function used identically by the coverage check and the
compute-time loader here, not two independent reimplementations.
"""
import os
import re
import gzip
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

DATA_ROOT = "/home/zz950/DATA"
HGNC_MAP_PATH = f"{DATA_ROOT}/1.Databases/HGNC_gene_id_mapping/processed/v0.1/hgnc_symbol_ensembl_map.tsv"

# GSE231559 sample -> paper classification, per Step 7's real, exact
# reconstruction (results/07_clim_external_data/GSE231559_inventory.tsv,
# PR #36) -- reused directly, not re-derived.
GSE231559_CLIM_GSMS = [
    "GSM7290761", "GSM7290767", "GSM7290775", "GSM7290778", "GSM7290779",
    "GSM7290781", "GSM7290782", "GSM7290783", "GSM7290785",
]  # 9 CLiM (liver tumor)
GSE231559_PRIMARY_GSMS = [
    "GSM7290763", "GSM7290769", "GSM7290772", "GSM7290773", "GSM7290774",
    "GSM7290777",
]  # 6 primary CRC (colon tumor)
GSE231559_NORMAL_GSMS = [
    "GSM7290760", "GSM7290762", "GSM7290764", "GSM7290765", "GSM7290766",
    "GSM7290768", "GSM7290770", "GSM7290771", "GSM7290776", "GSM7290780",
    "GSM7290784",
]  # 11 paired-normal
GSE231559_PRIMARY_POPULATION_GSMS = GSE231559_CLIM_GSMS + GSE231559_PRIMARY_GSMS  # 15

GSE285990_GSMS = [f"GSM8714{595+i}" for i in range(10)]  # P01_LM-P10_LM
GSE285990_PATIENTS = [f"P{i:02d}_LM" for i in range(1, 11)]

ENSEMBL_BARE_PATTERN = re.compile(r"^ENSG[0-9]{11}$")


def canonicalize_ensembl_ids(gene_ids):
    """Strip Ensembl version suffixes to bare IDs, assert no within-
    reference duplicate collisions. Returns a pd.Index of bare IDs,
    same length/order as input."""
    bare = pd.Index([str(g).split(".")[0] for g in gene_ids])
    if bare.has_duplicates:
        dupes = bare[bare.duplicated()].unique().tolist()
        raise ValueError(
            f"canonicalize_ensembl_ids: stripping version suffixes produces "
            f"{len(dupes)} duplicate bare-Ensembl-ID collisions "
            f"(e.g. {dupes[:5]}) -- refusing to silently corrupt the "
            f"gene-to-column mapping."
        )
    return bare


def map_symbols_to_ensembl(symbols, hgnc_map_path=HGNC_MAP_PATH):
    """Maps gene symbols to Ensembl IDs via this project's own existing
    HGNC table (built in Step 2, reused as-is). Drops any symbol with no
    unambiguous mapping (logged with counts, not silently dropped).
    Asserts no duplicate Ensembl IDs result from the mapping.

    Returns (kept_positions: np.ndarray of original-array positions to
    keep, ensembl_ids: pd.Index of the corresponding bare Ensembl IDs,
    stats: dict)."""
    hgnc = pd.read_csv(hgnc_map_path, sep="\t")
    sym2ens = {}
    for sym, ens in zip(hgnc["symbol"], hgnc["ensembl_id"]):
        if isinstance(ens, str) and ens:
            sym2ens.setdefault(sym, ens)  # first mapping wins; HGNC table
            # itself is the authority, not re-adjudicated here

    symbols = list(symbols)
    kept_positions, ensembl_ids = [], []
    n_unmapped = 0
    for i, s in enumerate(symbols):
        e = sym2ens.get(s)
        if e:
            kept_positions.append(i)
            ensembl_ids.append(e)
        else:
            n_unmapped += 1

    ensembl_ids = pd.Index(ensembl_ids)
    stats = {
        "n_input_symbols": len(symbols),
        "n_mapped": len(ensembl_ids),
        "n_unmapped": n_unmapped,
    }

    if ensembl_ids.has_duplicates:
        dupes = ensembl_ids[ensembl_ids.duplicated()].unique().tolist()
        raise ValueError(
            f"map_symbols_to_ensembl: symbol->Ensembl mapping produces "
            f"{len(dupes)} duplicate Ensembl-ID collisions (e.g. "
            f"{dupes[:5]}) -- two input symbols mapped onto the same gene. "
            f"Refusing to silently sum/average/pick-first."
        )

    return np.asarray(kept_positions), ensembl_ids, stats


def _read_mtx_triplet(barcodes_path, features_path, matrix_path):
    """Loads one 10x-style MTX triplet. Returns (gene_ids: list[str],
    barcodes: list[str], matrix: scipy.sparse genes x cells CSR)."""
    with gzip.open(features_path, "rt") as f:
        gene_ids = [line.split("\t")[0] for line in f]
    with gzip.open(barcodes_path, "rt") as f:
        barcodes = [line.strip() for line in f]
    mat = sio.mmread(matrix_path).tocsr()
    return gene_ids, barcodes, mat


def load_gse231559_population(gsm_list, population_name, extract_dir=None):
    """Loads a GSE231559 sample subset (the 15-sample primary population
    or the 11-sample normal population) into the crc_gut_scoring_core
    internal contract.

    Round-1/round-2 review fix (design doc): GSE231559 splits into two
    gene-reference batches (33,538-gene SC10/SC143 batch; 36,601-gene
    SC173/SC216 batch). Locked order: canonicalize each sample's native
    axis to bare Ensembl FIRST, assert no within-sample duplicate
    collisions, THEN intersect the canonical axes across all samples in
    this population -- not the reverse order (an earlier design draft
    had this backwards, corrected in round 2 review)."""
    if extract_dir is None:
        extract_dir = f"{DATA_ROOT}/scRNAseq/GSE231559/raw/extracted"
    files = os.listdir(extract_dir)

    per_sample = {}
    for gsm in gsm_list:
        matches = [f for f in files if f.startswith(gsm + "_")]
        b = next(f for f in matches if f.endswith("barcodes.tsv.gz"))
        feat = next(f for f in matches if f.endswith("features.tsv.gz"))
        mtx = next(f for f in matches if f.endswith("matrix.mtx.gz"))
        gene_ids, barcodes, mat = _read_mtx_triplet(
            os.path.join(extract_dir, b), os.path.join(extract_dir, feat),
            os.path.join(extract_dir, mtx),
        )
        canonical = canonicalize_ensembl_ids(gene_ids)
        per_sample[gsm] = {
            "canonical_genes": canonical, "matrix": mat,
            "barcodes": [f"{gsm}_{bc}" for bc in barcodes],
        }
        print(f"  [{population_name}] {gsm}: {mat.shape[1]} cells x "
              f"{mat.shape[0]} genes (canonical axis)", flush=True)

    # Intersect canonical axes across all samples in THIS population
    # (not just across the two reference batches globally -- a
    # population might draw from only one batch, in which case the
    # intersection is trivially that batch's own full axis).
    common_genes = None
    for gsm, d in per_sample.items():
        s = set(d["canonical_genes"])
        common_genes = s if common_genes is None else (common_genes & s)
    common_genes = pd.Index(sorted(common_genes))
    print(f"  [{population_name}] intersected gene axis: {len(common_genes)} "
          f"genes across {len(gsm_list)} samples", flush=True)

    # Subset + concatenate onto the common axis.
    mats, all_barcodes = [], []
    for gsm in gsm_list:
        d = per_sample[gsm]
        gene_pos = pd.Index(d["canonical_genes"]).get_indexer(common_genes)
        sub = d["matrix"][gene_pos, :]  # genes(common) x cells
        mats.append(sub)
        all_barcodes.extend(d["barcodes"])
    full = sp.hstack(mats, format="csr").T.tocsr()  # cells x genes

    adata = ad.AnnData(
        X=full,
        obs=pd.DataFrame(index=all_barcodes),
        var=pd.DataFrame(index=common_genes),
    )
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print(f"  [{population_name}] final: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)
    return adata


def load_gse285990_population(extract_dir=None):
    """Loads all 10 GSE285990 P01_LM-P10_LM samples -- uniform 37,487-gene
    reference (confirmed in Step 7 inventory, re-verified in Step 8
    design), no cross-sample reference split unlike GSE231559, but still
    canonicalized via the same shared function for consistency."""
    if extract_dir is None:
        extract_dir = f"{DATA_ROOT}/scRNAseq/GSE285990/raw"

    per_sample = {}
    for gsm, patient in zip(GSE285990_GSMS, GSE285990_PATIENTS):
        b = f"{extract_dir}/{gsm}_{patient}_barcodes.tsv.gz"
        feat = f"{extract_dir}/{gsm}_{patient}_features.tsv.gz"
        mtx = f"{extract_dir}/{gsm}_{patient}_matrix.mtx.gz"
        gene_ids, barcodes, mat = _read_mtx_triplet(b, feat, mtx)
        canonical = canonicalize_ensembl_ids(gene_ids)
        per_sample[gsm] = {
            "canonical_genes": canonical, "matrix": mat,
            "barcodes": [f"{gsm}_{bc}" for bc in barcodes],
        }
        print(f"  [GSE285990] {gsm} ({patient}): {mat.shape[1]} cells x "
              f"{mat.shape[0]} genes", flush=True)

    common_genes = None
    for gsm, d in per_sample.items():
        s = set(d["canonical_genes"])
        common_genes = s if common_genes is None else (common_genes & s)
    common_genes = pd.Index(sorted(common_genes))
    print(f"  [GSE285990] intersected gene axis: {len(common_genes)} genes", flush=True)

    mats, all_barcodes = [], []
    for gsm in GSE285990_GSMS:
        d = per_sample[gsm]
        gene_pos = pd.Index(d["canonical_genes"]).get_indexer(common_genes)
        mats.append(d["matrix"][gene_pos, :])
        all_barcodes.extend(d["barcodes"])
    full = sp.hstack(mats, format="csr").T.tocsr()

    adata = ad.AnnData(
        X=full,
        obs=pd.DataFrame(index=all_barcodes),
        var=pd.DataFrame(index=common_genes),
    )
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print(f"  [GSE285990] final: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)
    return adata


def load_gse225857_nonimmune_population(raw_dir=None):
    """Loads GSE225857's non-immune fraction only (GSM7058755, 41,892
    cells) -- the ONLY GSE225857 population in scope per the design's
    locked cell/sample-scoring-population contract; the immune fraction
    (GSM7058754) is never loaded here.

    Gene axis: symbols, mapped to Ensembl via this project's own HGNC
    table (map_symbols_to_ensembl), not the Ensembl-native path used for
    GSE231559/GSE285990."""
    if raw_dir is None:
        raw_dir = f"{DATA_ROOT}/scRNAseq/GSE225857/raw"
    counts_path = f"{raw_dir}/GSM7058755_non_immune_counts.txt.gz"

    print("  [GSE225857 non-immune] reading counts matrix "
          "(genes x cells TXT)...", flush=True)
    df = pd.read_csv(counts_path, sep="\t", index_col=0)
    symbols = df.index.tolist()
    cell_barcodes = [f"GSM7058755_{c}" for c in df.columns]
    print(f"  [GSE225857 non-immune] raw matrix: {df.shape[1]} cells x "
          f"{df.shape[0]} genes (symbols)", flush=True)

    kept_pos, ensembl_ids, stats = map_symbols_to_ensembl(symbols)
    print(f"  [GSE225857 non-immune] symbol->Ensembl mapping: "
          f"{stats['n_mapped']}/{stats['n_input_symbols']} mapped, "
          f"{stats['n_unmapped']} dropped (no unambiguous mapping)", flush=True)

    mat = sp.csr_matrix(df.values[kept_pos, :])  # genes(mapped) x cells

    adata = ad.AnnData(
        X=mat.T.tocsr(),  # cells x genes
        obs=pd.DataFrame(index=cell_barcodes),
        var=pd.DataFrame(index=ensembl_ids),
    )
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print(f"  [GSE225857 non-immune] final: {adata.n_obs} cells x "
          f"{adata.n_vars} genes", flush=True)
    return adata, stats


def coverage_check(adata, population_name, panels, load_panel_ensembl_ids_fn):
    """Coverage check reusing the exact same canonical gene axis the
    adata object was already built on (adata.var_names) -- per the
    design's locked contract, the SAME canonicalization the compute-time
    loader used, not an independent reimplementation."""
    atlas_genes = set(adata.var_names)
    rows = []
    for panel in panels:
        panel_genes_raw = load_panel_ensembl_ids_fn(panel)
        panel_genes = canonicalize_ensembl_ids(panel_genes_raw)
        n_testable = sum(1 for g in panel_genes if g in atlas_genes)
        rows.append({
            "population": population_name, "panel": panel,
            "n_panel_genes": len(panel_genes), "n_testable": n_testable,
        })
    return pd.DataFrame(rows)
