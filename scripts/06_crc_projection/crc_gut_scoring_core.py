#!/usr/bin/env python3
"""
Step 6 gut re-anchor: shared scoring machinery for the primary compute
(PR #26, docs/STEP6_GUT_SCORING_COMPUTE_DESIGN.md, APPROVE after 2 review
rounds). Imported by both the convergence-check driver and the full
665,473-cell driver so the two runs share byte-identical scoring logic.

Normalization note (not specified in the design doc -- standard scanpy
convention, applied explicitly here so it's auditable): `layers['counts']`
is the only confirmed-raw-integer layer in this atlas (`X`'s status is
undocumented -- `has_raw=False`, `X_integer_fraction` came back `None` in
the Phase-I inventory, i.e. not confirmed usable). `adata.X` is therefore
always rebuilt here as `normalize_total(target_sum=1e4)` + `log1p` of
`layers['counts']` before any `score_genes` call -- the standard
scanpy/Seurat-equivalent input `score_genes` expects.

Null-calibration method (per the approved design): for a real gene panel,
bin ALL atlas genes into `N_BINS` expression-detectability strata
(fraction of cells with `counts > 0` in the object being scored -- for the
convergence check this is the 20,000-cell subset itself, for the full run
it's all 665,473 cells, so the strata always match the population actually
being scored). Each of `N_PERM` null draws samples one random gene per
real-signature gene from the same stratum (without replacement within a
single draw; the same gene may appear across different draws or already
be in the real signature -- not excluded, same convention as the Step 4a/
PR#23-24 precedent's `matched_permutation_null`). `scanpy.tl.score_genes`
is called once for the real signature and once per null draw -- both use
the exact same function, so the real score and null scores share one
consistent construction. Primary common-scale value = each cell's
empirical percentile among its own `N_PERM` null draws; secondary
sensitivity value = null-calibrated z-score
`(observed - mean(null)) / std(null)`.
"""
import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

REPO = "/home/zz950/TWEAKR-OncoPlacental"
ATLAS_H5AD = "/home/zz950/DATA/scRNAseq/meta_study/CRC_single_cell_atlas_2025/adata_nmf.h5ad"
GENE_SET_DIR = f"{REPO}/results/06_crc_projection/gut_scoring"

# 13-panel inventory (per docs/STEP6_GUT_SCORING_COMPUTE_DESIGN.md)
REVCSC_PANELS = [
    "revCSC_primary27_full", "revCSC_primary27_minus_CLU",
    "revCSC_primary27_minus_ASS1", "revCSC_primary27_minus_CLU_ASS1",
    "revCSC_extended28_full", "revCSC_extended28_minus_CLU",
    "revCSC_extended28_minus_ASS1", "revCSC_extended28_minus_CLU_ASS1",
]
DFP_PANELS = ["D_Gut-shared", "F_Gut-specific", "F_Colon-specific",
              "F_SI-specific", "P_Gut-specific"]
ALL_PANELS = REVCSC_PANELS + DFP_PANELS

# revCSC(primary,extended) <-> D/F/P comparison pairs, per the locked
# overlap-exclusion contract (docs/STEP6_GUT_SCORING_COMPUTE_DESIGN.md's
# "revCSC scoring inventory" table).
COMPARISON_PAIRS = [
    ("revCSC_primary27_full", "D_Gut-shared"),
    ("revCSC_primary27_full", "P_Gut-specific"),
    ("revCSC_primary27_minus_CLU", "F_Colon-specific"),
    ("revCSC_primary27_minus_ASS1", "F_SI-specific"),
    ("revCSC_primary27_minus_CLU_ASS1", "F_Gut-specific"),
    ("revCSC_extended28_full", "D_Gut-shared"),
    ("revCSC_extended28_full", "P_Gut-specific"),
    ("revCSC_extended28_minus_CLU", "F_Colon-specific"),
    ("revCSC_extended28_minus_ASS1", "F_SI-specific"),
    ("revCSC_extended28_minus_CLU_ASS1", "F_Gut-specific"),
]

N_BINS = 20
SEED = 20260815


def load_panel_ensembl_ids(panel):
    path = f"{GENE_SET_DIR}/{panel}.ensembl.txt"
    with open(path) as f:
        return [l.strip() for l in f if l.strip()]


def load_atlas(n_cells_subset=None, seed=SEED):
    """Loads the full atlas, optionally stratified-subsetting to
    n_cells_subset cells proportional by study_id (used by the
    convergence check; full run passes n_cells_subset=None)."""
    print(f"Loading {ATLAS_H5AD} ...", flush=True)
    t0 = time.time()
    adata = ad.read_h5ad(ATLAS_H5AD)
    print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes in {time.time()-t0:.1f}s", flush=True)

    if n_cells_subset is not None:
        rng = np.random.default_rng(seed)
        study = adata.obs["study_id"].astype(str)
        frac = n_cells_subset / adata.n_obs
        keep_idx = []
        for s, idx in study.groupby(study).groups.items():
            idx = np.asarray(idx)
            n_take = max(1, int(round(len(idx) * frac)))
            n_take = min(n_take, len(idx))
            take = rng.choice(idx, size=n_take, replace=False)
            keep_idx.extend(take.tolist())
        keep_idx = np.asarray(sorted(set(keep_idx)))
        # positional indexing -- obs.index values are labels, not positions
        pos = adata.obs.index.get_indexer(keep_idx)
        adata = adata[pos].copy()
        print(f"Stratified subset (by study_id, seed={seed}): {adata.n_obs} cells "
              f"(target {n_cells_subset})", flush=True)

    # Normalization (see module docstring): rebuild X from layers['counts'],
    # the only confirmed-raw layer in this atlas.
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def compute_detectability(adata):
    """Fraction of cells with counts>0, per gene, from layers['counts'] --
    computed on whatever object is passed (subset or full), so the strata
    always match the population actually being scored."""
    counts = adata.layers["counts"]
    detected = (counts > 0)
    frac = np.asarray(detected.mean(axis=0)).ravel()
    return pd.Series(frac, index=adata.var_names)


def bin_detectability(detectability, n_bins=N_BINS):
    try:
        bins, _ = pd.qcut(detectability, q=n_bins, labels=False, retbins=True, duplicates="drop")
    except ValueError:
        bins, _ = pd.qcut(detectability, q=min(n_bins, detectability.nunique()),
                           labels=False, retbins=True, duplicates="drop")
    return pd.Series(bins, index=detectability.index)


def testable_genes(panel_ensembl_ids, atlas_var_names_set):
    return [g for g in panel_ensembl_ids if g in atlas_var_names_set]


def draw_null_gene_set(real_genes_testable, bin_labels, bin_pools, rng):
    """One null draw: for each real gene, one random gene from the same
    detectability bin, without replacement within this draw."""
    real_bins = bin_labels.loc[real_genes_testable]
    draw = []
    used = set()
    for b, genes_in_bin in real_bins.groupby(real_bins).groups.items():
        n_needed = len(genes_in_bin)
        pool = [g for g in bin_pools.get(b, []) if g not in used]
        if len(pool) >= n_needed:
            picked = rng.choice(pool, size=n_needed, replace=False)
        else:
            # bin exhausted (rare, only for very large panels in a small
            # bin) -- fall back to sampling with replacement from the pool,
            # reported via a warning rather than silently failing
            print(f"  WARNING: detectability bin {b} exhausted "
                  f"({len(pool)} available < {n_needed} needed), sampling with replacement", flush=True)
            picked = rng.choice(bin_pools.get(b, real_genes_testable), size=n_needed, replace=True)
        draw.extend(picked.tolist())
        used.update(picked.tolist())
    return draw


def score_panel(adata, panel_name, panel_ensembl_ids, bin_labels, n_perm, seed=SEED):
    """Returns (percentile: np.ndarray[n_cells], zscore: np.ndarray[n_cells],
    n_testable: int) for one panel."""
    atlas_var_set = set(adata.var_names)
    real_genes = testable_genes(panel_ensembl_ids, atlas_var_set)
    n_testable = len(real_genes)
    if n_testable == 0:
        raise ValueError(f"{panel_name}: 0 testable genes present in atlas")

    bin_pools = {b: idx.tolist() for b, idx in bin_labels.groupby(bin_labels).groups.items()}
    rng = np.random.default_rng(hash((seed, panel_name)) % (2**32))

    t0 = time.time()
    sc.tl.score_genes(adata, gene_list=real_genes, score_name="_tmp_real", use_raw=False)
    real_score = adata.obs["_tmp_real"].values.copy()
    print(f"    {panel_name}: real score computed in {time.time()-t0:.1f}s "
          f"(n_testable={n_testable})", flush=True)

    n_cells = adata.n_obs
    null_scores = np.zeros((n_cells, n_perm), dtype=np.float32)
    t0 = time.time()
    for i in range(n_perm):
        null_genes = draw_null_gene_set(real_genes, bin_labels, bin_pools, rng)
        sc.tl.score_genes(adata, gene_list=null_genes, score_name="_tmp_null", use_raw=False)
        null_scores[:, i] = adata.obs["_tmp_null"].values
        if (i + 1) % 20 == 0:
            elapsed = time.time() - t0
            print(f"    {panel_name}: {i+1}/{n_perm} null draws done "
                  f"({elapsed:.1f}s, {elapsed/(i+1):.2f}s/draw)", flush=True)
    print(f"    {panel_name}: all {n_perm} null draws done in {time.time()-t0:.1f}s", flush=True)

    n_less = (null_scores < real_score[:, None]).sum(axis=1)
    n_equal = (null_scores == real_score[:, None]).sum(axis=1)
    percentile = (n_less + 0.5 * n_equal) / n_perm * 100.0

    null_mean = null_scores.mean(axis=1)
    null_std = null_scores.std(axis=1, ddof=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        zscore = (real_score - null_mean) / null_std
    zscore = np.where(null_std > 0, zscore, np.nan)

    for c in ["_tmp_real", "_tmp_null"]:
        if c in adata.obs.columns:
            del adata.obs[c]

    return percentile.astype(np.float32), zscore.astype(np.float32), n_testable


def score_all_panels(adata, n_perm, panels=None, seed=SEED, checkpoint_dir=None):
    """Runs score_panel for every panel in `panels` (default ALL_PANELS),
    returns a DataFrame indexed by adata.obs_names with
    <panel>_percentile / <panel>_zscore columns. If checkpoint_dir is
    given, each panel's result is written immediately and skipped on
    re-run if already present (job-restart safety for a job this long)."""
    panels = panels or ALL_PANELS
    detectability = compute_detectability(adata)
    bin_labels = bin_detectability(detectability, N_BINS)
    print(f"Detectability strata: {N_BINS} bins over {len(detectability)} genes", flush=True)

    out = pd.DataFrame(index=adata.obs_names)
    n_testable_report = {}
    for panel in panels:
        ckpt_path = f"{checkpoint_dir}/{panel}.n{n_perm}.parquet" if checkpoint_dir else None
        if ckpt_path and os.path.exists(ckpt_path):
            print(f"  {panel}: checkpoint found, loading {ckpt_path}", flush=True)
            df = pd.read_parquet(ckpt_path)
            out[f"{panel}_percentile"] = df[f"{panel}_percentile"].values
            out[f"{panel}_zscore"] = df[f"{panel}_zscore"].values
            n_testable_report[panel] = int(df.attrs.get("n_testable", -1)) if hasattr(df, "attrs") else -1
            continue

        panel_ens = load_panel_ensembl_ids(panel)
        percentile, zscore, n_testable = score_panel(adata, panel, panel_ens, bin_labels, n_perm, seed)
        out[f"{panel}_percentile"] = percentile
        out[f"{panel}_zscore"] = zscore
        n_testable_report[panel] = n_testable

        if ckpt_path:
            os.makedirs(checkpoint_dir, exist_ok=True)
            df = out[[f"{panel}_percentile", f"{panel}_zscore"]].copy()
            df.to_parquet(ckpt_path)
            print(f"  {panel}: checkpointed to {ckpt_path}", flush=True)

    return out, n_testable_report
