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


def load_atlas(n_cells_subset=None, seed=SEED, to_csc=False):
    """Loads the full atlas, optionally stratified-subsetting to
    n_cells_subset cells proportional by study_id (used by the
    convergence check; full run passes n_cells_subset=None).

    `to_csc=True` converts `adata.X` to CSC (compressed sparse column)
    format after normalization -- profiled empirically (timing probe,
    not guessed): `layers['counts']` (and the X rebuilt from it) is CSR
    (row-oriented) by default, and repeated column-slicing (exactly what
    score_genes_fast does every call: `x[:, gene_pos]`) on a CSR matrix
    is a well-known scipy inefficiency, confirmed here at real scale --
    14.4s/call on CSR vs. 0.9s/call on CSC for a 27-gene panel, 3.2s/call
    on CSC for a 2,173-gene panel, at the full 665,473 cells. The
    one-time conversion cost (~80s) is negligible against the ~1300+
    calls a full 13-panel run makes. Only used for the full run's fast
    path -- the convergence check's real scanpy.tl.score_genes calls are
    left on the default CSR (matches how score_genes is normally used;
    changing that path's format isn't necessary since its total runtime
    was already acceptable and it must stay the untouched reference
    implementation)."""
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

    if to_csc:
        t0 = time.time()
        adata.X = adata.X.tocsc()
        print(f"Converted X to CSC (for fast repeated column-slicing) in "
              f"{time.time()-t0:.1f}s", flush=True)

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


SCORE_GENES_CTRL_SIZE = 50
SCORE_GENES_N_BINS = 25  # scanpy.tl.score_genes's own default -- NOT our
                          # N_BINS=20 detectability strata; these are two
                          # independent binning steps kept deliberately
                          # distinct (see score_genes_fast docstring).


def _nan_means_dense_or_sparse(x, axis):
    from scipy import sparse
    if sparse.issparse(x):
        return np.asarray(x.mean(axis=axis)).ravel()
    return np.nanmean(x, axis=axis)


def precompute_score_genes_bins(adata, gene_pool=None, n_bins=SCORE_GENES_N_BINS):
    """Precomputes scanpy.tl.score_genes's own internal control-gene
    expression-rank binning ONCE per scored population, so it can be
    reused across many score_genes-equivalent calls instead of being
    silently recomputed from scratch inside every call (the actual
    bottleneck profiled empirically: cost scales with n_cells x n_genes
    and is otherwise IDENTICAL for a 27-gene and a 2,192-gene panel,
    confirmed via a timing probe -- 1.6s/9.6s/21.5s per call at
    20k/100k/300k cells, extrapolating to ~45-55s/call at the full
    665,473 cells, i.e. an infeasible ~18h+ for the full 13-panel x 101-
    draw job with the naive per-call approach).

    Exact replica of scanpy 1.11's `_score_genes_bins` (verified by
    reading `scanpy/tools/_score_genes.py` directly on Argos, argos-codex
    env) -- `obs_avg` = mean expression per gene_pool gene, `obs_cut` =
    rank-based bin assignment with `n_items = round(len(obs_avg)/(n_bins-1))`.
    Returns (obs_avg, obs_cut, gene_pool_index)."""
    var_names = adata.var_names
    gene_pool_idx = pd.Index(gene_pool) if gene_pool is not None else var_names
    x = adata.X
    obs_avg = pd.Series(_nan_means_dense_or_sparse(x, axis=0), index=var_names)
    obs_avg = obs_avg.loc[gene_pool_idx]
    obs_avg = obs_avg[np.isfinite(obs_avg)]
    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items
    return obs_avg, obs_cut


def score_genes_fast(adata, gene_list, obs_cut, ctrl_size=SCORE_GENES_CTRL_SIZE, random_state=0,
                      return_control_genes=False):
    """Faithful, faster reimplementation of scanpy.tl.score_genes reusing
    a precomputed `obs_cut` (see precompute_score_genes_bins) instead of
    recomputing the expression-rank binning every call. Reproduces the
    same control-gene-selection algorithm and the same `random_state=0`
    reseed-per-call behavior scanpy's own default uses (so repeated calls
    with the same gene_list are deterministic, matching score_genes'
    behavior exactly) -- validated to reproduce scanpy.tl.score_genes'
    output to floating-point precision on a held-out check before use in
    the full-scale run (see validate_score_genes_fast.py).

    Returns a 1D np.ndarray of per-cell scores (gene_list mean expression
    minus control-gene-pool mean expression)."""
    gene_list = pd.Index(gene_list).intersection(adata.var_names)
    if random_state is not None:
        np.random.seed(random_state)

    control_genes = pd.Index([], dtype="string")
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes = obs_cut[obs_cut == cut].index
        if ctrl_size < len(r_genes):
            r_genes = r_genes.to_series().sample(ctrl_size).index
        r_genes = r_genes.difference(gene_list)
        control_genes = control_genes.union(r_genes)

    var_idx = adata.var_names
    gene_pos = var_idx.get_indexer(gene_list)
    ctrl_pos = var_idx.get_indexer(control_genes)
    x = adata.X
    mean_genes = _nan_means_dense_or_sparse(x[:, gene_pos], axis=1)
    mean_ctrl = _nan_means_dense_or_sparse(x[:, ctrl_pos], axis=1)
    score = np.asarray(mean_genes - mean_ctrl).ravel()
    if return_control_genes:
        return score, control_genes
    return score


def score_panel_fast(adata, panel_name, panel_ensembl_ids, bin_labels, obs_cut, n_perm, seed=SEED):
    """Same contract as score_panel (percentile, zscore, n_testable) but
    using score_genes_fast throughout -- for the full 665,473-cell run,
    where the naive per-call scanpy.tl.score_genes cost is infeasible.
    Numerically validated against score_panel's real-scanpy output before
    being trusted (see validate_score_genes_fast.py)."""
    atlas_var_set = set(adata.var_names)
    real_genes = testable_genes(panel_ensembl_ids, atlas_var_set)
    n_testable = len(real_genes)
    if n_testable == 0:
        raise ValueError(f"{panel_name}: 0 testable genes present in atlas")

    bin_pools = {b: idx.tolist() for b, idx in bin_labels.groupby(bin_labels).groups.items()}
    rng = np.random.default_rng(hash((seed, panel_name)) % (2**32))

    t0 = time.time()
    real_score = score_genes_fast(adata, real_genes, obs_cut)
    print(f"    {panel_name}: real score computed in {time.time()-t0:.1f}s "
          f"(n_testable={n_testable})", flush=True)

    n_cells = adata.n_obs
    null_scores = np.zeros((n_cells, n_perm), dtype=np.float32)
    t0 = time.time()
    for i in range(n_perm):
        null_genes = draw_null_gene_set(real_genes, bin_labels, bin_pools, rng)
        null_scores[:, i] = score_genes_fast(adata, null_genes, obs_cut)
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

    return percentile.astype(np.float32), zscore.astype(np.float32), n_testable


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


def score_all_panels(adata, n_perm, panels=None, seed=SEED, checkpoint_dir=None, fast=False):
    """Runs score_panel (or score_panel_fast if fast=True -- see that
    function's docstring for why/when) for every panel in `panels`
    (default ALL_PANELS), returns a DataFrame indexed by adata.obs_names
    with <panel>_percentile / <panel>_zscore columns. If checkpoint_dir is
    given, each panel's result is written immediately and skipped on
    re-run if already present (job-restart safety for a job this long).

    `fast=True` uses score_genes_fast (precomputed control-gene binning,
    reused across all calls) instead of repeated real
    scanpy.tl.score_genes calls -- numerically validated equivalent
    (validate_score_genes_fast.py), used for the full 665,473-cell run
    where the naive per-call approach is empirically ~18h+ (profiled),
    infeasible. The convergence check (20,000 cells) uses fast=False
    (the real function) since its whole purpose is validating the
    statistical design with the reference implementation."""
    panels = panels or ALL_PANELS
    detectability = compute_detectability(adata)
    bin_labels = bin_detectability(detectability, N_BINS)
    print(f"Detectability strata: {N_BINS} bins over {len(detectability)} genes", flush=True)

    obs_cut = None
    if fast:
        t0 = time.time()
        _, obs_cut = precompute_score_genes_bins(adata)
        print(f"Precomputed score_genes control-gene bins in {time.time()-t0:.1f}s "
              f"(reused across all {len(panels)} panels x {n_perm+1} calls each)", flush=True)

    out = pd.DataFrame(index=adata.obs_names)
    n_testable_report = {}
    for panel in panels:
        ckpt_path = f"{checkpoint_dir}/{panel}.n{n_perm}.parquet" if checkpoint_dir else None
        n_testable_path = f"{checkpoint_dir}/{panel}.n{n_perm}.n_testable.txt" if checkpoint_dir else None
        if ckpt_path and os.path.exists(ckpt_path):
            print(f"  {panel}: checkpoint found, loading {ckpt_path}", flush=True)
            df = pd.read_parquet(ckpt_path)
            out[f"{panel}_percentile"] = df[f"{panel}_percentile"].values
            out[f"{panel}_zscore"] = df[f"{panel}_zscore"].values
            # BUG FIX: DataFrame.attrs does not survive a to_parquet/
            # read_parquet round-trip, so n_testable was previously always
            # lost (silently reported as -1) for any checkpoint-reused
            # panel. Fixed by writing n_testable to a small sidecar file
            # instead of relying on .attrs.
            if n_testable_path and os.path.exists(n_testable_path):
                with open(n_testable_path) as f:
                    n_testable_report[panel] = int(f.read().strip())
            else:
                n_testable_report[panel] = -1
                print(f"  WARNING: {panel} checkpoint has no sidecar n_testable "
                      f"file (pre-fix checkpoint) -- reporting -1, re-derive "
                      f"manually if needed", flush=True)
            continue

        panel_ens = load_panel_ensembl_ids(panel)
        if fast:
            percentile, zscore, n_testable = score_panel_fast(
                adata, panel, panel_ens, bin_labels, obs_cut, n_perm, seed)
        else:
            percentile, zscore, n_testable = score_panel(adata, panel, panel_ens, bin_labels, n_perm, seed)
        out[f"{panel}_percentile"] = percentile
        out[f"{panel}_zscore"] = zscore
        n_testable_report[panel] = n_testable

        if ckpt_path:
            os.makedirs(checkpoint_dir, exist_ok=True)
            df = out[[f"{panel}_percentile", f"{panel}_zscore"]].copy()
            df.to_parquet(ckpt_path)
            with open(n_testable_path, "w") as f:
                f.write(str(n_testable))
            print(f"  {panel}: checkpointed to {ckpt_path}", flush=True)

    return out, n_testable_report
