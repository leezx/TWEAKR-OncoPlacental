#!/usr/bin/env python3
"""
Step 6 dataset extension: shared loader + coverage-check + provenance-audit
machinery for HTAN_CRC_progressive_plasticity + CRLM_NMP_ATLAS (per
docs/STEP6_DATASET_EXTENSION_DESIGN.md, PR #33, APPROVE after 3 review
rounds). Imported by dataset_extension_scoring.py. Reuses
crc_gut_scoring_core.py's scoring machinery unchanged -- this module only
builds a correctly-axis-aligned `adata` object matching the primary
atlas's internal contract, so score_all_panels() etc. work unmodified.

Pre-compute inspection (done interactively on Argos before writing this
script, not assumed): both extension datasets' `adata.raw.var_names` and
`adata.var_names` are byte-identical bare-Ensembl-ID arrays (confirmed
`(a.raw.var_names == a.var_names).all() == True` for both
HTAN_CRC_progressive_plasticity and CRLM_NMP_ATLAS), so the
version-suffix/canonicalization branch below is a real assertion, not a
guess that happens to be skipped -- it still runs and would raise if that
assumption were ever violated on a re-run against updated source data.
`adata.raw.X` is confirmed 100% integer-valued (raw counts); `adata.X` is
already normalized (non-integer, max ~8-9 consistent with log1p output)
-- confirmed directly, not inferred from the inventory doc's shorthand.
"""
import os
import re
import time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

REPO = "/home/zz950/TWEAKR-OncoPlacental"
DATA_ROOT = "/home/zz950/DATA/CRC-Atlas/phase2/03_data/raw"

EXTENSION_DATASETS = {
    "HTAN_CRC_progressive_plasticity": f"{DATA_ROOT}/HTAN_CRC_progressive_plasticity/epithelial.h5ad",
    "CRLM_NMP_ATLAS": f"{DATA_ROOT}/CRLM_NMP_ATLAS/crlm_nmp_atlas.h5ad",
}

# Primary meta-atlas path, reused only for the provenance-overlap audit
# (crc_gut_scoring_core.ATLAS_H5AD is the same path -- not re-imported to
# avoid pulling in that module's full 665,473-cell load path here).
META_ATLAS_H5AD = "/home/zz950/DATA/scRNAseq/meta_study/CRC_single_cell_atlas_2025/adata_nmf.h5ad"


def strip_ensembl_version(gene_ids):
    """'ENSG00000141510.4' -> 'ENSG00000141510'; no-op if already bare."""
    return pd.Index([g.split(".")[0] for g in gene_ids])


def load_extension_dataset(path, dataset_name, to_csc=True):
    """Loads an extension dataset and standardizes it into the SAME
    internal contract crc_gut_scoring_core.load_atlas() produces for the
    primary atlas: layers['counts'] holds raw integer counts, X holds
    normalize_total(1e4)+log1p of those same counts, and var_names is the
    canonical bare-Ensembl-ID gene axis every downstream function
    (compute_detectability, testable_genes, score_genes_fast) indexes
    against.

    Round-2-review-locked contract (docs/STEP6_DATASET_EXTENSION_DESIGN.md
    section 1/2): raw.X is NOT assumed to align positionally with
    adata.var_names just because the dimensions match. Explicit steps:
      (a) assert a 1:1 correspondence between adata.raw.var_names and
          adata.var_names (equal length + set equality, not just count);
      (b) if either axis carries version suffixes, canonicalize BOTH to
          stripped bare-Ensembl IDs and assert no duplicate collisions
          on either axis after stripping;
      (c) only then assign the (axis-verified) adata.raw.X into
          adata.layers['counts'];
      (d) rebuild X via normalize_total(target_sum=1e4) + log1p of that
          same layers['counts'] -- identical recipe to the primary
          atlas's loader, so the null-calibration machinery treats this
          object exactly like the primary atlas's.
    """
    print(f"Loading {dataset_name} from {path} ...", flush=True)
    adata = ad.read_h5ad(path)
    print(f"  Loaded {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

    if adata.raw is None:
        raise ValueError(
            f"{dataset_name}: adata.raw is None -- the confirmed-raw-counts "
            f"layer (.raw.X, verified at inventory time) is missing. "
            f"Refusing to guess a substitute count source."
        )

    raw_var_names = pd.Index(adata.raw.var_names)
    work_var_names = pd.Index(adata.var_names)

    # (a) 1:1 correspondence check -- not assumed from equal length alone.
    if len(raw_var_names) != len(work_var_names):
        raise ValueError(
            f"{dataset_name}: raw.var_names (n={len(raw_var_names)}) and "
            f"var_names (n={len(work_var_names)}) have different lengths -- "
            f"positional alignment cannot be assumed. Refusing to proceed."
        )

    byte_identical = raw_var_names.equals(work_var_names)
    print(f"  raw.var_names byte-identical to var_names: {byte_identical}", flush=True)

    if byte_identical:
        # Round-1 review correction: the byte-identical branch previously
        # accepted the axes as "no suffix stripping needed" without ever
        # actually checking they are bare Ensembl IDs -- byte-identical
        # is not the same claim as bare-ID, and the code silently assumed
        # the latter. Fixed to assert it directly (not inferred from
        # coverage numbers looking reasonable): every ID must match the
        # bare `ENSG` + 11-digit pattern, no version suffix, no other
        # format. Confirmed to hold for both extension datasets before
        # this fix was written (checked directly, not assumed), but the
        # assertion itself is what makes that a proven fact of this run,
        # not an inference from indirect evidence.
        bare_ensembl_pattern = re.compile(r"^ENSG[0-9]{11}$")
        non_bare = [v for v in work_var_names if not bare_ensembl_pattern.match(v)]
        if non_bare:
            raise ValueError(
                f"{dataset_name}: var_names is byte-identical to raw.var_names, "
                f"but {len(non_bare)} entries are not bare Ensembl IDs "
                f"(e.g. {non_bare[:5]}) -- refusing to assume this is a "
                f"simple no-suffix case without checking the ID format itself."
            )
        canonical_var_names = work_var_names
    else:
        # (b) versioned-ID branch: canonicalize BOTH axes, assert no
        # duplicate collisions on either.
        raw_stripped = strip_ensembl_version(raw_var_names)
        work_stripped = strip_ensembl_version(work_var_names)
        if raw_stripped.has_duplicates:
            dupes = raw_stripped[raw_stripped.duplicated()].unique().tolist()
            raise ValueError(
                f"{dataset_name}: stripping version suffixes from raw.var_names "
                f"produces {len(dupes)} duplicate bare-Ensembl-ID collisions "
                f"(e.g. {dupes[:5]}) -- refusing to silently corrupt the "
                f"gene-to-column mapping."
            )
        if work_stripped.has_duplicates:
            dupes = work_stripped[work_stripped.duplicated()].unique().tolist()
            raise ValueError(
                f"{dataset_name}: stripping version suffixes from var_names "
                f"produces {len(dupes)} duplicate bare-Ensembl-ID collisions "
                f"(e.g. {dupes[:5]}) -- refusing to silently corrupt the "
                f"gene-to-column mapping."
            )
        if not raw_stripped.equals(work_stripped):
            raise ValueError(
                f"{dataset_name}: stripped raw.var_names and stripped var_names "
                f"are still not identical after version-suffix removal -- the "
                f"two axes are not a simple versioned/unversioned pair of the "
                f"same gene order. Refusing to assume positional alignment."
            )
        canonical_var_names = raw_stripped
        adata.var_names = raw_stripped  # canonicalize the working axis too

    # (c) only now assign the axis-verified raw counts.
    adata.layers["counts"] = adata.raw.X.copy()
    adata.var_names = canonical_var_names

    counts_dense_sample = adata.layers["counts"][:1000]
    counts_dense_sample = (counts_dense_sample.toarray()
                            if hasattr(counts_dense_sample, "toarray")
                            else np.asarray(counts_dense_sample))
    frac_integer = np.mean(np.mod(counts_dense_sample, 1) == 0)
    print(f"  layers['counts'] integer fraction (1000-cell sample): {frac_integer:.4f}", flush=True)
    if frac_integer < 0.999:
        raise ValueError(
            f"{dataset_name}: layers['counts'] (from raw.X) is not confirmed "
            f"integer-valued (fraction={frac_integer:.4f}) -- refusing to "
            f"treat this as raw counts."
        )

    # (d) rebuild X -- identical recipe to crc_gut_scoring_core.load_atlas().
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    if to_csc:
        # Same fix already established for the primary atlas
        # (crc_gut_scoring_core.load_atlas): both extension datasets'
        # raw.X/X are CSR by default (confirmed, not assumed), and
        # repeated column-slicing (score_genes_fast's access pattern) on
        # CSR is a known severe scipy inefficiency. One-time conversion
        # cost is negligible against a 13-panel x 500-perm run even at
        # these datasets' smaller scale.
        t0 = time.time()
        adata.X = adata.X.tocsc()
        print(f"  Converted X to CSC in {time.time()-t0:.1f}s", flush=True)

    return adata


def canonical_marker_sentinel_check(adata, dataset_name, group_col, markers=None):
    """Biological sentinel check (structural assertions in
    load_extension_dataset() are the real safeguard; this is a
    supplementary check, not a substitute -- per round-1 review). Reports
    mean log-normalized expression of a few canonical epithelial
    (EPCAM/KRT8/KRT18/KRT19, expected high in malignant + normal
    epithelial cells) and immune (PTPRC, expected high in TME cells, low
    in epithelial cells) markers, per group in `group_col`. Markers are
    keyed directly by their stable, well-known GRCh38 Ensembl gene ID
    (hardcoded here, not resolved via any symbol-mapping file this
    project's other Ensembl-ID-primary panels already avoid depending
    on) -- checked for presence in `adata.var_names` directly, and
    reported as unresolvable (not silently skipped) if absent."""
    if markers is None:
        markers = {
            "EPCAM": "ENSG00000119888",
            "KRT8": "ENSG00000170421",
            "KRT18": "ENSG00000111057",
            "KRT19": "ENSG00000171345",
            "PTPRC": "ENSG00000081237",
        }

    rows = []
    for m, ens in markers.items():
        if ens not in adata.var_names:
            rows.append({"dataset": dataset_name, "marker": m, "ensembl_id": ens,
                         "status": "not_in_dataset"})
            continue
        gi = adata.var_names.get_loc(ens)
        expr = adata.X[:, gi]
        expr = expr.toarray().ravel() if hasattr(expr, "toarray") else np.asarray(expr).ravel()
        for grp, mask in adata.obs.groupby(group_col).groups.items():
            idx = adata.obs.index.get_indexer(mask)
            rows.append({
                "dataset": dataset_name, "marker": m, "ensembl_id": ens,
                "group": grp, "n_cells": len(idx),
                "mean_log1p_expr": float(np.mean(expr[idx])),
                "status": "ok",
            })
    return pd.DataFrame(rows)


def coverage_check(adata, dataset_name, panels, panel_ensembl_ids_fn):
    """Per-panel n_testable = |panel ∩ adata.var_names| against the
    canonicalized working axis (round-2 review correction -- not
    raw.var_names in isolation)."""
    var_set = set(adata.var_names)
    rows = []
    for panel in panels:
        ens_ids = set(panel_ensembl_ids_fn(panel))
        n_testable = len(ens_ids & var_set)
        rows.append({"dataset": dataset_name, "panel": panel,
                     "n_panel_genes": len(ens_ids), "n_testable": n_testable})
    return pd.DataFrame(rows)


def htan_provenance_overlap_audit():
    """Study-provenance overlap audit (docs/STEP6_DATASET_EXTENSION_DESIGN.md
    section 4, required before any 'replication'/'external validation'
    language). Cross-references HTAN_CRC_progressive_plasticity's donor/
    patient identifiers against every identifier column in the
    CRC_single_cell_atlas_2025 meta-atlas's obs (dataset/sample_id/
    patient_id/donor_id/*_accession), across all 36 constituent studies
    -- not just the one ('HTAPP_HTAN') whose name suggests HTAN
    provenance, since a real overlap could in principle appear under any
    study's identifiers."""
    htan_path = EXTENSION_DATASETS["HTAN_CRC_progressive_plasticity"]
    htan = ad.read_h5ad(htan_path, backed="r")
    htan_donor_ids = set(htan.obs["donor_id"].astype(str).unique())
    htan_patient_ids = set(htan.obs["Patient"].astype(str).unique())
    htan_sample_ids = set(htan.obs["Sample ID"].astype(str).unique())
    # donor_id carries the real cross-cohort-comparable identifier prefix
    # (e.g. "HTA8_6001"); Patient/Sample ID are dataset-internal short
    # codes ("103", "KG103M") unlikely to appear verbatim elsewhere, but
    # checked anyway for completeness.
    htan_id_prefixes = {i.split("_")[0] for i in htan_donor_ids if "_" in i}

    meta = ad.read_h5ad(META_ATLAS_H5AD, backed="r")
    id_cols = ["dataset", "sample_id", "patient_id", "donor_id",
               "NCBI_BioProject_accession", "SRA_sample_accession",
               "GEO_sample_accession", "ENA_sample_accession",
               "synapse_sample_accession"]
    id_cols = [c for c in id_cols if c in meta.obs.columns]

    rows = []
    direct_hits = set()
    for col in id_cols:
        vals = meta.obs[col].astype(str)
        exact_overlap = set(vals.unique()) & (htan_donor_ids | htan_patient_ids | htan_sample_ids)
        prefix_hits = 0
        for pfx in htan_id_prefixes:
            prefix_hits += int(vals.str.contains(pfx, na=False, regex=False).sum())
        rows.append({
            "meta_atlas_column": col,
            "n_exact_id_overlap": len(exact_overlap),
            "exact_overlap_examples": list(exact_overlap)[:5],
            "n_htan_id_prefix_hits": prefix_hits,
        })
        direct_hits |= exact_overlap

    studies_with_htan_prefix_hits = []
    for pfx in htan_id_prefixes:
        for col in ["donor_id", "patient_id", "sample_id"]:
            if col not in meta.obs.columns:
                continue
            hit_mask = meta.obs[col].astype(str).str.contains(pfx, na=False, regex=False)
            if hit_mask.any():
                studies_with_htan_prefix_hits.extend(
                    meta.obs.loc[hit_mask, "study_id"].astype(str).unique().tolist()
                )
    studies_with_htan_prefix_hits = sorted(set(studies_with_htan_prefix_hits))

    summary = pd.DataFrame(rows)
    n_total_exact_overlap = len(direct_hits)
    audit_conclusion = (
        "NO_OVERLAP_FOUND" if n_total_exact_overlap == 0 and not studies_with_htan_prefix_hits
        else "OVERLAP_FOUND_REVIEW_REQUIRED"
    )
    return summary, audit_conclusion, studies_with_htan_prefix_hits


def htan_name_similarity_supplementary_check():
    """Round-1 review correction: the results write-up asserted a
    specific narrative (the meta-atlas's name-similar `HTAPP_HTAN` study
    is actually a different HTAN sub-cohort, `HTA1_`-prefixed via Pelka
    et al. 2021's Synapse project, vs. our dataset's `HTA8_`-prefixed
    identifiers) that was checked interactively during development but
    never committed as auditable script output -- an unsupported claim
    by this project's own standing discipline. This function makes that
    specific check a real, reproducible, committed artifact: for every
    meta-atlas study whose name contains 'HTA' (a superset net of just
    'HTAPP_HTAN', so it doesn't silently miss a differently-named
    HTAN-affiliated study), reports its cell count and a sample of its
    donor/patient/sample identifiers, so the specific narrative claim is
    directly checkable against committed output rather than asserted."""
    meta = ad.read_h5ad(META_ATLAS_H5AD, backed="r")
    study_ids = meta.obs["study_id"].astype(str)
    candidate_studies = sorted(study_ids[study_ids.str.contains("HTA", na=False, regex=False)].unique())
    rows = []
    for study in candidate_studies:
        sub = meta.obs[study_ids == study]
        for col in ["donor_id", "patient_id", "sample_id"]:
            if col not in sub.columns:
                continue
            examples = sorted(sub[col].astype(str).unique().tolist())[:5]
            rows.append({
                "meta_atlas_study_id": study, "n_cells": len(sub),
                "identifier_column": col, "example_identifiers": examples,
            })
    return pd.DataFrame(rows)
