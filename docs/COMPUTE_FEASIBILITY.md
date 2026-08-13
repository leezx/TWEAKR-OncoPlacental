# Compute Feasibility — Mac mini vs. Argos

Standing reference. **Before starting any analysis step in this project, check this file first** to decide whether it runs locally on the Mac mini or needs to go to Argos. Update this file (new dataset, new analysis step, revised estimate) whenever the answer changes.

## Hardware baseline (measured, not assumed)

### Mac mini (local)

- Nominal 16GB, actual physical RAM: **17.18GB** (`sysctl hw.memsize`)
- **Real-world available headroom is much less than nominal.** Checked live on 2026-08-13: `top -l 1` showed only ~109MB "unused" with 15GB already in use (other apps running, 6.3GB in the memory compressor). This isn't a one-off — assume ordinary daily use (browser, Obsidian, etc.) already consumes several GB before any analysis starts.
- **Working rule: budget ~6-8GB peak for any single local R/Python analysis process**, not the full 16-17GB. This is a safety margin, not a hard OS limit — exceeding it risks heavy swapping (macOS compresses/swaps before killing, so it degrades to unusably slow rather than crashing cleanly, which is often worse because it's not obvious it's failing).
- CPU-bound tasks (not just memory) are also limited — check core count before assuming a "small" job is fast locally.

### Argos (DFCI HPC, remote)

- SSH: `argos.dfci.harvard.edu` (config alias exists; `argos5` needs `argos-login` ProxyJump, not reachable from outside DFCI network/VPN)
- Scheduler: **SGE** (`qsub -pe pvm <slots>`), not SLURM/LSF — existing project tooling already wraps this: `SOFTWARES/bin/argos-submit-cluster`, `SOFTWARES/Argos-Server/env/submit_cluster.sh`, `SOFTWARES/scripts/mount.argos.sh` (sshfs mount for browsing results locally)
- Nodes (from `qhost`, 2026-08-13):

  | Host | CPUs | RAM |
  |---|---|---|
  | argos1-8 | 64 | 376.6GB each |
  | argos9-10 | 160 | 629.3GB each |

- **CPU-only — no GPU** (per earlier project notes; confirm again before any GPU-dependent step, e.g. deep-learning-based methods)
- 20-40x the Mac mini's usable memory per node. Effectively no memory ceiling for anything in this project's current data (largest single file is 20GB on disk).

## Decision rule of thumb

For a dataset being loaded as a **sparse matrix** (scRNA/scATAC count matrix — the common case here), estimate the in-memory floor before attempting a local load:

```
base_load_GB ≈ nnz × 8 bytes / 1e9        (float32 data + int32 index, per matrix)
             × 2   if the file also carries a duplicate "raw" counts layer alongside normalized X
```

Then add margin for:
- **Processing spikes**: normalization, `ScaleData`/`scale()` (can produce a *dense* matrix over selected features × all cells), PCA, clustering, and especially building neighbor graphs typically create additional copies — peak memory during analysis is commonly **2-4x the base load**, not just the base load itself.
- **Multiple objects in one session**: loading dataset A, then also dataset B for a comparison, sums their footprints — don't reason about one file in isolation if the actual analysis needs several loaded at once.

If `base_load_GB × (2 to 4)` comfortably clears the ~6-8GB local budget → Mac mini is fine. If it's close to or exceeds it → Argos, no exceptions (better to submit a 2-minute SGE job than debug a thrashing Mac mini for an hour).

For **RDS files** (R serialized objects — the HDMA Seurat objects), the on-disk size is gzip-compressed and gives no reliable read on in-memory footprint without either loading it (risky for a large one) or inspecting a smaller file from the same generating pipeline as calibration (see HDMA section below).

## Per-dataset verdict

Computed from each file's actual `nnz` (non-zero count) via direct HDF5/MatrixMarket header inspection — not guessed from on-disk size, which is a poor proxy once compression and sparsity ratios differ between files.

### h5ad / mtx datasets (nnz known exactly)

| Dataset | On-disk | nnz | Has duplicate raw layer | Base load | Verdict |
|---|---|---|---|---|---|
| `Arutyunyan2023_MFI` primary_tissue (all donors, all cell states) | 12.38GB | 746,882,673 | ✅ yes (`raw.X`, same nnz as `X`) | **~11.95GB** | **Argos.** Base load alone is already past the local budget; default `sc.read_h5ad()` loads both `X` and `raw.X`. If a local pass is ever needed, load with `sc.read_h5ad(..., backed='r')` and subset before materializing, or strip `raw` on Argos first and ship a lighter file back. |
| `Arutyunyan2023_MFI` organoid PTO | 1.56GB | 92,964,929 | ✅ yes | ~1.49GB | Mac mini OK. |
| `Arutyunyan2023_MFI` organoid TSC | 0.79GB | 45,081,174 | ✅ yes | ~0.72GB | Mac mini OK. |
| `Arutyunyan2023_MFI` organoid Fig3 unstimulated | 2.01GB | 116,051,071 | ✅ yes | ~1.86GB | Mac mini OK. |
| `2026_human_maternal_fetal_Nature` scPlacenta_host | 3.06GB | 380,051,327 | ❌ no | ~3.04GB | Mac mini OK alone; avoid loading alongside `snRNA_raw_counts` in the same session. |
| `2026_human_maternal_fetal_Nature` snRNA_raw_counts | 4.66GB | 579,269,053 | ❌ no | ~4.63GB | Mac mini OK alone, but tight — this is already ~60-75% of the safe local budget before any processing spike. Treat any heavier downstream step (integration, clustering) on this file as Argos-first. |
| `VentoTormo_Nature_2018` decidua-v3 | 0.83GB | 169,413,359 | ❌ no | ~1.36GB | Mac mini OK. |
| `Greenbaum_NatMed_2024` RNA matrix | 2.6GB | 92,502,573 | ❌ no | ~0.74GB | Mac mini OK. |
| `Greenbaum_NatMed_2024` ATAC matrix | 4.9GB | 362,992,363 | ❌ no | ~2.90GB | Mac mini OK for loading the raw matrix. **Signac-style ATAC processing (peak calling, chromVAR, building the 444K-peak × cell graph) is a different, much heavier workload than just loading the matrix — reassess before running that specific step, don't assume the load-time verdict carries over.** |

### HDMA RDS Seurat objects — calibrated empirically 2026-08-13

On-disk sizes are gzip-compressed `.rds` files; unlike h5ad/mtx there's no header to inspect nnz without deserializing the whole object, so this was calibrated for real instead of guessed: loaded `Adrenal` and `Thyroid` locally (the 2 smallest, low-risk) via `readRDS()`, measured R process RSS before/after (`ps -o rss=`) and each object's `n_cells`/`n_genes`.

**Method note**: R's own `object.size()` is unreliable for S4/reference-semantics objects like Seurat's (it recursively over-counts, came out ~2x the actual RSS delta in both tests below) — the before/after process RSS is the trustworthy number for "will this crash my machine," and that's what's used here.

| Organ | On-disk | Cells × Genes | Measured RSS delta | disk/cell | RSS/cell |
|---|---|---|---|---|---|
| Adrenal | 650.5MB | 2,883 × 25,314 | **174MB** | 0.226 MB/cell | 0.060 MB/cell |
| Thyroid | 1,576.8MB | 9,299 × 26,163 | **675MB** | 0.170 MB/cell | 0.073 MB/cell |

Both assays include `RNA`+`SCT`+`decontX` layers plus `pca`/`umap` reductions and 4 neighbor graphs — i.e. these are fully-processed objects, not raw counts only, so the ratios above already include that overhead.

RSS/cell trended up (0.060→0.073) from the smaller to larger file, so extrapolating to the untested organs uses the higher value **×1.5 safety buffer** (0.109 MB/cell) rather than a flat average — deliberately conservative given only 2 calibration points and a machine where getting this wrong is costly (thrashing, not a clean crash). Cell counts for the untested organs aren't known directly, so first estimated from the on-disk/cell ratio (avg 0.198 MB/cell), then converted to an RSS estimate:

| Organ | On-disk | Est. cells | **Est. RSS (conservative)** | Verdict |
|---|---|---|---|---|
| Spleen | 7.41GB | ~37,500 | ~4.1GB | Local budget (~6-8GB) — plausible locally, but this is an *estimate*, not measured. Close other apps first, watch memory live (`top -l 1`), be ready to kill it. Argos if in any doubt or if anything else needs to be loaded at the same time. |
| Thymus | 10.62GB | ~53,700 | ~5.8GB | Same caveat — borderline, estimate only. |
| Skin | 14.36GB | ~72,700 | ~7.9GB | At the edge of the safe budget even on the conservative-but-optimistic estimate — **Argos**, don't risk it locally. |
| Liver | 17.31GB | ~87,600 | ~9.5GB | **Argos.** |
| StomachEsophagus | 19.92GB | ~100,800 | ~11.0GB | **Argos.** |

These 5 estimates are extrapolated, not measured — if any of them becomes analytically important, redo the Adrenal/Thyroid-style calibration properly **on Argos** (where a wrong guess costs nothing) rather than trusting the extrapolation for a load decision on the Mac mini.

### Small/negligible

| Dataset | Size | Verdict |
|---|---|---|
| `HPA_trophoblast_proteome` (`proteinatlas.tsv` + tier files) | 43MB | Mac mini, trivially fine. |
| `Arutyunyan2023_MFI_Visium` (8× Space Ranger outputs) | 2.2GB total, ~300MB/sample | Mac mini fine per-sample; Space Ranger output loads as a much smaller filtered matrix (thousands of spots, not hundreds of thousands of cells) plus an image — nowhere near the scRNA-scale numbers above. |

## How to run something on Argos

Existing tooling in this workspace (not written by this project, reuse it):

```bash
# one-time: configure SOFTWARES/Argos-Server/env/cluster.env (copy from cluster.env.example)
SOFTWARES/bin/argos-submit-cluster bash path/to/script.sh          # submit a job
SOFTWARES/scripts/mount.argos.sh                                     # sshfs-mount Argos locally to browse results
```

SSH is reachable directly for quick interactive checks (`ssh argos.dfci.harvard.edu`), but for actual compute, submit via `qsub`/the wrapper above rather than running heavy jobs on the login node.

## Standing rule for this project

Before starting any new analysis step (not just data loading — clustering, integration, signature scoring, anything that touches a full matrix):

1. Identify which dataset(s) and how many need to be in memory simultaneously.
2. Look them up (or compute fresh, if a new dataset) in the tables above.
3. If any required dataset is flagged Argos, or if the combined base load of everything needed at once pushes past ~6-8GB, run it on Argos — don't attempt it locally "just to see."
4. If a dataset isn't in this file yet, compute its `nnz`-based estimate the same way before deciding (h5ad: `h5py`, peek `X['data']` length; mtx: read the MatrixMarket header's declared nnz; RDS: no shortcut — calibrate with `readRDS()` + `ps -o rss=` on the smallest available file from that source first, the way `Adrenal`/`Thyroid` were done above, before trusting any extrapolation for a larger file in the same series).

## Live check before a local run

Quick pre-flight, not just a one-time read of this file — actual free memory fluctuates with whatever else is running:

```bash
top -l 1 | grep PhysMem     # macOS: how much is actually free right now
```
