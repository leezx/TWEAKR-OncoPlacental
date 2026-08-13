# Compute Feasibility — Mac mini vs. Argos

Standing reference. **Before starting any analysis step in this project, check this file first** to decide whether it runs locally on the Mac mini or needs to go to Argos. Update this file (new dataset, new analysis step, revised estimate) whenever the answer changes.

## Hardware baseline (measured, not assumed)

### Mac mini (local)

- Nominal 16GB, actual physical RAM: **17.18GB** (`sysctl hw.memsize`)
- **Real-world available headroom is much less than nominal.** Checked live on 2026-08-13: `top -l 1` showed only ~88-116MB "unused" across several checks, with 15GB already in use and 4.7-8.2GB cycling through the memory compressor. This isn't a one-off — assume ordinary daily use (browser, Obsidian, etc.) already consumes several GB before any analysis starts.
- **Working rule: budget ~6-8GB peak for any single local R/Python analysis process**, not the full 16-17GB. This is a safety margin, not a hard OS limit — exceeding it risks heavy swapping (macOS compresses/swaps before killing, so it degrades to unusably slow rather than crashing cleanly, which is often worse because it's not obvious it's failing).
- CPU-bound tasks (not just memory) are also limited — check core count before assuming a "small" job is fast locally.
- **⚠️ Live process RSS is not reliable evidence of memory safety on this machine** — see the measurement-methodology section below. Use the deterministic dtype/nnz-based accounting instead.

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

## Measurement methodology — revised 2026-08-13 after review

The first version of this doc used a single flat formula (`nnz × 8 bytes`, assuming float32 data + int32 indices for every file) and live process-RSS readings as the safety check. Both turned out to have real gaps, caught by review and then confirmed empirically rather than just patched by argument:

**1. The flat 8-bytes/nnz assumption doesn't hold across file formats.** Verified by actually loading files and inspecting the resulting array dtypes/nbytes directly (not estimating):

| File | Loaded via | Assumed dtype | **Actual dtype** | Formula estimate | **Actual measured** |
|---|---|---|---|---|---|
| `2026_human_maternal_fetal_Nature` snRNA_raw_counts (h5ad) | `anndata.read_h5ad` | float32+int32 | float32+int32 (confirmed) | 4.63GB | **4.67GB** (sparse payload + obs/var, matches formula within 1%) |
| `Greenbaum_NatMed_2024` RNA matrix (mtx) | `scipy.io.mmread` | float32+int32 | **float64** data | ~0.74GB | **1.11GB** (CSR, after `.tocsr()`) |
| `Greenbaum_NatMed_2024` ATAC matrix (mtx) | `scipy.io.mmread` | float32+int32 | **int64** data+indices | ~2.90GB | **4.36GB** (CSR) / **5.81GB** (COO, as initially loaded — `mmread` returns COO, which needs a full row *and* col index array, not one index + a compact indptr) |

h5ad files store their own dtype in the file and it happened to be float32 here, so the original formula was accurate *for those specific files* — but that was luck, not something the formula checked for. Generic MatrixMarket (`.mtx`) readers make no such promise: `scipy.io.mmread` upcasts integer-valued entries to `int64` and real-valued entries to `float64` by default, and returns COO (not CSR/CSC) until you explicitly convert it. Both of those independently increase memory ~1.5-2x over the naive float32/int32/CSR assumption. **This is exactly the "dtype upcasting" and "format" gap flagged in review** — now measured, not just acknowledged.

**Corrected rule: don't assume a dtype — check what the actual loading path produces**, either by inspecting the file's own declared dtype (h5py: `f['X']['data'].dtype`) or, for formats that don't declare one (mtx), by doing a quick empirical load-and-inspect (`m.dtype`, `m.data.nbytes` etc.) before trusting an estimate, the same way this revision did.

**2. Live process RSS underreports on this specific machine.** Because the Mac mini is chronically under memory pressure (see hardware section), macOS aggressively compresses resident pages in the background. Measured directly: loading the 4.67GB (true, nbytes-verified) `snRNA_raw_counts` object showed only a **994MB-1.8GB** `ps`-reported RSS delta across repeated runs — a large, inconsistent undercount, not measurement noise at the margins. The **deterministic nbytes accounting is now the trusted number**, not live RSS. This matters practically: a "looks fine, RSS is low" reading during an idle/background-compressed moment is not evidence of safety once real computation (e.g. normalization touching every element) forces the OS to decompress everything at once.

Net effect: base-load estimates for h5ad files with a known/verified dtype are validated as accurate. The two Greenbaum mtx-derived estimates and the HDMA RDS estimates (see below — calibrated via the now-suspect RSS method) are revised upward / flagged with wider uncertainty accordingly.

## Decision rule of thumb

For a dataset being loaded as a **sparse matrix**, estimate the in-memory floor before attempting a local load:

```
base_load_GB ≈ nnz × bytes_per_nonzero / 1e9
  where bytes_per_nonzero = (data_dtype_size + index_dtype_size), typically:
    8  bytes  if data is float32/int32 and matrix is CSR/CSC (verify the actual dtype — don't assume)
    12 bytes  if data is float64 (8B) + int32 index (4B), CSR/CSC
    16 bytes  if data is float64/int64 (8B) + int64 index (8B), OR the matrix is still in COO form (needs both row+col index arrays)
             × 2   if the file also carries a duplicate "raw" counts layer alongside normalized X
```

Then add margin for:
- **Processing spikes**: normalization, `ScaleData`/`scale()` (can produce a *dense* matrix over selected features × all cells), PCA, clustering, and especially building neighbor graphs typically create additional copies — peak memory during analysis is commonly **2-4x the base load**, not just the base load itself.
- **Multiple objects in one session**: loading dataset A, then also dataset B for a comparison, sums their footprints — don't reason about one file in isolation if the actual analysis needs several loaded at once.

If `base_load_GB × (2 to 4)` comfortably clears the ~6-8GB local budget → Mac mini is fine. If it's close to or exceeds it → Argos, no exceptions (better to submit a 2-minute SGE job than debug a thrashing Mac mini for an hour). **Verify the estimate came from a checked dtype, not an assumed one — see the methodology section above.**

## Per-dataset verdict

### h5ad datasets — dtype-verified via `h5py`/`anndata`, formula confirmed accurate

| Dataset | On-disk | nnz | Verified dtype | Has duplicate raw layer | Base load | Verdict |
|---|---|---|---|---|---|---|
| `Arutyunyan2023_MFI` primary_tissue (all donors, all cell states) | 12.38GB | 746,882,673 | float32+int32 | ✅ yes (`raw.X`, same nnz as `X`) | **~11.95GB** | **Argos.** Base load alone is already past the local budget; default `sc.read_h5ad()` loads both `X` and `raw.X`. If a local pass is ever needed, load with `sc.read_h5ad(..., backed='r')` and subset before materializing, or strip `raw` on Argos first and ship a lighter file back. |
| `Arutyunyan2023_MFI` organoid PTO | 1.56GB | 92,964,929 | float32+int32 | ✅ yes | ~1.49GB | Mac mini OK. |
| `Arutyunyan2023_MFI` organoid TSC | 0.79GB | 45,081,174 | float32+int32 | ✅ yes | ~0.72GB (measured: 754MB incl. obs/var/obsm/obsp) | Mac mini OK — empirically verified, see calibration note below. |
| `Arutyunyan2023_MFI` organoid Fig3 unstimulated | 2.01GB | 116,051,071 | float32+int32 | ✅ yes | ~1.86GB | Mac mini OK. |
| `2026_human_maternal_fetal_Nature` scPlacenta_host | 3.06GB | 380,051,327 | float32+int32 | ❌ no | ~3.04GB | Mac mini OK alone; avoid loading alongside `snRNA_raw_counts` in the same session. |
| `2026_human_maternal_fetal_Nature` snRNA_raw_counts | 4.66GB | 579,269,053 | float32+int32 | ❌ no | ~4.63GB (**measured: 4.67GB** incl. obs/var — formula confirmed accurate) | **Borderline, treat as Argos-first for anything beyond a quick load-and-inspect.** ~4.7GB already eats most of the 6-8GB safety budget before any processing spike (2-4x multiplier would push it to 9-19GB). |
| `VentoTormo_Nature_2018` decidua-v3 | 0.83GB | 169,413,359 | float32+int32 | ❌ no | ~1.36GB (measured: 645MB RSS delta — likely an undercount per the methodology note, but even the formula ceiling is small) | Mac mini OK. |

### mtx datasets — dtype-verified via actual load, revised upward from the original estimate

| Dataset | On-disk | nnz | Actual dtype | **Measured base load (CSR)** | Verdict |
|---|---|---|---|---|---|
| `Greenbaum_NatMed_2024` RNA matrix | 2.6GB | 92,502,573 | float64 (not float32) | **1.11GB** | Mac mini OK — revised up from the original 0.74GB estimate but still comfortably within budget. |
| `Greenbaum_NatMed_2024` ATAC matrix | 4.9GB | 362,992,363 | int64 (not float32/int32) | **4.36GB** (CSR) / up to **5.81GB** if a downstream tool keeps it in COO form | **Revised from "Mac mini OK" to borderline/Argos-first.** At 4.36GB this alone is already most of the safe local budget; loaded together with the RNA matrix (1.11GB) for any joint RNA+ATAC step, or run through Signac-style processing (peak calling, chromVAR, building the 444K-peak × cell graph — a categorically heavier workload than just loading the matrix), **default to Argos.** |

### HDMA RDS Seurat objects — calibrated empirically, but via the now-suspect RSS method; treat with wider uncertainty than the h5ad numbers above

Loaded `Adrenal` and `Thyroid` locally via `readRDS()`, measured R process RSS before/after (`ps -o rss=`). **Caveat added in this revision**: the same RSS-underreporting-under-memory-pressure effect confirmed for Python (see methodology section) plausibly affects these R measurements too — there's no RDS-native equivalent to `nbytes`-based verification the way sparse-matrix arrays have, so unlike the h5ad numbers above, **these have not been cross-checked against a dtype-verified ground truth and should be read as a lower bound, not a confident estimate.**

| Organ | On-disk | Cells × Genes | Measured RSS delta (lower-bound only) | disk/cell | RSS/cell |
|---|---|---|---|---|---|
| Adrenal | 650.5MB | 2,883 × 25,314 | **174MB** | 0.226 MB/cell | 0.060 MB/cell |
| Thyroid | 1,576.8MB | 9,299 × 26,163 | **675MB** | 0.170 MB/cell | 0.073 MB/cell |

Both objects include `RNA`+`SCT`+`decontX` assays plus `pca`/`umap` reductions and 4 neighbor graphs. Given the RSS-reliability caveat above, both are still small enough in absolute terms (sub-1GB) to trust as "Mac mini OK" even if the true figure is 2-3x higher than measured — but that safety margin runs out fast for the larger, untested organs:

| Organ | On-disk | Est. cells | Est. RSS (extrapolated, ×1.5 buffer on trend) | **Revised verdict** |
|---|---|---|---|---|
| Spleen | 7.41GB | ~37,500 | ~4.1GB | **Argos-first**, not "plausible locally" as the previous revision said — given RSS is now known to be an unreliable local safety signal on this machine, an estimate that's already in the 4-6GB range (before considering the possible 2-3x RSS undercount seen elsewhere) is too close to the budget to trust with a live-memory-only check. |
| Thymus | 10.62GB | ~53,700 | ~5.8GB | **Argos-first**, same reasoning. |
| Skin | 14.36GB | ~72,700 | ~7.9GB | **Argos.** |
| Liver | 17.31GB | ~87,600 | ~9.5GB | **Argos.** |
| StomachEsophagus | 19.92GB | ~100,800 | ~11.0GB | **Argos.** |

If any of Spleen/Thymus becomes analytically important, redo the calibration **on Argos** with a dtype/nbytes-verified method (e.g. `lobstr::obj_size()` cross-checked against `Matrix::nnzero()` on the sparse slots, not just process RSS) rather than trusting either the Mac mini extrapolation above or a live-RSS-only check locally.

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
4. If a dataset isn't in this file yet, **verify the actual dtype the loading path produces before estimating** — h5ad/mtx headers or a quick empirical load-and-inspect (`X.dtype`, `X.data.nbytes`) — don't assume float32/int32. For RDS: no shortcut, calibrate with `readRDS()` + `ps -o rss=` on the smallest available file first, and treat the result as a lower bound given the RSS-reliability caveat above.
5. **Don't trust a live "looks fine" RSS reading as sufficient evidence on this machine** — it's been shown to underreport under memory pressure. Prefer the deterministic nnz/dtype-based estimate; use live RSS only as an early-warning trip-wire during a run (if it's already high, stop), not as proof something fits.

## Live check before a local run

Quick pre-flight, not just a one-time read of this file — actual free memory fluctuates with whatever else is running:

```bash
top -l 1 | grep PhysMem     # macOS: how much is actually free right now
```
