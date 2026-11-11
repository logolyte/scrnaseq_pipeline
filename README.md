# scRNA-seq pipeline

Slurm wrapper that takes 10x FASTQs to a quality controlled, ambient-corrected, doublet-free AnnData, optionally with cell types:

1. **STARsolo** — [nf-core/scrnaseq](https://nf-co.re/scrnaseq) 4.2.0 (`starsolo_pipeline/`)
2. **Samplesheet** — `make_scdownstream_samplesheet.py` builds the scdownstream input from step 1
3. **scdownstream** — [nf-core/scdownstream](https://github.com/nf-core/scdownstream), pinned to commit `0048567` on `dev` (no releases yet), with the CellBender and scDblFinder containers replaced (`downstream_processing_pipeline/`)
4. **CellAssign** (optional) — `cellassign.py` assigns cell types from a marker gene list with [CellAssign](https://docs.scvi-tools.org/en/stable/user_guide/models/cellassign.html) (scvi-tools)
5. **Finalize** — `finalize_adata.py` assembles the final AnnData: scdownstream's output with log1p in `X`, simpler key names, doublet scores, the neighbor graph, cell types and sample metadata
6. **Velocyto** (optional) — a second nf-core/scrnaseq run producing spliced/unspliced counts (`velocyto_pipeline/`)

## What each step does

| Step | Settings |
|---|---|
| STARsolo | `GeneFull_Ex50pAS` features (exons + introns, for nuclei), `EmptyDrops_CR` cell calling, Cell Ranger UMI dedup (`1MM_CR`) and multi-gene UMI filtering, CellBender and RDS conversion skipped |
| Samplesheet | Passes only the **raw** matrix per sample, so scdownstream does the cell calling itself; chooses CellBender's expected cells and total droplets per sample (`cellbender_priors.py`, see Notes) |
| scdownstream | CellBender empty-droplet removal (GPU) → decontX ambient correction (empty droplets as background) → QC (only `min_genes` 200 and `min_counts_cell` 200, see Notes) → scDblFinder doublet removal (on the decontX-corrected counts) → log1p, PCA on the 4000 genes with the highest binomial deviance, kNN graph on 30 PCs, UMAP and Leiden; downstream analysis runs on the ambient-corrected counts |
| CellAssign | Trained on the marker genes from `MARKER_GENES` only, using the counts scdownstream analysed (`ambient_corrected_decontx`), rounded stochastically to integers, with size factors from all genes; CPU. Fitted 5 times from different seeds, keeping the fit with the best likelihood. Each cell gets its most probable type (see Notes) |
| Finalize | `X` = log1p (layers `raw_counts` and `ambient_corrected_decontx` kept). Renames `<key>-<res>_leiden` → `leiden_<res>`, `<key>-<res>:entropy` → `leiden_<res>_sample_entropy`, `X_<key>_umap` → `X_umap`, `X_pca_pca` → `X_pca`, `deviance:*` → `*`. Adds `doublet_score` (scDblFinder), `ambient_fraction_removed`, `total_counts_raw`, `n_genes_by_counts_raw`, the CellAssign `cell_type`, `cell_type_probability` and `obsm['cellassign_probabilities']`, the scdownstream kNN graph (`obsp`, `uns['neighbors']`), sample metadata from the samplesheet and provenance (`uns['pipeline']`). Drops the duplicates `n_counts`, `n_genes` and the placeholders `donor`, `batch`, `condition`, `label` and var `intersection` (see Notes) |
| Velocyto | `Gene Velocyto` features, same STAR settings; the h5ads carry `spliced`/`unspliced`/`ambiguous` layers |

## Setup

1. Edit `pipeline.conf`: `PROJECT_ROOT`, `SBATCH_ACCOUNT`, `CONTAINER_DIR`, `SAMPLESHEET`, `PROTOCOL`, `FASTA_PATH`, `GTF_PATH`, `STAR_INDEX_PATH`, `RUN_VELOCYTO`, `MARKER_GENES`. The STAR index must have been built from that FASTA/GTF. Leave `MARKER_GENES` empty to skip cell type assignment.
2. Run `bash setup.sh` once (in an interactive job if your cluster doesn't allow builds on the login node). It creates the conda environment and builds the pipeline's own containers. It skips anything that already exists.
3. Fill in the samplesheet (`SAMPLESHEET`, default `samplesheet.csv` next to `pipeline.conf`; a relative path is relative to the directory of `pipeline.conf`): `sample,fastq_1,fastq_2`, one row per lane, with lanes of one sample sharing the `sample` name. Leave out `expected_cells`: if it's set, nf-core/scrnaseq adds its own `--soloCellFilter CellRanger2.2`.
   - **Sample metadata:** add any other columns (e.g. `condition,age_weeks,time_days`). nf-core/scrnaseq ignores them, and `finalize_adata.py` copies them into `obs` for every cell of the sample. Lanes of one sample must have the same values. Columns where every value is a number become numeric, the others categorical, with the categories in samplesheet order. The columns nf-core/scrnaseq reads (`fastq_barcode`, `expected_cells`, `seq_center`, `sample_type`, `feature_type`) are not copied, and a name that already exists in `obs` (e.g. `total_counts`) stops the job.
4. Check the cluster settings in the `slurm.config` files: partition `smp`, 64 CPUs / 72 GB for `process_high`. The CellBender GPU job uses the `gpu` queue and an H100 slice (`--gpus=h100_1g.12gb:1`), hardcoded in `downstream_processing_pipeline/slurm.config`.

### Environments and containers

| Path | Contents |
|---|---|
| `envs/scrnaseq.yml` | Conda environment `scrnaseq`, top-level packages: Nextflow, nf-core, STAR and the Python packages the pipeline's own scripts use, including scvi-tools with CPU-only PyTorch for CellAssign. The scripts activate it by name. |
| `envs/scrnaseq.lock.yml` | The same environment with exact versions of every package. `setup.sh` uses it, and falls back to `scrnaseq.yml` if it no longer resolves. |
| `containers/cellbender/` | CellBender 0.4.0 on PyTorch 2.4.1/CUDA 12.1. Replaces the scdownstream module's GPU image (CellBender 0.3.2), which doesn't run on H100 GPUs. |
| `containers/scdblfinder/` | scDblFinder 1.24.10. Replaces the module's 1.24.0, which calls 2–3× too few doublets with xgboost ≥ 3.1 ([nf-core/scdownstream#382](https://github.com/nf-core/scdownstream/issues/382)). |

`setup.sh` builds each `containers/<name>/<image>.def` into `$CONTAINER_DIR/<image>.sif`, and `downstream_processing_pipeline/slurm.config` points the two processes at those images.
- **Lock files:** each recipe installs exact versions from its `lock.txt`. If a recipe has no `lock.txt` yet, the first build writes one from the image; commit it so later builds are identical.
- **Rebuilding:** delete the `.sif` and rerun `setup.sh`.
- **Updating the conda environment:** `setup.sh` skips an existing `scrnaseq` environment. After `envs/` changes (e.g. the CellAssign packages were added), run `conda env update -n scrnaseq -f envs/scrnaseq.lock.yml --prune`.
- **No `apptainer build --fakeroot`:** if your cluster doesn't allow it, build the `.def` elsewhere and copy the `.sif` into `$CONTAINER_DIR`.

## Running

Run these on the login node, from the repo root:

```bash
bash run_pipeline.sh [path/to/pipeline.conf]   # all steps, chained with Slurm dependencies
bash run_downstream.sh [path/to/pipeline.conf] # scdownstream only (samplesheet must exist)
```

CellAssign (if `MARKER_GENES` is set) starts after scdownstream, and `finalize_adata.py` after CellAssign (or after scdownstream without `MARKER_GENES`). To rerun CellAssign on its own, e.g. with another marker list or a probability cut-off, and then rebuild the final h5ad:

```bash
(
  set -a; source pipeline.conf; set +a
  export PIPELINE_CONF="$PWD/pipeline.conf"
  conda activate scrnaseq
  sbatch --chdir=. --output=slurm_logs/cellassign_%j.out --error=slurm_logs/cellassign_%j.err \
         cellassign.py [--markers marker_gene_sets/other.yaml] [--min-probability 0.9] [--outdir ...]
  # once it has finished (or chain it with --dependency=afterok:<job id>):
  sbatch --chdir=. --output=slurm_logs/finalize_%j.out --error=slurm_logs/finalize_%j.err \
         finalize_adata.py [--cellassign-dir ...] [--output ...]
)
```

`python cellassign.py --help` and `python finalize_adata.py --help` list all options. `finalize_adata.py` reads only the CellAssign outputs (`cell_types.tsv`, `run_info.json`), so rebuilding the final h5ad, e.g. after adding samplesheet columns, does not refit CellAssign. It stops if the CellAssign cells differ from the scdownstream cells, i.e. CellAssign is from an older scdownstream run.

The velocyto run (if `RUN_VELOCYTO="true"`) is submitted independently of the other jobs. To run it on its own:

```bash
(
  set -a; source pipeline.conf; set +a
  export PIPELINE_CONF="$PWD/pipeline.conf"
  sbatch --chdir=velocyto_pipeline \
         --output="$PWD/slurm_logs/velocyto_%j.out" --error="$PWD/slurm_logs/velocyto_%j.err" \
         velocyto_pipeline/velocyto_nfcore.sh
)
```

All runs use `-resume`, so resubmitting only redoes changed or failed tasks. Logs go to `slurm_logs/`.

## Outputs

All under `PROJECT_ROOT`:

| Path | Contents |
|---|---|
| `results/starsolo/star/<sample>/` | BAM, `Solo.out` matrices |
| `results/starsolo/star/mtx_conversions/` | `<sample>/<sample>_{raw,filtered}_matrix.h5ad`, `combined_*.h5ad` |
| `results/final/<name>.h5ad` | **The final AnnData** (see Finalize above) |
| `results/scdownstream/09_finalized/` | scdownstream's AnnData and provenance; `finalize_adata.py` reads it and the intermediates in `02_quality_control/` and `05_cluster_dimred/` (`save_intermediates` must stay on) |
| `results/cellassign/` | `cell_types.tsv` (per cell: sample, cell type, its probability, one probability column per type), `cell_type_counts.tsv`, `run_info.json` (settings and cell types, becomes `uns['cellassign']`), `restarts.tsv` (one row per fit), `training_history.png`, and the marker list used (`markers.yaml`, `marker_matrix.tsv`) |
| `results/velocyto/star/mtx_conversions/` | Same layout as starsolo, with velocity layers |
| `work/{starsolo,scdownstream,velocyto}/` | Nextflow work dirs (delete when done) |

## Marker gene sets

`marker_gene_sets/` holds the marker lists for CellAssign; `mouse_glomerulus.yaml` has the mouse glomerular and kidney cell types from the Lmx1b project. Each list is a YAML file with one cell type per line and its marker genes, named as in the scdownstream h5ad (gene symbols, e.g. `Nphs1` for mouse, `NPHS1` for human):

```yaml
Podocyte: [Nphs1, Nphs2, Podxl, Synpo]
Endothelial: [Pecam1, Emcn, Kdr, Flt1]
Immune: [Ptprc, Cd74, C1qa]
```

Name new lists by species and tissue (e.g. `mouse_glomerulus.yaml`) and point `MARKER_GENES` at one. A gene can be a marker for more than one type. The job logs markers that are missing from the data. It stops if a cell type has no markers left, or if more than half are missing (usually the wrong species or gene IDs).

## Notes

- **Final h5ad, counts and QC columns:** `raw_counts` holds the STARsolo counts of the cells CellBender called (CellBender only chooses the barcodes; its denoised counts are not used). scdownstream's QC columns (`total_counts`, `n_genes_by_counts`, `pct_counts_*`) are computed on the decontX counts; `total_counts_raw` and `n_genes_by_counts_raw` are on `raw_counts`.
  - `ambient_fraction_removed` = 1 − decontX counts / raw counts per cell. scdownstream drops decontX's own `contamination` estimate, so this is the share decontX actually removed, which is close to that estimate but not identical.
  - `doublet_score` is scDblFinder's score from the per-sample files in `02_quality_control/doublet_detection/scdblfinder/`. The cells called doublets are already removed, so the score shows how close the remaining cells came to the cut-off.
  - `obsp`/`uns['neighbors']` is the graph scdownstream computed `X_umap` and the Leiden clusters from (15 neighbors on 30 PCs of `X_pca`), copied from `05_cluster_dimred/pca/neighbors/`. Recomputing it with other genes or PCs gives a different graph (see `X_pca`: PCA of the 4000 genes with the highest binomial deviance, `var['highly_deviant']`).
  - `donor`, `batch`, `condition` and `label` are dropped because scdownstream fills them with the sample name or `Unknown`; give them as samplesheet columns to have them.

- **CellAssign fits can merge cell types:** depending on its random start, a fit can end in a local optimum where two cell types are merged into one and the other stays (nearly) empty. `cellassign.py` therefore fits `--n-restarts` times (default 5) and keeps the fit with the lowest negative log-likelihood. `restarts.tsv` and the job log show how far each fit agrees with the best one. If no other fit agrees on at least 95% of cells, the log warns that a better optimum may have been missed. Then check `cell_type_counts.tsv` for an expected type that is empty, and rerun with more restarts.
- **CellAssign has no "other" class:** every cell is assigned to one of the listed types, so cells of types missing from the marker list are forced into the closest one. List every type you expect in the tissue, and check the probabilities: `--min-probability` (e.g. 0.9) labels cells below it `Unassigned`. All probabilities are saved, so you can also set a cut-off later in the analysis.
- **Velocyto counts:** the velocyto run's `X` counts exons only (STAR requires `Gene` for `Velocyto`, and `Gene` can't be combined with `GeneFull_Ex50pAS` here). Use only its `spliced`/`unspliced` layers. Take them from the per-sample **raw** h5ads and match them to the scdownstream cells by sample and barcode. Correct the layers for ambient RNA before running scVelo. In snRNA-seq, ambient RNA is mostly spliced.
- **No `--twopassMode Basic`:** a STAR bug breaks it with `GeneFull_Ex50pAS` ([STAR#1533](https://github.com/alexdobin/STAR/issues/1533)). `validation.ignoreParams` lets the `star_feature` value through nf-core's parameter check.
- **CellBender settings per sample:** `make_scdownstream_samplesheet.py` calls `cellbender_priors.py` for each sample and writes `downstream_processing_pipeline/cellbender_priors.tsv`. `slurm.config` reads `--expected-cells` and `--total-droplets-included` from that file by sample name, because scdownstream only passes `expected_cells` from the samplesheet to CellBender. The run stops if a sample is missing from the file.
  - `expected_cells` = 0.75 × the STARsolo EmptyDrops cells. This keeps borderline EmptyDrops calls out of CellBender's cell prior; CellBender still calls more cells when the data support them.
  - `total_droplets` = the number of barcodes with at least as many UMI as the empty-droplet peak. The peak is the tallest peak of the log-UMI histogram of non-EmptyDrops barcodes with ≥ 10 UMI, which avoids the low plateau of barcode-swapped barcodes. This puts CellBender's empty prior on the empty droplets.
  - Without a usable peak (none, or below 20 UMI), `total_droplets` falls back to EmptyDrops cells + 20,000.
  - **Check before running scdownstream:** the samplesheet job's log and the file's `warnings` column. The codes are explained at the top of `cellbender_priors.py`. `POOR_SEPARATION` and `EMPTYDROPS_CALLS_EMPTIES` mean no prior settings will give clean cell calls, so inspect the sample (damage, ambient RNA, depth). `LARGE_TOTAL` means CellBender will be slow and may not fit on the GPU.
  - The thresholds were tuned on high-ambient snRNA-seq (mouse kidney). In low-ambient data, e.g. single-cell, expect more fallbacks; check the first run.
- **Per-sample `use_ambient_corrected_counts`:** the samplesheet builder writes this column for every sample. If it's empty, scdownstream silently turns ambient-corrected counts off.
- **QC is deliberately lenient:** the pipeline only removes empty droplets and debris, and stricter filtering is left to each analysis. Both filters are computed on the decontX-corrected counts.
  - `min_counts_cell` 200 removes droplets with fewer than 200 UMI left after ambient correction. These are mostly ambient RNA, such as debris.
  - `min_genes` 200 effectively counts raw genes, because decontX's fractional counts keep every detected gene above zero. It has no effect when the empty droplets already have more than 200 genes (high-ambient nuclei), but it removes near-empty droplets in single-cell data.
- **No mitochondrial filter (single-cell data):** the pipeline sets no mito cut-off, which is fine for nuclei, where mito reads are close to zero after ambient correction. Single-cell output is not QC-complete until a mito filter has been applied in the analysis.
  - Dying or broken cells keep their mitochondrial transcripts but lose cytoplasmic mRNA, so they have a high mito fraction and pass both gene and count floors.
  - In kidney, healthy proximal tubule and other metabolically active cells can have 20–50% mito, so choose the cut-off per tissue, ideally per cell type, or use a per-sample adaptive cut-off (e.g. `pct_counts_mt_nmads`), and check which cells it removes.
