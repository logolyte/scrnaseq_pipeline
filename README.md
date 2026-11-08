# scRNA-seq pipeline

Slurm wrapper that takes 10x FASTQs to a QC'd, ambient-corrected, doublet-free AnnData:

1. **STARsolo** — [nf-core/scrnaseq](https://nf-co.re/scrnaseq) 4.2.0 (`starsolo_pipeline/`)
2. **Samplesheet** — `make_scdownstream_samplesheet.py` builds the scdownstream input from step 1
3. **scdownstream** — fork [logolyte/scdownstream](https://github.com/logolyte/scdownstream), branch `scdblfinder-first` (`downstream_processing_pipeline/`)
4. **Velocyto** (optional) — a second nf-core/scrnaseq run producing spliced/unspliced counts (`velocyto_pipeline/`)

## What each step does

| Step | Settings |
|---|---|
| STARsolo | `GeneFull_Ex50pAS` features (exons + introns, for nuclei), `EmptyDrops_CR` cell calling, Cell Ranger UMI dedup (`1MM_CR`) and multi-gene UMI filtering, CellBender and RDS conversion skipped |
| Samplesheet | Passes only the **raw** matrix per sample, so scdownstream does the cell calling itself; `expected_cells` = number of STARsolo EmptyDrops cells |
| scdownstream | CellBender empty-droplet removal (GPU) → scDblFinder doublet removal → decontX ambient correction (empty droplets as background) → QC (only `min_genes` 200 and `min_counts_cell` 200, see Notes) → log1p, PCA (30 PCs); downstream analysis runs on the ambient-corrected counts |
| Velocyto | `Gene Velocyto` features, same STAR settings; the h5ads carry `spliced`/`unspliced`/`ambiguous` layers |

## Setup

1. Edit `pipeline.conf`: `PROJECT_ROOT`, `SBATCH_ACCOUNT`, `CONTAINER_DIR`, `PROTOCOL`, `FASTA_PATH`, `GTF_PATH`, `STAR_INDEX_PATH`, `RUN_VELOCYTO`. The STAR index must have been built from that FASTA/GTF.
2. Run `bash setup.sh` once (in an interactive job if your cluster doesn't allow builds on the login node). It creates the conda environment and builds the pipeline's own containers. It skips anything that already exists.
3. Fill in `starsolo_pipeline/samplesheet.csv`: `sample,fastq_1,fastq_2`, one row per lane, with lanes of one sample sharing the `sample` name. Leave out `expected_cells`: if it's set, nf-core/scrnaseq adds its own `--soloCellFilter CellRanger2.2`.
4. Check the cluster settings in the `slurm.config` files: partition `smp`, 64 CPUs / 72 GB for `process_high`. The CellBender GPU job uses the `gpu` queue and an H100 slice (`--gpus=h100_1g.12gb:1`), hardcoded in `downstream_processing_pipeline/slurm.config`.

### Environments and containers

| Path | Contents |
|---|---|
| `envs/scrnaseq.yml` | Conda environment `scrnaseq`, top-level packages: Nextflow, nf-core, STAR and the Python packages the pipeline's own scripts use. The scripts activate it by name. |
| `envs/scrnaseq.lock.yml` | The same environment with exact versions of every package. `setup.sh` uses it, and falls back to `scrnaseq.yml` if it no longer resolves. |
| `containers/cellbender/` | CellBender 0.4.0 on PyTorch 2.4.1/CUDA 12.1. Replaces the scdownstream module's GPU image (CellBender 0.3.2), which doesn't run on H100 GPUs. |
| `containers/scdblfinder/` | scDblFinder 1.24.10. Replaces the module's 1.24.0, which calls 2–3× too few doublets with xgboost ≥ 3.1 ([nf-core/scdownstream#382](https://github.com/nf-core/scdownstream/issues/382)). |

`setup.sh` builds each `containers/<name>/<image>.def` into `$CONTAINER_DIR/<image>.sif`, and `downstream_processing_pipeline/slurm.config` points the two processes at those images.
- **Lock files:** each recipe installs exact versions from its `lock.txt`. If a recipe has no `lock.txt` yet, the first build writes one from the image; commit it so later builds are identical.
- **Rebuilding:** delete the `.sif` and rerun `setup.sh`.
- **No `apptainer build --fakeroot`:** if your cluster doesn't allow it, build the `.def` elsewhere and copy the `.sif` into `$CONTAINER_DIR`.

## Running

Run these on the login node, from the repo root:

```bash
bash run_pipeline.sh [path/to/pipeline.conf]   # all steps, chained with Slurm dependencies
bash run_downstream.sh [path/to/pipeline.conf] # scdownstream only (samplesheet must exist)
```

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
| `results/scdownstream/09_finalized/` | Final AnnData and provenance |
| `results/velocyto/star/mtx_conversions/` | Same layout as starsolo, with velocity layers |
| `work/{starsolo,scdownstream,velocyto}/` | Nextflow work dirs (delete when done) |

## Notes

- **Velocyto counts:** the velocyto run's `X` counts exons only (STAR requires `Gene` for `Velocyto`, and `Gene` can't be combined with `GeneFull_Ex50pAS` here). Use only its `spliced`/`unspliced` layers. Take them from the per-sample **raw** h5ads and match them to the scdownstream cells by sample and barcode. Correct the layers for ambient RNA before running scVelo. In snRNA-seq, ambient RNA is mostly spliced.
- **No `--twopassMode Basic`:** a STAR bug breaks it with `GeneFull_Ex50pAS` ([STAR#1533](https://github.com/alexdobin/STAR/issues/1533)). `validation.ignoreParams` lets the `star_feature` value through nf-core's parameter check.
- **Per-sample `use_ambient_corrected_counts`:** the samplesheet builder writes this column for every sample. If it's empty, scdownstream silently turns ambient-corrected counts off.
- **QC is deliberately lenient:** the pipeline only removes empty droplets and debris, and stricter filtering is left to each analysis. Both filters are computed on the decontX-corrected counts.
  - `min_counts_cell` 200 removes droplets with fewer than 200 UMI left after ambient correction. These are mostly ambient RNA, such as debris.
  - `min_genes` 200 effectively counts raw genes, because decontX's fractional counts keep every detected gene above zero. It has no effect when the empty droplets already have more than 200 genes (high-ambient nuclei), but it removes near-empty droplets in single-cell data.
- **No mitochondrial filter (single-cell data):** the pipeline sets no mito cut-off, which is fine for nuclei, where mito reads are close to zero after ambient correction. Single-cell output is not QC-complete until a mito filter has been applied in the analysis.
  - Dying or broken cells keep their mitochondrial transcripts but lose cytoplasmic mRNA, so they have a high mito fraction and pass both gene and count floors.
  - In kidney, healthy proximal tubule and other metabolically active cells can have 20–50% mito, so choose the cut-off per tissue, ideally per cell type, or use a per-sample adaptive cut-off (e.g. `pct_counts_mt_nmads`), and check which cells it removes.
