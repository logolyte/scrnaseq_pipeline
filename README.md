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
| scdownstream | CellBender empty-droplet removal (GPU) → scDblFinder doublet removal → decontX ambient correction (empty droplets as background) → QC (only `min_genes` 200) → log1p, PCA (30 PCs); downstream analysis runs on the ambient-corrected counts |
| Velocyto | `Gene Velocyto` features, same STAR settings; the h5ads carry `spliced`/`unspliced`/`ambiguous` layers |

## Setup

1. Create the conda environment: `conda env create -f environment.yml` (exact versions) or `environment_without_dependencies.yml` (top-level packages only). The scripts activate an environment named `scrnaseq`.
2. Edit `pipeline.conf`: `PROJECT_ROOT`, `SBATCH_ACCOUNT`, `PROTOCOL`, `FASTA_PATH`, `GTF_PATH`, `STAR_INDEX_PATH`, `RUN_VELOCYTO`. The STAR index must have been built from that FASTA/GTF.
3. Fill in `starsolo_pipeline/samplesheet.csv`: `sample,fastq_1,fastq_2`, one row per lane, with lanes of one sample sharing the `sample` name. Leave out `expected_cells`: if it's set, nf-core/scrnaseq adds its own `--soloCellFilter CellRanger2.2`.
4. Check the cluster settings in the `slurm.config` files: partition `smp`, 64 CPUs / 72 GB for `process_high`. The CellBender GPU job uses the `gpu` queue, the account `ag-beyer` and a local container (`cellbender-0.4.0-gpu.sif`), all hardcoded in `downstream_processing_pipeline/slurm.config`.

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
