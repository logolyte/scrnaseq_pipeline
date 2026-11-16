# Preprocessing

Takes one samplesheet of FASTQs to a preprocessed AnnData with cell types, as a chain of Slurm jobs.
Works on in-house data and on the samplesheets of `public_data/`.

```
samplesheet.csv                      one row per sequencing run: sample, fastq_1, fastq_2, protocol, metadata
   │
   │  1_curate: you choose the samples to preprocess and whether they are sc or sn
   ▼
<outdir>/samplesheet.csv             the chosen runs and the metadata columns of config.yaml
<outdir>/config.yaml                 a copy of config.yaml with your answers; every later step reads it
   │
   │  2_scrnaseq: STARsolo
   ▼
results/starsolo/                    raw and EmptyDrops-filtered matrix per sample
   │
   │  3_scdownstream: CellBender priors, then CellBender → decontX → QC → scDblFinder → PCA, UMAP, Leiden
   ▼
results/scdownstream/                one h5ad of all samples
   │
   │  4_cellassign: cell types from marker genes (and 4_velocity_decontx with run_velocyto)
   ▼
results/cellassign/                  cell type of every cell
   │
   │  5_finalize
   ▼
results/final/<name>.h5ad            log1p in X, raw and decontX counts, cell types, sample metadata
```

Steps 2–5 are submitted as Slurm jobs, each starting when the ones it needs have finished, and `run_pipeline.py` returns.
Start the next dataset right away with another outdir.

## 1. Set up

Once, from the repository root:

1. Fill in the TODOs of `preprocessing/config.yaml`.
2. Create the conda environment `scrnaseq` and build our container images into `<cache_dir>/containers/`:
   ```bash
   bash setup.sh <cache_dir of config.yaml>
   ```

| Image | Why our own |
|---|---|
| `containers/cellbender/` | CellBender 0.4.0: the module's image (0.3.2) does not run on H100 GPUs |
| `containers/scdblfinder/` | scDblFinder 1.24.10: the module's 1.24.0 calls 2–3× too few doublets with xgboost ≥ 3.1 ([nf-core/scdownstream#382](https://github.com/nf-core/scdownstream/issues/382)) |
| `containers/decontx/` | the image of scdownstream's decontX module, for `4_velocity_decontx/` |

## 2. Write the samplesheet

One row per sequencing run; runs of one sample share `sample` and must have the same protocol and metadata. `protocol` is the scrnaseq `--protocol`, e.g. `10XV3`. FASTQ paths are absolute or relative to the samplesheet folder:

```
sample,fastq_1,fastq_2,protocol,sample_genotype
WT_1,/data/WT_1_L001_R1.fastq.gz,/data/WT_1_L001_R2.fastq.gz,10XV3,WT
WT_1,/data/WT_1_L002_R1.fastq.gz,/data/WT_1_L002_R2.fastq.gz,10XV3,WT
KO_1,/data/KO_1_L001_R1.fastq.gz,/data/KO_1_L001_R2.fastq.gz,10XV3,KO
```

`public_data/organize_fetchngs.py` writes this format for GEO studies. The metadata columns listed in `metadata_columns` of `config.yaml` go to obs of the final h5ad.

## 3. Start a run

On the login node, with `scrnaseq` active:

```bash
python preprocessing/run_pipeline.py <samplesheet.csv> <outdir>
```

It shows the samples with their metadata and asks:

| question | answer | meaning |
|---|---|---|
| Samples | `1,3,5`, `1-5` | these samples |
| | Enter | all samples |
| Take GSM1, GSM2? | Enter | yes |
| | anything else | choose again |
| sc or sn | `sc` / `sn` | cells / nuclei: the settings of `modes` in `config.yaml` (STAR features, QC) |

Then it submits the jobs and prints their IDs. Logs are in `<outdir>/slurm_logs/`.

After a failure, rerun the same command: it does not ask again (the copy of the config exists), and the Nextflow steps continue with `-resume`.
Editing `preprocessing/config.yaml` does not change a run already started, but the code of the repository is read by each job when it starts: do not change it while runs are queued.

## 4. Results

```
<outdir>/
├── config.yaml, samplesheet.csv     from 1_curate
├── scrnaseq/, velocyto/             launch folders of Nextflow: .nextflow, params.yaml, work/ (delete when done)
├── scdownstream/                    launch folder: params.yaml, samplesheet.csv, cellbender_priors.tsv, work/
├── slurm_logs/
└── results/
    ├── starsolo/                    BAM, Solo.out, mtx_conversions/<sample>/<sample>_{raw,filtered}_matrix.h5ad
    ├── scdownstream/                09_finalized/ and the intermediates 5_finalize reads
    ├── cellassign/                  cell_types.tsv, markers.yaml
    ├── velocyto/, velocyto_decontx/ with run_velocyto
    └── final/<name>.h5ad            the final AnnData
```

## What each step does

| Step | Settings |
|---|---|
| 2_scrnaseq | STARsolo with `star_feature` of the mode: `Gene` (sc) or `GeneFull_Ex50pAS` (sn, exons + introns); `EmptyDrops_CR` cell calling, Cell Ranger UMI dedup (`1MM_CR`) and multi-gene UMI filtering |
| 3_scdownstream | `prepare_scdownstream.py` gives scdownstream only the **raw** matrix per sample, so CellBender calls the cells, with expected cells and total droplets per sample from `cellbender_priors.py`; QC thresholds = lenient defaults + the `qc` of the mode. scdownstream: CellBender (GPU) → decontX (empty droplets as background) → QC → scDblFinder → log1p, PCA on the 4000 genes with the highest binomial deviance, kNN graph on 30 PCs, UMAP, Leiden; all on the decontX counts |
| 4_cellassign | CellAssign trained on the marker genes only, on the decontX counts rounded stochastically to integers, size factors from all genes; 5 fits from seeds 0–4, the one with the lowest negative log-likelihood is kept |
| 4_velocity_decontx | with `run_velocyto`: a second scrnaseq run with `Gene Velocyto`; its spliced/unspliced/ambiguous counts are corrected with one decontX run per sample on the three layers stacked |
| 5_finalize | `X` = log1p; layers `raw_counts` and `ambient_corrected_decontx`; adds doublet score, ambient fraction removed, raw count metrics, cell types, kNN graph, sample metadata and the run settings (`uns['pipeline']`); simpler key names (see the docstring of `finalize_adata.py`) |

## Files

```
preprocessing/
├── config.yaml               reference, Slurm, sc/sn settings, metadata columns, markers
├── run_pipeline.py           starts one run
├── 1_curate/                 curate_samples.py
├── 2_scrnaseq/               scrnaseq.config, params.yaml
├── 3_scdownstream/           prepare_scdownstream.py, cellbender_priors.py, scdownstream.config, params.yaml
├── 4_cellassign/             cellassign.py, marker_gene_sets/
├── 4_velocity_decontx/       velocity_decontx.py, velocity_decontx.R
└── 5_finalize/               finalize_adata.py
```

The number of a step folder is its order; the two steps 4 run in parallel.
`params.yaml` of a step holds the nf-core pipeline parameters of every run; `run_pipeline.py` adds the values of the run.

## Notes

- **Counts:** `raw_counts` holds the STARsolo counts of the cells CellBender called; CellBender only chooses the barcodes. scdownstream's QC columns (`total_counts`, `pct_counts_*`) are on the decontX counts, `total_counts_raw` and `n_genes_by_counts_raw` on `raw_counts`.
- **CellBender priors:** `expected_cells` = 0.75 × the EmptyDrops cells, which keeps borderline EmptyDrops calls out of the cell prior. `total_droplets` = the barcodes with at least as many UMI as the empty-droplet peak (the tallest peak of the log-UMI histogram of non-EmptyDrops barcodes with ≥ 10 UMI); without a usable peak, EmptyDrops cells + 20,000.
  Check `slurm_logs/prepare_scdownstream_*.log` and the `warnings` column of `scdownstream/cellbender_priors.tsv` before trusting the cell calls; the codes are explained at the top of `cellbender_priors.py`. The thresholds were tuned on high-ambient mouse kidney nuclei.
- **QC of sn is lenient:** `min_counts_cell` and `min_genes` 200 remove empty droplets and debris; stricter filtering is left to each analysis. `min_genes` effectively counts raw genes, because decontX's fractional counts keep every detected gene above zero.
- **Mitochondrial filter:** not needed for nuclei, where mito reads are close to zero after ambient correction. Cells need one (`max_mito_percentage` or `pct_counts_mt_nmads` in `modes.sc.qc`); in kidney, healthy proximal tubule cells can have 20–50% mito.
- **CellAssign has no "other" class:** cells of types missing from the marker list are forced into the closest listed type. List every type you expect and check `cell_type_probability`. A fit can merge two cell types; the job log warns when the 5 fits disagree.
- **Velocyto counts:** the velocyto run's `X` counts exons only (`Velocyto` needs `Gene`); only its spliced/unspliced/ambiguous layers are used, matched to the scdownstream cells by sample and barcode and to the genes by symbol.
- **No `--twopassMode Basic`:** a STAR bug breaks it with `GeneFull_Ex50pAS` ([STAR#1533](https://github.com/alexdobin/STAR/issues/1533)).
