# Preprocessing

Turns one samplesheet of FASTQs into a preprocessed AnnData with cell types.
The steps run as a chain of Slurm jobs.
It works on in-house data and on the samplesheets written by `public_data/`.

```
samplesheet.csv                      one row per sequencing run: sample, fastq_1, fastq_2, protocol, metadata
   │
   │  1_curate: you choose the samples and the mode (sc or sn)
   ▼
<outdir>/samplesheet.csv             the runs of the chosen samples
<outdir>/config.yaml                 a copy of config.yaml with your choices; every later step reads it
   │
   │  2_scrnaseq: STARsolo counts the reads per gene and cell barcode
   ▼
results/starsolo/                    one count matrix per sample
   │
   │  3_scdownstream: finds the cells, removes ambient RNA and doublets, then PCA, UMAP and clusters
   ▼
results/scdownstream/                one h5ad with all samples
   │
   │  4_cellassign: cell types from marker genes
   ▼
results/cellassign/                  the cell type of every cell
   │
   │  5_finalize: puts everything into one file
   ▼
results/final/merged.h5ad            the result
```

`run_pipeline.py` submits steps 2–5 as Slurm jobs and returns.
Each job starts when the jobs it needs have finished, so you can start the next dataset right away with another outdir.

## 1. Set up

Do this once, from the repository root:

1. Fill in the TODOs of `preprocessing/config.yaml`.
2. Create the conda environment `scrnaseq` and build our container image into `<cache_dir>/containers/`:
   ```bash
   bash setup.sh <cache_dir of config.yaml>
   ```

The only image we build ourselves is CellBender 0.4.0 (`containers/cellbender/`).
The image that comes with nf-core/scdownstream has CellBender 0.3.2, which does not run on H100 GPUs.

## 2. Write the samplesheet

Write one row per sequencing run.
Runs of the same sample have the same `sample` name and must have the same protocol and metadata.
`protocol` is the `--protocol` of nf-core/scrnaseq, e.g. `10XV3`.
FASTQ paths can be absolute or relative to the folder of the samplesheet.

```
sample,fastq_1,fastq_2,protocol,sample_genotype
WT_1,/data/WT_1_L001_R1.fastq.gz,/data/WT_1_L001_R2.fastq.gz,10XV3,WT
WT_1,/data/WT_1_L002_R1.fastq.gz,/data/WT_1_L002_R2.fastq.gz,10XV3,WT
KO_1,/data/KO_1_L001_R1.fastq.gz,/data/KO_1_L001_R2.fastq.gz,10XV3,KO
```

For GEO studies, `public_data/organize_fetchngs.py` writes this file for you.

Any other column is sample metadata, like `sample_genotype` in this example.
The columns you list in `metadata_columns` of `config.yaml` are copied to obs of the final h5ad.

## 3. Start a run

On the login node, with the `scrnaseq` environment active:

```bash
python preprocessing/run_pipeline.py <samplesheet.csv> <outdir>
```

The script shows the samples with their metadata and asks:

| question | answer | meaning |
|---|---|---|
| Samples | `1,3,5`, `1-5` | take these samples |
| | Enter | take all samples |
| Take GSM1, GSM2? | Enter | yes |
| | anything else | choose again |
| Mode (sc, sn) | `sc` / `sn` | cells or nuclei; each mode has its own settings under `modes` in `config.yaml` |

To skip the questions, add `--mode sc` or `--mode sn`. The run then takes all samples of the samplesheet.

The script then submits the jobs and prints their IDs. The log of each step is in `<outdir>/slurm_logs/`.

If a run fails, fix the problem and run the same command again.
It does not ask the questions again, and the Nextflow steps continue where they stopped.

Each run has its own copy of `config.yaml`, so later edits of `preprocessing/config.yaml` do not change a run that has already started.
The code is different: every job reads it when it starts, so do not change the repository while jobs are waiting in the queue.

## 4. Results

```
<outdir>/
├── README.txt                       what was run, by whom and on which samples
├── config.yaml, samplesheet.csv     the settings and samples of this run
├── slurm_logs/                      one log per step
├── scrnaseq/, scdownstream/         Nextflow launch folders with work/ (delete them when the run is done)
└── results/
    ├── starsolo/                    STARsolo output and the count matrices of each sample
    ├── scdownstream/                QC report, plots and tables of nf-core/scdownstream
    ├── cellassign/                  cell_types.tsv: the cell type and its probability for every cell
    └── final/merged.h5ad            the final AnnData
```

The final h5ad contains:

| Where | What |
|---|---|
| `X` | normalised counts (log1p) |
| `layers['raw_counts']` | the counts from STARsolo |
| `layers['ambient_corrected_decontx']` | the counts after decontX removed ambient RNA; all analysis steps use these |
| `obs` | QC metrics, `doublet_score`, `ambient_fraction_removed`, the clusters `leiden_0.25`, `leiden_0.5` and `leiden_1.0`, `cell_type`, `cell_type_probability`, `protocol` (the 10x chemistry) and your metadata columns |
| `obsm` | `X_pca`, `X_umap` and `cellassign_probabilities` (the probability of every cell type) |
| `obsp` | the neighbour graph behind the UMAP and the clusters |
| `uns['pipeline']` | the settings and software versions of the run |

## 5. Copy to storage

When a run has finished, copy its logs, reports, tables and the final h5ad from the server to your storage.
The large intermediate files (BAM, count matrices, h5/h5ad/pkl) stay in the outdir on the server until you delete it.
Of `work/`, only the logs of each task (`.command.*`, `.exitcode`) are copied.

```bash
rsync -avPh --prune-empty-dirs \
    --exclude='/*/.nextflow/' \
    --include='/*/work/??/*/.command.*' \
    --include='/*/work/??/*/.exitcode' \
    --include='/results/final/*.h5ad' \
    --include='*/' \
    --exclude='/*/work/**' \
    --exclude='*.bam' --exclude='*.tab' --exclude='*.gz' \
    --exclude='*.h5' --exclude='*.h5ad' --exclude='*.pkl' \
    username@server:/path/to/outdir/ /path/to/new/folder 
```

Add `-n --stats` to only list the files and their total size.

## What each step does

**2_scrnaseq.** STARsolo aligns the reads and counts them per gene.
Reads in introns are counted too (`GeneFull_Ex50pAS`, the `star_feature` of the mode in `config.yaml`).
Cell calling and UMI counting follow Cell Ranger.

**3_scdownstream.** First, `prepare_scdownstream.py` estimates for each sample how many cells and droplets CellBender should expect.
Then nf-core/scdownstream runs on each sample:

1. CellBender decides which barcodes are cells (on a GPU).
2. decontX removes ambient RNA.
3. Cells that fail the QC thresholds of `config.yaml` are removed.
4. scDblFinder finds doublets, and they are removed.

After that the samples are merged and normalised (log1p).
PCA uses the 4000 most informative genes (highest binomial deviance).
The neighbour graph, the UMAP and the Leiden clusters use the first 30 PCs.

**4_cellassign.** CellAssign gives every cell one of the cell types of the marker file (`marker_genes` in `config.yaml`).
It is fitted five times and the best fit is kept.

**4_velocity_decontx.** Runs only with `run_velocyto: true`.
A second STARsolo run counts spliced and unspliced reads, and decontX corrects them.

**5_finalize.** Collects everything in the final h5ad and gives scdownstream's keys shorter names.

## Files

```
preprocessing/
├── config.yaml               reference, Slurm, QC and mode settings, metadata columns, markers
├── run_pipeline.py           starts one run
├── 1_curate/                 curate_samples.py
├── 2_scrnaseq/               scrnaseq.config, params.yaml
├── 3_scdownstream/           prepare_scdownstream.py, cellbender_priors.py, scdownstream.config, params.yaml
├── 4_cellassign/             cellassign.py, marker_gene_sets/
├── 4_velocity_decontx/       velocity_decontx.py, velocity_decontx.R
└── 5_finalize/               finalize_adata.py
```

The number of a step folder is its order. The two steps 4 run at the same time.
The `params.yaml` of a step holds the nf-core parameters that are the same for every run.

## Good to know

**There are two kinds of counts.**
`raw_counts` are the counts from STARsolo. `ambient_corrected_decontx` are the same counts after decontX.
CellBender is only used to decide which barcodes are cells; its corrected counts are not kept.
The QC columns from scdownstream (`total_counts`, `pct_counts_mt`, …) are calculated on the decontX counts.
`total_counts_raw` and `n_genes_by_counts_raw` are calculated on the raw counts.

**Check the CellBender settings before you trust the number of cells.**
For every sample the pipeline chooses how many cells and droplets CellBender should expect:
75% of the cells STARsolo called, and all barcodes down to the peak of empty droplets.
Look at the `warnings` column of `<outdir>/scdownstream/cellbender_priors.tsv`.
The codes are explained at the top of `cellbender_priors.py`.
These rules were tuned on mouse kidney nuclei with a lot of ambient RNA, so look more carefully at other tissues.

**QC is mild on purpose.**
It removes cells with fewer than 200 genes or 200 counts, and cells with more than 5% (nuclei) or 15% (cells) mitochondrial counts.
Stricter filtering is left to your analysis.
You can change the thresholds under `qc` and `modes` in `config.yaml`.

**CellAssign has no "unknown" type.**
Every cell gets one of the cell types of the marker file, even if its real type is not listed.
So list every cell type you expect, and look at `cell_type_probability` before you trust a label.
If none of the other four fits agrees with the best one, `slurm_logs/cellassign_*.log` has a warning.

**Velocyto is off by default.**
It needs a second alignment of all samples.
With `run_velocyto: true` the final h5ad also has the layers `spliced`, `unspliced` and `ambiguous` (after decontX) and the same three with `_raw` (before decontX).
