# Public data from GEO

Downloads GEO studies with nf-core/fetchngs and writes one samplesheet per GSE for `preprocessing/run_pipeline.py`.
Only 10x Chromium 3' v2, v3 and v4 are supported now (`chemistries/tenx.py`).

```
ids.csv                              the GSEs to download
   │
   │  nf-core/fetchngs (step 2)
   ▼
FETCHNGS_OUT/
├── fastq/                           every FASTQ of every GSE together: <SRX>_<SRR>_<read>.fastq.gz
└── samplesheet/samplesheet.csv      one row per sequencing run, with its GSE and GSM
   │
   │  organize_fetchngs.py (step 3)
   ▼
DATASET/<GSE>/                       one folder per GSE
├── *.fastq.gz                       R1 and R2 of the chosen GSMs (index reads stay in FETCHNGS_OUT)
├── samplesheet.csv                  input of preprocessing/run_pipeline.py, with the sample metadata from GEO/SRA
└── README.txt                       the study, who sorted the files and when, the skipped GSMs
```

`organize_fetchngs.py` moves the files the preprocessing needs, with the sample metadata, into one folder per study,
following the data structure rules of AG Beyer.

## 1. List the studies

`ids.csv`: one GSE per line, no header, e.g.

```
GSE146912
GSE171722
```

## 2. Download with fetchngs

```bash
nextflow run nf-core/fetchngs -r 1.13.0 -profile apptainer -c fetchngs.config \
    --input ids.csv --outdir /path/to/FETCHNGS_OUT --download_method sratools
```

*Why sra-tools:* ENA (the default) gives only the reads SRA marks as biological; for 10x these are often just R2 (cDNA).
R1 (cell barcode + UMI) and the index reads are technical and come only from `fasterq-dump --include-technical`.

To add studies, append GSEs to the bottom of `ids.csv` and rerun the same command with `-resume` from the same folder.

## 3. Sort the FASTQs by GSE

`config.yaml` sets how many reads of each FASTQ are measured and the number of FASTQs read at the same time. Your email there is optional: it is sent to NCBI with each request, without it NCBI prints a warning.

In the `scrnaseq` environment (it needs biopython):

```bash
python organize_fetchngs.py --fetchngs /path/to/FETCHNGS_OUT --destination /path/to/DATASET
```

For each GSM it shows two tables and asks which files to take:

| answer | meaning |
|---|---|
| Enter | take the proposal |
| `tenx 1 2` | one library: file 1 = R1, file 2 = R2 |
| `tenx 1 2 3 4` | one library sequenced in two runs |
| `tenx 1 2, 3 4` | two libraries |
| `-` | skip this GSM |

Files table:
- `length`: median read length;
- `phred30_head` / `phred30_tail`: share of bases with Phred ≥ 30 in the first / last quarter of the read;
- `phred30_bp`: number of read positions where ≥ 70% of bases have Phred ≥ 30;
- `original_names`: file names as the authors uploaded them to SRA.

tenx table:
- `role`: R1, R2 or I (index, ≤ 10 bp);
- `protocol`: the 10x version whose whitelist matches best;
- `wl`: share of reads whose first 16 nt are in that whitelist; R1 is the read with the highest `wl`;
- `tso`: share of reads with the TSO right after the UMI, a sign of 10x 5';
- `library`: runs with the same cells, i.e. sharing at least half of their 1000 most frequent cell barcodes,
  are one library.

No proposal if `wl` < 0.75, if R1 is shorter than cell barcode + UMI, or if `tso` ≥ 0.1.

When a GSE is done, its R1/R2 move to `DATASET/<GSE>/` with `samplesheet.csv`: one row per run, `sample` = GSM
(GSM_LIB_1, GSM_LIB_2, … if a GSM has several libraries), `protocol`, the fetchngs columns and the BioSample attributes from SRA as `sample_*` columns
(e.g. `sample_genotype`). A GSE that already has a samplesheet is skipped.

Then run `preprocessing/run_pipeline.py DATASET/<GSE>/samplesheet.csv <outdir>`; choose the `sample_*` columns for the final h5ad in `metadata_columns` of `preprocessing/config.yaml`.

## Whitelists

Currently, `whitelists/` only contains copies of nf-core/scrnaseq 4.2.0 `assets/whitelist/` (10x v2, v3, v4), the files STARsolo uses.

## Adding a chemistry

A new chemistry is a module in `chemistries/`, registered in `CHEMISTRIES` in `organize_fetchngs.py`.
What the module must define is listed in the docstring of `organize_fetchngs.py`; `tenx.py` is the example.
Read structures of many chemistries: [scg_lib_structs](https://teichlab.github.io/scg_lib_structs/), [Booeshaghi et al., *Bioinformatics* 2024](https://doi.org/10.1093/bioinformatics/btae168).
