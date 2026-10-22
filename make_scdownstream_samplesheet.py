#!/usr/bin/env python3
#SBATCH --job-name=downstream_samplesheet
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --time=00:10:00
"""Build the nf-core/scdownstream samplesheet from the nf-core/scrnaseq run.

Submitted directly, no wrapper needed, since the #SBATCH directives above are
ordinary Python comments:

    sbatch --chdir=. --dependency=afterok:$STAR_JOB make_scdownstream_samplesheet.py

Without --starsolo-outdir, STARSOLO_OUTDIR is taken from pipeline.conf
($PIPELINE_CONF if set, else ./pipeline.conf).

Only the raw matrices are listed. scdownstream branches on whether 'filtered'
is set: when it is, the sample is treated as already filtered and CellBender
never runs empty-droplet removal. Giving it just 'unfiltered' is what makes
CellBender call the cells, using the per-sample expected_cells below.

Run from the project root; run_job.sh guarantees that with sbatch --chdir.

    pipeline.conf                                  <- read (STARSOLO_OUTDIR)
    starsolo_pipeline/nf-params.json               <- read (aligner, input)
    --starsolo-outdir                              <- read (STARsolo output)
    downstream_processing_pipeline/samplesheet.csv <- written
"""

import argparse
import csv
import glob
import gzip
import json
import os
import subprocess
import sys

DEFAULT_CONF = os.environ.get("PIPELINE_CONF", "pipeline.conf")
DEFAULT_PARAMS = "starsolo_pipeline/nf-params.json"
DEFAULT_OUTPUT = "downstream_processing_pipeline/samplesheet.csv"

# var["gene_symbol"] exists in every scrnaseq h5ad; naming it here makes
# scdownstream index on symbols rather than unversioned Ensembl IDs.
SYMBOL_COL = "gene_symbol"

# mtx_to_h5ad_star.py writes obs["sample"] per file, and one library is one
# mouse, so it serves as both the batch key and the biological replicate.
BATCH_COL = "sample"
DONOR_COL = "sample"

# Written explicitly per sample: on nf-core/scdownstream 083acc65 a missing
# value becomes [] in meta, which silently overrides the global
# --use_ambient_corrected_counts with false.
USE_AMBIENT_CORRECTED_COUNTS = "true"

# Only a 200-gene floor per nucleus; every other cell and gene filter in
# SCANPY_FILTER is set to its "off" value explicitly, so schema default
# changes cannot re-enable it.
QC = {
    "min_genes": 200,
    "min_counts_cell": 0,
    "min_cells": 0,
    "min_counts_gene": 0,
    "max_mito_percentage": 100,
    "min_ribo_percentage": 0,
    "max_hb_percentage": 100,
    "log1p_total_counts_nmads": 0,
    "log1p_n_genes_by_counts_nmads": 0,
    "pct_counts_in_top_20_genes_nmads": 0,
    "pct_counts_mt_nmads": 0,
}

COLUMNS = ["sample", "unfiltered", "symbol_col", "batch_col", "donor_col",
           "expected_cells", "use_ambient_corrected_counts"] + list(QC)


def read_starsolo_outdir(conf_path):
    """STARSOLO_OUTDIR from pipeline.conf, sourced by bash so ${PROJECT_ROOT} expands."""
    if not os.path.exists(conf_path):
        sys.exit(f"no --starsolo-outdir given and {conf_path} not found")
    result = subprocess.run(
        ["bash", "-c", 'source "$1" && printf %s "$STARSOLO_OUTDIR"', "_", conf_path],
        capture_output=True, text=True,
    )
    if result.returncode != 0 or not result.stdout:
        sys.exit(f"could not read STARSOLO_OUTDIR from {conf_path}: {result.stderr.strip()}")
    return result.stdout


def read_config(params_path):
    """Read aligner and samplesheet from the scrnaseq params; paths are relative to that file."""
    with open(params_path) as handle:
        config = json.load(handle)

    samplesheet = config.get("input", "samplesheet.csv")
    if not os.path.isabs(samplesheet):
        samplesheet = os.path.join(os.path.dirname(os.path.abspath(params_path)),
                                   samplesheet)

    return config.get("aligner", "star"), samplesheet


def read_samples(path):
    """Sample names in samplesheet order, de-duplicated across lanes."""
    with open(path) as handle:
        names = [row["sample"] for row in csv.DictReader(handle)]
    return list(dict.fromkeys(names))


def count_called_cells(outdir, aligner, sample):
    """Cells called by EmptyDrops, handed to CellBender as expected_cells."""
    pattern = os.path.join(outdir, aligner, sample, f"{sample}.Solo.out",
                           "Gene*", "filtered", "barcodes.tsv.gz")
    matches = glob.glob(pattern)
    if not matches:
        return ""
    with gzip.open(matches[0], "rt") as handle:
        return sum(1 for _ in handle)


def build_rows(outdir, aligner, samples):
    conversions = os.path.join(outdir, aligner, "mtx_conversions")
    rows = []

    for sample in samples:
        unfiltered = os.path.join(conversions, sample, f"{sample}_raw_matrix.h5ad")
        if not os.path.exists(unfiltered):
            sys.exit(f"missing raw matrix {unfiltered}; did the scrnaseq run finish?")

        cells = count_called_cells(outdir, aligner, sample)
        print(f"{sample}: expected_cells {cells or 'unknown'}", file=sys.stderr)

        row = {
            "sample": sample,
            "unfiltered": os.path.abspath(unfiltered),
            "symbol_col": SYMBOL_COL,
            "batch_col": BATCH_COL,
            "donor_col": DONOR_COL,
            "expected_cells": cells,
            "use_ambient_corrected_counts": USE_AMBIENT_CORRECTED_COUNTS,
        }
        row.update(QC)
        rows.append(row)

    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--starsolo-outdir",
                        help="outdir of the nf-core/scrnaseq run (default: STARSOLO_OUTDIR in --pipeline-conf)")
    parser.add_argument("--pipeline-conf", default=DEFAULT_CONF)
    parser.add_argument("--params", default=DEFAULT_PARAMS)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    outdir = args.starsolo_outdir or read_starsolo_outdir(args.pipeline_conf)
    outdir = outdir.rstrip("/")
    print(f"starsolo outdir: {outdir}", file=sys.stderr)
    aligner, samplesheet = read_config(args.params)
    rows = build_rows(outdir, aligner, read_samples(samplesheet))

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {args.output} with {len(rows)} samples", file=sys.stderr)


if __name__ == "__main__":
    main()
