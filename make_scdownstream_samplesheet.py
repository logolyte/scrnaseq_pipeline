#!/usr/bin/env python3
#SBATCH --job-name=downstream_samplesheet
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --time=00:10:00
#SBATCH --account=ag-beyer
"""Build the nf-core/scdownstream samplesheet from the nf-core/scrnaseq run.

Submitted directly, no wrapper needed, since the #SBATCH directives above are
ordinary Python comments:

    sbatch --chdir=. --dependency=afterok:$STAR_JOB make_downstream_samplesheet.py

Only the raw matrices are listed. scdownstream branches on whether 'filtered'
is set: when it is, the sample is treated as already filtered and CellBender
never runs empty-droplet removal. Giving it just 'unfiltered' is what makes
CellBender call the cells, using the per-sample expected_cells below.

Run from the project root; run_job.sh guarantees that with sbatch --chdir.

    starsolo_pipeline/nf-params.json               <- read
    downstream_processing_pipeline/samplesheet.csv <- written
"""

import argparse
import csv
import glob
import gzip
import json
import os
import sys

DEFAULT_PARAMS = "starsolo_pipeline/nf-params.json"
DEFAULT_OUTPUT = "downstream_processing_pipeline/samplesheet.csv"

# var["gene_symbol"] exists in every scrnaseq h5ad; naming it here makes
# scdownstream index on symbols rather than unversioned Ensembl IDs.
SYMBOL_COL = "gene_symbol"

# mtx_to_h5ad_star.py writes obs["sample"] per file, and one library is one
# mouse, so it serves as both the batch key and the biological replicate.
BATCH_COL = "sample"
DONOR_COL = "sample"

COLUMNS = ["sample", "unfiltered", "symbol_col", "batch_col", "donor_col",
           "expected_cells"]


def read_config(params_path):
    """Read the scrnaseq params; paths inside are relative to that file."""
    with open(params_path) as handle:
        config = json.load(handle)

    samplesheet = config.get("input", "samplesheet.csv")
    if not os.path.isabs(samplesheet):
        samplesheet = os.path.join(os.path.dirname(os.path.abspath(params_path)),
                                   samplesheet)

    return config["outdir"].rstrip("/"), config.get("aligner", "star"), samplesheet


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

        rows.append({
            "sample": sample,
            "unfiltered": os.path.abspath(unfiltered),
            "symbol_col": SYMBOL_COL,
            "batch_col": BATCH_COL,
            "donor_col": DONOR_COL,
            "expected_cells": cells,
        })

    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--params", default=DEFAULT_PARAMS)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    outdir, aligner, samplesheet = read_config(args.params)
    rows = build_rows(outdir, aligner, read_samples(samplesheet))

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {args.output} with {len(rows)} samples", file=sys.stderr)


if __name__ == "__main__":
    main()
