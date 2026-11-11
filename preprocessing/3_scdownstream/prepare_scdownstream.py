#!/usr/bin/env python3
"""Writes the scdownstream input of one run from its scrnaseq results: <outdir>/scdownstream/samplesheet.csv
and cellbender_priors.tsv with CellBender's --expected-cells and --total-droplets-included of each sample."""

import argparse
from pathlib import Path

import pandas as pd
import yaml

import cellbender_priors

# QC thresholds of every sample; the qc of the run's mode in config.yaml overrides them.
# Lenient on purpose: only empty droplets and debris are removed. 0 turns a *_nmads filter off.
QC_DEFAULTS = {
    "min_genes": 200,
    "min_counts_cell": 200,
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path, help="Output folder of the run")
    args = parser.parse_args()

    config = yaml.safe_load((args.outdir / "config.yaml").read_text())
    qc = QC_DEFAULTS | config["modes"][config["mode"]]["qc"]
    samples = pd.read_csv(args.outdir / "samplesheet.csv", dtype=str)["sample"].unique()
    matrices = args.outdir / "results" / "starsolo" / "star" / "mtx_conversions"

    samplesheet_rows = []
    priors_rows = []
    for sample in samples:
        raw_matrix = matrices / sample / f"{sample}_raw_matrix.h5ad"
        filtered_matrix = matrices / sample / f"{sample}_filtered_matrix.h5ad"
        priors = cellbender_priors.priors_for_sample(sample, raw_matrix, filtered_matrix)
        priors_rows.append(priors)
        print(
            f"{sample}: EmptyDrops cells {priors['emptydrops_cells']}, expected_cells {priors['expected_cells']}, "
            f"total_droplets {priors['total_droplets']} ({priors['method']})"
        )

        row = {
            "sample": sample,
            # only the raw matrix, so that CellBender calls the cells
            "unfiltered": raw_matrix,
            # index genes on symbols, not Ensembl IDs
            "symbol_col": "gene_symbol",
            # one library is one animal
            "batch_col": "sample",
            "donor_col": "sample",
            "expected_cells": priors["expected_cells"],
            # an empty value here would silently turn the global setting off for the sample
            "use_ambient_corrected_counts": "true",
        }
        samplesheet_rows.append(row | qc)

    scdownstream_dir = args.outdir / "scdownstream"
    pd.DataFrame(samplesheet_rows).to_csv(scdownstream_dir / "samplesheet.csv", index=False)
    # read by scdownstream.config, because scdownstream passes only expected_cells on to CellBender
    pd.DataFrame(priors_rows).to_csv(scdownstream_dir / "cellbender_priors.tsv", sep="\t", index=False)

    for priors in priors_rows:
        if priors["warnings"] != "none":
            print(f"CHECK {priors['sample']}: {priors['warnings']} (the codes are explained in cellbender_priors.py)")


if __name__ == "__main__":
    main()
