#!/usr/bin/env python3
"""Writes the scdownstream input of one run from its scrnaseq results: <outdir>/scdownstream/samplesheet.csv
and cellbender_priors.tsv with CellBender's --expected-cells and --total-droplets-included of each sample."""

import argparse
from pathlib import Path

import pandas as pd
import yaml

import cellbender_priors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path, help="Output folder of the run")
    args = parser.parse_args()

    config = yaml.safe_load((args.outdir / "config.yaml").read_text())
    qc = config["qc"] | config["modes"][config["mode"]]["qc"]
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
            # every library is its own batch
            "batch_col": "sample",
            # an empty value here would silently turn the global setting off for the sample
            "use_ambient_corrected_counts": "true",
        }
        samplesheet_rows.append(row | qc)

    scdownstream_dir = args.outdir / "scdownstream"
    pd.DataFrame(samplesheet_rows).to_csv(scdownstream_dir / "samplesheet.csv", index=False)
    # read by scdownstream.config, which passes the priors on to CellBender
    pd.DataFrame(priors_rows).to_csv(scdownstream_dir / "cellbender_priors.tsv", sep="\t", index=False)

    for priors in priors_rows:
        if priors["warnings"] != "none":
            print(
                f"CHECK {priors['sample']}: {priors['warnings']} (the codes are explained in cellbender_priors.py)"
            )


if __name__ == "__main__":
    main()
