#!/usr/bin/env python3
"""Corrects the velocyto spliced/unspliced/ambiguous counts of every sample for ambient RNA with decontX.
The three layers are stacked into one matrix (each gene once per layer) and decontX runs once on it, with the
cells CellBender called in scdownstream as cells and every other droplet as background.
Writes <sample>_cells.h5ad (stacked counts) and <sample>_decontx.h5ad (corrected, with obs['decontX_contamination'])
to <outdir>/results/velocyto_decontx/, where finalize_adata.py reads them."""

import argparse
import os
import subprocess
from pathlib import Path

import anndata
import h5py
import numpy as np
import pandas as pd
import scipy.sparse

LAYERS = ["spliced", "unspliced", "ambiguous"]
# the same R and celda as scdownstream's decontX module; built by setup.sh
CONTAINER = Path(os.environ["CONTAINER_DIR"]) / "decontx-celda-1.26.0.sif"
R_SCRIPT = Path(__file__).resolve().parent / "velocity_decontx.R"


def plain_strings(values):
    """Returns values as an index of plain Python strings, the only string encoding anndataR can read."""
    return pd.Index(np.asarray(values, dtype=object), dtype=object)


def write_counts(path, matrix, obs_names, var):
    """Writes a droplets x features count matrix as an h5ad that anndataR can read."""
    var = var.copy()
    var.index = plain_strings(var.index)
    for column in var.columns:
        var[column] = np.asarray(var[column], dtype=object)
    obs = pd.DataFrame(index=plain_strings(obs_names))
    adata = anndata.AnnData(scipy.sparse.csr_matrix(matrix, dtype=np.float32), obs=obs, var=var)
    adata.write_h5ad(path)


def correct_sample(sample, velocyto_path, cellbender_path, output_dir):
    """Stacks the velocyto layers of one sample, splits them into cells and background and runs decontX."""
    with h5py.File(cellbender_path, "r") as handle:
        cells = anndata.io.read_elem(handle["obs"]).index
    with h5py.File(velocyto_path, "r") as handle:
        droplets = anndata.io.read_elem(handle["obs"]).index
        var = anndata.io.read_elem(handle["var"])
        layers = []
        for layer in LAYERS:
            layers.append(scipy.sparse.csr_matrix(anndata.io.read_elem(handle["layers"][layer])))

    # get_indexer gives -1 for a missing barcode, which would silently take the last droplet
    cell_rows = droplets.get_indexer(cells)
    if (cell_rows < 0).any():
        raise ValueError(f"{sample}: {(cell_rows < 0).sum()} CellBender barcodes are not in {velocyto_path}")

    stacked = scipy.sparse.hstack(layers, format="csr")
    del layers
    features = pd.DataFrame({
        "gene_ids": np.tile(var.index.to_numpy(dtype=object), len(LAYERS)),
        "gene_symbol": np.tile(var["gene_symbol"].to_numpy(dtype=object), len(LAYERS)),
        "layer": np.repeat(np.array(LAYERS, dtype=object), len(var)),
    })
    features.index = features["gene_ids"] + "_" + features["layer"]

    is_cell = np.zeros(stacked.shape[0], dtype=bool)
    is_cell[cell_rows] = True
    droplet_counts = np.asarray(stacked.sum(axis=1)).ravel()
    background_rows = np.flatnonzero((droplet_counts > 0) & ~is_cell)
    print(f"{sample}: {len(cell_rows)} cells, {len(background_rows)} background droplets, {len(features)} features")

    cells_path = output_dir / f"{sample}_cells.h5ad"
    background_path = output_dir / f"{sample}_background.h5ad"
    output_path = output_dir / f"{sample}_decontx.h5ad"
    write_counts(cells_path, stacked[cell_rows], cells, features)
    write_counts(background_path, stacked[background_rows], droplets[background_rows], features[[]])
    del stacked

    binds = f"{output_dir},{R_SCRIPT.parent}"
    command = ["apptainer", "exec", "--bind", binds, CONTAINER, "Rscript", R_SCRIPT, cells_path, background_path, output_path]
    subprocess.run(command, check=True)
    background_path.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path, help="Output folder of the run")
    args = parser.parse_args()

    results = args.outdir / "results"
    output_dir = results / "velocyto_decontx"
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(args.outdir / "samplesheet.csv", dtype=str)["sample"].unique()
    for sample in samples:
        velocyto_path = results / "velocyto" / "star" / "mtx_conversions" / sample / f"{sample}_raw_matrix.h5ad"
        # the cells CellBender called (kept by save_intermediates)
        cellbender_path = results / "scdownstream" / "02_quality_control" / "empty_droplet_removal" / f"{sample}_subset.h5ad"
        correct_sample(sample, velocyto_path, cellbender_path, output_dir)


if __name__ == "__main__":
    main()
