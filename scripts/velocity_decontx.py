#!/usr/bin/env python3
#SBATCH --job-name=velocity_decontx
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
"""Correct the velocyto spliced/unspliced/ambiguous counts for ambient RNA with decontX.

Can be submitted directly with sbatch, from the project root:

    sbatch --chdir=. scripts/velocity_decontx.py

Defaults come from pipeline.conf ($PIPELINE_CONF if set, else ./pipeline.conf):
VELOCYTO_OUTDIR, DOWNSTREAM_OUTDIR, SAMPLESHEET, CONTAINER_DIR and
VELOCYTO_DECONTX_OUTDIR (the output).

Per sample, the three layers are stacked into one matrix (each gene three
times, as <gene id>_spliced, _unspliced and _ambiguous) and decontX runs once
on it, as in scdownstream's CELDA_DECONTX module:
  cells       the barcodes CellBender called (scdownstream's
              02_quality_control/empty_droplet_removal/<sample>_subset.h5ad), the
              same cells scdownstream's decontX corrected
  background  every other droplet with counts in the velocyto raw matrix
decontX clusters the cells itself. With one contamination estimate per cell and
one ambient profile over all three layers, the spliced counts are corrected more
than the unspliced ones where the ambient RNA is mostly spliced (as in nuclei).

decontX runs in the decontx container ($CONTAINER_DIR/decontx-celda-1.26.0.sif,
built by setup.sh from containers/decontx), the image scdownstream's module uses.

Outputs in VELOCYTO_DECONTX_OUTDIR, read by finalize_adata.py:
  <sample>_cells.h5ad    stacked counts of the CellBender cells; var has
                         gene_ids, gene_symbol and layer
  <sample>_decontx.h5ad  the same cells and features, decontX-corrected, with
                         obs['decontX_contamination']
  run_info.json          container, samples and inputs
"""

import argparse
import datetime
import json
import os
import subprocess
import sys

# sbatch runs a copy of this script from Slurm's spool directory, so look for
# pipeline_common.py next to the original and in scripts/ of the project root (cwd)
sys.path[:0] = [os.path.dirname(os.path.abspath(__file__)), os.path.join(os.getcwd(), "scripts")]
from pipeline_common import DEFAULT_CONF, conf_relative, read_conf, read_sample_columns  # noqa: E402

import anndata  # noqa: E402
import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse  # noqa: E402

LAYERS = ("spliced", "unspliced", "ambiguous")
CONTAINER = "decontx-celda-1.26.0.sif"
DEFAULT_PARAMS = "velocyto_pipeline/nf-params.json"


def find_script(name):
    """A file in scripts/, next to the original script or of the project root (cwd)."""
    for directory in (os.path.dirname(os.path.abspath(__file__)), os.path.join(os.getcwd(), "scripts")):
        path = os.path.join(directory, name)
        if os.path.exists(path):
            return path
    sys.exit(f"{name} not found in scripts/; submit from the project root")


def object_index(values):
    # anndataR (in the container) cannot read the nullable-string encoding that
    # anndata writes for pandas' default string dtype, so write plain strings
    return pd.Index(np.asarray(values, dtype=object), dtype=object)


def write_counts(path, matrix, obs_names, var):
    var = var.copy()
    var.index = object_index(var.index)
    for column in var.columns:
        var[column] = np.asarray(var[column], dtype=object)
    adata = anndata.AnnData(scipy.sparse.csr_matrix(matrix, dtype=np.float32),
                            obs=pd.DataFrame(index=object_index(obs_names)), var=var)
    adata.write_h5ad(path)


def read_velocyto(path):
    """Droplet names, genes and the three layers (droplets x genes, CSR) of a velocyto h5ad."""
    with h5py.File(path, "r") as handle:
        missing = [k for k in LAYERS if k not in handle.get("layers", {})]
        if missing:
            sys.exit(f"{path} has no layer(s) {', '.join(missing)}; was it made by the velocyto run?")
        obs_names = anndata.io.read_elem(handle["obs"]).index
        var = anndata.io.read_elem(handle["var"])
        layers = {k: scipy.sparse.csr_matrix(anndata.io.read_elem(handle["layers"][k])) for k in LAYERS}
    if "gene_symbol" not in var:
        sys.exit(f"{path} has no var['gene_symbol']")
    return obs_names, var, layers


def process_sample(sample, velocyto_path, cellbender_path, outdir, container):
    with h5py.File(cellbender_path, "r") as handle:
        cells = anndata.io.read_elem(handle["obs"]).index
    obs_names, var, layers = read_velocyto(velocyto_path)
    rows = obs_names.get_indexer(cells)
    if (rows < 0).any():
        sys.exit(f"{sample}: {(rows < 0).sum()} CellBender barcodes (e.g. {cells[rows < 0][0]}) "
                 f"are not in {velocyto_path}")

    stacked = scipy.sparse.hstack([layers[k] for k in LAYERS], format="csr")
    del layers
    features = pd.DataFrame({
        "gene_ids": np.tile(var.index.to_numpy(dtype=object), len(LAYERS)),
        "gene_symbol": np.tile(var["gene_symbol"].to_numpy(dtype=object), len(LAYERS)),
        "layer": np.repeat(np.array(LAYERS, dtype=object), len(var)),
    })
    features.index = features["gene_ids"] + "_" + features["layer"]
    is_cell = np.zeros(stacked.shape[0], dtype=bool)
    is_cell[rows] = True
    totals = np.asarray(stacked.sum(axis=1)).ravel()
    background = np.flatnonzero((totals > 0) & ~is_cell)
    print(f"{sample}: {len(rows)} cells, {len(background)} background droplets, "
          f"{len(features)} features", flush=True)

    cells_path = os.path.join(outdir, f"{sample}_cells.h5ad")
    background_path = os.path.join(outdir, f"{sample}_background.h5ad")
    output_path = os.path.join(outdir, f"{sample}_decontx.h5ad")
    write_counts(cells_path, stacked[rows], cells, features)
    write_counts(background_path, stacked[background], obs_names[background],
                 features[[]])
    del stacked

    script = find_script("velocity_decontx.R")
    command = ["apptainer", "exec", "--bind", f"{outdir},{os.path.dirname(script)}",
               container, "Rscript", script, cells_path, background_path, output_path]
    print(" ".join(command), flush=True)
    subprocess.run(command, check=True)
    os.remove(background_path)

    # The R side matches cells to background droplets by name, so check the names survived
    corrected = anndata.read_h5ad(output_path, backed="r")
    if not (corrected.obs_names.equals(pd.Index(cells))
            and corrected.var_names.equals(pd.Index(features.index))):
        sys.exit(f"{output_path}: cell or feature names differ from {cells_path}")
    contamination = corrected.obs["decontX_contamination"]
    print(f"{sample}: median contamination {contamination.median():.3f}", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--pipeline-conf", default=DEFAULT_CONF)
    parser.add_argument("--params", default=DEFAULT_PARAMS,
                        help="velocyto nf-params.json, for the aligner")
    parser.add_argument("--samples", nargs="+", help="only these samples (default: all)")
    parser.add_argument("--outdir", help="default: VELOCYTO_DECONTX_OUTDIR")
    parser.add_argument("--container", help=f"default: $CONTAINER_DIR/{CONTAINER}")
    args = parser.parse_args(argv)

    conf = read_conf(args.pipeline_conf, ["VELOCYTO_OUTDIR", "DOWNSTREAM_OUTDIR", "SAMPLESHEET",
                                          "CONTAINER_DIR", "VELOCYTO_DECONTX_OUTDIR"])
    outdir = args.outdir or conf["VELOCYTO_DECONTX_OUTDIR"]
    if not outdir:
        sys.exit("no output directory: set VELOCYTO_DECONTX_OUTDIR in pipeline.conf or give --outdir")
    container = args.container or os.path.join(conf["CONTAINER_DIR"], CONTAINER)
    if not os.path.exists(container):
        sys.exit(f"{container} not found; run setup.sh to build it")
    with open(args.params) as handle:
        aligner = json.load(handle).get("aligner", "star")
    samplesheet = conf_relative(conf["SAMPLESHEET"], args.pipeline_conf)
    samples = args.samples or list(read_sample_columns(samplesheet).index)

    os.makedirs(outdir, exist_ok=True)
    inputs = {}
    for sample in samples:
        velocyto_path = os.path.join(conf["VELOCYTO_OUTDIR"], aligner, "mtx_conversions",
                                     sample, f"{sample}_raw_matrix.h5ad")
        cellbender_path = os.path.join(conf["DOWNSTREAM_OUTDIR"], "02_quality_control",
                                       "empty_droplet_removal", f"{sample}_subset.h5ad")
        for path in (velocyto_path, cellbender_path):
            if not os.path.exists(path):
                sys.exit(f"{path} not found (the scdownstream intermediates need save_intermediates true)")
        inputs[sample] = {"velocyto": velocyto_path, "cellbender_cells": cellbender_path}

    for sample in samples:
        process_sample(sample, inputs[sample]["velocyto"], inputs[sample]["cellbender_cells"],
                       outdir, container)

    info_path = os.path.join(outdir, "run_info.json")
    info = {}
    if os.path.exists(info_path):  # keep the entries of samples not rerun
        with open(info_path) as handle:
            info = json.load(handle)
    info.setdefault("samples", {}).update(inputs)
    info.update({"container": os.path.abspath(container), "layers": list(LAYERS),
                 "updated": datetime.datetime.now().isoformat(timespec="seconds")})
    with open(info_path, "w") as handle:
        json.dump(info, handle, indent=2)
    print(f"wrote {outdir}")


if __name__ == "__main__":
    main()
