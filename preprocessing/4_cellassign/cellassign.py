#!/usr/bin/env python3
"""Assigns a cell type to every cell of scdownstream's finalized h5ad with CellAssign (scvi-tools),
from the marker genes of config.yaml, and writes <outdir>/results/cellassign/cell_types.tsv.
CellAssign has no 'other' class: cells of types missing from the marker file go to the closest listed type."""

import argparse
import os
import shutil
from pathlib import Path

# Slurm jobs inherit OMP_THREAD_LIMIT=1, which would keep PyTorch on one core; set before torch is imported
N_CPUS = os.environ["SLURM_CPUS_PER_TASK"]
os.environ["OMP_THREAD_LIMIT"] = N_CPUS
os.environ["OMP_NUM_THREADS"] = N_CPUS

import anndata  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse  # noqa: E402
import scvi  # noqa: E402
import torch  # noqa: E402
import yaml  # noqa: E402
from scvi.external import CellAssign  # noqa: E402

PREPROCESSING_DIR = Path(__file__).resolve().parents[1]
# the decontX counts scdownstream analysed
LAYER = "ambient_corrected_decontx"
# fits from seeds 0, 1, ...; the one with the lowest negative log-likelihood is kept
N_RESTARTS = 5
ROUNDING_SEED = 0
# fits that give at least this share of cells the same type as the best fit found the same optimum
AGREEMENT = 0.95


def finalized_h5ad(outdir):
    """Returns the path of the h5ad in scdownstream's 09_finalized folder (not the _cellxgene one)."""
    paths = []
    for path in (outdir / "results" / "scdownstream" / "09_finalized").glob("*.h5ad"):
        if not path.name.endswith("_cellxgene.h5ad"):
            paths.append(path)
    if len(paths) != 1:
        raise ValueError(f"expected one finalized h5ad, found {paths}")
    return paths[0]


def marker_matrix(markers, var_names):
    """Builds the gene x cell type matrix (1 = marker) of the markers that are in the data.
    Prints the missing markers; a cell type without any marker in the data is an error."""
    present_genes = set(var_names)
    marker_genes = []
    for cell_type, genes in markers.items():
        missing = [gene for gene in genes if gene not in present_genes]
        if missing:
            print(f"WARNING: {cell_type}: markers not in the data: {', '.join(missing)}")
        if len(missing) == len(genes):
            raise ValueError(f"no marker of {cell_type} is in the data; wrong species or gene names?")
        for gene in genes:
            if gene in present_genes and gene not in marker_genes:
                marker_genes.append(gene)

    matrix = pd.DataFrame(0, index=marker_genes, columns=list(markers))
    for cell_type, genes in markers.items():
        for gene in genes:
            if gene in present_genes:
                matrix.loc[gene, cell_type] = 1
    return matrix


def stochastic_round(counts):
    """Rounds the fractional decontX counts to the integers CellAssign needs, up with probability equal
    to the fractional part, which keeps the expected counts."""
    counts = scipy.sparse.csr_matrix(counts, dtype=np.float64, copy=True)
    floor = np.floor(counts.data)
    rng = np.random.default_rng(ROUNDING_SEED)
    counts.data = floor + (rng.random(counts.data.shape) < counts.data - floor)
    counts.eliminate_zeros()
    return counts.astype(np.float32)


@torch.no_grad()
def mean_nll(model, adata):
    """Returns the negative log-likelihood per cell, which compares fits better than the training loss:
    the loss can favour a fit that merges two cell types."""
    total = 0.0
    for tensors in model._make_data_loader(adata=adata, batch_size=4096):
        _, generative, _ = model.module(tensors)
        total += float(-torch.logsumexp(generative["p_x_c"], dim=1).sum())
    return total / adata.n_obs


def train_best(bdata, matrix):
    """Fits CellAssign N_RESTARTS times and returns the fit with the lowest negative log-likelihood.
    Prints how many cells each fit types as the best fit does."""
    models = []
    nlls = []
    for seed in range(N_RESTARTS):
        scvi.settings.seed = seed
        model = CellAssign(bdata, matrix)
        model.train(accelerator="cpu", enable_progress_bar=False)
        models.append(model)
        nlls.append(mean_nll(model, bdata))

    best = int(np.argmin(nlls))
    best_types = models[best].predict().idxmax(axis=1)
    n_agreeing = 0
    for seed, model in enumerate(models):
        agreement = (model.predict().idxmax(axis=1) == best_types).mean()
        print(f"seed {seed}: NLL {nlls[seed]:.4f}, same type as the best fit for {agreement:.1%} of cells")
        if agreement >= AGREEMENT:
            n_agreeing += 1
    # the best fit always agrees with itself
    if n_agreeing < 2:
        print(
            "WARNING: no other fit agrees with the best one, so a better optimum may have been missed. "
            "Check that no expected cell type is (nearly) empty."
        )
    return models[best]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path, help="Output folder of the run")
    args = parser.parse_args()

    config = yaml.safe_load((args.outdir / "config.yaml").read_text())
    markers_path = PREPROCESSING_DIR / config["marker_genes"]
    # {cell type: [genes]}, e.g. {"Podocyte": ["Nphs1", "Nphs2"]}
    markers = yaml.safe_load(markers_path.read_text())
    adata = anndata.read_h5ad(finalized_h5ad(args.outdir))
    matrix = marker_matrix(markers, adata.var_names)
    print(f"{len(matrix)} marker genes in the data")

    # size factors over all genes, not only the markers
    library_size = np.asarray(adata.layers[LAYER].sum(axis=1)).ravel()
    bdata = anndata.AnnData(
        X=stochastic_round(adata[:, matrix.index].layers[LAYER]),
        obs=pd.DataFrame({"size_factor": library_size / library_size.mean()}, index=adata.obs_names),
        var=pd.DataFrame(index=matrix.index),
    )
    CellAssign.setup_anndata(bdata, size_factor_key="size_factor")
    torch.set_num_threads(int(N_CPUS))
    model = train_best(bdata, matrix)

    probabilities = model.predict()
    probabilities.index = adata.obs_names
    probabilities.columns = probabilities.columns.astype(str)
    cell_types = pd.DataFrame(
        {
            "sample": adata.obs["sample"].to_numpy(),
            "cell_type": probabilities.idxmax(axis=1),
            "probability": probabilities.max(axis=1),
        },
        index=adata.obs_names,
    )
    cell_types = cell_types.join(probabilities)
    cell_types.index.name = "cell"
    print(pd.crosstab(cell_types["sample"], cell_types["cell_type"]).to_string())

    cellassign_dir = args.outdir / "results" / "cellassign"
    cellassign_dir.mkdir(parents=True, exist_ok=True)
    cell_types.to_csv(cellassign_dir / "cell_types.tsv", sep="\t", float_format="%.4g")
    shutil.copyfile(markers_path, cellassign_dir / "markers.yaml")


if __name__ == "__main__":
    main()
