#!/usr/bin/env python3
#SBATCH --job-name=cellassign
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
"""Assign cell types to the scdownstream output with CellAssign (scvi-tools).

Can be submitted directly with sbatch:

    sbatch --chdir=. cellassign.py

Defaults come from pipeline.conf ($PIPELINE_CONF if set, else ./pipeline.conf)
and downstream_processing_pipeline/nf-params.json:

    --markers   MARKER_GENES, relative to the project root (the working directory)
    --input     the one h5ad in $DOWNSTREAM_OUTDIR/09_finalized (not *_cellxgene.h5ad)
    --outdir    CELLASSIGN_OUTDIR
    --layer     ambient_corrected_<ambient_correction>, or raw_counts without
                ambient correction

The marker file is YAML, one cell type per line with its marker genes, named as
in the h5ad's var_names (gene symbols):

    Podocyte: [Nphs1, Nphs2, Podxl, Synpo]
    Endothelial: [Pecam1, Emcn, Kdr, Flt1]

Steps:
  1. Size factors from the counts layer over all genes (library size / mean).
  2. The layer is rounded stochastically to integers (decontX counts are
     fractional; CellAssign needs counts). Integer layers are left as they are.
  3. CellAssign is trained on the marker genes only, on CPU. A fit can end in a
     local optimum where two cell types merge, so it is fitted --n-restarts
     times from different seeds and the fit with the lowest negative
     log-likelihood is kept.
  4. Each cell gets the type with the highest probability. With
     --min-probability, cells whose highest probability is below it are
     labelled 'Unassigned'. CellAssign has no 'other' class, so cells of types
     missing from the marker file are forced into the closest listed type;
     check the probabilities.

Outputs in --outdir:
  cell_types.tsv          per cell: sample, cell_type, probability (the highest)
                          and one probability column per cell type
  cell_type_counts.tsv    cells per sample and cell type
  <input>_cellassign.h5ad the input with obs['cell_type'],
                          obs['cell_type_probability'],
                          obsm['cellassign_probabilities'] and uns['cellassign']
  restarts.tsv            NLL of each fit and how often it agrees with the best;
                          low agreement means the fit is unstable
  training_history.png    ELBO per epoch of the best fit; it should level off
  markers.yaml            copy of the marker file used
  marker_matrix.tsv       the gene x cell type matrix given to CellAssign
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys

# sbatch jobs inherit OMP_THREAD_LIMIT=1 on some clusters, which keeps PyTorch
# on one core; set it to the job's CPUs before torch is imported
N_CPUS = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
os.environ["OMP_THREAD_LIMIT"] = str(N_CPUS)
os.environ["OMP_NUM_THREADS"] = str(N_CPUS)

import anndata  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse  # noqa: E402
import torch  # noqa: E402
import yaml  # noqa: E402
from matplotlib import pyplot  # noqa: E402

import scvi  # noqa: E402
from scvi.external import CellAssign  # noqa: E402

DEFAULT_CONF = os.environ.get("PIPELINE_CONF", "pipeline.conf")
DEFAULT_PARAMS = "downstream_processing_pipeline/nf-params.json"
UNASSIGNED = "Unassigned"
# Fixed, not --seed: fits with different seeds (also in separate jobs) then see
# the same counts, so their NLLs can be compared
ROUNDING_SEED = 0
AGREEMENT = 0.95  # fits with at least this share of cells typed as in the best fit count as the same optimum


def read_conf(conf_path, names):
    """Variables from pipeline.conf, sourced by bash so ${PROJECT_ROOT} expands."""
    if not os.path.exists(conf_path):
        sys.exit(f"{conf_path} not found; give --pipeline-conf or the paths directly")
    script = 'source "$1" && shift && for name; do printf "%s\\0" "${!name-}"; done'
    result = subprocess.run(["bash", "-c", script, "_", conf_path, *names],
                            capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"could not read {conf_path}: {result.stderr.strip()}")
    return dict(zip(names, result.stdout.split("\0")))


def default_layer(params_path):
    """The counts layer scdownstream analysed, from its nf-params.json."""
    with open(params_path) as handle:
        params = json.load(handle)
    method = params.get("ambient_correction", "none")
    if params.get("use_ambient_corrected_counts") and method not in ("", "none"):
        return f"ambient_corrected_{method}"
    return "raw_counts"


def find_input(downstream_outdir):
    """The finalized h5ad of the scdownstream run."""
    pattern = os.path.join(downstream_outdir, "09_finalized", "*.h5ad")
    candidates = [p for p in sorted(glob.glob(pattern)) if not p.endswith("_cellxgene.h5ad")]
    if len(candidates) != 1:
        found = ", ".join(candidates) or "none"
        sys.exit(f"expected one h5ad matching {pattern} (found: {found}); give --input")
    return candidates[0]


def read_markers(path):
    """{cell type: [genes]} from a YAML file, checked for shape."""
    with open(path) as handle:
        markers = yaml.safe_load(handle)
    if not isinstance(markers, dict) or not markers:
        sys.exit(f"{path}: expected 'Cell type: [gene, ...]' lines")
    for cell_type, genes in markers.items():
        if not isinstance(genes, list) or not genes or not all(isinstance(g, str) for g in genes):
            sys.exit(f"{path}: '{cell_type}' needs a non-empty list of gene names")
        if cell_type == UNASSIGNED:
            sys.exit(f"{path}: '{UNASSIGNED}' is reserved for low-probability cells")
    return {str(cell_type): list(dict.fromkeys(genes)) for cell_type, genes in markers.items()}


def marker_matrix(markers, var_names):
    """Binary gene x cell type matrix of the markers found in var_names.

    Missing genes are reported; a cell type with no markers left is an error.
    """
    present = set(var_names)
    missing = {t: [g for g in genes if g not in present] for t, genes in markers.items()}
    n_missing = sum(len(g) for g in missing.values())
    n_total = sum(len(g) for g in markers.values())
    for cell_type, genes in missing.items():
        if genes:
            print(f"WARNING: {cell_type}: markers not in the data: {', '.join(genes)}",
                  file=sys.stderr)

    if n_missing > n_total / 2:
        # Most likely the wrong species or gene ID type
        lower = {v.lower(): v for v in var_names}
        case_hits = [g for genes in missing.values() for g in genes if g.lower() in lower]
        hint = (f"; {len(case_hits)} match with different capitalisation "
                f"(e.g. {case_hits[0]} -> {lower[case_hits[0].lower()]}), wrong species?"
                if case_hits else "")
        sys.exit(f"{n_missing} of {n_total} marker genes are not in the data (var_names look "
                 f"like {', '.join(var_names[:3])}){hint}")

    empty = [t for t, genes in markers.items() if len(genes) == len(missing[t])]
    if empty:
        sys.exit(f"no markers in the data for: {', '.join(empty)}")

    genes = sorted({g for genes in markers.values() for g in genes if g in present})
    matrix = pd.DataFrame(0, index=genes, columns=list(markers), dtype=int)
    for cell_type, type_genes in markers.items():
        matrix.loc[[g for g in type_genes if g in present], cell_type] = 1
    return matrix


def stochastic_round(X, seed):
    """Round each value up with probability equal to its fractional part (unbiased)."""
    X = scipy.sparse.csr_matrix(X, dtype=np.float64, copy=True)
    floor = np.floor(X.data)
    if np.array_equal(floor, X.data):
        return X.astype(np.float32)
    rng = np.random.default_rng(seed)
    X.data = floor + (rng.random(X.data.shape) < (X.data - floor))
    X.eliminate_zeros()
    return X.astype(np.float32)


@torch.no_grad()
def mean_nll(model, adata):
    """Negative log-likelihood per cell, -log sum_c p(x, c), comparable between fits.

    The training loss is not used for this: it is an EM-style objective that
    does not penalise wrong hard assignments, so a fit that merges cell types
    can score better on it.
    """
    total, n = 0.0, 0
    for tensors in model._make_data_loader(adata=adata, batch_size=4096):
        _, generative, _ = model.module(tensors)
        total += float(-torch.logsumexp(generative["p_x_c"], dim=1).sum())
        n += tensors["X"].shape[0]
    return total / n


def train_best(bdata, matrix, n_restarts, seed):
    """Fit CellAssign n_restarts times from different random starts; keep the lowest NLL.

    The fit depends on its random initialisation and can end in a local optimum
    where two cell types merge into one, so several fits are run and the one
    with the highest marginal likelihood is kept. This is the original CellAssign
    R package's num_runs option (default 1 there); scvi-tools has no equivalent.
    """
    fits, best_model = [], None
    for restart in range(n_restarts):
        scvi.settings.seed = seed + restart
        model = CellAssign(bdata, matrix)
        model.train(accelerator="cpu", enable_progress_bar=False)
        nll = mean_nll(model, bdata)
        print(f"fit {restart + 1}/{n_restarts} (seed {seed + restart}): NLL {nll:.4f}, "
              f"{len(model.history['elbo_train'])} epochs")
        fit = {"seed": seed + restart, "nll": nll,
               "labels": model.predict().idxmax(axis=1).to_numpy()}
        if not fits or nll < min(f["nll"] for f in fits):
            best_model = model
        fits.append(fit)

    best = min(fits, key=lambda f: f["nll"])
    restarts = pd.DataFrame({"seed": [f["seed"] for f in fits],
                             "nll": [f["nll"] for f in fits],
                             "agreement_with_best": [(f["labels"] == best["labels"]).mean()
                                                     for f in fits],
                             "best": [f is best for f in fits]})
    print(restarts.to_string(index=False))
    n_same = int((restarts["agreement_with_best"] >= AGREEMENT).sum())
    if n_restarts > 1 and n_same < 2:
        print(f"WARNING: no other fit agrees with the best one on >= {AGREEMENT:.0%} of cells; "
              f"the fits end in different optima, so a better one may have been missed. "
              f"Check that no expected cell type is (nearly) empty, and rerun with more "
              f"--n-restarts", file=sys.stderr)
    return best_model, restarts


def plot_history(model, path):
    figure, axis = pyplot.subplots(figsize=(6, 4), layout="constrained")
    for key, label in (("elbo_train", "train"), ("elbo_validation", "validation")):
        if key in model.history:
            axis.plot(model.history[key].index, model.history[key].iloc[:, 0], label=label)
    axis.set_xlabel("Epoch")
    axis.set_ylabel("ELBO")
    axis.legend()
    figure.savefig(path, dpi=150)
    pyplot.close(figure)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--pipeline-conf", default=DEFAULT_CONF)
    parser.add_argument("--params", default=DEFAULT_PARAMS,
                        help="scdownstream nf-params.json, for the default --layer")
    parser.add_argument("--markers", help="marker gene YAML (default: MARKER_GENES)")
    parser.add_argument("--input", help="h5ad (default: the scdownstream finalized h5ad)")
    parser.add_argument("--outdir", help="output directory (default: CELLASSIGN_OUTDIR)")
    parser.add_argument("--layer", help="counts layer (default: from --params)")
    parser.add_argument("--sample-key", default="sample", help="obs column with the sample")
    parser.add_argument("--min-probability", type=float, default=0.0,
                        help="label cells below this highest probability 'Unassigned' "
                             "(default: 0, every cell gets a type)")
    parser.add_argument("--n-restarts", type=int, default=3,
                        help="fits from different random starts; the lowest NLL is kept "
                             "(default: 3)")
    parser.add_argument("--seed", type=int, default=0, help="seed of the first fit; "
                        "fit i uses seed + i (default: 0)")
    args = parser.parse_args(argv)

    if not args.markers or not args.input or not args.outdir:
        conf = read_conf(args.pipeline_conf, ["MARKER_GENES", "DOWNSTREAM_OUTDIR",
                                              "CELLASSIGN_OUTDIR"])
        args.markers = args.markers or conf["MARKER_GENES"]
        args.input = args.input or find_input(conf["DOWNSTREAM_OUTDIR"])
        args.outdir = args.outdir or conf["CELLASSIGN_OUTDIR"]
    if not args.markers:
        sys.exit("no marker file: set MARKER_GENES in pipeline.conf or give --markers")
    if not args.outdir:
        sys.exit("no output directory: set CELLASSIGN_OUTDIR in pipeline.conf or give --outdir")
    args.layer = args.layer or default_layer(args.params)
    if args.n_restarts < 1:
        sys.exit("--n-restarts must be at least 1")
    if not 0 <= args.min_probability <= 1:
        sys.exit("--min-probability must be between 0 and 1")

    markers = read_markers(args.markers)
    scvi.settings.seed = args.seed
    torch.set_num_threads(N_CPUS)
    print(f"markers: {args.markers} ({len(markers)} cell types)")
    print(f"input:   {args.input}")
    print(f"layer:   {args.layer}")
    print(f"threads: {N_CPUS}")

    adata = anndata.read_h5ad(args.input)
    if args.layer not in adata.layers:
        sys.exit(f"layer '{args.layer}' not in {args.input} (layers: {', '.join(str(k) for k in adata.layers.keys() if k)})")
    if not adata.var_names.is_unique:
        sys.exit(f"{args.input}: var_names are not unique, so marker genes are ambiguous")
    matrix = marker_matrix(markers, adata.var_names)
    print(f"{len(matrix)} marker genes in the data")

    # Size factors over all genes, before subsetting to the markers
    counts = adata.layers[args.layer]
    library_size = np.asarray(counts.sum(axis=1)).ravel()
    if (library_size <= 0).any():
        sys.exit(f"{(library_size <= 0).sum()} cells have no counts in layer '{args.layer}'")

    bdata = anndata.AnnData(
        X=stochastic_round(adata[:, matrix.index].layers[args.layer], ROUNDING_SEED),
        obs=pd.DataFrame({"size_factor": library_size / library_size.mean()},
                         index=adata.obs_names),
        var=pd.DataFrame(index=matrix.index),
    )
    CellAssign.setup_anndata(bdata, size_factor_key="size_factor")
    model, restarts = train_best(bdata, matrix, args.n_restarts, args.seed)

    probabilities = model.predict()
    probabilities.index = adata.obs_names
    probabilities.columns = probabilities.columns.astype(str)
    best = probabilities.max(axis=1)
    cell_type = probabilities.idxmax(axis=1)
    categories = list(matrix.columns)
    if args.min_probability > 0:
        cell_type[best < args.min_probability] = UNASSIGNED
        categories.append(UNASSIGNED)
    cell_type = pd.Categorical(cell_type, categories=categories)

    os.makedirs(args.outdir, exist_ok=True)
    out = lambda name: os.path.join(args.outdir, name)  # noqa: E731

    table = pd.DataFrame({"cell_type": cell_type, "probability": best.to_numpy()},
                         index=adata.obs_names)
    if args.sample_key in adata.obs:
        table.insert(0, "sample", adata.obs[args.sample_key].to_numpy())
    table = table.join(probabilities)
    table.index.name = "cell"
    table.to_csv(out("cell_types.tsv"), sep="\t", float_format="%.4g")

    group = table["sample"] if "sample" in table else pd.Series("all", index=table.index)
    type_counts = pd.crosstab(group, table["cell_type"], dropna=False)
    type_counts.to_csv(out("cell_type_counts.tsv"), sep="\t")
    print(type_counts.to_string())

    restarts.to_csv(out("restarts.tsv"), sep="\t", index=False, float_format="%.6g")
    plot_history(model, out("training_history.png"))
    shutil.copyfile(args.markers, out("markers.yaml"))
    matrix.to_csv(out("marker_matrix.tsv"), sep="\t")

    adata.obs["cell_type"] = cell_type
    adata.obs["cell_type_probability"] = best.to_numpy()
    adata.obsm["cellassign_probabilities"] = probabilities
    adata.uns["cellassign"] = {
        "marker_file": os.path.abspath(args.markers),
        "markers": {t: [g for g in genes if g in matrix.index] for t, genes in markers.items()},
        "layer": args.layer,
        "min_probability": args.min_probability,
        "seed": int(restarts.loc[restarts["best"], "seed"].iloc[0]),
        "n_restarts": args.n_restarts,
        "scvi_tools_version": scvi.__version__,
    }
    stem = os.path.splitext(os.path.basename(args.input))[0]
    h5ad_path = out(f"{stem}_cellassign.h5ad")
    adata.write_h5ad(h5ad_path, compression="gzip")
    print(f"wrote {args.outdir}")


if __name__ == "__main__":
    main()
