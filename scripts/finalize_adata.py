#!/usr/bin/env python3
#SBATCH --job-name=finalize_adata
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48gb
#SBATCH --time=04:00:00
"""Assemble the final AnnData from the scdownstream output and the steps around it.

Can be submitted directly with sbatch, from the project root:

    sbatch --chdir=. scripts/finalize_adata.py

Defaults come from pipeline.conf ($PIPELINE_CONF if set, else ./pipeline.conf):

    --input              the one h5ad in $DOWNSTREAM_OUTDIR/09_finalized (not *_cellxgene.h5ad)
    --downstream-outdir  DOWNSTREAM_OUTDIR, for the scdownstream intermediates
    --samplesheet        SAMPLESHEET, relative to the directory of pipeline.conf
    --cellassign-dir     CELLASSIGN_OUTDIR, only if MARKER_GENES is set
    --output             $FINAL_OUTDIR/<input name>.h5ad

scdownstream must have run with save_intermediates true: the doublet scores and
the neighbor graph are only in its intermediate files.

Changes to the scdownstream h5ad:
  X        log1p-normalised counts (was the analysis counts; the 'log1p' layer is
           moved here). Layers 'raw_counts' (STARsolo counts of the cells CellBender
           called) and 'ambient_corrected_<method>' are kept.
  obs      added: doublet_score (scDblFinder, from the per-sample intermediates),
           ambient_fraction_removed (1 - corrected / raw counts per cell),
           total_counts_raw and n_genes_by_counts_raw (the QC columns from
           scdownstream are computed on the corrected counts), cell_type and
           cell_type_probability (CellAssign), and every samplesheet column that is
           neither an nf-core/scrnaseq column nor a QC column, as per-sample metadata.
           renamed: <key>-<res>_leiden -> leiden_<res>,
           <key>-<res>:entropy -> leiden_<res>_sample_entropy (key e.g. pca-global).
           dropped: n_counts and n_genes (duplicates of total_counts and
           n_genes_by_counts), and the scdownstream placeholders donor, batch,
           condition and label; a samplesheet column of the same name replaces them.
  var      renamed: deviance:<column> -> <column>; dropped: intersection.
  obsm     renamed: X_<key>_umap -> X_umap, X_pca_pca -> X_pca;
           added: cellassign_probabilities.
  obsp     connectivities and distances: the graph scdownstream computed its UMAP
           and Leiden clusters from (05_cluster_dimred/*/neighbors/).
  uns      colors renamed with their obs column; added: neighbors, cellassign
           (the run_info.json of the CellAssign run), qc_thresholds (per sample, the
           defaults with the samplesheet's QC columns applied) and pipeline (provenance).
"""

import argparse
import datetime
import glob
import json
import os
import re
import subprocess
import sys

# sbatch runs a copy of this script from Slurm's spool directory, so look for
# pipeline_common.py next to the original and in scripts/ of the project root (cwd)
sys.path[:0] = [os.path.dirname(os.path.abspath(__file__)), os.path.join(os.getcwd(), "scripts")]
from pipeline_common import (DEFAULT_CONF, SAMPLE_QC_COLUMNS, conf_relative,  # noqa: E402
                             find_finalized_h5ad, read_conf, read_sample_columns, sample_qc)

import anndata  # noqa: E402
import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

DROP_OBS = ["n_counts", "n_genes", "donor", "batch", "condition", "label"]
DROP_VAR = ["intersection"]
OBSM_RENAMES = {"X_pca_pca": "X_pca"}
UNASSIGNED = "Unassigned"  # as in cellassign.py


def read_h5_elem(path, key):
    """One element (e.g. 'obs', 'obsp') of an h5ad without loading the rest."""
    with h5py.File(path, "r") as handle:
        if key not in handle:
            sys.exit(f"{path} has no '{key}'")
        return anndata.io.read_elem(handle[key])


def embedding_key(adata):
    """The scdownstream key of the clustering, e.g. 'pca-global', from X_<key>_umap."""
    keys = [m.group(1) for k in adata.obsm if (m := re.fullmatch(r"X_(.+)_umap", k))]
    if len(keys) != 1:
        sys.exit(f"expected one obsm X_<key>_umap, found: {', '.join(keys) or 'none'}")
    return keys[0]


def rename_clustering(name, key):
    """<key>-<res>_leiden[...] -> leiden_<res>[...], <key>-<res>:entropy -> leiden_<res>_sample_entropy."""
    match = re.fullmatch(rf"{re.escape(key)}-([0-9.]+)(?:_leiden)?(.*)", name)
    if not match:
        return name
    resolution, rest = match.groups()
    if rest == ":entropy":
        # scdownstream's entropy of the 'batch' column (the sample) within each cluster
        return f"leiden_{resolution}_sample_entropy"
    return f"leiden_{resolution}{rest}"


def rename_keys(adata, key):
    adata.obs = adata.obs.rename(columns=lambda c: rename_clustering(c, key))
    for old in list(adata.uns):
        new = rename_clustering(old, key)
        if new != old:
            adata.uns[new] = adata.uns.pop(old)
    renames = dict(OBSM_RENAMES)
    renames.update({f"X_{key}_umap": "X_umap", f"X_{key}_tsne": "X_tsne"})
    for old, new in renames.items():
        if old in adata.obsm:
            adata.obsm[new] = adata.obsm.pop(old)
    adata.var = adata.var.rename(columns=lambda c: c.removeprefix("deviance:"))
    adata.var = adata.var.drop(columns=[c for c in DROP_VAR if c in adata.var])


def add_neighbors(adata, downstream_outdir, key):
    """Copy the kNN graph scdownstream computed the UMAP and Leiden clusters from."""
    pattern = os.path.join(downstream_outdir, "05_cluster_dimred", "*", "neighbors",
                           f"{key}_neighbors.h5ad")
    paths = sorted(glob.glob(pattern))
    if len(paths) != 1:
        sys.exit(f"expected one file matching {pattern} (found: {', '.join(paths) or 'none'}); "
                 f"scdownstream needs save_intermediates true")
    path = paths[0]
    graph_cells = read_h5_elem(path, "obs").index
    if not graph_cells.is_unique or set(graph_cells) != set(adata.obs_names):
        sys.exit(f"the cells in {path} differ from the finalized h5ad")
    order = graph_cells.get_indexer(adata.obs_names)
    obsp = read_h5_elem(path, "obsp")
    for name in ("connectivities", "distances"):
        adata.obsp[name] = obsp[name].tocsr()[order][:, order]
    params = dict(read_h5_elem(path, "uns/neighbors").get("params", {}))
    params["use_rep"] = "X_pca"  # renamed from scdownstream's name
    adata.uns["neighbors"] = {"connectivities_key": "connectivities",
                              "distances_key": "distances", "params": params}
    print(f"neighbors: {path} ({params})")


def doublet_scores(adata, downstream_outdir):
    """scDblFinder scores from the per-sample intermediates, which still hold the doublets."""
    pattern = os.path.join(downstream_outdir, "02_quality_control", "doublet_detection",
                           "scdblfinder", "*_scdblfinder.h5ad")
    paths = sorted(glob.glob(pattern))
    if not paths:
        sys.exit(f"no files match {pattern}; scdownstream needs save_intermediates true "
                 f"and doublet_detection scdblfinder")
    scores = []
    for path in paths:
        obs = read_h5_elem(path, "obs")
        if "scdblfinder_score" not in obs:
            sys.exit(f"{path} has no obs['scdblfinder_score']")
        scores.append(obs["scdblfinder_score"])
    scores = pd.concat(scores)
    if not scores.index.is_unique:
        sys.exit(f"cell names repeat across the files matching {pattern}")
    missing = adata.obs_names.difference(scores.index)
    if len(missing):
        sys.exit(f"{len(missing)} cells have no scDblFinder score (e.g. {missing[0]}); "
                 f"are the intermediates from the same scdownstream run?")
    print(f"doublet scores: {len(paths)} files in {os.path.dirname(pattern)}")
    return scores.reindex(adata.obs_names).to_numpy(dtype=np.float64)


def row_sums(matrix):
    return np.asarray(matrix.sum(axis=1), dtype=np.float64).ravel()


def add_count_metrics(adata):
    raw = adata.layers["raw_counts"]
    raw_total = row_sums(raw)
    adata.obs["total_counts_raw"] = raw_total
    adata.obs["n_genes_by_counts_raw"] = np.asarray((raw > 0).sum(axis=1)).ravel()
    layers = [k for k in adata.layers if str(k).startswith("ambient_corrected_")]
    if len(layers) > 1:
        sys.exit(f"more than one ambient-corrected layer: {', '.join(layers)}")
    if layers:
        if (raw_total <= 0).any():
            sys.exit(f"{(raw_total <= 0).sum()} cells have no raw counts")
        adata.obs["ambient_fraction_removed"] = 1 - row_sums(adata.layers[layers[0]]) / raw_total
        print(f"ambient_fraction_removed: from {layers[0]} vs raw_counts, median "
              f"{adata.obs['ambient_fraction_removed'].median():.3f}")


def read_sample_metadata(path, samples):
    """Per-sample metadata: the samplesheet columns that are neither nf-core/scrnaseq's nor QC."""
    metadata = read_sample_columns(path)
    metadata = metadata.drop(columns=[c for c in SAMPLE_QC_COLUMNS if c in metadata])
    columns = list(metadata.columns)
    missing = sorted(set(samples) - set(metadata.index))
    if missing:
        sys.exit(f"samples not in {path}: {', '.join(missing)}")
    for column in columns:
        values = metadata[column]
        numbers = pd.to_numeric(values, errors="coerce")
        if numbers.notna().sum() == values.notna().sum():
            metadata[column] = numbers
        else:  # categories in samplesheet order
            metadata[column] = pd.Categorical(values, categories=list(dict.fromkeys(values.dropna())))
    return metadata


def add_sample_metadata(adata, metadata):
    clashes = [c for c in metadata.columns if c in adata.obs]
    if clashes:
        sys.exit(f"samplesheet columns clash with obs columns: {', '.join(clashes)}; rename them")
    sample = adata.obs["sample"].astype(str)
    for column in metadata.columns:
        values = metadata[column].reindex(sample).to_numpy()
        if isinstance(metadata[column].dtype, pd.CategoricalDtype):
            values = pd.Categorical(values, categories=metadata[column].cat.categories)
        adata.obs[column] = values
    if len(metadata.columns):
        print(f"sample metadata: {', '.join(metadata.columns)}")


def add_cellassign(adata, cellassign_dir):
    """Cell types from the CellAssign outputs (cell_types.tsv and run_info.json)."""
    table_path = os.path.join(cellassign_dir, "cell_types.tsv")
    info_path = os.path.join(cellassign_dir, "run_info.json")
    for path in (table_path, info_path):
        if not os.path.exists(path):
            sys.exit(f"{path} not found; did cellassign.py finish? (--no-cellassign skips it)")
    with open(info_path) as handle:
        info = json.load(handle)
    table = pd.read_csv(table_path, sep="\t", index_col="cell", keep_default_na=False,
                        dtype={"cell_type": str})
    if set(table.index) != set(adata.obs_names):
        sys.exit(f"the cells in {table_path} differ from the finalized h5ad; rerun cellassign.py")
    table = table.reindex(adata.obs_names)
    types = list(info["cell_types"])
    categories = types + ([UNASSIGNED] if info["min_probability"] > 0 else [])
    unknown = set(table["cell_type"]) - set(categories)
    if unknown:
        sys.exit(f"{table_path}: cell types not in {info_path}: {', '.join(sorted(unknown))}")
    adata.obs["cell_type"] = pd.Categorical(table["cell_type"], categories=categories)
    adata.obs["cell_type_probability"] = table["probability"].to_numpy(dtype=np.float64)
    adata.obsm["cellassign_probabilities"] = table[types].astype(np.float64)
    adata.uns["cellassign"] = info
    print(f"cell types: {table_path}")


def provenance(input_path, downstream_outdir, samplesheet, cellassign_dir):
    def text(pattern):  # newest file matching pattern, or ''
        paths = sorted(glob.glob(os.path.join(downstream_outdir, "pipeline_info", pattern)),
                       key=os.path.getmtime)
        if not paths:
            return ""
        with open(paths[-1]) as handle:
            return handle.read()

    commit = subprocess.run(["git", "describe", "--always", "--dirty"],
                            capture_output=True, text=True)
    return {
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "pipeline_commit": commit.stdout.strip() if commit.returncode == 0 else "unknown",
        "scdownstream_h5ad": os.path.abspath(input_path),
        "samplesheet": os.path.abspath(samplesheet),
        "cellassign_dir": os.path.abspath(cellassign_dir) if cellassign_dir else "",
        # JSON/YAML text: version strings such as 1.10 survive unchanged
        "scdownstream_params": text("params_*.json"),
        "software_versions": text("*software*versions.yml"),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--pipeline-conf", default=DEFAULT_CONF)
    parser.add_argument("--input", help="scdownstream h5ad (default: its finalized h5ad)")
    parser.add_argument("--downstream-outdir", help="default: DOWNSTREAM_OUTDIR")
    parser.add_argument("--samplesheet", help="default: SAMPLESHEET")
    parser.add_argument("--cellassign-dir",
                        help="default: CELLASSIGN_OUTDIR if MARKER_GENES is set")
    parser.add_argument("--no-cellassign", action="store_true", help="leave out cell types")
    parser.add_argument("--output", help="default: $FINAL_OUTDIR/<input name>.h5ad")
    args = parser.parse_args(argv)

    conf = read_conf(args.pipeline_conf, ["DOWNSTREAM_OUTDIR", "SAMPLESHEET", "MARKER_GENES",
                                          "CELLASSIGN_OUTDIR", "FINAL_OUTDIR"])
    downstream_outdir = args.downstream_outdir or conf["DOWNSTREAM_OUTDIR"]
    input_path = args.input or find_finalized_h5ad(downstream_outdir)
    samplesheet = args.samplesheet or conf_relative(conf["SAMPLESHEET"], args.pipeline_conf)
    cellassign_dir = None
    if not args.no_cellassign:
        cellassign_dir = args.cellassign_dir or (conf["CELLASSIGN_OUTDIR"]
                                                 if conf["MARKER_GENES"] else None)
    output = args.output
    if not output:
        if not conf["FINAL_OUTDIR"]:
            sys.exit("no output: set FINAL_OUTDIR in pipeline.conf or give --output")
        output = os.path.join(conf["FINAL_OUTDIR"], os.path.basename(input_path))
    if not samplesheet or not os.path.exists(samplesheet):
        sys.exit(f"samplesheet {samplesheet!r} not found: set SAMPLESHEET or give --samplesheet")
    print(f"input:       {input_path}")
    print(f"samplesheet: {samplesheet}")
    print(f"cellassign:  {cellassign_dir or 'skipped'}")

    # Check the samplesheet before loading the whole h5ad
    samples = read_h5_elem(input_path, "obs")["sample"].astype(str).unique()
    metadata = read_sample_metadata(samplesheet, samples)
    qc = sample_qc(samplesheet, samples)  # also warns about misnamed QC columns

    adata = anndata.read_h5ad(input_path)
    if not adata.obs_names.is_unique:
        sys.exit(f"{input_path}: cell names are not unique")
    key = embedding_key(adata)

    adata.X = adata.layers.pop("log1p")
    adata.obs = adata.obs.drop(columns=[c for c in DROP_OBS if c in adata.obs])
    rename_keys(adata, key)
    add_neighbors(adata, downstream_outdir, key)
    adata.obs["doublet_score"] = doublet_scores(adata, downstream_outdir)
    add_count_metrics(adata)
    if cellassign_dir:
        add_cellassign(adata, cellassign_dir)
    add_sample_metadata(adata, metadata)
    # As in the current samplesheet; edited after the scdownstream run, it no longer matches
    adata.uns["qc_thresholds"] = pd.DataFrame.from_dict(qc, orient="index").rename_axis("sample")
    adata.uns["pipeline"] = provenance(input_path, downstream_outdir, samplesheet, cellassign_dir)

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    adata.write_h5ad(output, compression="gzip")
    print(adata)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
