#!/usr/bin/env python3
"""Assembles the final AnnData of one run from scdownstream's finalized h5ad and the steps around it,
and writes it to <outdir>/results/final/.

  X       log1p counts; layers raw_counts and ambient_corrected_decontx (with velocyto also spliced,
          unspliced, ambiguous and their _raw versions)
  obs     adds doublet_score, ambient_fraction_removed, total_counts_raw, n_genes_by_counts_raw, cell_type,
          cell_type_probability, velocity_contamination and the samplesheet metadata; drops scdownstream's
          duplicate and placeholder columns
  keys    scdownstream's names are shortened: leiden_<res>, X_umap, X_pca, var columns without 'deviance:'
  obsp    the kNN graph of the UMAP and Leiden clusters; uns['pipeline'] holds the settings and versions
"""

import argparse
import datetime
import os
import re
import subprocess
from pathlib import Path

import anndata
import h5py
import numpy as np
import pandas as pd
import scipy.sparse
import yaml

PREPROCESSING_DIR = Path(__file__).resolve().parents[1]
DROP_OBS = ["n_counts", "n_genes", "donor", "batch", "condition", "label"]
VELOCYTO_LAYERS = ["spliced", "unspliced", "ambiguous"]


def finalized_h5ad(outdir):
    """Returns the path of the h5ad in scdownstream's 09_finalized folder (not the _cellxgene one)."""
    paths = []
    for path in (outdir / "results" / "scdownstream" / "09_finalized").glob("*.h5ad"):
        if not path.name.endswith("_cellxgene.h5ad"):
            paths.append(path)
    if len(paths) != 1:
        raise ValueError(f"expected one finalized h5ad, found {paths}")
    return paths[0]


def read_h5ad_element(path, key):
    """Reads one element of an h5ad, e.g. 'obs' or 'obsp', without loading the rest."""
    with h5py.File(path, "r") as handle:
        return anndata.io.read_elem(handle[key])


def embedding_key(adata):
    """Returns scdownstream's key of the embedding, e.g. 'pca-global' from obsm 'X_pca-global_umap'."""
    keys = []
    for name in adata.obsm:
        match = re.fullmatch(r"X_(.+)_umap", name)
        if match:
            keys.append(match.group(1))
    if len(keys) != 1:
        raise ValueError(f"expected one obsm X_<key>_umap, found {keys}")
    return keys[0]


def rename_clustering(name, key):
    """Renames a clustering key of scdownstream: 'pca-global-0.5_leiden' -> 'leiden_0.5',
    'pca-global-0.5:entropy' -> 'leiden_0.5_sample_entropy'; other names stay."""
    match = re.fullmatch(rf"{re.escape(key)}-([0-9.]+)(?:_leiden)?(.*)", name)
    if not match:
        return name
    resolution, rest = match.groups()
    if rest == ":entropy":
        return f"leiden_{resolution}_sample_entropy"
    return f"leiden_{resolution}{rest}"


def rename_keys(adata, key):
    """Gives scdownstream's obs, uns, obsm and var keys simpler names."""
    adata.obs = adata.obs.rename(columns=lambda name: rename_clustering(name, key))
    for old_name in list(adata.uns):
        new_name = rename_clustering(old_name, key)
        if new_name != old_name:
            adata.uns[new_name] = adata.uns.pop(old_name)
    obsm_names = {"X_pca_pca": "X_pca", f"X_{key}_umap": "X_umap"}
    for old_name, new_name in obsm_names.items():
        adata.obsm[new_name] = adata.obsm.pop(old_name)
    adata.var = adata.var.rename(columns=lambda name: name.removeprefix("deviance:"))
    adata.var = adata.var.drop(columns=["intersection"])


def add_neighbors(adata, scdownstream_results, key):
    """Copies the kNN graph scdownstream computed the UMAP and Leiden clusters from."""
    paths = list((scdownstream_results / "05_cluster_dimred").glob(f"*/neighbors/{key}_neighbors.h5ad"))
    if len(paths) != 1:
        raise ValueError(f"expected one {key}_neighbors.h5ad, found {paths}")
    path = paths[0]
    graph_cells = read_h5ad_element(path, "obs").index
    if set(graph_cells) != set(adata.obs_names):
        raise ValueError(f"the cells in {path} differ from the finalized h5ad")
    order = graph_cells.get_indexer(adata.obs_names)
    obsp = read_h5ad_element(path, "obsp")
    for name in ["connectivities", "distances"]:
        adata.obsp[name] = obsp[name].tocsr()[order][:, order]
    params = dict(read_h5ad_element(path, "uns/neighbors")["params"])
    params["use_rep"] = "X_pca"
    adata.uns["neighbors"] = {"connectivities_key": "connectivities", "distances_key": "distances", "params": params}


def add_doublet_scores(adata, scdownstream_results):
    """Adds the scDblFinder score from the per-sample files, which still hold the removed doublets."""
    scores = []
    for path in (scdownstream_results / "02_quality_control" / "doublet_detection" / "scdblfinder").glob("*_scdblfinder.h5ad"):
        scores.append(read_h5ad_element(path, "obs")["scdblfinder_score"])
    scores = pd.concat(scores)
    missing = adata.obs_names.difference(scores.index)
    if len(missing):
        raise ValueError(f"{len(missing)} cells have no scDblFinder score")
    adata.obs["doublet_score"] = scores.reindex(adata.obs_names).to_numpy(dtype=np.float64)


def add_count_metrics(adata):
    """Adds the total counts and detected genes of the raw counts and the share decontX removed."""
    raw_counts = adata.layers["raw_counts"]
    raw_totals = np.asarray(raw_counts.sum(axis=1), dtype=np.float64).ravel()
    corrected_totals = np.asarray(adata.layers["ambient_corrected_decontx"].sum(axis=1), dtype=np.float64).ravel()
    adata.obs["total_counts_raw"] = raw_totals
    adata.obs["n_genes_by_counts_raw"] = np.asarray((raw_counts > 0).sum(axis=1)).ravel()
    adata.obs["ambient_fraction_removed"] = 1 - corrected_totals / raw_totals


def add_cell_types(adata, cell_types_path):
    """Adds the CellAssign cell type and its probability to obs and the probabilities of all types to obsm."""
    table = pd.read_csv(cell_types_path, sep="\t", index_col="cell", dtype={"cell_type": str})
    if set(table.index) != set(adata.obs_names):
        raise ValueError(f"the cells in {cell_types_path} differ from the finalized h5ad; rerun cellassign.py")
    table = table.reindex(adata.obs_names)
    # one column per cell type after 'probability', in marker file order
    cell_type_names = list(table.columns[table.columns.get_loc("probability") + 1 :])
    adata.obs["cell_type"] = pd.Categorical(table["cell_type"], categories=cell_type_names)
    adata.obs["cell_type_probability"] = table["probability"].to_numpy(dtype=np.float64)
    adata.obsm["cellassign_probabilities"] = table[cell_type_names].astype(np.float64)


def add_velocyto(adata, velocyto_dir):
    """Adds the velocyto layers, corrected and raw, from the velocity_decontx.py outputs.
    They come from a separate STAR run: cells are matched by sample and barcode, genes by symbol."""
    cell_samples = adata.obs["sample"].astype(str).to_numpy()
    blocks = {}
    for layer in VELOCYTO_LAYERS:
        blocks[layer] = []
        blocks[f"{layer}_raw"] = []
    positions = []
    contamination = np.full(adata.n_obs, np.nan)
    for sample in pd.unique(cell_samples):
        raw = anndata.read_h5ad(velocyto_dir / f"{sample}_cells.h5ad")
        corrected = anndata.read_h5ad(velocyto_dir / f"{sample}_decontx.h5ad")
        cell_positions = np.flatnonzero(cell_samples == sample)
        # scdownstream names cells <sample>_<barcode>
        barcodes = adata.obs_names[cell_positions].str.removeprefix(f"{sample}_")
        # get_indexer gives -1 for a missing name, which would silently take the last row
        rows = raw.obs_names.get_indexer(barcodes)
        if (rows < 0).any():
            raise ValueError(f"{sample}: {(rows < 0).sum()} cells are not in {velocyto_dir}/{sample}_cells.h5ad")
        raw_counts = raw.X.tocsr()[rows]
        corrected_counts = corrected.X.tocsr()[rows]
        for layer in VELOCYTO_LAYERS:
            layer_columns = np.flatnonzero(raw.var["layer"].to_numpy() == layer)
            symbols = pd.Index(raw.var["gene_symbol"].to_numpy()[layer_columns])
            columns = symbols.get_indexer(adata.var_names)
            if (columns < 0).any():
                raise ValueError(f"{(columns < 0).sum()} genes have no velocyto counts; were both runs made with the same GTF?")
            blocks[layer].append(corrected_counts[:, layer_columns[columns]])
            blocks[f"{layer}_raw"].append(raw_counts[:, layer_columns[columns]])
        contamination[cell_positions] = corrected.obs["decontX_contamination"].to_numpy()[rows]
        positions.append(cell_positions)

    # puts the cells, stacked sample by sample, back in the order of adata
    order = np.argsort(np.concatenate(positions))
    for name, matrices in blocks.items():
        adata.layers[name] = scipy.sparse.vstack(matrices, format="csr")[order].astype(np.float32)
    adata.obs["velocity_contamination"] = contamination


def add_sample_metadata(adata, samplesheet_path):
    """Copies the metadata columns of the samplesheet to obs of every cell of their sample.
    A column of numbers becomes numeric, any other categorical with the categories in samplesheet order."""
    samplesheet = pd.read_csv(samplesheet_path, dtype=str)
    metadata = samplesheet.drop(columns=["fastq_1", "fastq_2"]).drop_duplicates("sample").set_index("sample")
    cell_samples = adata.obs["sample"].astype(str)
    for column in metadata.columns:
        if column in adata.obs:
            raise ValueError(f"samplesheet column {column} is already in obs; rename it")
        values = metadata[column].reindex(cell_samples).to_numpy()
        numbers = pd.to_numeric(values, errors="coerce")
        if np.isnan(numbers).sum() == pd.isna(values).sum():
            adata.obs[column] = numbers
        else:
            adata.obs[column] = pd.Categorical(values, categories=metadata[column].dropna().unique())


def pipeline_info(outdir, input_path):
    """Returns the settings and software versions of the run for uns['pipeline']."""
    scdownstream_info = outdir / "results" / "scdownstream" / "pipeline_info"
    commit = subprocess.run(
        ["git", "describe", "--always", "--dirty"], cwd=PREPROCESSING_DIR, capture_output=True, text=True, check=True
    )
    return {
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "pipeline_commit": commit.stdout.strip(),
        "scdownstream_h5ad": str(input_path),
        "config": (outdir / "config.yaml").read_text(),
        "scdownstream_samplesheet": (outdir / "scdownstream" / "samplesheet.csv").read_text(),
        "scdownstream_params": max(scdownstream_info.glob("params_*.json"), key=os.path.getmtime).read_text(),
        "software_versions": max(scdownstream_info.glob("*software*versions.yml"), key=os.path.getmtime).read_text(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir", type=Path, help="Output folder of the run")
    args = parser.parse_args()

    config = yaml.safe_load((args.outdir / "config.yaml").read_text())
    results = args.outdir / "results"
    input_path = finalized_h5ad(args.outdir)
    adata = anndata.read_h5ad(input_path)
    key = embedding_key(adata)

    adata.X = adata.layers.pop("log1p")
    adata.obs = adata.obs.drop(columns=DROP_OBS)
    rename_keys(adata, key)
    add_neighbors(adata, results / "scdownstream", key)
    add_doublet_scores(adata, results / "scdownstream")
    add_count_metrics(adata)
    if config["marker_genes"]:
        add_cell_types(adata, results / "cellassign" / "cell_types.tsv")
    if config["run_velocyto"]:
        add_velocyto(adata, results / "velocyto_decontx")
    add_sample_metadata(adata, args.outdir / "samplesheet.csv")
    adata.uns["pipeline"] = pipeline_info(args.outdir, input_path)

    output_dir = results / "final"
    output_dir.mkdir(exist_ok=True)
    adata.write_h5ad(output_dir / input_path.name, compression="gzip")
    print(adata)


if __name__ == "__main__":
    main()
