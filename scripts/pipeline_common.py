"""Helpers shared by the pipeline's Python scripts (pipeline.conf and scdownstream paths).

The scripts are submitted with sbatch, which runs a copy from Slurm's spool
directory, so they find this module in scripts/ of the project root (the working
directory) as well as in their own directory.
"""

import glob
import math
import os
import subprocess
import sys

import pandas as pd

DEFAULT_CONF = os.environ.get("PIPELINE_CONF", "pipeline.conf")

# Columns nf-core/scrnaseq 4.2.0 reads from the samplesheet (assets/schema_input.json);
# the others are per-sample QC thresholds (SAMPLE_QC_COLUMNS) or sample metadata
SCRNASEQ_COLUMNS = {"sample", "fastq_1", "fastq_2", "fastq_barcode", "expected_cells",
                    "seq_center", "sample_type", "feature_type"}

# scdownstream's per-sample QC thresholds (its samplesheet columns) and the pipeline's
# defaults. Lenient: only empty droplets and debris are removed, stricter filtering is
# left to each analysis. A samplesheet column of the same name overrides a cell filter
# for its sample; 0 turns a *_nmads filter off.
QC_DEFAULTS = {
    "min_genes": 200,         # removes near-empty droplets in samples with little ambient RNA
    "min_counts_cell": 200,   # >= 200 UMI left after decontX, removes ambient-like debris
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
# Gene filters stay at the default: per-sample values would remove genes from some samples only
GENE_FILTERS = {"min_cells", "min_counts_gene"}
SAMPLE_QC_COLUMNS = [c for c in QC_DEFAULTS if c not in GENE_FILTERS]
# Typed integer in scdownstream's input schema; the *_nmads columns are numbers
QC_INTEGER_COLUMNS = {"min_genes", "min_counts_cell", "min_cells", "min_counts_gene",
                      "max_mito_percentage", "min_ribo_percentage", "max_hb_percentage"}


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


def conf_relative(path, conf_path):
    """A path from pipeline.conf; relative paths are relative to the directory of that file."""
    if not path or os.path.isabs(path):
        return path
    return os.path.join(os.path.dirname(os.path.abspath(conf_path)), path)


def read_sample_columns(path):
    """The samplesheet columns nf-core/scrnaseq does not use, one row per sample.

    Values are strings, empty cells are missing. Lanes of one sample are separate
    rows; their non-empty values must agree.
    """
    sheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    if "sample" not in sheet:
        sys.exit(f"{path} has no 'sample' column")
    columns = [c for c in sheet.columns if c not in SCRNASEQ_COLUMNS]
    table = sheet[["sample", *columns]].replace("", pd.NA)
    n_values = table.groupby("sample")[columns].nunique()
    conflicts = [f"{s}: {c}" for c in columns for s in n_values.index[n_values[c] > 1]]
    if conflicts:
        sys.exit(f"{path}: rows of one sample have different values for {', '.join(conflicts)}")
    return table.groupby("sample", sort=False).first()  # first non-empty value per column


def warn_misnamed_qc(columns, path):
    """Warn about columns that look like QC filters but are not applied as one."""
    for column in columns:
        if column in SAMPLE_QC_COLUMNS:
            continue
        if column in GENE_FILTERS:
            reason = "a gene filter, which cannot be set per sample"
        elif column.startswith(("min_", "max_")) or column.endswith("_nmads"):
            reason = "not one of scdownstream's QC columns"
        else:
            continue
        print(f"WARNING: {path}: column '{column}' looks like a QC filter but is {reason}. "
              f"It is not applied as a filter and is kept as sample metadata. Per-sample QC "
              f"columns: {', '.join(SAMPLE_QC_COLUMNS)}", file=sys.stderr)


def sample_qc(path, samples):
    """QC thresholds per sample: QC_DEFAULTS, overridden by the samplesheet's QC columns."""
    table = read_sample_columns(path)
    warn_misnamed_qc(table.columns, path)
    qc, errors = {}, []
    for sample in samples:
        values = dict(QC_DEFAULTS)
        for column in SAMPLE_QC_COLUMNS:
            if column not in table or sample not in table.index or pd.isna(table.at[sample, column]):
                continue
            text = table.at[sample, column]
            try:
                value = float(text)
            except ValueError:
                value = math.nan
            if not math.isfinite(value) or value < 0:
                errors.append(f"{sample}: {column} = {text!r} is not a number >= 0")
            elif column.endswith("_percentage") and value > 100:
                errors.append(f"{sample}: {column} = {text!r} is above 100")
            elif column in QC_INTEGER_COLUMNS and not value.is_integer():
                errors.append(f"{sample}: {column} = {text!r} must be a whole number")
            else:
                values[column] = int(value) if column in QC_INTEGER_COLUMNS else value
        qc[sample] = values
    if errors:
        sys.exit(f"{path}: invalid QC values:\n  " + "\n  ".join(errors))
    return qc


def find_finalized_h5ad(downstream_outdir):
    """The finalized h5ad of the scdownstream run."""
    pattern = os.path.join(downstream_outdir, "09_finalized", "*.h5ad")
    candidates = [p for p in sorted(glob.glob(pattern)) if not p.endswith("_cellxgene.h5ad")]
    if len(candidates) != 1:
        found = ", ".join(candidates) or "none"
        sys.exit(f"expected one h5ad matching {pattern} (found: {found}); give --input")
    return candidates[0]
