"""Helpers shared by the pipeline's Python scripts (pipeline.conf and scdownstream paths).

The scripts are submitted with sbatch, which runs a copy from Slurm's spool
directory, so they find this module through the project root (the working
directory) as well as their own directory.
"""

import glob
import os
import subprocess
import sys

DEFAULT_CONF = os.environ.get("PIPELINE_CONF", "pipeline.conf")


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


def find_finalized_h5ad(downstream_outdir):
    """The finalized h5ad of the scdownstream run."""
    pattern = os.path.join(downstream_outdir, "09_finalized", "*.h5ad")
    candidates = [p for p in sorted(glob.glob(pattern)) if not p.endswith("_cellxgene.h5ad")]
    if len(candidates) != 1:
        found = ", ".join(candidates) or "none"
        sys.exit(f"expected one h5ad matching {pattern} (found: {found}); give --input")
    return candidates[0]
