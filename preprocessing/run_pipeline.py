#!/usr/bin/env python3
"""Preprocesses one dataset: asks which samples to take, submits the steps to Slurm as a chain of jobs and returns.
Run it on the login node in the scrnaseq environment; rerun it with the same outdir after a failure."""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import yaml

PIPELINE_DIR = Path(__file__).resolve().parent
SCRNASEQ = "nf-core/scrnaseq -r 4.2.0"
# scdownstream has no releases yet
SCDOWNSTREAM = "nf-core/scdownstream -r 0048567589248aac5669ffd0151b8f6d15417804"
# the Nextflow jobs run only the driver; the tasks are separate jobs
SLURM_OPTIONS = {
    "scrnaseq": "--cpus-per-task=2 --mem=8G --time=72:00:00",
    "velocyto": "--cpus-per-task=2 --mem=8G --time=72:00:00",
    "prepare_scdownstream": "--cpus-per-task=1 --mem=8G --time=01:00:00",
    "scdownstream": "--cpus-per-task=2 --mem=8G --time=72:00:00",
    "cellassign": "--cpus-per-task=8 --mem=16G --time=24:00:00",
    "velocity_decontx": "--cpus-per-task=4 --mem=64G --time=24:00:00",
    "finalize": "--cpus-per-task=2 --mem=64G --time=04:00:00",
}


def write_params(launch_dir, step_params, run_params):
    """Writes launch_dir/params.yaml: the step's params of every run plus the values of this run."""
    params = yaml.safe_load(step_params.read_text()) | run_params
    launch_dir.mkdir(parents=True, exist_ok=True)
    (launch_dir / "params.yaml").write_text(yaml.safe_dump(params, sort_keys=False))


def nextflow_command(pipeline, configs):
    """Builds the command that runs a Nextflow pipeline from its launch folder."""
    config_options = " ".join(f"-c {config}" for config in configs)
    return f"nextflow run {pipeline} -profile singularity,gpu {config_options} -params-file params.yaml -resume"


def python_command(script, outdir):
    """Builds the command that runs one of our Python steps on the output folder."""
    return f"{sys.executable} {PIPELINE_DIR / script} {outdir}"


def submit(step, command, launch_dir, after, outdir):
    """Submits a step that starts once the jobs in `after` have succeeded and returns its job ID."""
    options = [
        "--parsable",
        "--partition=smp",
        f"--job-name={step}",
        f"--chdir={launch_dir}",
        f"--output={outdir / 'slurm_logs' / step}_%j.log",
    ]
    options += SLURM_OPTIONS[step].split()
    if after:
        options += ["--dependency=afterok:" + ":".join(after), "--kill-on-invalid-dep=yes"]
    result = subprocess.run(["sbatch", *options, "--wrap", command], capture_output=True, text=True, check=True)
    job_id = result.stdout.strip()
    print(f"{step}: {job_id}")
    return job_id


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samplesheet", type=Path, help="One row per run: sample, fastq_1, fastq_2, protocol, metadata")
    parser.add_argument("outdir", type=Path, help="Output folder of this run")
    args = parser.parse_args()
    outdir = args.outdir.resolve()

    # curate_samples.py writes config.yaml last, so an existing one means the samples are chosen
    if not (outdir / "config.yaml").exists():
        curate_script = PIPELINE_DIR / "1_curate" / "curate_samples.py"
        subprocess.run([sys.executable, curate_script, args.samplesheet, outdir], check=True)
    config = yaml.safe_load((outdir / "config.yaml").read_text())

    # the jobs inherit this environment; Nextflow and the .config files read these variables
    cache_dir = Path(config["cache_dir"])
    os.environ.update({
        "NXF_HOME": str(cache_dir / "nextflow"),
        "NXF_SINGULARITY_CACHEDIR": str(cache_dir / "singularity_cache"),
        "APPTAINER_CACHEDIR": str(cache_dir / "apptainer_cache"),
        "CONTAINER_DIR": str(cache_dir / "containers"),
        "NXF_OPTS": "-Xms1g -Xmx4g",
        "SBATCH_ACCOUNT": config["slurm_account"],
        "CELLBENDER_PRIORS": str(outdir / "scdownstream" / "cellbender_priors.tsv"),
    })
    for variable in ["NXF_HOME", "NXF_SINGULARITY_CACHEDIR", "APPTAINER_CACHEDIR"]:
        Path(os.environ[variable]).mkdir(parents=True, exist_ok=True)
    (outdir / "slurm_logs").mkdir(exist_ok=True)

    reference = config["reference"]
    scrnaseq_params = {
        "input": str(outdir / "samplesheet.csv"),
        "outdir": str(outdir / "results" / "starsolo"),
        "protocol": config["protocol"],
        "star_feature": config["modes"][config["mode"]]["star_feature"],
        "fasta": reference["fasta"],
        "gtf": reference["gtf"],
        "star_index": reference["star_index"],
    }
    scrnaseq_dir = PIPELINE_DIR / "2_scrnaseq"
    write_params(outdir / "scrnaseq", scrnaseq_dir / "params.yaml", scrnaseq_params)
    command = nextflow_command(SCRNASEQ, [scrnaseq_dir / "scrnaseq.config"])
    scrnaseq_job = submit("scrnaseq", command, outdir / "scrnaseq", [], outdir)

    command = python_command("3_scdownstream/prepare_scdownstream.py", outdir)
    prepare_job = submit("prepare_scdownstream", command, outdir, [scrnaseq_job], outdir)

    scdownstream_dir = PIPELINE_DIR / "3_scdownstream"
    scdownstream_params = {
        "input": str(outdir / "scdownstream" / "samplesheet.csv"),
        "outdir": str(outdir / "results" / "scdownstream"),
        "species": config["species"],
    }
    write_params(outdir / "scdownstream", scdownstream_dir / "params.yaml", scdownstream_params)
    command = nextflow_command(SCDOWNSTREAM, [scdownstream_dir / "scdownstream.config"])
    scdownstream_job = submit("scdownstream", command, outdir / "scdownstream", [prepare_job], outdir)

    finalize_after = [scdownstream_job]
    if config["marker_genes"]:
        command = python_command("4_cellassign/cellassign.py", outdir)
        finalize_after.append(submit("cellassign", command, outdir, [scdownstream_job], outdir))
    if config["run_velocyto"]:
        # the same scrnaseq run, counting spliced and unspliced reads
        velocyto_params = scrnaseq_params | {
            "outdir": str(outdir / "results" / "velocyto"),
            "star_feature": "Gene Velocyto",
        }
        write_params(outdir / "velocyto", scrnaseq_dir / "params.yaml", velocyto_params)
        configs = [scrnaseq_dir / "scrnaseq.config", scrnaseq_dir / "velocyto.config"]
        command = nextflow_command(SCRNASEQ, configs)
        velocyto_job = submit("velocyto", command, outdir / "velocyto", [], outdir)
        # decontX of the velocyto counts needs the cells CellBender called in scdownstream
        command = python_command("4_velocity_decontx/velocity_decontx.py", outdir)
        finalize_after.append(submit("velocity_decontx", command, outdir, [velocyto_job, scdownstream_job], outdir))

    command = python_command("5_finalize/finalize_adata.py", outdir)
    submit("finalize", command, outdir, finalize_after, outdir)


if __name__ == "__main__":
    main()
