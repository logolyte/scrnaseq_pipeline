#!/bin/bash -l
#SBATCH --job-name=nf_scdownstream
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=72:00:00

# initialize conda and activate the environment
conda activate scrnaseq

# fail on any error, and on a failure anywhere in a pipe
set -euo pipefail

# load the global config (path exported by run_job.sh)
set -a; source "${PIPELINE_CONF:?submit via run_job.sh or set PIPELINE_CONF}"; set +a
export NXF_OPTS='-Xms1g -Xmx4g'
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$APPTAINER_CACHEDIR" "$NXF_HOME"

# the images slurm.config points to; built by setup.sh
for image in cellbender-0.4.0-gpu.sif scdblfinder-1.24.10.sif; do
    if [[ ! -e "$CONTAINER_DIR/$image" ]]; then
        echo "$CONTAINER_DIR/$image not found; run setup.sh first" >&2
        exit 1
    fi
done

# execute the pipeline
nextflow run nf-core/scdownstream \
    -r 0048567589248aac5669ffd0151b8f6d15417804 \
    -profile singularity,gpu \
    -c slurm.config \
    -work-dir "$WORK_DIR/scdownstream" \
    -params-file nf-params.json \
    --outdir "$DOWNSTREAM_OUTDIR" \
    -resume
