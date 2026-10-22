#!/bin/bash -l
#SBATCH --job-name=nf_starsolo
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

# execute the pipeline
nextflow run nf-core/scrnaseq \
    -r 4.2.0 \
    -profile singularity,gpu \
    -c slurm.config \
    -work-dir "$WORK_DIR/starsolo" \
    -params-file nf-params.json \
    --outdir "$STARSOLO_OUTDIR" \
    --protocol "$PROTOCOL" \
    --fasta "$FASTA_PATH" \
    --gtf "$GTF_PATH" \
    --star_index "$STAR_INDEX_PATH" \
    -resume
