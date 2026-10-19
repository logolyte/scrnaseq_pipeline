#!/bin/bash -l
#SBATCH --job-name=nf_downstream
#SBATCH --output=slurm_logs/nf_downstream_%j.out
#SBATCH --error=slurm_logs/nf_downstream_%j.err
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=72:00:00
#SBATCH --account=ag-beyer

# initialize conda and activate the environment
conda activate scrnaseq

# fail on any error, and on a failure anywhere in a pipe
set -euo pipefail

# route cache directories to scratch storage
export NXF_SINGULARITY_CACHEDIR="/scratch/$USER/cache/singularity_cache"
export APPTAINER_CACHEDIR="/scratch/$USER/cache/apptainer_cache"
export NXF_HOME="/scratch/$USER/cache/nextflow"
export NXF_OPTS='-Xms1g -Xmx4g'
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$APPTAINER_CACHEDIR" "$NXF_HOME"

# execute the pipeline
nextflow run nf-core/scdownstream -r dev -profile singularity,gpu -work-dir /scratch/$USER/work/downstream -resume -params-file nf-params.json -c slurm.config
