#!/bin/bash
# Run on the login node: bash run_job.sh [path/to/pipeline.conf]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

PIPELINE_CONF="$(realpath "${1:-$ROOT/pipeline.conf}")"
export PIPELINE_CONF
set -a; source "$PIPELINE_CONF"; set +a
mkdir -p slurm_logs

submit() {   # submit <name> <chdir> <script> [sbatch options...]
    local name=$1 chdir=$2 script=$3
    shift 3 # Do not pass name, chdir and script do the submitted script
    sbatch --parsable \
           --chdir="$ROOT/$chdir" \
           --output="$ROOT/slurm_logs/${name}_%j.out" \
           --error="$ROOT/slurm_logs/${name}_%j.err" \
           "$@" "$ROOT/$script"
}

STAR_JOB=$(submit starsolo starsolo_pipeline starsolo_pipeline/starsolo_nfcore.sh)
echo "starsolo:     $STAR_JOB"

SHEET_JOB=$(submit samplesheet . make_scdownstream_samplesheet.py \
            --dependency=afterok:"$STAR_JOB" --kill-on-invalid-dep=yes)
echo "samplesheet:  $SHEET_JOB"

DOWN_JOB=$(submit scdownstream downstream_processing_pipeline \
           downstream_processing_pipeline/downstream_processing_nfcore.sh \
           --dependency=afterok:"$SHEET_JOB" --kill-on-invalid-dep=yes)
echo "scdownstream: $DOWN_JOB"
