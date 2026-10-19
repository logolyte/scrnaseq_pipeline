#!/bin/bash
# Run on the login node: bash run_job.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p slurm_logs

submit() {                       # submit <name> <chdir> <script>
    sbatch --parsable \
           --chdir="$ROOT/$2" \
           --output="$ROOT/slurm_logs/$1_%j.out" \
           --error="$ROOT/slurm_logs/$1_%j.err" \
           "${@:4}" \
           "$ROOT/$3"
}

STAR_JOB=$(submit starsolo starsolo_pipeline starsolo_pipeline/starsolo_nfcore.sh)
echo "starsolo:     $STAR_JOB"

SHEET_JOB=$(submit samplesheet . make_downstream_samplesheet.py \
            --dependency=afterok:"$STAR_JOB" --kill-on-invalid-dep=yes)
echo "samplesheet:  $SHEET_JOB"

DOWN_JOB=$(submit scdownstream downstream_processing_pipeline \
           downstream_processing_pipeline/downstream_processing_nfcore.sh \
           --dependency=afterok:"$SHEET_JOB" --kill-on-invalid-dep=yes)
echo "scdownstream: $DOWN_JOB"
