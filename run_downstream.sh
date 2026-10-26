#!/bin/bash
# Run on the login node: bash run_downstream.sh [path/to/pipeline.conf]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

PIPELINE_CONF="$(realpath "${1:-$ROOT/pipeline.conf}")"
export PIPELINE_CONF
set -a; source "$PIPELINE_CONF"; set +a
mkdir -p slurm_logs

sbatch --chdir="$ROOT/downstream_processing_pipeline" \
       --output="$ROOT/slurm_logs/scdownstream_%j.out" \
       --error="$ROOT/slurm_logs/scdownstream_%j.err" \
       "$ROOT/downstream_processing_pipeline/downstream_processing_nfcore.sh"
