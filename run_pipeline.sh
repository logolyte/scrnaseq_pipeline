#!/bin/bash
# Run on the login node: bash run_job.sh [path/to/pipeline.conf]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

PIPELINE_CONF="$(realpath "${1:-$ROOT/pipeline.conf}")"
export PIPELINE_CONF
set -a; source "$PIPELINE_CONF"; set +a
mkdir -p slurm_logs

# Check before submitting anything (MARKER_GENES is relative to this repo,
# SAMPLESHEET to the directory of pipeline.conf)
SAMPLESHEET_PATH="$(cd "$(dirname "$PIPELINE_CONF")" && realpath -m "${SAMPLESHEET:-}")"
if [[ -z "${SAMPLESHEET:-}" || ! -f "$SAMPLESHEET_PATH" ]]; then
    echo "SAMPLESHEET ${SAMPLESHEET:-(not set)} not found (looked for $SAMPLESHEET_PATH)" >&2
    exit 1
fi
if [[ "${RUN_VELOCYTO:-false}" == true && ! -f "${CONTAINER_DIR:-}/decontx-celda-1.26.0.sif" ]]; then
    echo "RUN_VELOCYTO is true but ${CONTAINER_DIR:-}/decontx-celda-1.26.0.sif is missing; run setup.sh" >&2
    exit 1
fi
if [[ -n "${MARKER_GENES:-}" && ! -f "$MARKER_GENES" ]]; then
    echo "MARKER_GENES file $MARKER_GENES not found" >&2
    exit 1
fi

# The samplesheet, cellassign and finalize jobs run Python directly and need the scrnaseq
# packages; sbatch passes this shell's environment on to them
set +u  # conda's activation scripts use unset variables
eval "$(conda shell.bash hook)"
conda activate scrnaseq
set -u

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

SHEET_JOB=$(submit samplesheet . scripts/make_scdownstream_samplesheet.py \
            --dependency=afterok:"$STAR_JOB" --kill-on-invalid-dep=yes)
echo "samplesheet:  $SHEET_JOB"

DOWN_JOB=$(submit scdownstream downstream_processing_pipeline \
           downstream_processing_pipeline/downstream_processing_nfcore.sh \
           --dependency=afterok:"$SHEET_JOB" --kill-on-invalid-dep=yes)
echo "scdownstream: $DOWN_JOB"

FINAL_DEPENDENCY=$DOWN_JOB
if [[ -n "${MARKER_GENES:-}" ]]; then
    CELLASSIGN_JOB=$(submit cellassign . scripts/cellassign.py \
                     --dependency=afterok:"$DOWN_JOB" --kill-on-invalid-dep=yes)
    echo "cellassign:   $CELLASSIGN_JOB"
    FINAL_DEPENDENCY=$CELLASSIGN_JOB
fi

# The velocyto run realigns all samples, independently of the other jobs; its
# ambient correction needs the CellBender cells from scdownstream
if [[ "${RUN_VELOCYTO:-false}" == true ]]; then
    VELO_JOB=$(submit velocyto velocyto_pipeline velocyto_pipeline/velocyto_nfcore.sh)
    echo "velocyto:     $VELO_JOB"
    VELO_DECONTX_JOB=$(submit velocity_decontx . scripts/velocity_decontx.py \
                       --dependency=afterok:"$VELO_JOB:$DOWN_JOB" --kill-on-invalid-dep=yes)
    echo "velo decontX: $VELO_DECONTX_JOB"
    FINAL_DEPENDENCY="$FINAL_DEPENDENCY:$VELO_DECONTX_JOB"
fi

FINAL_JOB=$(submit finalize . scripts/finalize_adata.py \
            --dependency=afterok:"$FINAL_DEPENDENCY" --kill-on-invalid-dep=yes)
echo "finalize:     $FINAL_JOB"
