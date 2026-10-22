PIPELINE_CONF="$(realpath "${1:-$ROOT/pipeline.conf}")"
export PIPELINE_CONF
set -a; source "$PIPELINE_CONF"; set +a
mkdir -p slurm_logs

sbatch --chdir="$PWD/downstream_processing_pipeline" \
       --output="$PWD/slurm_logs/scdownstream_%j.out" \
       --error="$PWD/slurm_logs/scdownstream_%j.err" \
       "$PWD/downstream_processing_pipeline/downstream_processing_nfcore.sh"
