#!/bin/bash
# One-time setup after cloning: bash setup.sh [path/to/pipeline.conf]
#
#  1. Creates the conda environment "scrnaseq" (Nextflow, nf-core, STAR, Python
#     packages for the pipeline's own scripts) from envs/scrnaseq.lock.yml, or
#     from envs/scrnaseq.yml if the exact versions no longer resolve.
#  2. Builds each container in containers/<name>/<image>.def into
#     $CONTAINER_DIR/<image>.sif, unless that image already exists.
#     If a recipe has no lock.txt yet, the exact package list of the new image
#     is written to it; commit it so later builds are identical.
#
# Safe to rerun: it skips what already exists. To rebuild an image, delete the
# .sif first. Building needs `apptainer build --fakeroot`; if your cluster does
# not allow it, build the .def elsewhere and copy the .sif into $CONTAINER_DIR.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

PIPELINE_CONF="$(realpath "${1:-$ROOT/pipeline.conf}")"
set -a; source "$PIPELINE_CONF"; set +a
: "${CONTAINER_DIR:?set CONTAINER_DIR in $PIPELINE_CONF}"

# 1. Conda environment
if ! command -v conda >/dev/null; then
    echo "conda not found on PATH" >&2
    exit 1
fi
if conda env list | awk '{print $1}' | grep -qx scrnaseq; then
    echo "conda environment scrnaseq exists; skipping"
else
    echo "creating conda environment scrnaseq from envs/scrnaseq.lock.yml"
    if ! conda env create -n scrnaseq -f envs/scrnaseq.lock.yml; then
        echo "exact versions failed; creating it from envs/scrnaseq.yml" >&2
        conda env remove -n scrnaseq -y >/dev/null 2>&1 || true
        conda env create -n scrnaseq -f envs/scrnaseq.yml
    fi
fi

# 2. Containers
mkdir -p "$CONTAINER_DIR"
for recipe in containers/*/; do
    defs=("$recipe"*.def)
    if [[ ${#defs[@]} -ne 1 || ! -e ${defs[0]} ]]; then
        echo "$recipe: expected exactly one .def file" >&2
        exit 1
    fi
    def=${defs[0]}
    image="$CONTAINER_DIR/$(basename "$def" .def).sif"
    if [[ -e "$image" ]]; then
        echo "$image exists; skipping"
        continue
    fi

    # Build from a copy of the recipe, with an empty lock.txt if there is none
    # yet, so the .def can always copy the same files
    stage=$(mktemp -d)
    trap 'rm -rf "$stage"' EXIT
    cp -r "$recipe". "$stage"/
    touch "$stage/lock.txt"
    echo "building $image from $def"
    rm -f "$image.tmp"
    # Clusters often bind host directories by default (e.g. SINGULARITY_BIND=/scratch).
    # %post runs in a writable container where those mount points don't exist, so
    # the build fails; it needs no host directories, so drop the binds for the build.
    (cd "$stage" && env -u SINGULARITY_BIND -u APPTAINER_BIND \
        apptainer build --fakeroot "$image.tmp" "$(basename "$def")")
    mv "$image.tmp" "$image"
    rm -rf "$stage"
    trap - EXIT

    if [[ ! -s "$recipe/lock.txt" ]]; then
        apptainer run --app lock "$image" > "$recipe/lock.txt"
        echo "wrote $recipe/lock.txt; commit it so later builds use the same versions"
    fi
done

echo "setup done"
