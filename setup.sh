#!/bin/bash
# One-time setup after cloning: bash setup.sh <cache_dir of preprocessing/config.yaml>
# Creates the conda environment "scrnaseq" from envs/scrnaseq.yml and builds every
# containers/<name>/<image>.def into <cache_dir>/containers/<image>.sif.
# Safe to rerun: it skips what already exists. To rebuild an image, delete its .sif first.
set -euo pipefail

# absolute path, because the builds below run inside each recipe folder
CONTAINER_DIR="$(realpath "$1")/containers"
cd "$(dirname "$0")"

if conda env list | awk '{print $1}' | grep -qx scrnaseq; then
    echo "conda environment scrnaseq exists; skipping"
else
    conda env create -n scrnaseq -f envs/scrnaseq.yml
fi

mkdir -p "$CONTAINER_DIR"
for def in containers/*/*.def; do
    recipe=$(dirname "$def")
    image="$CONTAINER_DIR/$(basename "$def" .def).sif"
    if [[ -e "$image" ]]; then
        echo "$image exists; skipping"
        continue
    fi
    # the .def copies lock.txt; an empty one means "no lock yet"
    touch "$recipe/lock.txt"
    echo "building $image from $def"
    # host directories bound by default (e.g. SINGULARITY_BIND) do not exist inside the build
    (cd "$recipe" && env -u SINGULARITY_BIND -u APPTAINER_BIND \
        apptainer build --fakeroot "$image.tmp" "$(basename "$def")")
    mv "$image.tmp" "$image"
    if [[ ! -s "$recipe/lock.txt" ]]; then
        apptainer run --app lock "$image" > "$recipe/lock.txt"
        echo "wrote $recipe/lock.txt; commit it so later builds use the same versions"
    fi
done

echo "setup done"
