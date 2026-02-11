#!/bin/bash
set -e

# Build Singularity images in stages.
# Each stage only rebuilds if its .sif doesn't already exist,
# or if you pass the stage number to force a rebuild.
#
# Usage:
#   ./build.sh          # build all missing stages
#   ./build.sh 2        # force rebuild from stage 2 onward
#   ./build.sh 3        # force rebuild stage 3 only

FORCE_FROM=${1:-99}  # default: don't force anything, only build missing

build_stage() {
    local num=$1
    local def=$2
    local sif=$3

    if [ "$FORCE_FROM" -le "$num" ] || [ ! -f "$sif" ]; then
        echo "=== Building stage $num: $def -> $sif ==="
        singularity build "$sif" "$def"
    else
        echo "=== Stage $num ($sif) already exists, skipping ==="
    fi
}

build_stage 1 stage1-base.def        stage1-base.sif
build_stage 2 stage2-r-install.def   stage2-r-install.sif
build_stage 3 stage3-r-packages.def  stage3-r-packages.sif
build_stage 4 stage4-r-propr-cuda.def     stage4-r-propr-cuda.sif
build_stage 4 stage4-r-propr-master.def   stage4-r-propr-master.sif

echo ""
echo "Done! Final image finished!"
