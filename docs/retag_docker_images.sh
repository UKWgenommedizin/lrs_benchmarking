#!/bin/bash
# ============================================================
# retag_docker_images.sh
# Re-tag schimar/lrs-* Docker images from 'latest' to pinned
# version tags.
#
# Usage: bash docs/retag_docker_images.sh <version>
#
# Example: bash docs/retag_docker_images.sh 1.0.0
#
# This tags all 6 images with the same version string.
# Run this on genmedbfx where the Docker images are stored.
# ============================================================

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <version-tag>"
    echo "Example: $0 1.0.0"
    echo ""
    echo "This tags all 6 images:"
    echo "  schimar/lrs-graphaligner:latest     -> schimar/lrs-graphaligner:<version>"
    echo "  schimar/lrs-parahat:latest           -> schimar/lrs-parahat:<version>"
    echo "  schimar/lrs-quicked:latest           -> schimar/lrs-quicked:<version>"
    echo "  schimar/lrs-vacmap:latest            -> schimar/lrs-vacmap:<version>"
    echo "  schimar/lrs-vg:latest                -> schimar/lrs-vg:<version>"
    echo "  schimar/lrs-minimap2-ntlink:latest   -> schimar/lrs-minimap2-ntlink:<version>"
    exit 1
fi

VERSION="$1"

IMAGES=(
    "schimar/lrs-graphaligner"
    "schimar/lrs-parahat"
    "schimar/lrs-quicked"
    "schimar/lrs-vacmap"
    "schimar/lrs-vg"
    "schimar/lrs-minimap2-ntlink"
)

echo "=== Tagging all images with version: $VERSION ==="
echo ""

for IMAGE in "${IMAGES[@]}"; do
    echo "  docker tag $IMAGE:latest $IMAGE:$VERSION"
    docker tag "$IMAGE:latest" "$IMAGE:$VERSION"
    echo "  docker push $IMAGE:$VERSION"
    docker push "$IMAGE:$VERSION" 2>&1 | tail -1
    echo "  Done: $IMAGE:$VERSION"
    echo ""
done

echo "=== All 6 images tagged and pushed ==="
echo ""
echo "Now update the .smk files with the new tag:"
echo ""
echo "  sed -i 's|:latest|:$VERSION|g' ont.read_mapping.*.smk pb.read_mapping.*.smk"
echo "  git add -A && git commit -m \"Pin Docker images to v$VERSION\" && git push"