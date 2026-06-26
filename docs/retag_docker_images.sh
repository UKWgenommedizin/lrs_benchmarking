#!/bin/bash
# ============================================================
# retag_docker_images.sh
# Re-tag schimar/lrs-* Docker images from 'latest' to pinned
# version tags.
#
# Usage: bash docs/retag_docker_images.sh
#
# This script inspects each image's metadata for a version
# label, then creates a versioned tag. If no version label
# is found, it asks you to provide one.
#
# Run this on genmedbfx where the Docker images are stored.
# ============================================================

set -euo pipefail

IMAGES=(
    "schimar/lrs-graphaligner"
    "schimar/lrs-parahat"
    "schimar/lrs-quicked"
    "schimar/lrs-vacmap"
    "schimar/lrs-vg"
    "schimar/lrs-minimap2-ntlink"
)

echo "=== Re-tagging schimar/lrs-* images from 'latest' to pinned versions ==="
echo ""

for IMAGE in "${IMAGES[@]}"; do
    echo "--- Processing $IMAGE ---"

    # Pull latest
    docker pull "$IMAGE:latest" 2>&1 | tail -1

    # Try to extract version from Docker label
    VERSION=$(docker inspect "$IMAGE:latest" \
        --format '{{index .Config.Labels "org.opencontainers.image.version"}}' 2>/dev/null || true)

    # Try alternative label
    if [ -z "$VERSION" ] || [ "$VERSION" = "<no value>" ]; then
        VERSION=$(docker inspect "$IMAGE:latest" \
            --format '{{index .Config.Labels "version"}}' 2>/dev/null || true)
    fi

    # If still no version, prompt
    if [ -z "$VERSION" ] || [ "$VERSION" = "<no value>" ]; then
        echo "  No version label found in image metadata."
        read -p "  Enter version tag for $IMAGE: " VERSION
    fi

    # Tag and push
    echo "  Tagging: $IMAGE:latest -> $IMAGE:$VERSION"
    docker tag "$IMAGE:latest" "$IMAGE:$VERSION"
    docker push "$IMAGE:$VERSION" 2>&1 | tail -1

    echo "  Done: $IMAGE:$VERSION"
    echo ""
done

echo "=== All images tagged ==="
echo ""
echo "To update the workflow files, replace 'latest' with the"
echo "version tags shown above in all ont/pb.read_mapping.*.smk files."
echo ""
echo "Example sed command for all files:"
echo "  sed -i 's|schimar/lrs-graphaligner:latest|schimar/lrs-graphaligner:<version>|g' *.smk"