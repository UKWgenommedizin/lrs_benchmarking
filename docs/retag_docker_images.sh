#!/bin/bash

set -euo pipefail

IMAGES=(
    "schimar/lrs-graphaligner"
    "schimar/lrs-parahat"
    "schimar/lrs-quicked"
    "schimar/lrs-vacmap"
    "schimar/lrs-vg"
    "schimar/lrs-minimap2-ntlink"
)

for IMAGE in "${IMAGES[@]}"; do
    echo "=== Processing $IMAGE ==="

    VERSION=$(docker image inspect "$IMAGE:latest" \
        --format '{{ index .Config.Labels "org.opencontainers.image.version" }}')

    if [[ -z "$VERSION" || "$VERSION" == "<no value>" ]]; then
        echo "ERROR: No version label found for $IMAGE"
        exit 1
    fi

    echo "Version: $VERSION"

    docker tag "$IMAGE:latest" "$IMAGE:$VERSION"
    docker push "$IMAGE:$VERSION"

    echo "Done: $IMAGE:$VERSION"
    echo
done

echo "Finished."
