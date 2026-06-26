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

docker pull schimar/lrs-graphaligner:v1.0.20 \
  schimar/lrs-fmalign2:v2.1.0 \
  schimar/lrs-parahat:v1.0.0-cuda \
  schimar/lrs-quicked:v0.1.0 \
  schimar/lrs-vacmap:v1.2.0 \
  schimar/lrs-vg:v1.73.0 \
  schimar/lrs-minimap2-ntlink:v2.28-1.3.9

docker rmi \
  schimar/lrs-graphaligner:latest schimar/lrs-graphaligner:24.1.2-0 \
  schimar/lrs-fmalign2:latest \
  schimar/lrs-parahat:latest \
  schimar/lrs-quicked:latest \
  schimar/lrs-vacmap:latest \
  schimar/lrs-vg:latest \
  schimar/lrs-minimap2-ntlink:latest