#!/usr/bin/env bash
# Build the aligner images used by alignment_analysis/run_metrics/ram_time/,
# re-run the workflows' tool check on each, and push them to Docker Hub
# (nicolasardila1). Needs `docker login -u nicolasardila1` before pushing.
#
#   bash containers/build_and_push_aligner_images.sh              # all images
#   bash containers/build_and_push_aligner_images.sh vg           # one image
#   bash containers/build_and_push_aligner_images.sh --no-push    # build + check only

set -euo pipefail

declare -A IMAGES=(
    [vg]="nicolasardila1/lrs-vg:v1.73.0"
    [vacmap]="nicolasardila1/lrs-vacmap:v1.2.0"
)

CONTAINERS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PUSH=1
TOOLS=()
for arg in "$@"; do
    case "$arg" in
        --no-push) PUSH=0 ;;
        *) [[ -n "${IMAGES[$arg]:-}" ]] || { echo "ERROR: unknown image '$arg' (known: ${!IMAGES[*]})"; exit 1; }
           TOOLS+=("$arg") ;;
    esac
done
[[ ${#TOOLS[@]} -gt 0 ]] || TOOLS=(vg vacmap)

docker info >/dev/null || { echo "ERROR: Docker is not available"; exit 1; }

for tool in "${TOOLS[@]}"; do
    image="${IMAGES[$tool]}"
    echo "=== $tool -> $image"

    docker build --pull -t "$image" "$CONTAINERS/$tool"

    # Same check as the *_run rules of the ram_time workflows.
    docker run --rm --entrypoint sh "$image" \
        -c 'command -v /usr/bin/time && command -v python3 && du -sb /tmp >/dev/null' \
        && echo "$image OK"

    if [[ $PUSH -eq 1 ]]; then
        docker push "$image"
    fi
done
