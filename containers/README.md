# Shared containers

Reserved for project-wide container definitions and lock files. Tool-specific assembly containers currently live under [`../assemblers/containers/`](../assemblers/containers/). Container images must be pinned according to [`../CONSTITUTION.md`](../CONSTITUTION.md).

## Aligner images (ram_time workflows)

All images used by [`../alignment_analysis/run_metrics/ram_time/`](../alignment_analysis/run_metrics/ram_time/README.md) are published under the `nicolasardila1` Docker Hub account. Each Dockerfile starts from the production image, so the aligner and samtools binaries are the ones used in the production mapping runs, and adds `/usr/bin/time`, `python3` and GNU `du` if they are missing.

| Directory | Image | Built from |
|---|---|---|
| [`vg/`](vg/Dockerfile) | `nicolasardila1/lrs-vg:v1.73.0` | `schimar/lrs-vg:v1.73.0` |
| [`vacmap/`](vacmap/Dockerfile) | `nicolasardila1/lrs-vacmap:v1.2.0` | `schimar/lrs-vacmap:v1.2.0` |

Build, check and push (after `docker login -u nicolasardila1`):

```bash
bash containers/build_and_push_aligner_images.sh            # all images
bash containers/build_and_push_aligner_images.sh vacmap     # one image
bash containers/build_and_push_aligner_images.sh --no-push  # build + check only
```

The base images are only pulled at build time; the workflows pull from `nicolasardila1` only.
