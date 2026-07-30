#!/usr/bin/env bash
set -u

echo "Checking Docker installation..."

if ! command -v docker >/dev/null 2>&1; then
    echo "WARNING: docker command is not available inside WSL."
    echo "Open Docker Desktop and enable WSL integration for Ubuntu."
    echo "Docker image checks were skipped."
    exit 0
fi

docker --version

echo
echo "Checking whether Docker daemon is running..."
if docker info >/dev/null 2>&1; then
    echo "Docker daemon is running."
else
    echo "WARNING: Docker command exists, but the Docker daemon is not running."
    echo "Open Docker Desktop, wait until it is fully started, and run this script again."
    exit 0
fi

echo
echo "Docker/container references found in this project:"
grep -RniE "docker run|docker://|container|image|storage-node|google/deepvariant|ensemblorg|nfcore" \
    *.smk header.smk config*.yaml 2>/dev/null || true

echo
echo "Image check finished. This script only identifies Docker images; it does not run variant calling."
