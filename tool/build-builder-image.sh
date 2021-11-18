#!/usr/bin/env sh

set -e

docker build -f Dockerfile_build \
    -t traviswheelerlab/fastavector-build:latest \
    $@ \
    .
