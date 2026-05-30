# Reproducible validation environment for apsis.
#
# Provides the Rust toolchain plus a Python environment with the figure and
# REBOUNDx-parity dependencies (rebound/reboundx, which do NOT build on
# Windows/MSVC — C99 VLAs in gr_full.c), so the validation harnesses and the
# figure pipeline reproduce on any host with Docker.
#
# SCOPE: this is the reproducible *validation / figures* environment. It is
# NOT the platform for the §3.5 cross-platform bit-exactness claim — that
# requires native Windows vs native Linux; an all-Linux container would mask
# it. Keep the two separate.
#
# Build:  docker build -f docker/Dockerfile -t apsis-validation .
# Use:    docker run --rm -it -v "$PWD:/apsis" apsis-validation
#           # then, inside the container:
#           make figures
#           cd validation/reboundx-parity/gr-mercury && python run.py
#
# See docker/README.md.

FROM ubuntu:24.04

# >= the workspace's rust-version (Cargo.toml: 1.89, edition 2024).
ARG RUST_VERSION=1.89.0

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    VIRTUAL_ENV=/opt/venv

# Build deps: build-essential (gcc — reboundx's C build), python venv/headers,
# curl (rustup), git, pkg-config; fontconfig for deterministic matplotlib text.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential curl git pkg-config \
        python3 python3-venv python3-dev \
        fontconfig \
    && rm -rf /var/lib/apt/lists/*

# Pinned Rust toolchain (minimal profile — no docs/clippy bloat; add via
# `rustup component add` if a target needs them).
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
        | sh -s -- -y --default-toolchain "${RUST_VERSION}" --profile minimal
ENV PATH=/root/.cargo/bin:${VIRTUAL_ENV}/bin:$PATH

# Bind-mounted repo is host-UID-owned; as root, git refuses it ("dubious
# ownership") and the build scripts' rev-parse silently yields an empty git_sha.
RUN git config --system --add safe.directory /apsis

# Python venv with the validation/figure deps. reboundx is compiled fresh
# against its co-installed rebound (--no-cache-dir): reusing a cached wheel
# bakes a stale librebound RPATH that fails to load. rebound 4.6.0 is pinned
# by reboundx 4.6.2 (no rebound-5 reboundx release yet).
RUN python3 -m venv "${VIRTUAL_ENV}"
COPY docker/requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -U pip \
    && pip install --no-cache-dir -r /tmp/requirements.txt \
    && python -c "import numpy, scipy, matplotlib, pandas, rebound, reboundx; \
print('deps OK:', rebound.__version__, reboundx.__version__)"

# The repo is bind-mounted at run time (-v "$PWD:/apsis"); nothing is copied
# into the image, so it stays a thin toolchain layer and the source is always
# the host's current checkout.
WORKDIR /apsis
CMD ["bash"]
