#!/bin/bash
set -euo pipefail

# Activate conda env
export PATH="/opt/conda/envs/fibertools/bin:$PATH"

# Set compiler env vars (conda gcc_linux-64)
export CC=x86_64-conda-linux-gnu-gcc
export CXX=x86_64-conda-linux-gnu-g++

# libtorch from libtorch conda package
export LIBTORCH=/opt/conda/envs/fibertools/lib
export LD_LIBRARY_PATH="${LIBTORCH}:${LD_LIBRARY_PATH:-}"
export LIBTORCH_CXX11_ABI=0

# Compile fibertools-rs with pytorch backend (supports Revio + all binding kits)
export CARGO_HOME=/opt/conda/.cargo
cargo install fibertools-rs --all-features --version ">=0.13.0"

# Copy ft binary to conda env PATH
cp ${CARGO_HOME}/bin/ft /opt/conda/envs/fibertools/bin/ft

# Verify
ft --version

# Clean cargo cache
rm -rf ${CARGO_HOME}/registry ${CARGO_HOME}/git /tmp/*