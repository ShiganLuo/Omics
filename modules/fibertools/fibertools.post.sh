set +euo pipefail
export PATH="/opt/conda/envs/fibertools/bin:$PATH"
export CC=x86_64-conda-linux-gnu-gcc
export CXX=x86_64-conda-linux-gnu-g++

# Source conda activation scripts
for f in /opt/conda/envs/fibertools/etc/conda/activate.d/*.sh; do
    source "$f" 2>/dev/null || true
done

# Install CPU-only PyTorch via pip (includes C++ headers)
pip install --no-cache-dir torch==2.9.0 --index-url https://download.pytorch.org/whl/cpu

# Set LIBTORCH from installed torch
export LIBTORCH=$(python3 -c "import os, torch; print(os.path.dirname(torch.__file__))")
export LD_LIBRARY_PATH="${LIBTORCH}/lib:${LD_LIBRARY_PATH:-}"
export LIBTORCH_CXX11_ABI=1
export CARGO_HOME=/opt/conda/.cargo
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=x86_64-conda-linux-gnu-gcc
export OPENSSL_DIR=/opt/conda/envs/fibertools
export OPENSSL_NO_VENDOR=1
set -e
cargo install fibertools-rs --all-features --version ">=0.13.0"
cp ${CARGO_HOME}/bin/ft /opt/conda/envs/fibertools/bin/ft
ft --version

# Persist LD_LIBRARY_PATH for Apptainer runtime (%environment auto-sources this)
mkdir -p /.singularity.d/env
echo "export LD_LIBRARY_PATH=\"${LIBTORCH}/lib:\${LD_LIBRARY_PATH:-}\"" >> /.singularity.d/env/99-fibertools.sh
rm -rf ${CARGO_HOME}/registry ${CARGO_HOME}/git