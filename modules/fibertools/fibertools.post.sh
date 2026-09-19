set +euo pipefail
export PATH="/opt/conda/envs/fibertools/bin:$PATH"
export CC=x86_64-conda-linux-gnu-gcc
export CXX=x86_64-conda-linux-gnu-g++

# Source conda activation scripts
for f in /opt/conda/envs/fibertools/etc/conda/activate.d/*.sh; do
    source "$f" 2>/dev/null || true
done

# Find torch include dir dynamically
TORCH_DIR=$(python3 -c "import torch; print(torch.utils.cmake_prefix_path)")
export LIBTORCH=$(python3 -c "import os, torch; print(os.path.dirname(torch.__file__))")
export LD_LIBRARY_PATH="${LIBTORCH}/lib:${LD_LIBRARY_PATH:-}"
export LIBTORCH_CXX11_ABI=0
export CARGO_HOME=/opt/conda/.cargo
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=x86_64-conda-linux-gnu-gcc
export OPENSSL_DIR=/opt/conda/envs/fibertools
export OPENSSL_NO_VENDOR=1
set -e
cargo install fibertools-rs --all-features --version ">=0.13.0"
cp ${CARGO_HOME}/bin/ft /opt/conda/envs/fibertools/bin/ft
ft --version
rm -rf ${CARGO_HOME}/registry ${CARGO_HOME}/git /tmp/*