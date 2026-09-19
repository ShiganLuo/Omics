export PATH="/opt/conda/envs/fibertools/bin:$PATH"
export CC=x86_64-conda-linux-gnu-gcc
export CXX=x86_64-conda-linux-gnu-g++
export LIBTORCH=/opt/conda/envs/fibertools/lib
export LD_LIBRARY_PATH="${LIBTORCH}:${LD_LIBRARY_PATH:-}"
export LIBTORCH_CXX11_ABI=0
export CARGO_HOME=/opt/conda/.cargo
cargo install fibertools-rs --all-features --version ">=0.13.0"
cp ${CARGO_HOME}/bin/ft /opt/conda/envs/fibertools/bin/ft
ft --version
rm -rf ${CARGO_HOME}/registry ${CARGO_HOME}/git /tmp/*