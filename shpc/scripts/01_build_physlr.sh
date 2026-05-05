#!/bin/bash
#SBATCH --job-name=physlr-build
#SBATCH --output=./logs/build.%j.out
#SBATCH --error=./logs/build.%j.err
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

# Build Physlr (Rust) on SHPC.
# Assumes Rust toolchain is available (install via rustup if needed).
#
# Usage: sbatch scripts/01_build_physlr.sh

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
SRCDIR="${BASEDIR}/physlr-next"
mkdir -p "${BASEDIR}/logs"

echo "=== Building Physlr ==="

# Check for Rust
if ! command -v cargo &>/dev/null; then
    echo "Rust not found. Installing via rustup..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source "$HOME/.cargo/env"
fi

# Clone or update source
if [ ! -d "${SRCDIR}" ]; then
    echo "ERROR: Physlr source not found at ${SRCDIR}"
    echo "Copy the physlr-next directory to ${SRCDIR} first."
    exit 1
fi

cd "${SRCDIR}"
cargo build --release

echo "Binary: ${SRCDIR}/target/release/physlr"
"${SRCDIR}/target/release/physlr" --version 2>/dev/null || "${SRCDIR}/target/release/physlr" --help | head -3
echo "=== Build complete ==="
