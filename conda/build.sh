#!/bin/bash
set -euo pipefail

# Build the physlr binary in release mode
cargo build --release

# Install the binary
mkdir -p "${PREFIX}/bin"
cp target/release/physlr "${PREFIX}/bin/physlr"

# Install helper scripts
cp scripts/plotpaf.py "${PREFIX}/bin/physlr-plotpaf"
chmod +x "${PREFIX}/bin/physlr-plotpaf"
