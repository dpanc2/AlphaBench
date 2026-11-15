#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
source_dir=data/eqtl/source
if [[ ! -d "$source_dir/.git" ]]; then
    mkdir -p data/eqtl
    git clone https://github.com/calengil/AlphaGenome-eQTL.git "$source_dir"
fi
git -C "$source_dir" checkout --detach 79c4cdf144a2ea21ec249b6e0d61902ce73e6ebb
