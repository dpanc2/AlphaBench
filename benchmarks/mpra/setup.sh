#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
source_dir=data/alphabench
if [[ ! -d "$source_dir/.git" ]]; then
    mkdir -p data
    git clone https://github.com/tsoiegor/alphabench.git "$source_dir"
fi
git -C "$source_dir" checkout --detach 58ed6e7a8d52cf88fdfea4125f5db21f723c5623
patch="$PWD/benchmarks/mpra/alphabench-splice-reference.patch"
if git -C "$source_dir" apply --check "$patch"; then
    git -C "$source_dir" apply "$patch"
else
    git -C "$source_dir" apply --reverse --check "$patch"
fi
python -m pip install -e "$source_dir"
