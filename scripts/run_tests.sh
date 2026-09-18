#!/usr/bin/env bash
set -euo pipefail

# Execute full pymol-visualize skill test suite using uv and pytest in headless OSMesa mode
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${ROOT_DIR}"

export PYOPENGL_PLATFORM=osmesa

echo "Running pymol-visualize skill test suite..."
uv run --python 3.12 --with pytest --with pymol-open-source-whl pytest -v "$@"
echo "All tests passed successfully."
