#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

TARGET="$(rustc -vV | sed -n 's/^host: //p')"

echo "Building sidecars for: $TARGET"

# Require an Apple signing identity
: "${APPLE_SIGNING_IDENTITY:?Apple signing identity is missing}"

echo "Available signing identities:"
security find-identity -v -p codesigning

python -c "import sys; print('Python:', sys.executable)"
python -c "import torch, torch_geometric, pyg_lib, rdkit"

libpyg="$(python -c "import importlib.util; s = importlib.util.find_spec('libpyg'); print(s.origin if s else '')")"

if [[ -z "$libpyg" || ! -f "$libpyg" ]]; then
  echo "Could not locate libpyg."
  python -m pip show pyg-lib
  exit 1
fi

echo "Bundling: $libpyg"

# Common PyInstaller arguments
COMMON_ARGS=(
  --onefile
  --noconfirm
  --clean
  --codesign-identity "$APPLE_SIGNING_IDENTITY"
)

# Build metabolite sidecar
python -m PyInstaller ../src/metabolite.py \
  "${COMMON_ARGS[@]}" \
  --name "metabolite-${TARGET}"

# Build metabolite database sidecar
python -m PyInstaller ../src/metabolite_for_db.py \
  "${COMMON_ARGS[@]}" \
  --name "metabolite_for_db-${TARGET}"

# Build CCS prediction sidecar
python -m PyInstaller ../ccs_prediction/model.py \
  "${COMMON_ARGS[@]}" \
  --console \
  --collect-all torch_geometric \
  --collect-all pyg_lib \
  --collect-all rdkit \
  --copy-metadata torch-geometric \
  --collect-binaries torch \
  --hidden-import pyg_lib \
  --hidden-import libpyg \
  --add-binary "${libpyg}:." \
  --add-data "../ccs_prediction/trained_weights.pt:." \
  --name "prediction_model-${TARGET}"

# Verify generated sidecars
for binary in metabolite metabolite_for_db prediction_model; do
  path="dist/${binary}-${TARGET}"

  test -f "$path"

  echo "Verifying: $path"

  codesign --verify --strict --verbose=2 "$path"
  codesign -dv --verbose=4 "$path" 2>&1
done

echo "All macOS sidecars built and signature checks passed."
