#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

TARGET="$(rustc -vV | sed -n 's/^host: //p')"

echo "Building sidecars for: $TARGET"
python -c "import sys; print('Python:', sys.executable)"
python -c "import torch, torch_geometric, pyg_lib, rdkit"

libpyg="$(python -c "import importlib.util; s = importlib.util.find_spec('libpyg'); print(s.origin if s else '')")"

if [[ -z "$libpyg" || ! -f "$libpyg" ]]; then
  echo "Could not locate libpyg."
  python -m pip show pyg-lib
  exit 1
fi

echo "Bundling: $libpyg"

python -m PyInstaller ../src/metabolite.py \
  --onefile --noconfirm \
  --name "metabolite-${TARGET}"

python -m PyInstaller ../src/metabolite_for_db.py \
  --onefile --noconfirm \
  --name "metabolite_for_db-${TARGET}"

python -m PyInstaller ../ccs_prediction/model.py \
  --onefile --noconfirm --clean --console \
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

for binary in metabolite metabolite_for_db prediction_model; do
  test -f "dist/${binary}-${TARGET}"
done