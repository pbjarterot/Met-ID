#!/bin/bash
pyinstaller ../src/metabolite.py --onefile -n metabolite-x86_64-apple-darwin
pyinstaller ../src/metabolite.py --onefile -n metabolite-aarch64-apple-darwin
pyinstaller ../src/metabolite_for_db.py --onefile -n metabolite_for_db-x86_64-apple-darwin
pyinstaller ../src/metabolite_for_db.py --onefile -n metabolite_for_db-aarch64-apple-darwin

libpyg="$(python -c "import importlib.util; s = importlib.util.find_spec('libpyg'); print(s.origin if s else '')")"

if [[ -z "$libpyg" || ! -f "$libpyg" ]]; then
  echo "Could not locate libpyg in this Python environment."
  python -m pip show pyg-lib
  exit 1
fi

echo "Bundling: $libpyg"

pyinstaller ../ccs_prediction/model.py \
    --onefile --noconfirm --clean --console \
    --collect-all torch_geometric --collect-all pyg_lib --collect-all rdkit --copy-metadata torch-geometric \
    --collect-binaries torch --copy-metadata torch-geometric --hidden-import pyg_lib --hidden-import libpyg --add-binary "${libpyg}:." \
    --add-data "../ccs_prediction/trained_weights.pt:." -n prediction_model-x86_64-pc-apple-darwin

pyinstaller ../ccs_prediction/model.py \
    --onefile --noconfirm --clean --console \
    --collect-all torch_geometric --collect-all pyg_lib --collect-all rdkit --copy-metadata torch-geometric \
    --collect-binaries torch --copy-metadata torch-geometric --hidden-import pyg_lib --hidden-import libpyg --add-binary "${libpyg}:." \
    --add-data "../ccs_prediction/trained_weights.pt:." -n prediction_model-aarch64-apple-darwin
