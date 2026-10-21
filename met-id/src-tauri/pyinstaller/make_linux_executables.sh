#!/bin/bash
pyinstaller ../src/metabolite.py --onefile -n metabolite-x86_64-unknown-linux-gnu
pyinstaller ../src/metabolite_for_db.py --onefile -n metabolite_for_db-x86_64-unknown-linux-gnu

pyinstaller ../ccs_prediction/model.py \
    --onefile --noconfirm --clean --console \
    --collect-all torch_geometric --collect-all pyg_lib --collect-all rdkit --copy-metadata torch-geometric \
    --collect-binaries torch --copy-metadata torch-geometric --hidden-import pyg_lib --hidden-import libpyg --add-binary "${libpyg}:." \
    --add-data "../ccs_prediction/trained_weights.pt:." -n prediction_model-x86_64-unknown-linux-gnu

