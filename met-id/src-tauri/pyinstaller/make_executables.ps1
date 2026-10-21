pyinstaller ../src/metabolite.py --onefile -n metabolite-x86_64-pc-windows-msvc.exe

pyinstaller ../src/metabolite_for_db.py --onefile -n metabolite_for_db-x86_64-pc-windows-msvc.exe

$libpyg = (& python -c "import importlib.util; s=importlib.util.find_spec('libpyg'); print(s.origin if s else '')").Trim()

if ([string]::IsNullOrWhiteSpace($libpyg) -or -not (Test-Path -LiteralPath $libpyg)) {
    throw "Could not locate libpyg.pyd"
}

Write-Host "Bundling: $libpyg"

pyinstaller ../ccs_prediction/model.py `
    --onefile --noconfirm --clean --console `
    --collect-all torch_geometric --collect-all pyg_lib --collect-all rdkit --copy-metadata torch-geometric `
    --collect-binaries torch --copy-metadata torch-geometric --hidden-import pyg_lib --hidden-import libpyg --add-binary "${libpyg}:." `
    --add-data "../ccs_prediction/trained_weights.pt:." -n prediction_model-x86_64-pc-windows-msvc.exe

