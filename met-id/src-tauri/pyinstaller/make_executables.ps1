
pyinstaller ../src/validation.py --onefile -n validation-x86_64-pc-windows-msvc.exe
pyinstaller ../src/metabolite.py --onefile -n metabolite-x86_64-pc-windows-msvc.exe
pyinstaller ../src/regression.py --onefile -n regression-x86_64-pc-windows-msvc.exe
pyinstaller ../src/metabolite_for_db.py --onefile -n metabolite_for_db-x86_64-pc-windows-msvc.exe
pyinstaller ../src/msms_to_db.py --onefile -n msms_to_db-x86_64-pc-windows-msvc.exe