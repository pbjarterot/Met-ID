pyinstaller ../src/metabolite.py --onefile -n metabolite-x86_64-pc-windows-msvc.exe
pyinstaller ../src/metabolite_for_db.py --onefile -n metabolite_for_db-x86_64-pc-windows-msvc.exe
pyinstaller ../src/matrix_for_db.py --onefile -n matrix_for_db-x86_64-pc-windows-msvc.exe