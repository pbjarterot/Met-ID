#!/bin/bash
pyinstaller ../src/metabolite.py --onefile -n metabolite-x86_64-apple-darwin
pyinstaller ../src/metabolite.py --onefile -n metabolite-aarch64-apple-darwin
pyinstaller ../src/metabolite_for_db.py --onefile -n metabolite_for_db-x86_64-apple-darwin
pyinstaller ../src/metabolite_for_db.py --onefile -n metabolite_for_db-aarch64-apple-darwin