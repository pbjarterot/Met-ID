#!/bin/bash
pyinstaller ../src/metabolite.py --onefile -n metabolite-x86_64-pc-linux-gnu
pyinstaller ../src/metabolite_for_db.py --onefile -n metabolite_for_db-x86_64-pc-linux-gnu