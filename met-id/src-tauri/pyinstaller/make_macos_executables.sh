#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Developer ID Application Certificate
DEVELOPER_ID= ${DEVELOPER_ID:-"Developer ID Application: Unknown"}
TEAM_ID= ${TEAM_ID:-"UnknownTeamID"}  

# Paths
OUTPUT_DIR="dist"
BUILD_DIR="build"

# Cleanup previous builds
if [ -d "$OUTPUT_DIR" ]; then
  rm -rf "$OUTPUT_DIR"
fi
if [ -d "$BUILD_DIR" ]; then
  rm -rf "$BUILD_DIR"
fi

# Ensure PyInstaller is installed
pip install --upgrade pyinstaller

# Build executables for both architectures
pyinstaller ../src/metabolite.py --onefile -n metabolite-x86_64-apple-darwin
pyinstaller ../src/metabolite.py --onefile -n metabolite-aarch64-apple-darwin
pyinstaller ../src/metabolite_for_db.py --onefile -n metabolite_for_db-x86_64-apple-darwin
pyinstaller ../src/metabolite_for_db.py --onefile -n metabolite_for_db-aarch64-apple-darwin

# Code signing and verification for each binary
echo "Signing binaries..."
for BINARY in $OUTPUT_DIR/*; do
  echo "Signing $BINARY..."
  codesign --deep --force --verify --verbose --sign "$DEVELOPER_ID" "$BINARY"
  codesign --verify --verbose "$BINARY"
  spctl -a -t exec -vv "$BINARY"
done

# Zip all binaries for distribution
echo "Creating zip files for distribution..."
cd $OUTPUT_DIR
zip -r metabolite-macos.zip metabolite-*
zip -r metabolite_for_db-macos.zip metabolite_for_db-*

echo "Build, signing, and packaging complete. Files available in 'dist/'"
