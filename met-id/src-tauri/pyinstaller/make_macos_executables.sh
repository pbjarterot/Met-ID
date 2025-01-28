#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Developer ID and Apple credentials
APPLE_DEVELOPER_ID="${APPLE_DEVELOPER_ID:-Developer ID Application: Your Name (APPLE_TEAM_ID)}"
APPLE_KEYCHAIN_PASSWORD="${APPLE_KEYCHAIN_PASSWORD:-}"
APPLE_ID="${APPLE_ID:-}"
APPLE_PASSWORD="${APPLE_PASSWORD:-}"

# Paths
OUTPUT_DIR="dist"
BUILD_DIR="build"

# Unlock the keychain
echo "Unlocking keychain..."
security unlock-keychain -p "$APPLE_KEYCHAIN_PASSWORD" ~/Library/Keychains/login.keychain-db
security list-keychains -s ~/Library/Keychains/login.keychain-db
security default-keychain -s ~/Library/Keychains/login.keychain-db

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
pyinstaller ../src/metabolite.py --onefile -n metabolite-aarch64-apple-darwin

# Code signing
echo "Signing binaries with Developer ID: $APPLE_DEVELOPER_ID..."
for BINARY in $OUTPUT_DIR/*; do
  echo "Signing $BINARY..."
  codesign --deep --force --verify --verbose --sign "$APPLE_DEVELOPER_ID" "$BINARY"
  codesign --verify --verbose "$BINARY"
done

# Notarization
echo "Submitting binary for notarization..."
for BINARY in $OUTPUT_DIR/*; do
  xcrun altool --notarize-app \
    --primary-bundle-id "com.farmbio.metid" \
    --username "$APPLE_ID" \
    --password "$APPLE_PASSWORD" \
    --file "$BINARY"

  echo "Stapling notarization ticket..."
  xcrun stapler staple "$BINARY"
done

# Zip all binaries for distribution
echo "Creating zip files for distribution..."
cd $OUTPUT_DIR
zip -r metabolite-macos.zip metabolite-*

echo "Build, signing, and packaging complete. Files available in 'dist/'"
