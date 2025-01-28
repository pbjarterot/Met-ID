#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Developer ID and Apple credentials
APPLE_DEVELOPER_ID="${APPLE_DEVELOPER_ID:-Developer ID Application: Your Name (TEAM_ID)}"
NOTARYTOOL_KEYCHAIN_PROFILE="AC_API"

# Paths
OUTPUT_DIR="dist"
BUILD_DIR="build"

# Unlock and set the keychain
echo "Unlocking keychain..."
security unlock-keychain -p "$APPLE_KEYCHAIN_PASSWORD" ~/Library/Keychains/build.keychain-db
security list-keychains -s ~/Library/Keychains/build.keychain-db
security default-keychain -s ~/Library/Keychains/build.keychain-db

# Cleanup previous builds
if [ -d "$OUTPUT_DIR" ]; then
  rm -rf "$OUTPUT_DIR"
fi
if [ -d "$BUILD_DIR" ]; then
  rm -rf "$BUILD_DIR"
fi

# Ensure PyInstaller is installed
pip install --upgrade pyinstaller

# Build executables
pyinstaller ../src/metabolite.py --onefile -n metabolite-aarch64-apple-darwin

# Sign binaries
echo "Signing binaries with Developer ID: $APPLE_DEVELOPER_ID..."
for BINARY in $OUTPUT_DIR/*; do
  echo "Signing $BINARY..."
  codesign --deep --force --verify --verbose --sign "$APPLE_DEVELOPER_ID" "$BINARY"
  codesign --verify --verbose "$BINARY"
done

# Create a zip file for notarization
echo "Creating zip file for notarization..."
for BINARY in $OUTPUT_DIR/*; do
  zip "${BINARY}.zip" "$BINARY"
done

# Notarization
echo "Submitting zip file for notarization..."
for ZIPFILE in $OUTPUT_DIR/*.zip; do
  echo "Submitting $ZIPFILE..."
  xcrun notarytool submit "$ZIPFILE" --keychain-profile "$NOTARYTOOL_KEYCHAIN_PROFILE" --wait
done

# Validate stapling
echo "Stapling notarization ticket..."
for BINARY in $OUTPUT_DIR/*; do
  if [[ "$BINARY" != *.zip ]]; then
    echo "Stapling $BINARY..."
    xcrun stapler staple "$BINARY"
    xcrun stapler validate "$BINARY"
  fi
done

echo "Build, signing, and packaging complete."
