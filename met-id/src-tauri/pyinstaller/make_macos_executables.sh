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
echo "Building the executable..."
pyinstaller ../src/metabolite.py --onefile -n metabolite-aarch64-apple-darwin

# Sign the main binary
echo "Signing the main binary with Developer ID: $APPLE_DEVELOPER_ID..."
codesign --force --verify --verbose --sign "$APPLE_DEVELOPER_ID" dist/metabolite-aarch64-apple-darwin

# Re-sign all embedded components
echo "Re-signing all embedded components..."
find dist/metabolite-aarch64-apple-darwin -type f -exec codesign --force --verify --verbose --sign "$APPLE_DEVELOPER_ID" {} \;

# Verify the entire binary
echo "Verifying the signed binary..."
codesign --verify --deep --verbose dist/metabolite-aarch64-apple-darwin

# Create a zip file for notarization
echo "Creating zip file for notarization..."
cd dist
zip metabolite-aarch64-apple-darwin.zip metabolite-aarch64-apple-darwin
cd ..

# Notarization
echo "Submitting zip file for notarization..."
xcrun notarytool submit dist/metabolite-aarch64-apple-darwin.zip --keychain-profile "$NOTARYTOOL_KEYCHAIN_PROFILE" --wait

echo "Verifying Notarization"
xcrun notarytool history --keychain-profile "$NOTARYTOOL_KEYCHAIN_PROFILE"

# Staple notarization ticket
echo "Stapling notarization ticket..."
xcrun stapler staple dist/metabolite-aarch64-apple-darwin

# Verify stapling
xcrun stapler validate dist/metabolite-aarch64-apple-darwin

echo "Build, signing, notarization, and stapling completed successfully."
