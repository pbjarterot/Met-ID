#!/bin/bash

set -e

# Developer ID and Apple credentials
APPLE_DEVELOPER_ID="${APPLE_DEVELOPER_ID:-Developer ID Application: Your Name (TEAM_ID)}"
NOTARYTOOL_KEYCHAIN_PROFILE="AC_API"

OUTPUT_DIR="dist"
BUILD_DIR="build"

echo "Unlocking keychain..."
security unlock-keychain -p "$APPLE_KEYCHAIN_PASSWORD" ~/Library/Keychains/build.keychain-db
security list-keychains -s ~/Library/Keychains/build.keychain-db
security default-keychain -s ~/Library/Keychains/build.keychain-db

# Cleanup
rm -rf "$OUTPUT_DIR" "$BUILD_DIR"

pip install --upgrade pyinstaller

echo "Building the executable..."
pyinstaller ../src/metabolite.py --onefile -n metabolite-aarch64-apple-darwin

echo "SHA256 before signing:"
shasum -a 256 dist/metabolite-aarch64-apple-darwin

echo "Signing the main binary with Hardened Runtime..."
codesign --force --verify --verbose --sign "$APPLE_DEVELOPER_ID" --options runtime dist/metabolite-aarch64-apple-darwin

echo "Verifying signed binary:"
codesign -d --verbose=4 dist/metabolite-aarch64-apple-darwin

echo "SHA256 after signing:"
shasum -a 256 dist/metabolite-aarch64-apple-darwin

echo "Creating .zip file for notarization..."
rm -f dist/metabolite-aarch64-apple-darwin.zip
ditto -c -k --keepParent dist/metabolite-aarch64-apple-darwin dist/metabolite-aarch64-apple-darwin.zip

echo "SHA256 of .zip before submission:"
shasum -a 256 dist/metabolite-aarch64-apple-darwin.zip

echo "Submitting for notarization..."
xcrun notarytool submit dist/metabolite-aarch64-apple-darwin.zip --keychain-profile "$NOTARYTOOL_KEYCHAIN_PROFILE" --wait

echo "Stapling notarization ticket..."
shasum -a 256 dist/metabolite-aarch64-apple-darwin
xcrun stapler staple dist/metabolite-aarch64-apple-darwin

xcrun stapler validate dist/metabolite-aarch64-apple-darwin
echo "Build, signing, notarization, and stapling completed successfully."
