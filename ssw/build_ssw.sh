#!/usr/bin/env bash
set -euo pipefail

# Build a local ssw_test binary under foldseek_new_3di/ssw
# Usage: ./build_ssw.sh
# This will clone the upstream Complete-Striped-Smith-Waterman-Library into tmp/ssw
# build it, and copy the resulting ssw_test into this directory as ./ssw_test

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP_DIR="$ROOT_DIR/tmp/ssw"
OUT_BIN="$ROOT_DIR/ssw_test"

echo "Working in: $ROOT_DIR"
mkdir -p "$ROOT_DIR/tmp"

if [ -f "$OUT_BIN" ]; then
  echo "Found existing $OUT_BIN — skipping build. Remove it to force rebuild."
  exit 0
fi

# If conda environment is active, export include/lib flags so the build can find zlib and other headers/libs
if [ -n "${CONDA_PREFIX:-}" ]; then
  echo "Detected CONDA_PREFIX=$CONDA_PREFIX — exporting CPPFLAGS/LDFLAGS to help the build find headers/libs"
  export CPPFLAGS="-I$CONDA_PREFIX/include ${CPPFLAGS:-}"
  export LDFLAGS="-L$CONDA_PREFIX/lib ${LDFLAGS:-}"
  # conservative extras
  export C_INCLUDE_PATH="$CONDA_PREFIX/include:${C_INCLUDE_PATH:-}"
  export LIBRARY_PATH="$CONDA_PREFIX/lib:${LIBRARY_PATH:-}"
else
  echo "CONDA_PREFIX not set — proceeding without extra include/lib flags. If build fails due to missing zlib, set CONDA_PREFIX or install system zlib-dev packages."
fi

if [ -d "$TMP_DIR" ]; then
  echo "Using existing clone at $TMP_DIR"
else
  echo "Cloning Complete-Striped-Smith-Waterman-Library into $TMP_DIR"
  git clone --depth 1 https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.git "$TMP_DIR"
fi

echo "Building ssw (this will run 'make' in $TMP_DIR/src)"
pushd "$TMP_DIR/src" > /dev/null
make
popd > /dev/null

if [ ! -f "$TMP_DIR/src/ssw_test" ]; then
  echo "Build failed: $TMP_DIR/src/ssw_test not found" >&2
  exit 1
fi

cp -v "$TMP_DIR/src/ssw_test" "$OUT_BIN"
chmod +x "$OUT_BIN"

echo "Built and installed: $OUT_BIN"

# Quick smoke test: print help line (non-failing)
echo "ssw_test --help (first lines):"
"$OUT_BIN" -h 2>&1 | sed -n '1,20p' || true

exit 0
