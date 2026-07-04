#!/bin/bash
# Create and populate a Python virtual environment for SIMPLE.
# Usage: ./setup-venv.sh [--no-pysimple] [--recreate]
#
# Options:
#   --no-pysimple   Skip building pysimple (Fortran-Python bindings).
#                   Use this if you only need Python test/plot dependencies.
#   --recreate      Delete an existing .venv before creating a fresh one.
#
# Prerequisites (system packages):
#   - python3 (>= 3.9) with venv module
#   - gfortran, cmake, ninja-build
#   - libnetcdff-dev (or netcdf-fortran-devel), liblapack-dev
#   See README.md for distro-specific install commands.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV_DIR="${SCRIPT_DIR}/.venv"
BUILD_PYSIMPLE=1
RECREATE=0

for arg in "$@"; do
    case "$arg" in
        --no-pysimple) BUILD_PYSIMPLE=0 ;;
        --recreate) RECREATE=1 ;;
        -h|--help)
            sed -n '2,/^$/s/^# //p' "$0"
            exit 0
            ;;
        *)
            echo "Unknown option: $arg" >&2
            exit 1
            ;;
    esac
done

if [ "$RECREATE" -eq 1 ] && [ -d "$VENV_DIR" ]; then
    echo "Removing existing virtual environment ${VENV_DIR} ..."
    rm -rf "$VENV_DIR"
fi

if [ -d "$VENV_DIR" ]; then
    echo "Reusing virtual environment in ${VENV_DIR} ..."
else
    echo "Creating virtual environment in ${VENV_DIR} ..."
    python3 -m venv "$VENV_DIR"
fi

# shellcheck disable=SC1091
source "${VENV_DIR}/bin/activate"

echo "Upgrading pip ..."
python -m pip install --upgrade pip

echo "Installing Python dependencies from requirements.txt ..."
python -m pip install --prefer-binary -r "${SCRIPT_DIR}/requirements.txt"

if [ "$BUILD_PYSIMPLE" -eq 1 ]; then
    echo "Building pysimple (Fortran-Python bindings) ..."
    # Build SIMPLE first if not already built
    if [ ! -f "${SCRIPT_DIR}/build/build.ninja" ]; then
        echo "  Configuring CMake ..."
        cmake -S "$SCRIPT_DIR" -B"${SCRIPT_DIR}/build" -GNinja \
            -DCMAKE_BUILD_TYPE=Release -DCMAKE_COLOR_DIAGNOSTICS=ON
    fi
    echo "  Building Fortran library ..."
    cmake --build "${SCRIPT_DIR}/build" --config Release

    echo "  Installing pysimple in editable mode ..."
    python -m pip install --no-build-isolation -e "$SCRIPT_DIR"
else
    echo "Skipping pysimple build (--no-pysimple)."
fi

echo ""
echo "Done. Activate with:"
echo "  source .venv/bin/activate"
