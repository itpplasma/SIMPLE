#!/bin/bash
CLONE_URL="https://github.com/itpplasma/SIMPLE.git"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CUR_VER=$(git -C "$SCRIPT_DIR" describe --tags --always --dirty)
REF_VER=${1:-"main"}

# Allow override of directories via environment variables
GOLDEN_RECORD_BASE_DIR=${GOLDEN_RECORD_BASE_DIR:-"$(pwd)/golden_record"}
GOLDEN_RECORD_REF_DIR=${GOLDEN_RECORD_REF_DIR:-""}
GOLDEN_RECORD_CUR_DIR=${GOLDEN_RECORD_CUR_DIR:-""}

# Set up deterministic directories
if [ -z "$GOLDEN_RECORD_REF_DIR" ]; then
    PROJECT_ROOT_REF="$GOLDEN_RECORD_BASE_DIR/simple_$REF_VER"
else
    PROJECT_ROOT_REF="$GOLDEN_RECORD_REF_DIR"
fi

if [ -z "$GOLDEN_RECORD_CUR_DIR" ]; then
    PROJECT_ROOT_CUR=$(cd "$SCRIPT_DIR/../.." && pwd)
else
    PROJECT_ROOT_CUR="$GOLDEN_RECORD_CUR_DIR"
fi

RUN_DIR_REF="$GOLDEN_RECORD_BASE_DIR/runs/run_$REF_VER"
RUN_DIR_CUR="$GOLDEN_RECORD_BASE_DIR/runs/run_$CUR_VER"
TEST_DATA_DIR="$GOLDEN_RECORD_BASE_DIR/test_data"

# Find test cases - they should be copied by CMake to the golden_record directory
if [ -d "$GOLDEN_RECORD_BASE_DIR/test_cases" ]; then
    TEST_CASES="$(cd "$GOLDEN_RECORD_BASE_DIR/test_cases" && find . -name simple.in -exec dirname {} \; | sed 's|^\./||' | sort)"
else
    # Fallback to original location if CMake hasn't copied them yet
    TEST_CASES="$(cd "$SCRIPT_DIR" && find . -name simple.in -exec dirname {} \; | sed 's|^\./||' | sort)"
fi

echo "Golden record base directory: $GOLDEN_RECORD_BASE_DIR"
echo "Reference version: $REF_VER"
echo "Current version: $CUR_VER"
echo "Project root (reference): $PROJECT_ROOT_REF"
echo "Project root (current): $PROJECT_ROOT_CUR"
echo "Test run directories:"
echo "  Current: $RUN_DIR_CUR"
echo "  Reference: $RUN_DIR_REF"

handle_success() {
    echo "Tests completed successfully."
    echo "Results preserved in: $GOLDEN_RECORD_BASE_DIR"
}

handle_failure() {
    local exit_code=$1
    echo "Tests failed with exit code: $exit_code"
    echo "Results preserved for debugging in: $GOLDEN_RECORD_BASE_DIR"
    exit $exit_code
}

main() {
    set -e  # Exit on any error
    
    # Create base directories
    mkdir -p "$GOLDEN_RECORD_BASE_DIR"
    mkdir -p "$TEST_DATA_DIR"
    
    # Check if we need to build reference version
    if [ ! -f "$PROJECT_ROOT_REF/build/simple.x" ]; then
        echo "Reference build not found, cloning and building..."
        clone "$REF_VER" "$PROJECT_ROOT_REF" || handle_failure $?
        build "$PROJECT_ROOT_REF" || handle_failure $?
    else
        echo "Using existing reference build at: $PROJECT_ROOT_REF/build/simple.x"
    fi
    
    # Check if current version needs building (it should already be built)
    if [ ! -f "$PROJECT_ROOT_CUR/build/simple.x" ]; then
        echo "Current build not found at: $PROJECT_ROOT_CUR/build/simple.x"
        handle_failure 1
    else
        echo "Using current build at: $PROJECT_ROOT_CUR/build/simple.x"
    fi
    
    # Use the new scripts to run tests and compare
    "$SCRIPT_DIR/run_golden_tests.sh" "$PROJECT_ROOT_REF" "$RUN_DIR_REF" "$TEST_DATA_DIR" || handle_failure $?
    "$SCRIPT_DIR/run_golden_tests.sh" "$PROJECT_ROOT_CUR" "$RUN_DIR_CUR" "$TEST_DATA_DIR" || handle_failure $?
    
    # Compare results
    "$SCRIPT_DIR/compare_golden_results.sh" "$RUN_DIR_REF" "$RUN_DIR_CUR"
    comparison_result=$?
    
    if [ $comparison_result -eq 0 ]; then
        handle_success
    else
        handle_failure $comparison_result
    fi
}


clone() {
    local VERSION="$1"
    local PROJECT_ROOT="$2"
    if [ ! -d "$PROJECT_ROOT" ]; then
        echo "Cloning SIMPLE version $VERSION"
        git clone --filter=blob:none --no-checkout "$CLONE_URL" "$PROJECT_ROOT"
    fi

    cd "$PROJECT_ROOT"
    git fetch --all --quiet
    git checkout "$VERSION" --quiet
}

build() {
    local PROJECT_ROOT="$1"
    echo "Building SIMPLE in $PROJECT_ROOT"
    cd $PROJECT_ROOT
    cmake -S . -Bbuild -GNinja -DCMAKE_BUILD_TYPE=Release -DENABLE_PYTHON_INTERFACE=OFF  > $PROJECT_ROOT/configure.log 2>&1
    cmake --build build --config Release  > $PROJECT_ROOT/build.log 2>&1
}


main "$@"
