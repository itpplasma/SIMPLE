#!/bin/bash
CLONE_URL="https://github.com/itpplasma/SIMPLE.git"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CUR_VER=$(git -C "$SCRIPT_DIR" describe --tags --always --dirty)
REF_VER=${1:-"main"}
TMP_DIR=${2:-"$(mktemp -d)"}
CLEANUP_ON_FAILURE=${3:-"false"}  # New parameter for cleanup behavior

if [ "$REF_VER" == "$CUR_VER" ]; then
    echo "Reference version and current version are the same. Exiting."
    exit 1
fi

PROJECT_ROOT_CUR=$(cd "$SCRIPT_DIR/../.." && pwd)
RUN_DIR_CUR="$TMP_DIR/run_$CUR_VER"

PROJECT_ROOT_REF="$TMP_DIR/simple_$REF_VER"
RUN_DIR_REF="$TMP_DIR/run_$REF_VER"

TEST_CASES="$(cd "$SCRIPT_DIR" && find . -name simple.in -exec dirname {} \; | sed 's|^\./||' | sort)"

echo "Temporary directory: $TMP_DIR"
echo "Reference version: $REF_VER"
echo "Current version: $CUR_VER"
echo "Project root (reference): $PROJECT_ROOT_REF"
echo "Project root (current): $PROJECT_ROOT_CUR"
echo "Test run directories:"
echo "  Current: $RUN_DIR_CUR"
echo "  Reference: $RUN_DIR_REF"

cleanup_temp_dir() {
    echo "Cleaning up temporary directory: $TMP_DIR"
    rm -rf "$TMP_DIR"
    echo "Cleanup completed."
}

handle_success() {
    echo "Tests completed successfully."
    if [ "$CLEANUP_ON_FAILURE" = "true" ]; then
        cleanup_temp_dir
    else
        echo "Temporary directory preserved at: $TMP_DIR"
        echo "To clean up manually: rm -rf '$TMP_DIR'"
    fi
}

handle_failure() {
    local exit_code=$1
    echo "Tests failed with exit code: $exit_code"
    if [ "$CLEANUP_ON_FAILURE" = "true" ]; then
        cleanup_temp_dir
    else
        echo "Temporary directory preserved for debugging at: $TMP_DIR"
        echo "To clean up manually: rm -rf '$TMP_DIR'"
    fi
    exit $exit_code
}

main() {
    if [ "$CLEANUP_ON_FAILURE" = "true" ]; then
        # For expected failure tests, don't exit on error immediately
        echo "Running in expected-failure mode with cleanup"
    else
        set -e  # Exit on any error to preserve temp directory for debugging
    fi
    
    build "$PROJECT_ROOT_CUR" || handle_failure $?
    clone "$REF_VER" "$PROJECT_ROOT_REF" || handle_failure $?
    build "$PROJECT_ROOT_REF" || handle_failure $?
    get_test_data || handle_failure $?
    run_cases "$PROJECT_ROOT_CUR" "$RUN_DIR_CUR" "$TEST_CASES" || handle_failure $?
    run_cases "$PROJECT_ROOT_REF" "$RUN_DIR_REF" "$TEST_CASES" || handle_failure $?
    
    run_comparisons "$RUN_DIR_REF" "$RUN_DIR_CUR" "$TEST_CASES"
    comparison_result=$?
    
    if [ $comparison_result -eq 0 ]; then
        handle_success
    else
        handle_failure $comparison_result
    fi
}

get_test_data() {
    cd "$TMP_DIR"
    # Fetch if not exists
    if [ -f "wout.nc" ]; then
        echo "wout.nc already exists. Skipping download."
        return
    fi
    echo "Fetching test data..."
    wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc > /dev/null 2>&1
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

run_cases() {
    local PROJECT_ROOT="$1"
    local RUN_DIR="$2"
    local TEST_CASES="$3"

    echo "Running test cases in $RUN_DIR"

    for CASE in $TEST_CASES; do
        echo "Running case: $CASE in $RUN_DIR/$CASE"
        mkdir -p "$RUN_DIR/$CASE"
        cp "$SCRIPT_DIR/$CASE/simple.in" "$RUN_DIR/$CASE/simple.in"
        cd "$RUN_DIR/$CASE"
        ln -sf "$TMP_DIR/wout.nc" wout.nc
        echo "  Working directory: $(pwd)"
        $PROJECT_ROOT/build/simple.x < simple.in > stdout.log 2>&1
        echo "  Generated files: $(ls -la | grep -E '\.(dat|log)$' | wc -l) output files"
    done

}

run_comparisons() {
    local REF_DIR="$1"
    local CUR_DIR="$2"
    local TEST_CASES="$3"

    for CASE in $TEST_CASES; do
        REF_FILE="$REF_DIR/$CASE/times_lost.dat"
        CUR_FILE="$CUR_DIR/$CASE/times_lost.dat"
        echo "Comparing $CASE case..."
        python "$SCRIPT_DIR/compare_files.py" "$REF_FILE" "$CUR_FILE"
    done
}

main "$@"
