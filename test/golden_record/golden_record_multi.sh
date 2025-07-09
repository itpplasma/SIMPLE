#!/bin/bash
# Golden record test script that runs working copy tests once
# and compares against multiple reference versions

CLONE_URL="https://github.com/itpplasma/SIMPLE.git"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CUR_VER=$(git -C "$SCRIPT_DIR" describe --tags --always --dirty)

# Parse arguments: reference versions and optional --expect-fail flags
REFS=()
EXPECT_FAIL=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --expect-fail)
            EXPECT_FAIL+=("true")
            shift
            ;;
        *)
            REFS+=("$1")
            EXPECT_FAIL+=("false")
            shift
            ;;
    esac
done

if [ ${#REFS[@]} -eq 0 ]; then
    echo "Usage: $0 <ref1> [--expect-fail] <ref2> [--expect-fail] <ref3> ..."
    echo "Example: $0 46f1f53 --expect-fail e25ca7e main"
    exit 1
fi

# Allow override of directories via environment variables
GOLDEN_RECORD_BASE_DIR=${GOLDEN_RECORD_BASE_DIR:-"$(pwd)/golden_record"}
PROJECT_ROOT_CUR=$(cd "$SCRIPT_DIR/../.." && pwd)

RUN_DIR_CUR="$GOLDEN_RECORD_BASE_DIR/runs/run_$CUR_VER"
TEST_DATA_DIR="$GOLDEN_RECORD_BASE_DIR/test_data"

# Find test cases
if [ -d "$GOLDEN_RECORD_BASE_DIR/test_cases" ]; then
    TEST_CASES="$(cd "$GOLDEN_RECORD_BASE_DIR/test_cases" && find . -name simple.in -exec dirname {} \; | sed 's|^\./||' | sort)"
else
    TEST_CASES="$(cd "$SCRIPT_DIR" && find . -name simple.in -exec dirname {} \; | sed 's|^\./||' | sort)"
fi

echo "Golden record multi-comparison test"
echo "==================================="
echo "Current version: $CUR_VER"
echo "Reference versions to compare: ${REFS[@]}"
echo "Expected failures: ${EXPECT_FAIL[@]}"
echo ""

# Function to clone and build a reference version
prepare_reference() {
    local ref_ver="$1"
    local project_root="$GOLDEN_RECORD_BASE_DIR/simple_$ref_ver"

    if [ ! -f "$project_root/build/simple.x" ]; then
        echo "Preparing reference version $ref_ver..."

        # Clone if needed
        if [ ! -d "$project_root" ]; then
            echo "  Cloning..."
            git clone --filter=blob:none --no-checkout "$CLONE_URL" "$project_root"
        fi

        cd "$project_root"
        git fetch --all --quiet
        git checkout "$ref_ver" --quiet

        # Build
        echo "  Building..."

        # Handle old project structure with libneo
        if [ -f "SRC/CMakeLists.txt" ] && grep -q "../libneo" "SRC/CMakeLists.txt" 2>/dev/null; then
            LIBNEO_PATH="/proj/plasma/CODE/ert/libneo"
            if [ -d "$LIBNEO_PATH" ]; then
                ln -sf "$LIBNEO_PATH" "../libneo" 2>/dev/null || true
            fi
        fi

        cmake -S . -Bbuild -GNinja -DCMAKE_BUILD_TYPE=Release -DENABLE_PYTHON_INTERFACE=OFF -DENABLE_GVEC=OFF > configure.log 2>&1
        if [ $? -ne 0 ]; then
            echo "  CMake configuration failed for $ref_ver. Check $project_root/configure.log"
            return 1
        fi

        cmake --build build --config Release > build.log 2>&1
        if [ $? -ne 0 ]; then
            echo "  Build failed for $ref_ver. Check $project_root/build.log"
            return 1
        fi

        echo "  Build complete for $ref_ver"
    else
        echo "Using existing build for $ref_ver"
    fi

    return 0
}

# Main execution
main() {
    local overall_status=0

    # Create base directories
    mkdir -p "$GOLDEN_RECORD_BASE_DIR"
    mkdir -p "$TEST_DATA_DIR"

    # Check current build
    if [ ! -f "$PROJECT_ROOT_CUR/build/simple.x" ]; then
        echo "ERROR: Current build not found at: $PROJECT_ROOT_CUR/build/simple.x"
        exit 1
    fi

    # Prepare all reference versions first
    echo "Step 1: Preparing reference versions"
    echo "------------------------------------"
    for ref in "${REFS[@]}"; do
        prepare_reference "$ref"
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to prepare reference version $ref"
            exit 1
        fi
    done

    # Run working copy tests ONCE
    echo ""
    echo "Step 2: Running tests on working copy"
    echo "-------------------------------------"
    "$SCRIPT_DIR/run_golden_tests.sh" "$PROJECT_ROOT_CUR" "$RUN_DIR_CUR" "$TEST_DATA_DIR"
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to run tests on working copy"
        exit 1
    fi

    # Compare against each reference version
    echo ""
    echo "Step 3: Comparing against reference versions"
    echo "--------------------------------------------"

    for i in "${!REFS[@]}"; do
        ref="${REFS[$i]}"
        expect_fail="${EXPECT_FAIL[$i]}"
        run_dir_ref="$GOLDEN_RECORD_BASE_DIR/runs/run_$ref"
        project_root_ref="$GOLDEN_RECORD_BASE_DIR/simple_$ref"

        echo ""
        echo "Comparing against $ref (expect failure: $expect_fail)..."

        # Run tests for reference version
        "$SCRIPT_DIR/run_golden_tests.sh" "$project_root_ref" "$run_dir_ref" "$TEST_DATA_DIR"
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to run tests for reference $ref"
            if [ "$expect_fail" = "false" ]; then
                overall_status=1
            fi
            continue
        fi

        # Compare results
        "$SCRIPT_DIR/compare_golden_results.sh" "$run_dir_ref" "$RUN_DIR_CUR"
        comparison_result=$?

        if [ $comparison_result -eq 0 ]; then
            echo "  ✓ Results match with $ref"
            if [ "$expect_fail" = "true" ]; then
                echo "  WARNING: Expected this comparison to fail, but it passed!"
                overall_status=1
            fi
        else
            echo "  ✗ Results differ from $ref"
            if [ "$expect_fail" = "false" ]; then
                echo "  ERROR: Unexpected difference!"
                overall_status=1
            else
                echo "  (This was expected to fail)"
            fi
        fi
    done

    # Summary
    echo ""
    echo "Summary"
    echo "======="
    echo "Working copy version: $CUR_VER"
    echo "Reference versions tested: ${REFS[@]}"
    echo "Results preserved in: $GOLDEN_RECORD_BASE_DIR"

    if [ $overall_status -eq 0 ]; then
        echo "Overall result: PASSED"
    else
        echo "Overall result: FAILED"
    fi

    exit $overall_status
}

main "$@"
