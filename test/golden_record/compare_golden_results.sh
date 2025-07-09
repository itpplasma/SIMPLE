#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage() {
    echo "Usage: $0 REFERENCE_DIR CURRENT_DIR [TEST_CASES]"
    echo ""
    echo "Compares golden record test results between two directories"
    echo ""
    echo "Arguments:"
    echo "  REFERENCE_DIR  Directory containing reference test results"
    echo "  CURRENT_DIR    Directory containing current test results to compare"
    echo "  TEST_CASES     Optional space-separated list of test cases to compare"
    echo "                 (default: all test cases found in script directory)"
    echo ""
    echo "Example:"
    echo "  $0 /tmp/run_ref /tmp/run_current"
    echo "  $0 /tmp/run_ref /tmp/run_current \"boozer canonical\""
    exit 1
}

if [ $# -lt 2 ]; then
    usage
fi

REFERENCE_DIR="$1"
CURRENT_DIR="$2"
TEST_CASES_ARG="$3"

if [ ! -d "$REFERENCE_DIR" ]; then
    echo "Error: Reference directory not found: $REFERENCE_DIR"
    exit 1
fi

if [ ! -d "$CURRENT_DIR" ]; then
    echo "Error: Current directory not found: $CURRENT_DIR"
    exit 1
fi

# If test cases not specified, find all available test cases
if [ -z "$TEST_CASES_ARG" ]; then
    TEST_CASES="$(cd "$SCRIPT_DIR" && find . -name simple.in -exec dirname {} \; | sed 's|^\./||' | sort)"
else
    TEST_CASES="$TEST_CASES_ARG"
fi

echo "Reference directory: $REFERENCE_DIR"
echo "Current directory: $CURRENT_DIR"
echo "Test cases to compare: $(echo "$TEST_CASES" | wc -w)"
echo ""

compare_cases() {
    local total_cases=0
    local passed_cases=0
    local failed_cases=0
    local missing_cases=0
    
    for CASE in $TEST_CASES; do
        total_cases=$((total_cases + 1))
        
        echo "Comparing $CASE case..."
        
        # Check if this is the classifier_fast case with multiple files
        if [ "$CASE" = "classifier_fast" ]; then
            # List of files to compare for classifier_fast (excluding simple.in and wout.nc)
            # Note: fort.* files are excluded due to non-deterministic ordering in parallel execution
            CLASSIFIER_FILES="avg_inverse_t_lost.dat class_parts.dat confined_fraction.dat healaxis.dat start.dat times_lost.dat"
            
            # Run multi-file comparison
            python "$SCRIPT_DIR/compare_files_multi.py" "$REFERENCE_DIR/$CASE" "$CURRENT_DIR/$CASE" --files $CLASSIFIER_FILES
            
            if [ $? -eq 0 ]; then
                echo "  ✓ PASSED"
                passed_cases=$((passed_cases + 1))
            else
                echo "  ✗ FAILED"
                failed_cases=$((failed_cases + 1))
            fi
        else
            # Original single-file comparison for other test cases
            REF_FILE="$REFERENCE_DIR/$CASE/times_lost.dat"
            CUR_FILE="$CURRENT_DIR/$CASE/times_lost.dat"
            
            if [ ! -f "$REF_FILE" ]; then
                echo "  ✗ Reference file missing: $REF_FILE"
                missing_cases=$((missing_cases + 1))
                continue
            fi
            
            if [ ! -f "$CUR_FILE" ]; then
                echo "  ✗ Current file missing: $CUR_FILE"
                missing_cases=$((missing_cases + 1))
                continue
            fi
            
            # Run comparison
            python "$SCRIPT_DIR/compare_files.py" "$REF_FILE" "$CUR_FILE"
            
            if [ $? -eq 0 ]; then
                echo "  ✓ PASSED"
                passed_cases=$((passed_cases + 1))
            else
                echo "  ✗ FAILED"
                failed_cases=$((failed_cases + 1))
            fi
        fi
        echo ""
    done
    
    # Summary
    echo "========================================="
    echo "Summary:"
    echo "  Total test cases: $total_cases"
    echo "  Passed: $passed_cases"
    echo "  Failed: $failed_cases"
    echo "  Missing files: $missing_cases"
    echo "========================================="
    
    # Return non-zero if any tests failed or had missing files
    if [ $failed_cases -gt 0 ] || [ $missing_cases -gt 0 ]; then
        return 1
    else
        return 0
    fi
}

# Main execution
compare_cases
exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo "All tests passed!"
else
    echo "Some tests failed or had missing files."
fi

exit $exit_code