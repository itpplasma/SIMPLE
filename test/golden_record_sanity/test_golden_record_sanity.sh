#!/bin/bash
# Golden record sanity test script
# Tests the golden record comparison system itself

set -e

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
GOLDEN_RECORD_DIR=$(cd "$SCRIPT_DIR/.." && pwd)/golden_record

# Colors for output (disabled if NO_COLOR is set)
if [ -n "$NO_COLOR" ]; then
    RED=''
    GREEN=''
    YELLOW=''
    NC=''
else
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    NC='\033[0m' # No Color
fi

echo "Golden Record Sanity Tests"
echo "=========================="
echo ""

# Test counters
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Function to run a test
run_test() {
    local test_name="$1"
    local ref_dir="$2"
    local cur_dir="$3"
    local expected_result="$4"  # "pass" or "fail"
    local test_cases="$5"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    echo -n "Test: $test_name ... "
    
    # Run the comparison
    "$GOLDEN_RECORD_DIR/compare_golden_results.sh" "$ref_dir" "$cur_dir" "$test_cases" > /tmp/golden_test_output.txt 2>&1 </dev/null
    result=$?
    
    if [ "$expected_result" = "pass" ]; then
        if [ $result -eq 0 ]; then
            echo -e "${GREEN}PASSED${NC}"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            echo -e "${RED}FAILED${NC} (expected pass but got fail)"
            echo "Output:"
            cat /tmp/golden_test_output.txt | sed 's/^/  /'
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    else  # expected_result = "fail"
        if [ $result -ne 0 ]; then
            echo -e "${GREEN}PASSED${NC} (correctly detected differences)"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            echo -e "${RED}FAILED${NC} (expected fail but got pass)"
            echo "Output:"
            cat /tmp/golden_test_output.txt | sed 's/^/  /'
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    fi
}

# Test 1: Matching simple_test case
run_test "simple_test matching" \
    "$SCRIPT_DIR/reference" \
    "$SCRIPT_DIR/matching" \
    "pass" \
    "simple_test"

# Test 2: Non-matching simple_test case
run_test "simple_test non-matching" \
    "$SCRIPT_DIR/reference" \
    "$SCRIPT_DIR/non_matching" \
    "fail" \
    "simple_test"

# Test 3: Matching classifier_fast case (all files)
run_test "classifier_fast matching (all files)" \
    "$SCRIPT_DIR/reference" \
    "$SCRIPT_DIR/matching" \
    "pass" \
    "classifier_fast"

# Test 4: Non-matching classifier_fast case (multiple files differ)
run_test "classifier_fast non-matching (multiple files)" \
    "$SCRIPT_DIR/reference" \
    "$SCRIPT_DIR/non_matching" \
    "fail" \
    "classifier_fast"

# Test 5: Both test cases at once - matching
run_test "all test cases matching" \
    "$SCRIPT_DIR/reference" \
    "$SCRIPT_DIR/matching" \
    "pass" \
    "simple_test classifier_fast"

# Test 6: Both test cases at once - non-matching
run_test "all test cases non-matching" \
    "$SCRIPT_DIR/reference" \
    "$SCRIPT_DIR/non_matching" \
    "fail" \
    "simple_test classifier_fast"

# Test 7: Test the multi-file comparison directly
echo ""
echo "Direct multi-file comparison tests:"
echo "-----------------------------------"

# Test matching files
echo -n "Test: Direct comparison matching ... "
python "$GOLDEN_RECORD_DIR/compare_files_multi.py" \
    "$SCRIPT_DIR/reference/classifier_fast" \
    "$SCRIPT_DIR/matching/classifier_fast" \
    --files times_lost.dat confined_fraction.dat fort.10022 > /tmp/direct_test.txt 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}PASSED${NC}"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    echo -e "${RED}FAILED${NC}"
    cat /tmp/direct_test.txt | sed 's/^/  /'
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi
TOTAL_TESTS=$((TOTAL_TESTS + 1))

# Test non-matching files
echo -n "Test: Direct comparison non-matching ... "
python "$GOLDEN_RECORD_DIR/compare_files_multi.py" \
    "$SCRIPT_DIR/reference/classifier_fast" \
    "$SCRIPT_DIR/non_matching/classifier_fast" \
    --files times_lost.dat confined_fraction.dat fort.10022 > /tmp/direct_test.txt 2>&1
if [ $? -ne 0 ]; then
    echo -e "${GREEN}PASSED${NC} (correctly detected differences)"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    echo -e "${RED}FAILED${NC} (should have detected differences)"
    cat /tmp/direct_test.txt | sed 's/^/  /'
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi
TOTAL_TESTS=$((TOTAL_TESTS + 1))

# Summary
echo ""
echo "================================"
echo "Summary:"
echo "  Total tests: $TOTAL_TESTS"
echo -e "  Passed: ${GREEN}$PASSED_TESTS${NC}"
echo -e "  Failed: ${RED}$FAILED_TESTS${NC}"
echo "================================"

# Clean up
rm -f /tmp/golden_test_output.txt /tmp/direct_test.txt

# Exit with appropriate code
if [ $FAILED_TESTS -gt 0 ]; then
    echo -e "\n${RED}Some tests failed!${NC}"
    exit 1
else
    echo -e "\n${GREEN}All tests passed!${NC}"
    exit 0
fi