#!/bin/bash
# Simplified golden record sanity test

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
GOLDEN_RECORD_DIR=$(cd "$SCRIPT_DIR/.." && pwd)/golden_record

echo "Golden Record Sanity Tests"
echo "=========================="
echo ""

# Test 1: Single file comparison - matching
echo -n "Test 1: Single file comparison (matching)... "
python "$GOLDEN_RECORD_DIR/compare_files.py" \
    "$SCRIPT_DIR/reference/simple_test/times_lost.dat" \
    "$SCRIPT_DIR/matching/simple_test/times_lost.dat" >/dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "PASSED"
else
    echo "FAILED"
fi

# Test 2: Single file comparison - non-matching
echo -n "Test 2: Single file comparison (non-matching)... "
python "$GOLDEN_RECORD_DIR/compare_files.py" \
    "$SCRIPT_DIR/reference/simple_test/times_lost.dat" \
    "$SCRIPT_DIR/non_matching/simple_test/times_lost.dat" >/dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "PASSED (correctly detected difference)"
else
    echo "FAILED (should have detected difference)"
fi

# Test 3: Multi-file comparison - matching
echo -n "Test 3: Multi-file comparison (matching)... "
python "$GOLDEN_RECORD_DIR/compare_files_multi.py" \
    "$SCRIPT_DIR/reference/classifier_fast" \
    "$SCRIPT_DIR/matching/classifier_fast" \
    --files times_lost.dat confined_fraction.dat fort.10022 >/dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "PASSED"
else
    echo "FAILED"
fi

# Test 4: Multi-file comparison - non-matching
echo -n "Test 4: Multi-file comparison (non-matching)... "
python "$GOLDEN_RECORD_DIR/compare_files_multi.py" \
    "$SCRIPT_DIR/reference/classifier_fast" \
    "$SCRIPT_DIR/non_matching/classifier_fast" \
    --files times_lost.dat confined_fraction.dat fort.10022 >/dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "PASSED (correctly detected differences)"
else
    echo "FAILED (should have detected differences)"
fi

# Test 5: Full comparison script - simple_test matching
echo -n "Test 5: Full comparison script (simple_test matching)... "
"$GOLDEN_RECORD_DIR/compare_golden_results.sh" \
    "$SCRIPT_DIR/reference" \
    "$SCRIPT_DIR/matching" \
    "simple_test" >/dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "PASSED"
else
    echo "FAILED"
fi

# Test 6: Full comparison script - classifier_fast non-matching
echo -n "Test 6: Full comparison script (classifier_fast non-matching)... "
"$GOLDEN_RECORD_DIR/compare_golden_results.sh" \
    "$SCRIPT_DIR/reference" \
    "$SCRIPT_DIR/non_matching" \
    "classifier_fast" >/dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "PASSED (correctly detected differences)"
else
    echo "FAILED (should have detected differences)"
fi

echo ""
echo "All tests completed!"