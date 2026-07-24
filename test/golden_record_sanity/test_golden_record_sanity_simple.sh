#!/bin/bash
# Simplified golden record sanity test

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
GOLDEN_RECORD_DIR=$(cd "$SCRIPT_DIR/.." && pwd)/golden_record

echo "Golden Record Sanity Tests"
echo "=========================="
echo ""
status=0

# Test 1: Single file comparison - matching
echo -n "Test 1: Single file comparison (matching)... "
python "$GOLDEN_RECORD_DIR/compare_files.py" \
    "$SCRIPT_DIR/reference/simple_test/times_lost.dat" \
    "$SCRIPT_DIR/matching/simple_test/times_lost.dat" >/dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "PASSED"
else
    echo "FAILED"
    status=1
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
    status=1
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
    status=1
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
    status=1
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
    status=1
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
    status=1
fi

echo -n "Test 7: Paired NaNs compare equal... "
python "$GOLDEN_RECORD_DIR/compare_files.py" \
    "$SCRIPT_DIR/reference/nan_test/times_lost.dat" \
    "$SCRIPT_DIR/matching/nan_test/times_lost.dat" >/dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "PASSED"
else
    echo "FAILED"
    status=1
fi

echo -n "Test 8: One-sided NaN compares unequal... "
python "$GOLDEN_RECORD_DIR/compare_files.py" \
    "$SCRIPT_DIR/reference/nan_test/times_lost.dat" \
    "$SCRIPT_DIR/non_matching/nan_test/times_lost.dat" >/dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "PASSED (correctly detected difference)"
else
    echo "FAILED (should have detected difference)"
    status=1
fi

echo -n "Test 9: Multi-file paired NaNs compare equal... "
python "$GOLDEN_RECORD_DIR/compare_files_multi.py" \
    "$SCRIPT_DIR/reference/nan_test" \
    "$SCRIPT_DIR/matching/nan_test" \
    --files times_lost.dat >/dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "PASSED"
else
    echo "FAILED"
    status=1
fi

echo -n "Test 10: Multi-file one-sided NaN compares unequal... "
python "$GOLDEN_RECORD_DIR/compare_files_multi.py" \
    "$SCRIPT_DIR/reference/nan_test" \
    "$SCRIPT_DIR/non_matching/nan_test" \
    --files times_lost.dat >/dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "PASSED (correctly detected difference)"
else
    echo "FAILED (should have detected difference)"
    status=1
fi

SANITY_TMP=$(mktemp -d)
trap 'rm -rf "$SANITY_TMP"' EXIT
for variant in reference matching non_matching; do
    mkdir -p "$SANITY_TMP/$variant/albert_coils"
    printf '1.0 2.0\n' > \
        "$SANITY_TMP/$variant/albert_coils/albert_coils_diagnostic.dat"
    printf '3.0 4.0\n' > \
        "$SANITY_TMP/$variant/albert_coils/albert_intermediate.dat"
    printf '1 0 1.0e-4 0 0\n' > \
        "$SANITY_TMP/$variant/albert_coils/orbit_exit_code.dat"
    printf '&config\ntrace_time = 1.0e-4\n/\n' > \
        "$SANITY_TMP/$variant/albert_coils/simple.in"
    printf '1.0\n' > \
        "$SANITY_TMP/$variant/albert_coils/runtime_seconds.txt"
done
printf '1 1.0e-4 0.5 0.3 0.1 0.3 0 0 1 0\n' > \
    "$SANITY_TMP/reference/albert_coils/times_lost.dat"
cp "$SANITY_TMP/reference/albert_coils/times_lost.dat" \
    "$SANITY_TMP/matching/albert_coils/times_lost.dat"
printf '1 1.0e-4 0.6 0.3 0.1 0.3 0 0 1 0\n' > \
    "$SANITY_TMP/non_matching/albert_coils/times_lost.dat"

echo -n "Test 11: Albert diagnostics and orbit results matching... "
if "$GOLDEN_RECORD_DIR/compare_golden_results.sh" \
        "$SANITY_TMP/reference" "$SANITY_TMP/matching" \
        "albert_coils"; then
    echo "PASSED"
else
    echo "FAILED"
    status=1
fi

echo -n "Test 12: Albert orbit difference is detected... "
if "$GOLDEN_RECORD_DIR/compare_golden_results.sh" \
        "$SANITY_TMP/reference" "$SANITY_TMP/non_matching" \
        "albert_coils"; then
    echo "FAILED (diagnostics hid an orbit difference)"
    status=1
else
    echo "PASSED (correctly detected orbit difference)"
fi

echo ""
echo "All tests completed!"
exit "$status"
