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

GOLDEN_RECORD_RTOL="${GOLDEN_RECORD_RTOL:-1e-7}"
GOLDEN_RECORD_ATOL="${GOLDEN_RECORD_ATOL:-1e-12}"
GOLDEN_RECORD_MAX_SLOWDOWN="${GOLDEN_RECORD_MAX_SLOWDOWN:-1.2}"
GOLDEN_RECORD_MIN_REF_RUNTIME_S="${GOLDEN_RECORD_MIN_REF_RUNTIME_S:-1.0}"
GOLDEN_RECORD_MAX_ABS_SLOWDOWN_S="${GOLDEN_RECORD_MAX_ABS_SLOWDOWN_S:-0.2}"

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

# GOLDEN_RECORD_SKIP_CASES: space-separated list of case names whose ref/cur
# divergence is expected and accepted (e.g. an intentional bug fix changes
# the output). Each entry must be the exact case directory name. Set it in
# CI workflow or the calling Makefile target rather than here, so the skip
# stays explicit per branch.
SKIP_CASES="${GOLDEN_RECORD_SKIP_CASES:-}"

case_is_skipped() {
    local case_name="$1"
    [ -z "$SKIP_CASES" ] && return 1
    for s in $SKIP_CASES; do
        [ "$s" = "$case_name" ] && return 0
    done
    return 1
}

compare_cases() {
    local total_cases=0
    local passed_cases=0
    local failed_cases=0
    local missing_cases=0
    local skipped_cases=0
    local perf_failed_cases=0

    for CASE in $TEST_CASES; do
        total_cases=$((total_cases + 1))
        recovery_transition=0

        if case_is_skipped "$CASE"; then
            echo "Comparing $CASE case..."
            echo "  ○ SKIPPED (in GOLDEN_RECORD_SKIP_CASES)"
            skipped_cases=$((skipped_cases + 1))
            echo ""
            continue
        fi

        echo "Comparing $CASE case..."

        # Check if reference has this test case (new tests may not exist in reference)
        if [ ! -d "$REFERENCE_DIR/$CASE" ]; then
            echo "  ○ SKIPPED (not present in reference version)"
            skipped_cases=$((skipped_cases + 1))
            echo ""
            continue
        fi

        # Check if this is albert_coils case with diagnostic file
        if [ "$CASE" = "albert_coils" ]; then
            # Albert coordinates with coils field: compare diagnostic files
            # With deterministic FP builds, results should be bit-identical
            REF_DIAG="$REFERENCE_DIR/$CASE/albert_coils_diagnostic.dat"
            CUR_DIAG="$CURRENT_DIR/$CASE/albert_coils_diagnostic.dat"
            REF_INTER="$REFERENCE_DIR/$CASE/albert_intermediate.dat"
            CUR_INTER="$CURRENT_DIR/$CASE/albert_intermediate.dat"

            if [ -f "$REF_DIAG" ] && [ -f "$CUR_DIAG" ]; then
                echo "  (comparing Albert coordinate diagnostic)"
                python "$SCRIPT_DIR/compare_files_multi.py" \
                    "$REFERENCE_DIR/$CASE" "$CURRENT_DIR/$CASE" \
                    --files albert_coils_diagnostic.dat \
                    --rtol "$GOLDEN_RECORD_RTOL" --atol "$GOLDEN_RECORD_ATOL"
                result=$?
                if [ -f "$REF_INTER" ] && [ -f "$CUR_INTER" ]; then
                    echo "  (also comparing intermediate values)"
                    python "$SCRIPT_DIR/compare_files_multi.py" \
                        "$REFERENCE_DIR/$CASE" "$CURRENT_DIR/$CASE" \
                        --files albert_intermediate.dat \
                        --rtol "$GOLDEN_RECORD_RTOL" --atol "$GOLDEN_RECORD_ATOL"
                    inter_result=$?
                    if [ $inter_result -ne 0 ]; then
                        echo "  FAILED: Intermediate values differ"
                        result=$inter_result
                    fi
                fi
            elif [ ! -f "$REF_DIAG" ] && [ -f "$CUR_DIAG" ]; then
                echo "  (reference lacks diagnostic - SKIPPING albert_coils)"
                result=0
            else
                echo "  (no diagnostic, comparing times_lost.dat)"
                python "$SCRIPT_DIR/compare_files_multi.py" \
                    "$REFERENCE_DIR/$CASE" "$CURRENT_DIR/$CASE" \
                    --files times_lost.dat \
                    --rtol "$GOLDEN_RECORD_RTOL" --atol "$GOLDEN_RECORD_ATOL"
                result=$?
            fi

            # Diagnostics used to overwrite these files, so this case only
            # compared field diagnostics and silently ignored orbit changes.
            # They now run in an isolated directory; enforce the orbit result
            # in addition to the diagnostic contract.
            REF_EXIT="$REFERENCE_DIR/$CASE/orbit_exit_code.dat"
            CUR_EXIT="$CURRENT_DIR/$CASE/orbit_exit_code.dat"
            echo "  (also comparing orbit results)"
            if [ -f "$REF_EXIT" ] && [ -f "$CUR_EXIT" ]; then
                python "$SCRIPT_DIR/compare_orbit_results.py" \
                    "$REFERENCE_DIR/$CASE" "$CURRENT_DIR/$CASE" \
                    --rtol "$GOLDEN_RECORD_RTOL" --atol "$GOLDEN_RECORD_ATOL"
                orbit_result=$?
                if [ $orbit_result -eq 3 ]; then
                    recovery_transition=1
                    orbit_result=0
                fi
            else
                GOLDEN_RECORD_RTOL="$GOLDEN_RECORD_RTOL" \
                    GOLDEN_RECORD_ATOL="$GOLDEN_RECORD_ATOL" \
                    python "$SCRIPT_DIR/compare_files.py" \
                    "$REFERENCE_DIR/$CASE/times_lost.dat" \
                    "$CURRENT_DIR/$CASE/times_lost.dat"
                orbit_result=$?
            fi
            if [ $orbit_result -ne 0 ]; then
                echo "  FAILED: Orbit results differ"
                result=$orbit_result
            fi

            if [ $result -eq 0 ]; then
                echo "  ✓ PASSED"
                passed_cases=$((passed_cases + 1))
            else
                echo "  ✗ FAILED"
                failed_cases=$((failed_cases + 1))
            fi
        # Check if this is a classifier case with multiple files
        elif [ "$CASE" = "classifier_fast" ] || [ "$CASE" = "classifier_combined" ]; then
            # List of files to compare for the classifier cases (excluding
            # simple.in and wout.nc). fort.* are excluded because of
            # non-deterministic ordering in parallel execution.
            # avg_inverse_t_lost.dat is excluded because the file is only
            # written when at least one particle is actually lost; with the
            # fast_class regression fix the small classifier_fast test loses
            # zero particles, so neither ref nor cur write the file and the
            # comparator's "missing" check spuriously fails.
            # healaxis.dat is excluded because the axis-healing diagnostic is no
            # longer written (dropped in the near-axis healing work / #331), so it
            # is always missing and the comparator spuriously fails.
            CLASSIFIER_FILES="class_parts.dat confined_fraction.dat start.dat times_lost.dat"
            
            # Run multi-file comparison
            python "$SCRIPT_DIR/compare_files_multi.py" \
                "$REFERENCE_DIR/$CASE" "$CURRENT_DIR/$CASE" \
                --files $CLASSIFIER_FILES \
                --rtol "$GOLDEN_RECORD_RTOL" --atol "$GOLDEN_RECORD_ATOL"
            
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
            
            REF_EXIT="$REFERENCE_DIR/$CASE/orbit_exit_code.dat"
            CUR_EXIT="$CURRENT_DIR/$CASE/orbit_exit_code.dat"
            if [ -f "$REF_EXIT" ] && [ -f "$CUR_EXIT" ]; then
                python "$SCRIPT_DIR/compare_orbit_results.py" \
                    "$REFERENCE_DIR/$CASE" "$CURRENT_DIR/$CASE" \
                    --rtol "$GOLDEN_RECORD_RTOL" --atol "$GOLDEN_RECORD_ATOL"
                result=$?
                if [ $result -eq 3 ]; then
                    recovery_transition=1
                    result=0
                fi
            else
                GOLDEN_RECORD_RTOL="$GOLDEN_RECORD_RTOL" \
                    GOLDEN_RECORD_ATOL="$GOLDEN_RECORD_ATOL" \
                    python "$SCRIPT_DIR/compare_files.py" "$REF_FILE" "$CUR_FILE"
                result=$?
            fi

            if [ $result -eq 0 ]; then
                echo "  ✓ PASSED"
                passed_cases=$((passed_cases + 1))
            else
                echo "  ✗ FAILED"
                failed_cases=$((failed_cases + 1))
            fi
        fi

        # Performance comparison (guard against >20% slowdown by default)
        ref_time_file="$REFERENCE_DIR/$CASE/runtime_seconds.txt"
        cur_time_file="$CURRENT_DIR/$CASE/runtime_seconds.txt"
        if [ -f "$ref_time_file" ] && [ -f "$cur_time_file" ]; then
            python - <<PY
import os
ref = float(open("${ref_time_file}", "r", encoding="utf-8").read().strip())
cur = float(open("${cur_time_file}", "r", encoding="utf-8").read().strip())
max_slow = float("${GOLDEN_RECORD_MAX_SLOWDOWN}")
min_ref = float("${GOLDEN_RECORD_MIN_REF_RUNTIME_S}")
max_abs = float("${GOLDEN_RECORD_MAX_ABS_SLOWDOWN_S}")
recovery_transition = int("${recovery_transition}")
ratio = cur / ref if ref > 0.0 else float("inf")
delta = (ratio - 1.0) * 100.0
abs_delta = cur - ref
print(f"  perf: ref={ref:.3f}s cur={cur:.3f}s ratio={ratio:.3f} (delta={delta:+.1f}%)")
if recovery_transition:
    import numpy as np
    ref_exit = np.loadtxt("${REFERENCE_DIR}/${CASE}/orbit_exit_code.dat")
    cur_exit = np.loadtxt("${CURRENT_DIR}/${CASE}/orbit_exit_code.dat")
    if ref_exit.ndim == 1:
        ref_exit = ref_exit.reshape(1, -1)
    if cur_exit.ndim == 1:
        cur_exit = cur_exit.reshape(1, -1)
    ref_times = np.loadtxt("${REFERENCE_DIR}/${CASE}/times_lost.dat")
    cur_times = np.loadtxt("${CURRENT_DIR}/${CASE}/times_lost.dat")
    if ref_times.ndim == 1:
        ref_times = ref_times.reshape(1, -1)
    if cur_times.ndim == 1:
        cur_times = cur_times.reshape(1, -1)
    ref_traced = np.isin(ref_exit[:, 1], (0, 1, 2)) & np.isfinite(ref_times[:, 1])
    cur_traced = np.isin(cur_exit[:, 1], (0, 1, 2)) & np.isfinite(cur_times[:, 1])
    ref_marker_time = float(np.sum(np.maximum(ref_times[ref_traced, 1], 0.0)))
    cur_marker_time = float(np.sum(np.maximum(cur_times[cur_traced, 1], 0.0)))
    marker_time_scale = cur_marker_time / max(ref_marker_time, 1.0e-300)
    ref_resolved = int(np.count_nonzero(np.isin(ref_exit[:, 1], (0, 1, 2, 3))))
    cur_resolved = int(np.count_nonzero(np.isin(cur_exit[:, 1], (0, 1, 2, 3))))
    resolved_scale = cur_resolved / max(ref_resolved, 1)
    workload_scale = min(marker_time_scale, resolved_scale)
    max_slow *= max(workload_scale, 1.0)
    print(
        "  perf: recovery-adjusted limit="
        f"{max_slow:.3f} (work scale={workload_scale:.3f}; "
        f"marker-time={marker_time_scale:.3f}, resolved={resolved_scale:.3f})"
    )
if ref < min_ref:
    print(f"  perf: using absolute guard (ref<{min_ref:.3f}s): delta={abs_delta:+.3f}s (limit={max_abs:.3f}s)")
    if abs_delta > max_abs:
        raise SystemExit(1)
else:
    if ratio > max_slow:
        raise SystemExit(1)
PY
            perf_rc=$?
            if [ $perf_rc -ne 0 ]; then
                if [ -n "$ref_time_file" ]; then
                    echo "  ✗ FAILED (performance regression)"
                else
                    echo "  ✗ FAILED (performance regression: ratio > ${GOLDEN_RECORD_MAX_SLOWDOWN})"
                fi
                perf_failed_cases=$((perf_failed_cases + 1))
                failed_cases=$((failed_cases + 1))
            fi
        else
            echo "  ⚠ Missing runtime_seconds.txt for performance comparison"
            missing_cases=$((missing_cases + 1))
        fi

        echo ""
    done
    
    # Summary
    echo "========================================="
    echo "Summary:"
    echo "  Total test cases: $total_cases"
    echo "  Passed: $passed_cases"
    echo "  Failed: $failed_cases"
    echo "  Skipped (new tests): $skipped_cases"
    echo "  Missing files: $missing_cases"
    echo "  Performance regressions: $perf_failed_cases"
    echo "========================================="
    
    # Return non-zero if any tests failed or had missing files
    if [ $failed_cases -gt 0 ] || [ $missing_cases -gt 0 ]; then
        return 1
    else
        return 0
    fi
}

compare_sympl_tokamak() {
    local case_dir="sympl_tokamak"
    local ref_dir="$REFERENCE_DIR/$case_dir"
    local cur_dir="$CURRENT_DIR/$case_dir"

    if [ ! -d "$ref_dir" ] || [ ! -d "$cur_dir" ]; then
        echo "Skipping symplectic tokamak comparison (missing outputs)."
        return 0
    fi

    echo "Comparing symplectic tokamak outputs..."
    local files="tokamak_testfield_euler1.dat tokamak_testfield_euler2.dat tokamak_testfield_midpoint.dat tokamak_testfield_gauss4.dat"
    python "$SCRIPT_DIR/compare_files_multi.py" "$ref_dir" "$cur_dir" \
        --files $files --rtol "$GOLDEN_RECORD_RTOL" --atol "$GOLDEN_RECORD_ATOL"
    return $?
}

# Main execution
compare_cases
exit_code=$?

sympl_exit=0
compare_sympl_tokamak || sympl_exit=$?

final_exit=$exit_code
if [ $sympl_exit -ne 0 ]; then
    final_exit=$sympl_exit
fi

if [ $final_exit -eq 0 ]; then
    echo "All tests passed!"
else
    echo "Some tests failed or had missing files."
fi

exit $final_exit
