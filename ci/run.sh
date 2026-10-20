#!/bin/bash
# Run one CI tier inside a bundle made by ci/bundle.sh: ci/run.sh TIER
#
#   check  report the environment and verify the offline mirrors
#   smoke  build + smoke set (libneo reverse-dependency gate)
#   fast   build + fast suite + Python package install
#   full   fast + slow suite, with the regression suite (golden record)
#          built and run concurrently in a second copy of the tree; a third
#          tree reports fast-only coverage without deciding the gate status
#
# Expects to run from the bundle's src/ directory. SCI_CPUS (default nproc)
# sets the parallelism; OMP_PROC_BIND/OMP_PLACES stay unset so concurrent
# tests do not pin their threads onto the same cores.
set -uo pipefail

tier=${1:?usage: ci/run.sh check|smoke|fast|full}
src=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
root=$(dirname "$src")
cpus=${SCI_CPUS:-$(nproc)}
[ -f "$root/bundle.env" ] && . "$root/bundle.env"
cd "$src" || exit 2

clamp() { local v=$1; [ "$v" -lt "$2" ] && v=$2; [ "$v" -gt "$3" ] && v=$3; echo "$v"; }
stamp() { echo "== $(date -u +%T) $*"; }

# --- offline sources -----------------------------------------------------------
# Redirect every forge URL the build, tests and pip use to the bundled mirrors.
export GIT_CONFIG_GLOBAL=$root/gitconfig
: > "$GIT_CONFIG_GLOBAL"
git config --global uploadpack.allowFilter true
git config --global uploadpack.allowAnySHA1InWant true
git config --global advice.detachedHead false
for m in "$root"/mirrors/*/*.git; do
    slug=${m#"$root"/mirrors/}
    slug=${slug%.git}
    git config --global --add "url.file://$m.insteadOf" "https://github.com/$slug.git"
    git config --global --add "url.file://$m.insteadOf" "https://github.com/$slug"
    git config --global --add "url.file://$m.insteadOf" "git@github.com:$slug.git"
done

# wget is used by the golden-record runner for its test data; serve it from
# the bundle and refuse anything else.
mkdir -p "$root/bin"
cat > "$root/bin/wget" <<EOF
#!/bin/bash
out=""; url=""
while [ \$# -gt 0 ]; do
    case \$1 in -O) out=\$2; shift ;; -*) ;; *) url=\$1 ;; esac
    shift
done
name=\$(awk -F'\t' -v u="\$url" '\$1 == u { print \$2 }' "$root/data/urls.tsv")
[ -n "\$name" ] && [ -n "\$out" ] || { echo "offline wget: no bundled copy of \$url" >&2; exit 1; }
cp "$root/data/\$name" "\$out"
EOF
chmod +x "$root/bin/wget"
export PATH=$root/bin:$PATH
export PIP_NO_INDEX=1 PIP_FIND_LINKS=$root/pip
export MPLBACKEND=Agg
unset DISPLAY

unset OMP_PROC_BIND OMP_PLACES
export OMP_NUM_THREADS=2
export CMAKE_BUILD_PARALLEL_LEVEL=$cpus
make_args=()
pip_args=()
[ -n "${SCI_LIBNEO_REF:-}" ] && make_args+=("LIBNEO_REF=$SCI_LIBNEO_REF")

stamp "tier $tier, commit ${SCI_COMMIT:-unknown} on ${SCI_HEAD_REF:-?}, $cpus CPUs, $(hostname -s)"

python_env() {
    stamp "python venv + requirements"
    local env_dir=${1:-$root/venv}
    python3 -m venv "$env_dir" && . "$env_dir/bin/activate" &&
        python -m pip install -q --upgrade pip &&
        python -m pip install -q -r requirements.txt
}

fast_tests() {
    export CTEST_PARALLEL_LEVEL=$1
    stamp "build-deterministic" && make build-deterministic "${make_args[@]}" &&
        stamp "test-fast (ctest -j$CTEST_PARALLEL_LEVEL)" && make test-fast "${make_args[@]}"
}

install_package() {
    stamp "pip install -e .[imas]" &&
        python -m pip install -q -e '.[imas]' --no-build-isolation \
            -Ccmake.define.SIMPLE_DETERMINISTIC_FP=ON "${pip_args[@]}"
}

regression_suite() {
    local jobs=$1 threads=$2 golden
    golden=$(cmake -DDEPENDENCY=libneo_golden -P cmake/read_dependency_version.cmake)
    export CTEST_PARALLEL_LEVEL=$jobs GOLDEN_RECORD_OMP_THREADS=$threads
    export GOLDEN_RECORD_PERF_GATE=0 GOLDEN_LIBNEO_REF=$golden
    # The golden reference runner uses Release; match its architecture flags.
    stamp "build-deterministic-nopy (libneo $golden)" &&
        make build-deterministic-nopy CONFIG=Release \
            FLAGS="-DLIBNEO_REF=$golden -DENABLE_COVERAGE=OFF -DSIMPLE_TEST_MEMORY_GB=$reg_mem" || return
    # Start the longest serial test first; ctest orders by recorded cost.
    mkdir -p build/test/Testing/Temporary
    printf '%s\n' 'test_fo_gauss2_resolution_large_matched_policy 1 1200' \
        'golden_record_albert_coils 1 420' 'golden_record_meiss_coils 1 210' \
        > build/test/Testing/Temporary/CTestCostData.txt
    stamp "regression (ctest -j$jobs, golden $threads threads)"
    (cd build && ctest --test-dir test --resource-spec-file "$PWD/ctest-resources.json" --output-on-failure -L regression)
}

coverage_report() {
    local bytes
    local search_paths=(build)
    [ ! -d _skbuild ] || search_paths+=(_skbuild)
    # GCC 12 cannot resolve the concept-only target's pathless generated
    # includes; that target is outside the production src/app filter.
    gcovr --root . --filter 'src/|app/' --include-internal-functions \
        --gcov-executable /usr/bin/gcov -j "$cov_cpus" \
        --gcov-exclude-directory '.*/CMakeFiles/fo_gauss2_concept.dir$' \
        --json coverage.json "${search_paths[@]}" || return
    gcovr --root . --filter 'src/|app/' --include-internal-functions \
        --add-tracefile coverage.json \
        --xml-pretty --xml "$result/coverage.xml" \
        --html-details "$result/coverage.html" --print-summary > "$result/coverage-summary.txt" || return
    # Leave headroom for the status files added after reporting. If detailed
    # HTML exceeds the result budget, retain the Pages-compatible overview.
    bytes=$( (cd "$result" && tar -cf - coverage*) | zstd -q -T1 | wc -c) || return
    if [ "$bytes" -gt $((63 * 1024 * 1024)) ]; then
        rm -f "$result"/coverage*.html "$result"/coverage*.css
        gcovr --root . --filter 'src/|app/' --include-internal-functions \
            --add-tracefile coverage.json --html "$result/coverage.html" || return
        bytes=$( (cd "$result" && tar -cf - coverage*) | zstd -q -T1 | wc -c) || return
    fi
    stamp "coverage result: $bytes bytes (zstd), transport limit 67108864"
    cat "$result/coverage-summary.txt"
    [ "$bytes" -le $((63 * 1024 * 1024)) ] || return 4
}

coverage_suite() {
    cd "$root/src-coverage" || exit
    export CMAKE_BUILD_PARALLEL_LEVEL=$cov_cpus
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
    # Separate environments keep editable Python installs from crossing builds.
    python_env "$root/venv-coverage" &&
        python -m pip install -q -r ci/coverage-requirements.txt || exit
    make_args+=("CONFIG=Fast" "FLAGS=-DENABLE_COVERAGE=ON -DSIMPLE_TEST_MEMORY_GB=$cov_mem")
    pip_args+=("-Ccmake.define.ENABLE_COVERAGE=ON" "-Ccmake.build-type=Fast")
    fast_tests "$cov_cpus"; rc_tests=$?
    install_package; rc_package=$?
    coverage_report; rc_report=$?
    if [ "$rc_report" -ne 0 ]; then
        # Never let partial/oversized coverage make the transport fail the gate.
        rm -f "$result"/coverage.*
    fi
    printf 'Coverage fast tests: exit %s\nCoverage Python install: exit %s\nCoverage gcovr: exit %s\n' \
        "$rc_tests" "$rc_package" "$rc_report" > "$result/coverage-status.txt"
    [ "$rc_tests" -eq 0 ] || exit "$rc_tests"
    [ "$rc_package" -eq 0 ] || exit "$rc_package"
    exit "$rc_report"
}

correctness_suite() {
    export CMAKE_BUILD_PARALLEL_LEVEL=$gate_cpus
    python_env || exit
    python -m pip install -q 'numpy==2.2.6' 'scipy==1.16.3' vtk \
        'git+https://github.com/itpplasma/map2disc.git@14da2f094ab2505a639b9cd04af57f5c0dbfaa54' || exit
    make_args+=("CONFIG=Fast" "FLAGS=-DENABLE_COVERAGE=OFF -DSIMPLE_ENABLE_CGAL=ON -DFETCHCONTENT_SOURCE_DIR_CGAL=$root/cgal -DSIMPLE_TEST_MEMORY_GB=$fast_mem")
    pip_args+=("-Ccmake.define.ENABLE_COVERAGE=OFF" "-Ccmake.build-type=Fast")
    # Regression retains Release to match the golden reference's build contract.
    # Both gate trees and coverage get disjoint CPU/memory scheduling budgets.
    export CMAKE_BUILD_PARALLEL_LEVEL=$fast_cpus
    (cd "$root/src-reg" && CMAKE_BUILD_PARALLEL_LEVEL=$reg_cpus regression_suite "$reg_cpus" 2) \
        > "$root/regression.log" 2>&1 &
    reg_pid=$!
    stamp "regression started in background (log follows at the end)"
    fast_tests "$fast_cpus"; rc_fast=$?
    install_package && stamp "test-slow (ctest -j$fast_cpus)" &&
        (cd build && ctest --test-dir test -N -L slow -R '^test_collisional_stl_wall$' |
            grep 'Test.*test_collisional_stl_wall$') && make test-slow "${make_args[@]}"
    rc_slow=$?
    stamp "fast/slow finished: exit $rc_fast/$rc_slow; waiting for regression"
    wait "$reg_pid"; rc_reg=$?
    stamp "regression log"
    cat "$root/regression.log"
    printf 'Gate fast tests: exit %s\nGate slow tests / Python install: exit %s\nGate regression: exit %s\n' \
        "$rc_fast" "$rc_slow" "$rc_reg" > "$result/gate-status.txt"
    [ "$rc_fast" -eq 0 ] || exit "$rc_fast"
    [ "$rc_slow" -eq 0 ] || exit "$rc_slow"
    exit "$rc_reg"
}

case $tier in
check)
    git describe --tags --always --dirty
    for m in "$root"/mirrors/*/*.git; do
        printf '%s %s\n' "${m#"$root"/mirrors/}" "$(git -C "$m" rev-parse --short HEAD)"
    done
    git ls-remote --heads https://github.com/itpplasma/libneo.git main
    nproc; free -g | head -2
    ;;
smoke)
    python_env || exit
    CTEST_PARALLEL_LEVEL=$(clamp $((cpus / 2)) 2 16)
    export CTEST_PARALLEL_LEVEL
    stamp "build + smoke" && make test-smoke "${make_args[@]}"
    ;;
fast)
    python_env || exit
    fast_tests "$(clamp $((cpus / 2)) 2 16)" && install_package
    ;;
full)
    # The full profile is unchanged: one allocation, with coverage reporting
    # alongside the authoritative uninstrumented correctness gate.
    [ "$cpus" -ge 3 ] || { echo "full needs at least 3 CPUs" >&2; exit 2; }
    cov_cpus=$(clamp $((cpus / 4)) 1 16)
    gate_cpus=$((cpus - cov_cpus))
    reg_cpus=$((gate_cpus / 2))
    fast_cpus=$((gate_cpus - reg_cpus))
    # Reserve 25% of the allocation for builds, OS and unreserved processes.
    test_mem=$(( ${SCI_MEM_GB:-192} * 3 / 4 ))
    cov_mem=$((test_mem * cov_cpus / cpus))
    reg_mem=$(((test_mem - cov_mem) / 2))
    fast_mem=$((test_mem - cov_mem - reg_mem))
    result=$src/ci-out
    mkdir -p "$result" "$root/cgal" &&
        tar -xf "$root/data/CGAL-6.0.1.tar.xz" -C "$root/cgal" --strip-components=1 || exit
    # Copy before building so counters, generated files and Python installs
    # remain isolated. Every CTest uses its generated Slurm/cgroup-aware spec.
    cp -a "$src" "$root/src-reg" && cp -a "$src" "$root/src-coverage" || exit
    stamp "split: gate $gate_cpus CPUs (fast/slow $fast_cpus, regression $reg_cpus), coverage $cov_cpus CPUs"
    stamp "CTest memory: fast/slow $fast_mem GiB, regression $reg_mem GiB, coverage $cov_mem GiB"
    (coverage_suite) > "$root/coverage.log" 2>&1 &
    cov_pid=$!
    (correctness_suite); rc_gate=$?
    stamp "correctness gate finished: exit $rc_gate; waiting for coverage"
    wait "$cov_pid"; rc_cov=$?
    stamp "coverage log (reporting only)"
    cat "$root/coverage.log"
    printf 'Correctness gate: exit %s (authoritative)\nCoverage worker: exit %s (reporting only)\n' \
        "$rc_gate" "$rc_cov" > "$result/full-status.txt"
    cat "$result/full-status.txt"
    if [ "$rc_cov" -ne 0 ]; then
        stamp "WARNING: coverage failed; correctness gate alone decides the check"
    fi
    (exit "$rc_gate")
    ;;
*)
    echo "unknown tier: $tier" >&2
    exit 2
    ;;
esac
rc=$?
stamp "tier $tier exit $rc"
exit "$rc"
