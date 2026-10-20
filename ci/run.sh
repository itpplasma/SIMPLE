#!/bin/bash
# Run one CI tier inside a bundle made by ci/bundle.sh: ci/run.sh TIER
#
#   check  report the environment and verify the offline mirrors
#   smoke  build + smoke set (libneo reverse-dependency gate)
#   fast   build + fast suite + Python package install
#   full   fast + slow suite, with the regression suite (golden record)
#          built and run concurrently in a second copy of the tree
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

unset OMP_PROC_BIND OMP_PLACES
export OMP_NUM_THREADS=2
export CMAKE_BUILD_PARALLEL_LEVEL=$cpus
make_args=()
[ -n "${SCI_LIBNEO_REF:-}" ] && make_args+=("LIBNEO_REF=$SCI_LIBNEO_REF")

stamp "tier $tier, commit ${SCI_COMMIT:-unknown} on ${SCI_HEAD_REF:-?}, $cpus CPUs, $(hostname -s)"

python_env() {
    stamp "python venv + requirements"
    python3 -m venv "$root/venv" && . "$root/venv/bin/activate" &&
        python -m pip install -q --upgrade pip &&
        python -m pip install -q -r requirements.txt
}

fast_suite() {
    export CTEST_PARALLEL_LEVEL=$1
    stamp "build-deterministic" && make build-deterministic "${make_args[@]}" &&
        stamp "test-fast (ctest -j$CTEST_PARALLEL_LEVEL)" && make test-fast &&
        stamp "pip install -e .[imas]" &&
        python -m pip install -q -e '.[imas]' --no-build-isolation \
            -Ccmake.define.SIMPLE_DETERMINISTIC_FP=ON
}

regression_suite() {
    local jobs=$1 threads=$2 golden
    golden=$(cmake -DDEPENDENCY=libneo_golden -P cmake/read_dependency_version.cmake)
    export CTEST_PARALLEL_LEVEL=$jobs GOLDEN_RECORD_OMP_THREADS=$threads
    export GOLDEN_RECORD_PERF_GATE=0 GOLDEN_LIBNEO_REF=$golden
    stamp "build-deterministic-nopy (libneo $golden)" &&
        make build-deterministic-nopy FLAGS="-DLIBNEO_REF=$golden" || return
    # Start the longest serial test first; ctest orders by recorded cost.
    mkdir -p build/test/Testing/Temporary
    printf '%s\n' 'test_fo_gauss2_resolution_large_matched_policy 1 1200' \
        'golden_record_albert_coils 1 420' 'golden_record_meiss_coils 1 210' \
        > build/test/Testing/Temporary/CTestCostData.txt
    stamp "regression (ctest -j$jobs, golden $threads threads)"
    (cd build && ctest --test-dir test --output-on-failure -L regression)
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
    fast_suite "$(clamp $((cpus / 2)) 2 16)"
    ;;
full)
    python_env || exit
    # Enable the unchanged collisional wall test in the merge-queue slow set.
    # CGAL's release archive is supplied offline by ci/bundle.sh.
    mkdir -p "$root/cgal" &&
        tar -xf "$root/data/CGAL-6.0.1.tar.xz" -C "$root/cgal" --strip-components=1 || exit
    python -m pip install -q 'numpy==2.2.6' 'scipy==1.16.3' vtk \
        'git+https://github.com/itpplasma/map2disc.git@14da2f094ab2505a639b9cd04af57f5c0dbfaa54' || exit
    make_args+=("FLAGS=-DSIMPLE_ENABLE_CGAL=ON -DFETCHCONTENT_SOURCE_DIR_CGAL=$root/cgal")
    # Second tree for the regression build (it needs its own build/ and a
    # different libneo); started first because it holds the critical path.
    cp -a "$src" "$root/src-reg"
    ja=$(clamp $((cpus / 4)) 2 16)
    jb=$(clamp $((cpus / 12)) 2 8)
    tb=$(clamp $(((cpus / 2) / jb)) 2 16)
    (cd "$root/src-reg" && regression_suite "$jb" "$tb") > "$root/regression.log" 2>&1 &
    reg_pid=$!
    stamp "regression started in background (log follows at the end)"
    fast_suite "$ja" && stamp "test-slow (ctest -j$ja)" &&
        (cd build && ctest --test-dir test -N -L slow -R '^test_collisional_stl_wall$' |
            grep 'Test.*test_collisional_stl_wall$') && make test-slow
    rc_a=$?
    stamp "fast/slow finished: exit $rc_a; waiting for regression"
    wait "$reg_pid"; rc_b=$?
    stamp "regression log"
    cat "$root/regression.log"
    stamp "regression finished: exit $rc_b"
    [ "$rc_a" -eq 0 ] && [ "$rc_b" -eq 0 ]
    ;;
*)
    echo "unknown tier: $tier" >&2
    exit 2
    ;;
esac
rc=$?
stamp "tier $tier exit $rc"
exit "$rc"
