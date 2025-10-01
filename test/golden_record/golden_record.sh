#!/bin/bash
set -euo pipefail

# Golden record comparison against a single reference (default: main)

CLONE_URL="git@github.com:itpplasma/SIMPLE.git"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT_CUR_DEFAULT=$(cd "$SCRIPT_DIR/../.." && pwd)
REF_SPEC=${1:-main}

GOLDEN_RECORD_BASE_DIR=${GOLDEN_RECORD_BASE_DIR:-"$(pwd)/golden_record"}
PROJECT_ROOT_CUR=${GOLDEN_RECORD_CUR_DIR:-$PROJECT_ROOT_CUR_DEFAULT}
PROJECT_ROOT_REF=${GOLDEN_RECORD_REF_DIR:-"$GOLDEN_RECORD_BASE_DIR/simple_$REF_SPEC"}

RUN_DIR_CUR="$GOLDEN_RECORD_BASE_DIR/runs/run_current"
RUN_DIR_REF="$GOLDEN_RECORD_BASE_DIR/runs/run_reference"
TEST_DATA_DIR="$GOLDEN_RECORD_BASE_DIR/test_data"

log() {
    echo "[golden] $*"
}

setup_libneo() {
    local project_root="$1"
    local base_dir
    base_dir="$(dirname "$project_root")"
    local libneo_url="git@github.com:itpplasma/libneo.git"

    local simple_ts
    simple_ts=$(git -C "$project_root" log -1 --format=%ct)

    local tags_dir="$base_dir/libneo_tags"
    if [ ! -d "$tags_dir/.git" ]; then
        git clone --quiet --no-checkout "$libneo_url" "$tags_dir"
    else
        git -C "$tags_dir" fetch --tags --quiet
    fi

    local libneo_tag
    libneo_tag=$(git -C "$tags_dir" for-each-ref --sort=-creatordate \
        --format '%(refname:short) %(creatordate:unix)' refs/tags | while read -r tag ts; do
            if [ -n "$ts" ] && [ "$ts" -le "$simple_ts" ] 2>/dev/null; then
                echo "$tag"
                break
            fi
        done)

    if [ -z "$libneo_tag" ]; then
        log "No libneo tag older than commit timestamp; falling back to oldest"
        libneo_tag=$(git -C "$tags_dir" for-each-ref --sort=creatordate \
            --format '%(refname:short)' refs/tags | head -n 1)
    fi

    log "Using libneo tag: $libneo_tag"
    local libneo_dir="$base_dir/libneo_$libneo_tag"
    if [ ! -d "$libneo_dir/.git" ]; then
        git clone --quiet --branch "$libneo_tag" --depth 1 "$libneo_url" "$libneo_dir"
    fi

    ln -sfn "$libneo_dir" "$base_dir/libneo"
}

prepare_reference() {
    mkdir -p "$GOLDEN_RECORD_BASE_DIR"
    mkdir -p "$(dirname "$PROJECT_ROOT_REF")"

    if [ ! -d "$PROJECT_ROOT_REF/.git" ]; then
        log "Cloning reference tree for $REF_SPEC"
        git clone --filter=blob:none --no-checkout "$CLONE_URL" "$PROJECT_ROOT_REF"
    fi

    cd "$PROJECT_ROOT_REF" >/dev/null
    git remote set-url origin "$CLONE_URL" >/dev/null 2>&1 || true
    git fetch --all --tags --prune --quiet

    local target_ref
    if [ "$REF_SPEC" = "main" ]; then
        target_ref="origin/main"
        git checkout --force -B main "$target_ref" --quiet
    else
        target_ref="$REF_SPEC"
        git checkout "$target_ref" --quiet
    fi
    git reset --hard "$target_ref" --quiet

    local build_dir="$PROJECT_ROOT_REF/build"
    local marker_file="$build_dir/.golden_reference_sha"
    local current_sha
    current_sha=$(git rev-parse HEAD)

    if [ -f "$marker_file" ] && [ -f "$build_dir/simple.x" ] && \
            [ "$(cat "$marker_file")" = "$current_sha" ]; then
        log "Reference build @${current_sha} already prepared"
        cd "$SCRIPT_DIR"
        return
    fi

    log "Building reference SIMPLE at $current_sha"
    rm -rf "$build_dir"
    git clean -fdx -q

    if [ -f "SRC/CMakeLists.txt" ] && grep -q "../libneo" "SRC/CMakeLists.txt" 2>/dev/null; then
        setup_libneo "$PROJECT_ROOT_REF"
    fi

    cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_PYTHON_INTERFACE=OFF -DENABLE_GVEC=OFF > configure.log 2>&1
    if [ $? -ne 0 ]; then
        log "CMake configuration failed; see $PROJECT_ROOT_REF/configure.log"
        exit 1
    fi

    cmake --build build --config Release > build.log 2>&1
    if [ $? -ne 0 ]; then
        log "Build failed; see $PROJECT_ROOT_REF/build.log"
        exit 1
    fi

    echo "$current_sha" > "$marker_file"
    cd "$SCRIPT_DIR"
}

run_suite() {
    local project_root="$1"
    local run_dir="$2"

    "$SCRIPT_DIR/run_golden_tests.sh" "$project_root" "$run_dir" "$TEST_DATA_DIR"
}

main() {
    mkdir -p "$GOLDEN_RECORD_BASE_DIR" "$TEST_DATA_DIR"

    if [ ! -f "$PROJECT_ROOT_CUR/build/simple.x" ]; then
        log "Current build missing at $PROJECT_ROOT_CUR/build/simple.x"
        exit 1
    fi

    prepare_reference

    local ref_sha
    ref_sha=$(git -C "$PROJECT_ROOT_REF" rev-parse HEAD)
    local cur_sha
    cur_sha=$(git -C "$PROJECT_ROOT_CUR" rev-parse HEAD 2>/dev/null || echo "working-tree")

    log "Reference: $REF_SPEC @ $ref_sha"
    log "Current:   $cur_sha"

    run_suite "$PROJECT_ROOT_REF" "$RUN_DIR_REF"
    run_suite "$PROJECT_ROOT_CUR" "$RUN_DIR_CUR"

    "$SCRIPT_DIR/compare_golden_results.sh" "$RUN_DIR_REF" "$RUN_DIR_CUR"
}

main "$@"
