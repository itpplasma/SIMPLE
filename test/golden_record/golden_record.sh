#!/bin/bash
set -euo pipefail

# Golden record comparison against a single reference (default: main)

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT_CUR_DEFAULT=$(cd "$SCRIPT_DIR/../.." && pwd)
REF_SPEC=${1:-main}

REMOTE_URL=$(git -C "$PROJECT_ROOT_CUR_DEFAULT" remote get-url origin 2>/dev/null || echo "git@github.com:itpplasma/SIMPLE.git")
if [[ "$REMOTE_URL" == https://github.com/* ]]; then
    REPO_PREFIX="https://github.com/itpplasma"
else
    REPO_PREFIX="git@github.com:itpplasma"
fi
LIBNEO_URL="${REPO_PREFIX%/}/libneo.git"

GOLDEN_RECORD_BASE_DIR=${GOLDEN_RECORD_BASE_DIR:-"$(pwd)/golden_record"}
PROJECT_ROOT_CUR=${GOLDEN_RECORD_CUR_DIR:-$PROJECT_ROOT_CUR_DEFAULT}
PROJECT_ROOT_REF=${GOLDEN_RECORD_REF_DIR:-"$GOLDEN_RECORD_BASE_DIR/reference_$REF_SPEC"}

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

    local simple_ts
    simple_ts=$(git -C "$project_root" log -1 --format=%ct)

    local tags_dir="$base_dir/libneo_tags"
    if [ ! -d "$tags_dir/.git" ]; then
        git clone --quiet --no-checkout "$LIBNEO_URL" "$tags_dir"
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
        git clone --quiet --branch "$libneo_tag" --depth 1 "$LIBNEO_URL" "$libneo_dir"
    fi

    ln -sfn "$libneo_dir" "$base_dir/libneo"
}

resolve_reference() {
    git -C "$PROJECT_ROOT_CUR_DEFAULT" fetch --tags --quiet origin || true

    if git -C "$PROJECT_ROOT_CUR_DEFAULT" rev-parse --verify --quiet "$REF_SPEC^{commit}" >/dev/null; then
        git -C "$PROJECT_ROOT_CUR_DEFAULT" rev-parse "$REF_SPEC"
        return
    fi

    if git -C "$PROJECT_ROOT_CUR_DEFAULT" rev-parse --verify --quiet "origin/$REF_SPEC^{commit}" >/dev/null; then
        git -C "$PROJECT_ROOT_CUR_DEFAULT" rev-parse "origin/$REF_SPEC"
        return
    fi

    log "Reference '$REF_SPEC' not found locally or on origin"
    exit 1
}

prepare_reference() {
    mkdir -p "$GOLDEN_RECORD_BASE_DIR"

    local resolved_sha
    resolved_sha=$(resolve_reference)

    if [ ! -e "$PROJECT_ROOT_REF/.git" ]; then
        rm -rf "$PROJECT_ROOT_REF"
        log "Creating reference worktree for $REF_SPEC"
        git -C "$PROJECT_ROOT_CUR_DEFAULT" worktree add --force --detach "$PROJECT_ROOT_REF" "$resolved_sha" >/dev/null
    else
        log "Updating reference worktree to $resolved_sha"
        git -C "$PROJECT_ROOT_REF" reset --hard "$resolved_sha" --quiet
    fi

    local build_dir="$PROJECT_ROOT_REF/build"
    local marker_file="$build_dir/.golden_reference_sha"

    if [ -f "$marker_file" ] && [ -f "$build_dir/simple.x" ] && \
            [ "$(cat "$marker_file")" = "$resolved_sha" ]; then
        log "Reference build @$resolved_sha already prepared"
        return
    fi

    log "Building reference SIMPLE at $resolved_sha"
    rm -rf "$build_dir"

    if [ -f "$PROJECT_ROOT_REF/SRC/CMakeLists.txt" ] && \
            grep -q "../libneo" "$PROJECT_ROOT_REF/SRC/CMakeLists.txt" 2>/dev/null; then
        setup_libneo "$PROJECT_ROOT_REF"
    fi

    (
        export LIBNEO_BRANCH="main"

        if ! cmake -S "$PROJECT_ROOT_REF" -B "$build_dir" -G Ninja \
                -DCMAKE_BUILD_TYPE=Release -DENABLE_PYTHON_INTERFACE=OFF -DENABLE_GVEC=OFF \
                > "$PROJECT_ROOT_REF/configure.log" 2>&1; then
            log "CMake configuration failed; see $PROJECT_ROOT_REF/configure.log"
            tail -n 50 "$PROJECT_ROOT_REF/configure.log"
            exit 1
        fi

        if ! cmake --build "$build_dir" --config Release > "$PROJECT_ROOT_REF/build.log" 2>&1; then
            log "Build failed; see $PROJECT_ROOT_REF/build.log"
            tail -n 50 "$PROJECT_ROOT_REF/build.log"
            exit 1
        fi
    ) || exit 1

    echo "$resolved_sha" > "$marker_file"
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
