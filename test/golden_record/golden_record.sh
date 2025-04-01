#!/bin/bash
CLONE_URL="https://github.com/itpplasma/SIMPLE.git"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CUR_VER=$(git -C "$SCRIPT_DIR" describe --tags --always --dirty)
REF_VER=${1:-"main"}
TMP_DIR=${2:-"$(mktemp -d)"}

if [ "$REF_VER" == "$CUR_VER" ]; then
    echo "Reference version and current version are the same. Exiting."
    exit 1
fi

PROJECT_ROOT_CUR=$(cd "$SCRIPT_DIR/../.." && pwd)
RUN_DIR_CUR="$TMP_DIR/run_$CUR_VER"

PROJECT_ROOT_REF="$TMP_DIR/simple_$REF_VER"
RUN_DIR_REF="$TMP_DIR/run_$REF_VER"

TEST_CASES="$(cd "$SCRIPT_DIR" && ls simple*.in)"

echo "Temporary directory: $TMP_DIR"
echo "Reference version: $REF_VER"
echo "Current version: $CUR_VER"
echo "Project root (reference): $PROJECT_ROOT_REF"
echo "Project root (current): $PROJECT_ROOT_CUR"

main() {
    build "$PROJECT_ROOT_CUR"
    clone "$REF_VER" "$PROJECT_ROOT_REF"
    build "$PROJECT_ROOT_REF"
    get_test_data
    run_cases "$PROJECT_ROOT_CUR" "$RUN_DIR_CUR" "$TEST_CASES"
    run_cases "$PROJECT_ROOT_REF" "$RUN_DIR_REF" "$TEST_CASES"
    run_comparisons "$RUN_DIR_REF" "$RUN_DIR_CUR" "$TEST_CASES"
}

get_test_data() {
    cd "$TMP_DIR"
    # Fetch if not exists
    if [ -f "wout.nc" ]; then
        echo "wout.nc already exists. Skipping download."
        return
    fi
    echo "Fetching test data..."
    wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc > /dev/null 2>&1
}

clone() {
    local VERSION="$1"
    local PROJECT_ROOT="$2"
    if [ ! -d "$PROJECT_ROOT" ]; then
        echo "Cloning SIMPLE version $VERSION"
        git clone --filter=blob:none --no-checkout "$CLONE_URL" "$PROJECT_ROOT"
    fi

    cd "$PROJECT_ROOT"
    git fetch --all --quiet
    git checkout "$VERSION" --quiet
}

build() {
    local PROJECT_ROOT="$1"
    echo "Building SIMPLE in $PROJECT_ROOT"
    cd $PROJECT_ROOT
    cmake -S . -Bbuild -GNinja -DCMAKE_BUILD_TYPE=Release -DENABLE_PYTHON_INTERFACE=OFF  > $PROJECT_ROOT/configure.log 2>&1
    cmake --build build --config Release  > $PROJECT_ROOT/build.log 2>&1
}

run_cases() {
    local PROJECT_ROOT="$1"
    local RUN_DIR="$2"
    local TEST_CASES="$3"

    echo "Running test cases in $RUN_DIR"

    for INPUT in $TEST_CASES; do
        CASE="${INPUT%.in}"
        echo "Running case: $CASE"
        mkdir -p "$RUN_DIR/$CASE"
        cp "$SCRIPT_DIR/$INPUT" "$RUN_DIR/$CASE/simple.in"
        cd "$RUN_DIR/$CASE"
        ln -sf "$TMP_DIR/wout.nc" wout.nc
        $PROJECT_ROOT/build/simple.x < simple.in > stdout.log 2>&1
    done

}

run_comparisons() {
    local REF_DIR="$1"
    local CUR_DIR="$2"
    local TEST_CASES="$3"

    for INPUT in $TEST_CASES; do
        CASE="${INPUT%.in}"
        REF_FILE="$REF_DIR/$CASE/times_lost.dat"
        CUR_FILE="$CUR_DIR/$CASE/times_lost.dat"
        python "$SCRIPT_DIR/compare_files.py" "$REF_FILE" "$CUR_FILE"
    done
}

main "$@"
