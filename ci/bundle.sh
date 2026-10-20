#!/bin/bash
# Pack an offline CI bundle of the checked-out commit: ci/bundle.sh OUT.tar.zst
#
# The bundle holds everything ci/run.sh needs without network access to the
# code forges: src/ (this commit, checked out on SCI_HEAD_REF so dependency
# branch matching works), bare git mirrors of SIMPLE (main + this commit, for
# the golden-record reference) and of every FetchContent dependency, and the
# downloaded test data. Python packages still come from PyPI.
#
# Environment:
#   SCI_HEAD_REF    branch name to check out on (default: current branch)
#   SCI_LIBNEO_REF  libneo ref for the build (default: branch matching)
#   SCI_TIER        full bundles include the CGAL headers (default: full)
#   SCI_CACHE       mirror and data cache (default: ~/.cache/sci-bundle)
set -euo pipefail

out=$(realpath -m "${1:?usage: ci/bundle.sh OUT.tar.zst}")
src=$(git rev-parse --show-toplevel)
cache=${SCI_CACHE:-$HOME/.cache/sci-bundle}
head_ref=${SCI_HEAD_REF:-$(git -C "$src" rev-parse --abbrev-ref HEAD)}
[ "$head_ref" = HEAD ] && head_ref=ci-head
git check-ref-format --branch "$head_ref" > /dev/null

deps="itpplasma/libneo itpplasma/map2disc lazy-fortran/fortnum lazy-fortran/fortplot lazy-fortran/fortio"
if [ "${SCI_TIER:-full}" = full ]; then
    deps="$deps dbstein/pybie2d dbstein/function_generator"
fi

# url<TAB>file name in test/test_data (shared with the golden-record runner).
data_urls="https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc	wout.nc
https://raw.githubusercontent.com/hiddenSymmetries/simsopt/master/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc	wout_qh.nc
https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_c09r00_fixedBoundary_0.5T_vacuum_ns201.nc	wout_ncsx.nc
https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00	coils.c09r00"

stage=$(mktemp -d "${TMPDIR:-/tmp}/sci-bundle.XXXXXX")
trap 'rm -rf "$stage"' EXIT
mkdir -p "$cache/mirrors" "$cache/data" "$stage/mirrors" "$stage/data"

# Bare mirror of branches and tags only (no pull-request refs).
update_mirror() {
    local slug=$1 dir=$cache/mirrors/$1.git
    if [ ! -d "$dir" ]; then
        mkdir -p "$(dirname "$dir")"
        git init -q --bare "$dir"
        git -C "$dir" remote add origin "https://github.com/$slug.git"
        git -C "$dir" config remote.origin.fetch '+refs/heads/*:refs/heads/*'
        git -C "$dir" config --add remote.origin.fetch '+refs/tags/*:refs/tags/*'
    fi
    git -C "$dir" fetch -q --prune origin
    # map2disc's two Git dependencies use master; our libraries use main.
    case $slug in
        dbstein/*) git -C "$dir" symbolic-ref HEAD refs/heads/master ;;
        *) git -C "$dir" symbolic-ref HEAD refs/heads/main ;;
    esac
}

for slug in $deps; do
    update_mirror "$slug"
    mkdir -p "$stage/mirrors/${slug%/*}"
    cp -a "$cache/mirrors/$slug.git" "$stage/mirrors/$slug.git"
done

# SIMPLE: main (golden-record reference), this commit as $head_ref, tags.
git -C "$src" fetch -q origin "+refs/heads/main:refs/remotes/origin/main"
simple=$stage/mirrors/itpplasma/SIMPLE.git
git init -q --bare "$simple"
git -C "$src" push -q "$simple" "HEAD:refs/heads/$head_ref" \
    "refs/remotes/origin/main:refs/heads/main" "+refs/tags/*:refs/tags/*"
git -C "$simple" symbolic-ref HEAD refs/heads/main
git clone -q --branch "$head_ref" "$simple" "$stage/src"
git -C "$stage/src" remote set-url origin https://github.com/itpplasma/SIMPLE.git

while IFS=$'\t' read -r url name; do
    key=$(printf '%s' "$url" | sha256sum | cut -c1-16)
    if [ ! -s "$cache/data/$key" ]; then
        curl -fsSL --retry 3 -o "$cache/data/$key.part" "$url"
        mv "$cache/data/$key.part" "$cache/data/$key"
    fi
    cp "$cache/data/$key" "$stage/data/$name"
    ln "$stage/data/$name" "$stage/src/test/test_data/$name"
    printf '%s\t%s\n' "$url" "$name" >> "$stage/data/urls.tsv"
done <<< "$data_urls"

# Compute nodes cannot fetch release archives. Keep CGAL out of fast bundles.
if [ "${SCI_TIER:-full}" = full ]; then
    cgal=CGAL-6.0.1.tar.xz
    if [ ! -s "$cache/data/$cgal" ]; then
        curl -fsSL --retry 3 -o "$cache/data/$cgal.part" \
            "https://github.com/CGAL/cgal/releases/download/v6.0.1/$cgal"
        mv "$cache/data/$cgal.part" "$cache/data/$cgal"
    fi
    cp "$cache/data/$cgal" "$stage/data/$cgal"
fi

{
    printf 'SCI_COMMIT=%s\n' "$(git -C "$src" rev-parse HEAD)"
    printf 'SCI_HEAD_REF=%s\n' "$head_ref"
    printf 'SCI_LIBNEO_REF=%s\n' "${SCI_LIBNEO_REF:-}"
} > "$stage/bundle.env"

tar -C "$stage" -cf - . | zstd -q -1 -T0 -o "$out" -f
echo "bundle: $out ($(du -h "$out" | cut -f1), commit $(git -C "$src" rev-parse --short HEAD) on $head_ref)"
