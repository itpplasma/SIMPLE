#!/bin/bash
# Pack an offline CI bundle of the checked-out commit: ci/bundle.sh OUT.tar.zst
#
# The bundle holds everything ci/run.sh needs without network access to the
# code forges: src/ (this commit, checked out on SCI_HEAD_REF so dependency
# branch matching works), bare git mirrors of SIMPLE (main + this commit, for
# the golden-record reference) and of every FetchContent dependency, and the
# downloaded test data and a pip wheel mirror for the cluster's Python 3.11.
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

# url<TAB>file name in test/test_data<TAB>SHA-256 (shared with the golden-record
# runner and test/tests/CMakeLists.txt, which pin the same files).
data_urls="https://raw.githubusercontent.com/hiddenSymmetries/simsopt/912bdcc729d95433518c71df82b8cc45287908b3/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc	wout.nc	333ab8cef64c4b4f406e76209266fc3bf4d7976e9c28c162983c2120e112e771
https://raw.githubusercontent.com/hiddenSymmetries/simsopt/912bdcc729d95433518c71df82b8cc45287908b3/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc	wout_qh.nc	249c397d67357949b2ab4690145224736d8388f5332f23c67c9cb0bf2aa9931e
https://raw.githubusercontent.com/hiddenSymmetries/simsopt/912bdcc729d95433518c71df82b8cc45287908b3/tests/test_files/wout_c09r00_fixedBoundary_0.5T_vacuum_ns201.nc	wout_ncsx.nc	36926d748099bc7abc138d82bdd4b447d2579f08f50189f27d441cf96c1b51e0
https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00	coils.c09r00	3c429da06f4c062887a497a16e2d2bd10f0ecb0b8858c252698631f3853da428"

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
# A run on main itself (push or dispatch) tests this commit as main; two
# refspecs for refs/heads/main in one push would be rejected.
simple=$stage/mirrors/itpplasma/SIMPLE.git
git init -q --bare "$simple"
refspecs=("HEAD:refs/heads/$head_ref" "+refs/tags/*:refs/tags/*")
if [ "$head_ref" != main ]; then
    git -C "$src" fetch -q origin "+refs/heads/main:refs/remotes/origin/main"
    refspecs+=("refs/remotes/origin/main:refs/heads/main")
fi
git -C "$src" push -q "$simple" "${refspecs[@]}"
git -C "$simple" symbolic-ref HEAD refs/heads/main
git clone -q --branch "$head_ref" "$simple" "$stage/src"
git -C "$stage/src" remote set-url origin https://github.com/itpplasma/SIMPLE.git

while IFS=$'\t' read -r url name sha; do
    key=${sha:0:16}
    if ! echo "$sha  $cache/data/$key" | sha256sum -c --status 2>/dev/null; then
        curl -fsSL --retry 5 --retry-all-errors --retry-delay 10 \
            -o "$cache/data/$key.part" "$url"
        echo "$sha  $cache/data/$key.part" | sha256sum -c --quiet
        mv "$cache/data/$key.part" "$cache/data/$key"
    fi
    cp "$cache/data/$key" "$stage/data/$name"
    ln "$stage/data/$name" "$stage/src/test/test_data/$name"
    printf '%s\t%s\t%s\n' "$url" "$name" "$sha" >> "$stage/data/urls.tsv"
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

# Download on the connected runner for Debian 12's CPython 3.11/glibc 2.36.
# The trusted build dispatcher supplies python3/pip; accept compatible
# manylinux wheels, never a wheel requiring a newer cluster glibc.
mkdir -p "$stage/pip"
pip_args=(--only-binary=:all: --python-version 3.11 --implementation cp --abi cp311)
for glibc_minor in {17..36}; do
    pip_args+=(--platform "manylinux_2_${glibc_minor}_x86_64")
done
pip_args+=(--platform manylinux2014_x86_64)
# GVEC is source-only. Download its sdist separately, then include its
# declared dependencies in the cross-platform wheel resolution.
python3 -m pip download --dest "$stage/pip" --no-deps 'gvec==1.4.1'
sed '/^gvec==/d' "$src/requirements.txt" > "$stage/pip-requirements.txt"
python3 - "$stage" <<'PYTHON'
from pathlib import Path
import sys
import tarfile
import tomllib

stage = Path(sys.argv[1])
with tarfile.open(stage / "pip/gvec-1.4.1.tar.gz") as archive:
    metadata = tomllib.loads(archive.extractfile("gvec-1.4.1/pyproject.toml").read().decode())
with (stage / "pip-requirements.txt").open("a") as requirements:
    requirements.write("\n" + "\n".join(metadata["project"]["dependencies"]) + "\n")
PYTHON
# GVEC also requests setuptools-scm dynamically from its metadata backend.
requirements=(-r "$stage/pip-requirements.txt" pip setuptools setuptools-scm wheel)
if [ "${SCI_TIER:-full}" = full ]; then
    requirements+=(-r "$src/ci/coverage-requirements.txt" 'numpy==2.2.6' 'scipy==1.16.3' vtk \
        hatchling hatch-vcs mpmath sympy numexpr numba)
fi
python3 -m pip download --dest "$stage/pip" "${pip_args[@]}" "${requirements[@]}"

{
    printf 'SCI_COMMIT=%s\n' "$(git -C "$src" rev-parse HEAD)"
    printf 'SCI_HEAD_REF=%s\n' "$head_ref"
    printf 'SCI_LIBNEO_REF=%s\n' "${SCI_LIBNEO_REF:-}"
} > "$stage/bundle.env"

tar -C "$stage" -cf - . | zstd -q -1 -T0 -o "$out" -f
echo "bundle: $out ($(du -h "$out" | cut -f1), commit $(git -C "$src" rev-parse --short HEAD) on $head_ref)"
