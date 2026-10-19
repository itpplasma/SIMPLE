#!/usr/bin/env bash
# Reproduce the L1 provider-jet identities. Run from the SIMPLE root:
#   bash concept/fo_gc_derivations/B_provider_jets/run_provider_jets.sh
# Needs FortSym (FORTSYM_BUILD, default ../fortsym/build) and read access to
# the pinned sources below. Writes results_symbolic.json in this directory.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../../.." && pwd)
fortsym_build=${FORTSYM_BUILD:-$root/../fortsym/build}
fo_gauss2=${FO_GAUSS2_SOURCE:-$root/concept/fo_symplectic_proof/spectre_fo_gauss2.f90}
libneo=${LIBNEO_SRC:-$root/build/_deps/libneo-src}
mkdir -p "$here/build/src"
cp "$fo_gauss2" "$here/build/src/spectre_fo_gauss2.f90"
cp "$libneo/src/interpolate/batch_interpolate_3d.f90" "$here/build/src/"
cp "$libneo/src/spectre/spectre_basis.f90" "$here/build/src/"
( cd "$here/build/src" && sha256sum -c --quiet - ) <<'SUMS'
6890146f4539968c2f4f660830837be5f555503eb714b860751a6478c9b5f1ad  spectre_fo_gauss2.f90
a249da9ca99fc217f621f498ad727d1c7dd3cc6ff02fd7ea079e6b43ab8adfc9  batch_interpolate_3d.f90
21c3cefd7c6a518a9090e9a98fc2547a14fe7521c46479ff6e8e8a5f2d406d7d  spectre_basis.f90
SUMS
cd "$here/build"
${FC:-gfortran} -O0 -I "$fortsym_build/include" "$here/provider_src.f90" \
    "$here/gen_provider_jets.f90" -L "$fortsym_build/lib" \
    -Wl,-rpath,"$fortsym_build/lib" -lfortsym -o gen_provider_jets
cd "$here"
build/gen_provider_jets build/src/batch_interpolate_3d.f90 \
    build/src/spectre_fo_gauss2.f90 build/src/spectre_basis.f90
