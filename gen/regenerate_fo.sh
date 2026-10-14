#!/usr/bin/env bash
set -euo pipefail
# Run from SIMPLE root. FortSym's build is supplied explicitly; no CAS is
# required for normal SIMPLE builds. This only runs the consumer generator.
fortsym_build=${FORTSYM_BUILD:-../fortsym/build}
mkdir -p build/fo_generator
${FC:-gfortran} -I "$fortsym_build/include" gen/gen_fo_hamiltonian.f90 \
    -L "$fortsym_build/lib" -Wl,-rpath,"$fortsym_build/lib" -lfortsym \
    -o build/fo_generator/gen_fo_hamiltonian
build/fo_generator/gen_fo_hamiltonian
