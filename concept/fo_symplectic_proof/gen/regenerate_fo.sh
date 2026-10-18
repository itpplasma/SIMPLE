#!/usr/bin/env bash
set -euo pipefail
# Run from SIMPLE root. FortSym's build is supplied explicitly; no CAS is
# required for normal SIMPLE builds. This only runs the concept generator.
# The generator writes <out>/src/generated/* and <out>/concept/fo_symplectic_proof/*.json;
# without FO_GENERATOR_OUTPUT the leaves are copied into the concept tree.
proof=concept/fo_symplectic_proof
fortsym_build=${FORTSYM_BUILD:-../fortsym/build}
mkdir -p build/fo_generator
${FC:-gfortran} -I "$fortsym_build/include" $proof/gen/gen_fo_hamiltonian.f90 \
    -L "$fortsym_build/lib" -Wl,-rpath,"$fortsym_build/lib" -lfortsym \
    -o build/fo_generator/gen_fo_hamiltonian
if [ -n "${FO_GENERATOR_OUTPUT:-}" ]; then
    build/fo_generator/gen_fo_hamiltonian
    exit 0
fi
out=build/fo_generator/out
mkdir -p "$out/src/generated" "$out/$proof"
FO_GENERATOR_OUTPUT="$out" build/fo_generator/gen_fo_hamiltonian
cp "$out"/src/generated/* "$proof/generated/"
cp "$out/$proof/symbolic.json" "$out/$proof/monitor.json" "$proof/"
