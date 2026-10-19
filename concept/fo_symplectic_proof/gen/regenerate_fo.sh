#!/usr/bin/env bash
set -euo pipefail
# Run from SIMPLE root. FortSym's build is supplied explicitly; no CAS is
# required for normal SIMPLE builds. This only runs the concept generator.
# The array facade uses fprettify 0.3.7, matching changed-line lint.
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
    python3 "$proof/gen/format_first_facade.py" \
        "$FO_GENERATOR_OUTPUT/src/generated/fo_kernel.f90"
    exit 0
fi
out=build/fo_generator/out
mkdir -p "$out/src/generated" "$out/$proof"
FO_GENERATOR_OUTPUT="$out" build/fo_generator/gen_fo_hamiltonian
python3 "$proof/gen/format_first_facade.py" "$out/src/generated/fo_kernel.f90"
cp "$out"/src/generated/* "$proof/generated/"
cp "$out/$proof/symbolic.json" "$out/$proof/monitor.json" "$proof/"
