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
for generator in geometry hessian phases poly poly_axis joint; do
    ${FC:-gfortran} -I "$fortsym_build/include" "$proof/gen/gen_fo_${generator}.f90" \
        -L "$fortsym_build/lib" -Wl,-rpath,"$fortsym_build/lib" -lfortsym \
        -o "build/fo_generator/gen_fo_${generator}"
done
finish_leaves() {
    local directory=$1
    FO_GENERATOR_OUTPUT="$directory" build/fo_generator/gen_fo_geometry
    FO_GENERATOR_OUTPUT="$directory" build/fo_generator/gen_fo_hessian
    FO_GENERATOR_OUTPUT="$directory" build/fo_generator/gen_fo_phases
    FO_GENERATOR_OUTPUT="$directory" build/fo_generator/gen_fo_poly
    FO_GENERATOR_OUTPUT="$directory" build/fo_generator/gen_fo_poly_axis
    python3 "$proof/gen/axis_facades.py" "$directory"
    FO_GENERATOR_OUTPUT="$directory" build/fo_generator/gen_fo_joint
    python3 "$proof/gen/joint_facades.py" "$directory"
    fprettify -c .fprettify "$directory/fo_joint_kernel.f90"
    python3 "$proof/gen/wrap_poly.py" "$directory"
    python3 "$proof/gen/array_facades.py" "$directory"
    fprettify -c .fprettify "$directory/fo_kernel.f90"
}
if [ -n "${FO_GENERATOR_OUTPUT:-}" ]; then
    build/fo_generator/gen_fo_hamiltonian
    python3 "$proof/gen/format_first_facade.py" \
        "$FO_GENERATOR_OUTPUT/src/generated/fo_kernel.f90"
    finish_leaves "$FO_GENERATOR_OUTPUT/src/generated"
    exit 0
fi
out=build/fo_generator/out
mkdir -p "$out/src/generated" "$out/$proof"
FO_GENERATOR_OUTPUT="$out" build/fo_generator/gen_fo_hamiltonian
python3 "$proof/gen/format_first_facade.py" "$out/src/generated/fo_kernel.f90"
finish_leaves "$out/src/generated"
cp "$out"/src/generated/* "$proof/generated/"
cp "$out/$proof/symbolic.json" "$out/$proof/monitor.json" "$proof/"
