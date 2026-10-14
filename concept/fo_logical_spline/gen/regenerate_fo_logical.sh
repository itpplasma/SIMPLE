#!/usr/bin/env bash
set -euo pipefail
# Run from the SIMPLE root. Regenerates src/generated/fo_logical_rhs.inc from
# gen/gen_fo_logical.f90. FortSym is needed only here, never for SIMPLE builds.
fortsym_build=${FORTSYM_BUILD:-../fortsym/build}
out=${FO_GENERATOR_OUTPUT:-.}
mkdir -p build/fo_generator "$out/src/generated" "$out/build/fo_generator"
${FC:-gfortran} -I "$fortsym_build/include" gen/gen_fo_logical.f90 \
    -L "$fortsym_build/lib" -Wl,-rpath,"$fortsym_build/lib" -lfortsym \
    -o build/fo_generator/gen_fo_logical
FO_GENERATOR_OUTPUT="$out" build/fo_generator/gen_fo_logical
