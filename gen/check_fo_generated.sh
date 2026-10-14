#!/usr/bin/env bash
set -euo pipefail
sha256sum -c concept/fo_symplectic_proof/reference.sha256
output=build/fo_generation_check
mkdir -p "$output/src/generated" "$output/concept/fo_symplectic_proof"
FO_GENERATOR_OUTPUT="$output" bash gen/regenerate_fo.sh
for leaf in fo_kernel.f90 fo_scalar.inc fo_monitor.inc fo_product.inc fo_radial.inc fo_metric_monitor.inc; do
    diff -u "src/generated/$leaf" "$output/src/generated/$leaf"
done
