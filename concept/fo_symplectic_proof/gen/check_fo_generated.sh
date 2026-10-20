#!/usr/bin/env bash
set -euo pipefail
proof=concept/fo_symplectic_proof
sha256sum -c $proof/reference.sha256
output=build/fo_generation_check
mkdir -p "$output/src/generated" "$output/$proof"
FO_GENERATOR_OUTPUT="$output" bash $proof/gen/regenerate_fo.sh
for leaf in fo_kernel.f90 fo_scalar.inc fo_scalar_first.inc fo_monitor.inc fo_product.inc fo_radial.inc fo_metric_monitor.inc \
    fo_geometry_2.inc fo_geometry_3.inc fo_potential_1.inc fo_potential_2.inc \
    fo_factor_velocity.inc fo_factor_hessian.inc fo_factor_flow.inc fo_hessian_proof.json; do
    diff -u "$proof/generated/$leaf" "$output/src/generated/$leaf"
done
