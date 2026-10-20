#!/usr/bin/env bash
set -euo pipefail
proof=concept/fo_symplectic_proof
sha256sum -c $proof/reference.sha256
output=build/fo_generation_check
mkdir -p "$output/src/generated" "$output/$proof"
FO_GENERATOR_OUTPUT="$output" bash $proof/gen/regenerate_fo.sh
for leaf in fo_joint_kernel.f90 fo_joint_proof.json \
    fo_poly_axis_1.inc fo_poly_axis_2.inc fo_poly_axis_proof.json fo_poly_coeff.inc fo_poly_point.inc fo_poly_horner_1.inc fo_poly_horner_2.inc \
    fo_poly_angular_1.inc fo_poly_angular_2.inc fo_poly_proof.json fo_kernel.f90 fo_scalar.inc fo_scalar_first.inc fo_monitor.inc fo_product.inc fo_radial.inc fo_metric_monitor.inc \
    fo_geometry_2.inc fo_geometry_3.inc fo_potential_1.inc fo_potential_2.inc \
    fo_radial_first.inc fo_metric_monitor_first.inc \
    fo_phase_sum.inc fo_phase_difference.inc fo_phase_proof.json \
    fo_phase_seed.inc fo_phase_orders.inc \
    fo_factor_velocity.inc fo_factor_hessian.inc fo_factor_flow.inc fo_hessian_proof.json; do
    diff -u "$proof/generated/$leaf" "$output/src/generated/$leaf"
done
