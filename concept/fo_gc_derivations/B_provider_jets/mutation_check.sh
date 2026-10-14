#!/usr/bin/env bash
# Non-vacuity control: each mutation of a cited source line must make the
# generator fail (an exact identity refuted or a transcription anchor broken).
# Requires a prior run_provider_jets.sh (uses build/src and build/gen_provider_jets).
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
muts=(
 "batch_interpolate_3d.f90|2171|s/k1\*c/(k1 + 1)*c/"
 "batch_interpolate_3d.f90|2216|s/k2\*coeff_23/2*k2*coeff_23/"
 "batch_interpolate_3d.f90|1835|s/kk\*(kk - 1)/kk*(kk + 1)/"
 "batch_interpolate_3d.f90|1805|s/b1(k3, 3)/b0(k3, 3)/"
 "spectre_fo_gauss2.f90|120|s/2.0_dp\*s(0)/3.0_dp*s(0)/"
 "spectre_fo_gauss2.f90|117|s/-2.0_dp\/t/-3.0_dp\/t/"
 "spectre_fo_gauss2.f90|336|s/d2(:, j, k)/d2(:, j, l)/"
 "spectre_basis.f90|43|s/2.0_dp\*cheby(ll - 1, 0)/3.0_dp*cheby(ll - 1, 0)/"
 "spectre_basis.f90|289|s/+ rmi\*rni/- rmi*rni/"
 "spectre_basis.f90|93|s/(m + 1)\*m\*(m - 1)/(m + 1)*m*(m + 1)/"
)
status=0
for m in "${muts[@]}"; do
  IFS='|' read -r file line expr <<<"$m"
  rm -rf build/mut && cp -r build/src build/mut
  sed -i "${line}${expr}" "build/mut/$file"
  if cmp -s "build/src/$file" "build/mut/$file"; then echo "MUTATION NOT APPLIED: $m"; status=1; continue; fi
  out=$(build/gen_provider_jets build/mut/batch_interpolate_3d.f90 build/mut/spectre_fo_gauss2.f90 build/mut/spectre_basis.f90 2>&1)
  rc=$?
  nf=$(grep -c "FAILED:" <<<"$out")
  anchor=$(grep -c "anchor mismatch" <<<"$out")
  if [ $rc -ne 0 ]; then echo "detected rc=$rc failed_identities=$nf anchor=$anchor : $file:$line $expr"
  else echo "NOT DETECTED: $m"; status=1; fi
done
rm -rf build/mut
exit $status
