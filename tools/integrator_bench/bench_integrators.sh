#!/bin/bash
# Sweep all SIMPLE integrators on a small case and report completion, confined
# fraction, inner-solver non-convergence, and wall time. Boy-scout check that
# every integrator runs cleanly + a quick performance benchmark. Run one case at
# a time so it does not compete with production.
#
# Usage: bench_integrators.sh <simple.x> <wout.nc> <base_simple.in> [outdir]
set -u
BIN=${1:?simple.x}
WOUT=${2:?wout.nc}
BASE=${3:?base simple.in}
OUT=${4:-/tmp/integ_bench}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}

# integmode -> name
declare -A NAME=( [0]=RK45 [1]=Euler1 [2]=Euler2 [3]=midpoint \
                  [4]=Gauss1 [5]=Gauss2 [6]=Gauss3 [7]=Gauss4 [15]=Lobatto3 )

printf '%-9s %-9s %-6s %-10s %-14s %-9s\n' integmode name rc conf_frac newton_maxit wall_s
for im in 0 1 2 3 4 5 6 7 15; do
  d="$OUT/im${im}"
  mkdir -p "$d"
  cp "$WOUT" "$d/wout.nc"
  cp "$BASE" "$d/simple.in"
  sed -i "s/^integmode = .*/integmode = ${im}/" "$d/simple.in"
  grep -q '^integmode' "$d/simple.in" || echo "integmode = ${im}" >> "$d/simple.in"
  cd "$d"
  t0=$(date +%s.%N)
  "$BIN" > log 2>&1
  rc=$?
  t1=$(date +%s.%N)
  wall=$(awk "BEGIN{printf \"%.2f\", $t1-$t0}")
  cf=$(tail -1 confined_fraction.dat 2>/dev/null | awk '{print $2}')
  mx=$(grep -oE '(newton1_maxit|newton2_maxit|rk_gauss_maxit|rk_lobatto_maxit|fixpoint_maxit)=[0-9]+' log | \
        awk -F= '{s+=$2} END{print s+0}')
  printf '%-9s %-9s %-6s %-10s %-14s %-9s\n' "$im" "${NAME[$im]:-?}" "$rc" "${cf:-NA}" "${mx:-0}" "$wall"
done
