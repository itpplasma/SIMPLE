#!/usr/bin/env bash
# Build the item-D switch study (no SIMPLE build; libneo/fortio snapshot from
# /home/ert/code/SIMPLE/build/_deps, libneo d30174fb, copied into build/neo).
set -euo pipefail
D="$(cd "$(dirname "$0")/.." && pwd)"
B="$D/build"
mkdir -p "$B/mod"
FC=${FC:-gfortran}
FFLAGS="-O2 -fPIC -Wall -Wextra -fimplicit-none -ffree-line-length-none"
cd "$B"
[ gc_switch_generated.o -nt "$D/generated/gc_switch_generated.f90" ] || \
  $FC $FFLAGS -O0 -Wno-unused -J mod -c "$D/generated/gc_switch_generated.f90" -o gc_switch_generated.o
$FC $FFLAGS -J mod -I mod -c "$D/src/gc_fo_switch.f90" -o gc_fo_switch.o
$FC $FFLAGS -J mod -I mod -I neo/include -c "$D/src/spectre_jets.f90" -o spectre_jets.o
$FC $FFLAGS -J mod -I mod -I neo/include -c "$D/tests/switch_study.f90" -o switch_study.o
$FC -O2 switch_study.o spectre_jets.o gc_fo_switch.o gc_switch_generated.o -L neo/lib \
  -lspectre_field -lspectre_support -lmagnetic_field_base -lcoordinates -lfield_base \
  -linterpolate -lneo_util -lfortio -lhdf5_fortran -lhdf5 -llapack -lblas -o switch_study
echo built "$B/switch_study"
