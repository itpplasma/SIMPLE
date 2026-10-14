#!/usr/bin/env bash
# Check the GC symplecticity theorems with the pinned toolchain (run from SIMPLE root).
# Fails on any error, warning (including sorry) or non-standard axiom.
set -euo pipefail
out=$(lake +leanprover/lean4:v4.16.0 -d concept/gc_symplectic_proof/lean \
    env lean concept/gc_symplectic_proof/lean/GcSymplectic.lean 2>&1)
printf '%s\n' "$out"
if printf '%s\n' "$out" | grep -Eiq 'error|warning|sorry'; then
    echo "FAIL: Lean reported errors, warnings or sorry"; exit 1
fi
expected=22
standard=$(printf '%s\n' "$out" | grep -c \
    "depends on axioms: \[propext, Classical.choice, Quot.sound\]$" || true)
total=$(printf '%s\n' "$out" | grep -c "depends on axioms" || true)
if [ "$standard" -ne "$expected" ] || [ "$total" -ne "$expected" ]; then
    echo "FAIL: expected $expected theorems on standard axioms, got $standard/$total"
    exit 1
fi
echo "PASS: $expected theorems, standard axioms only"
