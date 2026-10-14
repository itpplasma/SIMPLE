#!/usr/bin/env bash
set -euo pipefail
# Explicit toolchain: elan selects from invocation cwd, not Lake's -d argument.
lake +leanprover/lean4:v4.16.0 -d concept/fo_symplectic_proof/lean \
    env lean concept/fo_symplectic_proof/lean/Gauss2.lean
