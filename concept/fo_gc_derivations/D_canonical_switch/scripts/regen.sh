#!/usr/bin/env bash
# Verify regeneration in the native build; preserve committed output for review.
set -euo pipefail
switch_root="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$switch_root"
exec make test CONFIG=Fast TEST='^test_switch_regeneration$' VERBOSE=1 "$@"
