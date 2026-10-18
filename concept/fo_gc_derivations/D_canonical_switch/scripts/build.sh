#!/usr/bin/env bash
# Build through SIMPLE's owning CMake/Ninja authority.
set -euo pipefail
study_root="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$study_root"
exec make CONFIG=Fast "$@"
