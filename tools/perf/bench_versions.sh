#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  tools/perf/bench_versions.sh [options]

Benchmarks non-SoA CPU tracing across multiple git refs/tags and field modes.
It creates temporary git worktrees under /tmp, applies print-suppression patches
to historical tags (only in /tmp worktrees), builds with `make clean && make`,
and runs `build/simple.x bench.in` with 32-thread OpenMP by default.

Options:
  --refs <csv>          Git refs to benchmark (default: HEAD,main,v1.5.1,v1.4.2,v1.3.1)
  --modes <csv>         Modes (default: boozer,canflux,meiss,albert)
  --threads <n>         OpenMP threads (default: 32)
  --particles <n>       ntestpart (default: 1024)
  --trace-time <x>      trace_time (default: 1d-2)
  --npoiper <n>         npoiper (default: 100)
  --npoiper2 <n>        npoiper2 (default: 128)
  --ntimstep <n>        ntimstep (default: 10000)
  --nper <n>            nper (default: 1000)
  --integmode <n>       integmode (default: 1)
  --relerr <x>          relerr (default: 1d-13)
  --wout <path>         VMEC file (default: ./wout.nc from repo root)
  --keep-worktrees      Do not delete /tmp worktrees created by this run
  -h|--help             Show help

Output:
  Prints an ASCII table to stdout and writes per-run logs under /tmp.
EOF
}

ROOT="$(git rev-parse --show-toplevel)"

REFS_CSV="HEAD,main,v1.5.1,v1.4.2,v1.3.1"
MODES_CSV="boozer,canflux,meiss,albert"

THREADS=32
PARTICLES=1024
TRACE_TIME="1d-2"
NPOIPER=100
NPOIPER2=128
NTIMSTEP=10000
NPER=1000
INTEGMODE=1
RELERR="1d-13"
WOUT="${ROOT}/wout.nc"
KEEP_WORKTREES=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --refs) REFS_CSV="$2"; shift 2 ;;
    --modes) MODES_CSV="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --particles) PARTICLES="$2"; shift 2 ;;
    --trace-time) TRACE_TIME="$2"; shift 2 ;;
    --npoiper) NPOIPER="$2"; shift 2 ;;
    --npoiper2) NPOIPER2="$2"; shift 2 ;;
    --ntimstep) NTIMSTEP="$2"; shift 2 ;;
    --nper) NPER="$2"; shift 2 ;;
    --integmode) INTEGMODE="$2"; shift 2 ;;
    --relerr) RELERR="$2"; shift 2 ;;
    --wout) WOUT="$2"; shift 2 ;;
    --keep-worktrees) KEEP_WORKTREES=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ ! -f "$WOUT" ]]; then
  echo "Missing VMEC file: $WOUT" >&2
  exit 1
fi

declare -a CREATED_WORKTREES=()

cleanup() {
  if [[ "$KEEP_WORKTREES" -eq 1 ]]; then
    return
  fi
  for d in "${CREATED_WORKTREES[@]}"; do
    git -C "$ROOT" worktree remove --force "$d" >/dev/null 2>&1 || true
  done
}
trap cleanup EXIT

ref_slug() {
  echo "$1" | tr '/:@' '____' | tr -cs 'A-Za-z0-9_.-' '_' | sed 's/^_//;s/_$//'
}

mode_field_type() {
  case "$1" in
    canflux) echo 0 ;;
    boozer) echo 2 ;;
    meiss) echo 3 ;;
    albert) echo 4 ;;
    *) echo "Unknown mode: $1" >&2; exit 2 ;;
  esac
}

worktree_dir_for_ref() {
  local ref="$1"
  echo "/tmp/SIMPLE_bench_$(ref_slug "$ref")"
}

worktree_ensure() {
  local ref="$1"
  local dir
  dir="$(worktree_dir_for_ref "$ref")"
  if [[ -d "$dir" ]]; then
    echo "$dir"
    return
  fi
  git -C "$ROOT" worktree add --detach "$dir" "$ref" >/dev/null
  CREATED_WORKTREES+=("$dir")
  echo "$dir"
}

patch_file_python() {
  local file_path="$1"
  local tag="$2"
  python - "$file_path" "$tag" <<'PY'
import re
import sys

path = sys.argv[1]
tag = sys.argv[2]

with open(path, "r", encoding="utf-8", newline="") as f:
    text = f.read()

orig = text

def replace_once(pattern, repl, *, flags=0):
    global text
    new, n = re.subn(pattern, repl, text, count=1, flags=flags)
    if n == 1:
        text = new
        return True
    return False

def ensure_insert_after(anchor_pat, insert_text, *, flags=0):
    global text
    if insert_text in text:
        return
    m = re.search(anchor_pat, text, flags)
    if not m:
        raise RuntimeError(f"{path}: anchor not found: {anchor_pat!r}")
    i = m.end()
    text = text[:i] + insert_text + text[i:]

def ensure_print_disabled(msg_prefix):
    global text
    pat = rf"(^[ \t]*)print[ \t]*\\*,[ \t]*'{re.escape(msg_prefix)}(.*)$"
    repl = r"\\1if (.false.) print *,'" + msg_prefix + r"\\2"
    if re.search(rf"(^[ \t]*)if[ \t]*\\(\\.false\\.\\)[ \t]*print[ \t]*\\*,[ \t]*'{re.escape(msg_prefix)}", text, flags=re.M):
        return
    if replace_once(pat, repl, flags=re.M):
        return
    # Fallback for formatting variations
    needle = "print *,'" + msg_prefix
    guarded = "if (.false.) print *,'" + msg_prefix
    if needle in text and guarded not in text:
        text = text.replace(needle, guarded, 1)
        return
    raise RuntimeError(f"{path}: expected print line for: {msg_prefix}")

def patch_check_orbit_type_common():
    ensure_print_disabled("detect_oneline: nfp > nfp_max, nturns not reached")
    ensure_print_disabled("detect_oneline: stochastic, iper < ipermin")
    ensure_print_disabled("detect_oneline: return period")

def patch_get_cancoord_print_progress():
    global text
    # Disable progress printing entirely to avoid I/O dominating perf comparisons.
    if "subroutine print_progress" not in text:
        return
    if re.search(r"subroutine[ \t]+print_progress\([^\)]*\)[\s\S]*?^[ \t]*return[ \t]*$", text, flags=re.M | re.I):
        return
    # Insert return as first executable statement (right before the first write).
    if not replace_once(
        r"(^[ \t]*write[ \t]*\(\*,)",
        "      return\n\\g<0>",
        flags=re.M,
    ):
        raise RuntimeError(f"{path}: could not insert return before write(*) in print_progress")

def patch_field_can_meiss_print_progress():
    global text
    if "subroutine print_progress(i_phi)" not in text:
        return
    if re.search(r"subroutine[ \t]+print_progress\(i_phi\)[\s\S]*?^[ \t]*return[ \t]*$", text, flags=re.M | re.I):
        return
    if not replace_once(
        r"(^[ \t]*!\$omp[ \t]+critical[ \t]*$)",
        "    return\n\\g<0>",
        flags=re.M,
    ):
        raise RuntimeError(f"{path}: could not insert return before !$omp critical in print_progress(i_phi)")

def patch_field_splined_print_progress():
    global text
    if "subroutine print_progress(i, n)" not in text:
        return
    if re.search(r"subroutine[ \t]+print_progress\(i, n\)[\s\S]*?^[ \t]*return[ \t]*$", text, flags=re.M | re.I):
        return
    if not replace_once(
        r"(^[ \t]*!\$omp[ \t]+critical[ \t]*$)",
        "        return\n\\g<0>",
        flags=re.M,
    ):
        raise RuntimeError(f"{path}: could not insert return before !$omp critical in print_progress(i, n)")

def patch_field_coils_print_progress():
    global text
    if "subroutine print_progress(i)" not in text:
        return
    if re.search(r"subroutine[ \t]+print_progress\(i\)[\s\S]*?^[ \t]*return[ \t]*$", text, flags=re.M | re.I):
        return
    if not replace_once(
        r"(^[ \t]*!\$omp[ \t]+critical[ \t]*$)",
        "        return\n\\g<0>",
        flags=re.M,
    ):
        raise RuntimeError(f"{path}: could not insert return before !$omp critical in print_progress(i)")

def patch_v131_get_canonical_coordinates():
    global text
    # Disable omp-critical progress printing blocks.
    if "if (.false.) then" in text and "integrate ODE:" in text:
        return
    text = text.replace(
        "!$omp critical\n    i_ctr = i_ctr + 1\n    print *,'integrate ODE: ',i_ctr,' of ',n_theta_c\n!$omp end critical\n",
        "if (.false.) then\n!$omp critical\n    i_ctr = i_ctr + 1\n    print *,'integrate ODE: ',i_ctr,' of ',n_theta_c\n!$omp end critical\nend if\n",
    )
    text = text.replace(
        "!$omp critical\n    i_ctr = i_ctr + 1\n    print *,'compute components: ',i_ctr,' of ',n_theta_c\n!$omp end critical\n",
        "if (.false.) then\n!$omp critical\n    i_ctr = i_ctr + 1\n    print *,'compute components: ',i_ctr,' of ',n_theta_c\n!$omp end critical\nend if\n",
    )

def patch_v131_simple_progress():
    global text
    if "kpart = kpart+1" not in text:
        return
    if "if (.false.) then" in text:
        return
    # Disable the progress counter+prints (critical section).
    text = text.replace("!$omp critical\n    kpart = kpart+1\n",
                        "if (.false.) then\n    !$omp critical\n    kpart = kpart+1\n")
    text = text.replace("!$omp end critical\n    call trace_orbit(norb, i)\n",
                        "    !$omp end critical\nend if\n    call trace_orbit(norb, i)\n")

if tag in ("v1.4.2", "v1.5.1"):
    if path.endswith("check_orbit_type.f90"):
        patch_check_orbit_type_common()
    elif path.endswith("get_canonical_coordinates.F90"):
        patch_get_cancoord_print_progress()
    elif path.endswith("field_can_meiss.f90"):
        patch_field_can_meiss_print_progress()
    elif path.endswith("field_splined.f90"):
        patch_field_splined_print_progress()
    elif path.endswith("field_coils.f90"):
        patch_field_coils_print_progress()
elif tag == "v1.3.1":
    if path.endswith("check_orbit_type.f90"):
        patch_check_orbit_type_common()
    elif path.endswith("get_canonical_coordinates.f90"):
        patch_v131_get_canonical_coordinates()
    elif path.endswith("simple.f90"):
        patch_v131_simple_progress()
else:
    raise RuntimeError(f"Unsupported tag for patching: {tag}")

if text == orig:
    # Idempotent: allow no-op only if tag was already patched.
    pass

with open(path, "w", encoding="utf-8", newline="") as f:
    f.write(text)
PY
}

patch_suppress_prints_for_tag() {
  local ref="$1"
  local dir="$2"

  case "$ref" in
    v1.3.1)
      patch_file_python "${dir}/SRC/check_orbit_type.f90" "v1.3.1"
      patch_file_python "${dir}/SRC/get_canonical_coordinates.f90" "v1.3.1"
      patch_file_python "${dir}/SRC/simple.f90" "v1.3.1"
      ;;
    v1.4.2)
      patch_file_python "${dir}/src/check_orbit_type.f90" "v1.4.2"
      patch_file_python "${dir}/src/get_canonical_coordinates.F90" "v1.4.2"
      patch_file_python "${dir}/src/field/field_can_meiss.f90" "v1.4.2"
      patch_file_python "${dir}/src/field/field_coils.f90" "v1.4.2"
      ;;
    v1.5.1)
      patch_file_python "${dir}/src/check_orbit_type.f90" "v1.5.1"
      patch_file_python "${dir}/src/get_canonical_coordinates.F90" "v1.5.1"
      patch_file_python "${dir}/src/field/field_can_meiss.f90" "v1.5.1"
      patch_file_python "${dir}/src/field/field_splined.f90" "v1.5.1"
      ;;
    *)
      return
      ;;
  esac
}

write_bench_in_v131() {
  local out="$1"
  local field_type="$2"
  cat >"$out" <<EOF
&config
notrace_passing = 0
nper = ${NPER}
npoiper = ${NPOIPER}
ntimstep = ${NTIMSTEP}
ntestpart = ${PARTICLES}
trace_time = ${TRACE_TIME}
sbeg = 0.6d0
phibeg = 0.d0
thetabeg = 0.d0
contr_pp = -1.0d0
facE_al = 1.0
npoiper2 = ${NPOIPER2}
n_e = 2
n_d = 4
netcdffile = 'wout.nc'
ns_s = 5
ns_tp = 5
multharm = 5
isw_field_type = ${field_type}
startmode = 1
integmode = ${INTEGMODE}
relerr = ${RELERR}
tcut = -1d0
debug = .False.
class_plot = .False.
cut_in_per = 0d0
fast_class = .False.
vmec_B_scale = 1.0d0
vmec_RZ_scale = 1.0d0
swcoll = .False.
deterministic = .True.
old_axis_healing = .True.
old_axis_healing_boundary = .True.
/
EOF
}

write_bench_in_v142() {
  local out="$1"
  local field_type="$2"
  cat >"$out" <<EOF
&config
notrace_passing = 0
nper = ${NPER}
npoiper = ${NPOIPER}
ntimstep = ${NTIMSTEP}
ntestpart = ${PARTICLES}
trace_time = ${TRACE_TIME}
num_surf = 1
sbeg = 0.6d0
phibeg = 0.d0
thetabeg = 0.d0
contr_pp = -1.0d0
facE_al = 1.0
npoiper2 = ${NPOIPER2}
n_e = 2
n_d = 4
netcdffile = 'wout.nc'
ns_s = 5
ns_tp = 5
multharm = 5
isw_field_type = ${field_type}
generate_start_only = .False.
startmode = 1
grid_density = 0d0
special_ants_file = .False.
integmode = ${INTEGMODE}
relerr = ${RELERR}
tcut = -1d0
debug = .False.
class_plot = .False.
cut_in_per = 0d0
fast_class = .False.
vmec_B_scale = 1.0d0
vmec_RZ_scale = 1.0d0
swcoll = .False.
deterministic = .True.
old_axis_healing = .True.
old_axis_healing_boundary = .True.
batch_size = 2000000000
ran_seed = 12345
reuse_batch = .False.
field_input = ''
output_error = .False.
output_orbits_macrostep = .False.
/
EOF
}

write_bench_in_v151plus() {
  local out="$1"
  local field_type="$2"
  cat >"$out" <<EOF
&config
notrace_passing = 0
nper = ${NPER}
npoiper = ${NPOIPER}
ntimstep = ${NTIMSTEP}
ntestpart = ${PARTICLES}
trace_time = ${TRACE_TIME}
num_surf = 1
sbeg = 0.6d0
phibeg = 0.d0
thetabeg = 0.d0
contr_pp = -1.0d0
facE_al = 1.0
npoiper2 = ${NPOIPER2}
n_e = 2
n_d = 4
netcdffile = 'wout.nc'
ns_s = 5
ns_tp = 5
multharm = 5
isw_field_type = ${field_type}
generate_start_only = .False.
startmode = 1
grid_density = 0d0
special_ants_file = .False.
integmode = ${INTEGMODE}
relerr = ${RELERR}
tcut = -1d0
nturns = 1
debug = .False.
class_plot = .False.
cut_in_per = 0d0
fast_class = .False.
vmec_B_scale = 1.0d0
vmec_RZ_scale = 1.0d0
swcoll = .False.
deterministic = .True.
old_axis_healing = .True.
old_axis_healing_boundary = .True.
batch_size = 2000000000
ran_seed = 12345
reuse_batch = .False.
field_input = ''
coord_input = ''
integ_coords = -1000
output_error = .False.
output_orbits_macrostep = .False.
/
EOF
}

write_bench_in_for_dir() {
  local dir="$1"
  local out="$2"
  local field_type="$3"

  if [[ -d "${dir}/SRC" ]]; then
    write_bench_in_v131 "$out" "$field_type"
    return
  fi

  if rg -n "coord_input" -S "${dir}/src/params.f90" >/dev/null 2>&1; then
    write_bench_in_v151plus "$out" "$field_type"
    return
  fi

  write_bench_in_v142 "$out" "$field_type"
}

bench_exe_for_dir() {
  local dir="$1"
  if [[ -x "${dir}/build/simple.x" ]]; then
    echo "${dir}/build/simple.x"
    return
  fi
  if [[ -x "${dir}/build/simple" ]]; then
    echo "${dir}/build/simple"
    return
  fi
  echo "Missing build executable under ${dir}/build" >&2
  exit 1
}

supports_mode_for_ref() {
  local ref="$1"
  local mode="$2"
  case "$ref" in
    v1.3.1)
      [[ "$mode" == "boozer" || "$mode" == "canflux" ]]
      ;;
    *)
      return 0
      ;;
  esac
}

build_dir() {
  local dir="$1"
  local ref="$2"
  local build_log
  build_log="/tmp/simple_build_$(ref_slug "$ref").log"
  (cd "$dir" && make clean && make -j32) 2>&1 | tee "$build_log" >/dev/null
}

extract_metrics() {
  local log_file="$1"
  local trace_time_num="$2"
  local times_file="$3"

  local wall
  wall="$(awk '/^WALL[[:space:]]/ {print $2}' "$log_file" | tail -n 1)"
  local trace
  trace="$(awk '/INFO: Parallel particle tracing completed/ {print $(NF-4)}' "$log_file" | tail -n 1)"
  if [[ -z "${trace:-}" ]]; then
    trace="-"
  fi

  local n neg lost end
  read -r n neg lost end < <(
    awk -v T="$trace_time_num" 'BEGIN{n=0;neg=0;lost=0;end=0}
      NF>0{
        n++; t=$2;
        if (t < 0) {neg++}
        else if (t < T-1e-12) {lost++}
        else {end++}
      }
      END{print n, neg, lost, end}' "$times_file"
  )

  printf "%s %s %s %s %s %s" "$wall" "$trace" "$n" "$neg" "$lost" "$end"
}

now_tag="$(date +%Y%m%d_%H%M%S)"
trace_time_num="${TRACE_TIME//d/e}"
trace_time_num="${trace_time_num//D/e}"

echo "CPU: $(lscpu | awk -F: '/Model name/ {gsub(/^[ \t]+/,"",$2); print $2; exit}')"
echo "OpenMP: threads=${THREADS} OMP_PROC_BIND=spread OMP_PLACES=cores OMP_DYNAMIC=false"
echo "Case: ntestpart=${PARTICLES} trace_time=${TRACE_TIME} nper=${NPER} npoiper=${NPOIPER} npoiper2=${NPOIPER2} ntimstep=${NTIMSTEP} integmode=${INTEGMODE} relerr=${RELERR}"
echo

printf "%-12s %-8s %10s %10s %8s %8s %8s %8s\n" "REF" "MODE" "WALL(s)" "TRACE(s)" "N" "NEG" "LOST" "END"
printf "%-12s %-8s %10s %10s %8s %8s %8s %8s\n" "------------" "--------" "--------" "--------" "--------" "--------" "--------" "--------"

IFS=',' read -r -a REFS <<<"$REFS_CSV"
IFS=',' read -r -a MODES <<<"$MODES_CSV"

for ref in "${REFS[@]}"; do
  dir="$(worktree_ensure "$ref")"
  patch_suppress_prints_for_tag "$ref" "$dir"
  build_dir "$dir" "$ref"
  exe="$(bench_exe_for_dir "$dir")"

  for mode in "${MODES[@]}"; do
    if ! supports_mode_for_ref "$ref" "$mode"; then
      continue
    fi
    field_type="$(mode_field_type "$mode")"
    run_dir="/tmp/simple_bench_${now_tag}_$(ref_slug "$ref")_${mode}"
    mkdir -p "$run_dir"
    ln -sf "$WOUT" "$run_dir/wout.nc"
    write_bench_in_for_dir "$dir" "$run_dir/bench.in" "$field_type"

    log_file="$run_dir/run.log"
    (
      cd "$run_dir"
      export OMP_NUM_THREADS="$THREADS"
      export OMP_PROC_BIND=spread
      export OMP_PLACES=cores
      export OMP_DYNAMIC=false
      export LC_ALL=C
      TIMEFORMAT='WALL %3R'
      time "$exe" bench.in
    ) 2>&1 | tee "$log_file" >/dev/null

    metrics="$(extract_metrics "$log_file" "$trace_time_num" "$run_dir/times_lost.dat")"
    wall="$(awk '{print $1}' <<<"$metrics")"
    trace="$(awk '{print $2}' <<<"$metrics")"
    n="$(awk '{print $3}' <<<"$metrics")"
    neg="$(awk '{print $4}' <<<"$metrics")"
    lost="$(awk '{print $5}' <<<"$metrics")"
    end="$(awk '{print $6}' <<<"$metrics")"

    printf "%-12s %-8s %10s %10s %8s %8s %8s %8s\n" "$ref" "$mode" "$wall" "$trace" "$n" "$neg" "$lost" "$end"
  done
done
