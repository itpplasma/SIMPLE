# Non-SoA CPU performance (OpenMP) cross-version comparison

This note documents one focused performance check: **non-SoA** orbit tracing on CPU with OpenMP, comparing the current `perf/restore-v131` branch against older release tags where the same field mode exists.

## Benchmark environment

- CPU: AMD Ryzen 9 5950X (32 hardware threads)
- OpenMP:
  - `OMP_NUM_THREADS=32`
  - `OMP_PROC_BIND=spread`
  - `OMP_PLACES=cores`
  - `OMP_DYNAMIC=false`
- Build: CMake via `make` (`make clean && make -j32`), Release configuration

## Benchmark case (all runs)

These parameters were used to keep setup overhead small and to emphasize tracing throughput:

- `ntestpart = 1024`
- `trace_time = 1d-2`
- `ntimstep = 10000`
- `nper = 1000`
- `npoiper = 100`
- `npoiper2 = 128`
- `integmode = 1` (Euler1)
- `relerr = 1d-13`

## Output metrics

- `WALL(s)`: wall clock time from the shell `time` wrapper (`WALL <seconds>`)
- `TRACE(s)`: SIMPLE internal timing line `INFO: Parallel particle tracing completed ...`
- `NEG/LOST/END`: derived from `times_lost.dat` second column (`t_lost`)
  - `NEG`: `t_lost < 0`
  - `LOST`: `0 <= t_lost < trace_time`
  - `END`: `t_lost >= trace_time` (survived until end of trace)

## Current branch numbers (perf/restore-v131)

Mode mapping:
- CANFLUX = `isw_field_type = 0`
- BOOZER  = `isw_field_type = 2`
- MEISS   = `isw_field_type = 3`
- ALBERT  = `isw_field_type = 4`

Measured with the benchmark case above:

| MODE   | WALL(s) | TRACE(s) | N    | NEG | LOST | END |
|--------|---------|----------|------|-----|------|-----|
| BOOZER | 14.450  | 9.530    | 1024 | 497 | 108  | 419 |
| CANFLUX| 30.064  | 24.799   | 1024 | 497 | 105  | 422 |
| MEISS  | 46.233  | 38.118   | 1024 | 497 | 108  | 419 |
| ALBERT | 34.089  | 23.556   | 1024 | 497 | 164  | 363 |

## Cross-version comparison (same benchmark case)

Numbers below are from a single run of:

`tools/perf/bench_versions.sh --refs HEAD,v1.5.1,v1.4.2,v1.3.1 --modes canflux,boozer,meiss,albert --threads 32 --particles 1024 --trace-time 1d-2 --npoiper 100 --npoiper2 128 --integmode 1`

Important:
- Older tags printed progress unconditionally. For apples-to-apples timing, the script **patches those worktrees in `/tmp`** to disable progress output (prints are replaced by constant `.false.` guards or `return` in progress helpers).
- `TRACE(s)` is `-` for `v1.3.1` because that version does not emit the modern internal timing line.

| REF    | MODE    | WALL(s) | TRACE(s) | N    | NEG | LOST | END |
|--------|---------|---------|----------|------|-----|------|-----|
| HEAD   | CANFLUX | 30.064  | 24.799   | 1024 | 497 | 105  | 422 |
| HEAD   | BOOZER  | 14.450  | 9.530    | 1024 | 497 | 108  | 419 |
| HEAD   | MEISS   | 46.233  | 38.118   | 1024 | 497 | 108  | 419 |
| HEAD   | ALBERT  | 34.089  | 23.556   | 1024 | 497 | 164  | 363 |
| v1.5.1 | CANFLUX | 28.665  | 23.539   | 1024 | 474 | 128  | 422 |
| v1.5.1 | BOOZER  | 18.696  | 13.954   | 1024 | 474 | 126  | 424 |
| v1.5.1 | MEISS   | 48.254  | 39.767   | 1024 | 474 | 126  | 424 |
| v1.5.1 | ALBERT  | 34.051  | 23.479   | 1024 | 474 | 181  | 369 |
| v1.4.2 | CANFLUX | 29.024  | 23.840   | 1024 | 474 | 128  | 422 |
| v1.4.2 | BOOZER  | 18.792  | 14.081   | 1024 | 474 | 126  | 424 |
| v1.4.2 | MEISS   | 46.948  | 38.781   | 1024 | 474 | 126  | 424 |
| v1.4.2 | ALBERT  | 34.014  | 23.371   | 1024 | 474 | 180  | 370 |
| v1.3.1 | CANFLUX | 25.433  | -        | 1024 | 474 | 127  | 423 |
| v1.3.1 | BOOZER  | 16.050  | -        | 1024 | 474 | 127  | 423 |

## Automation

Run the reproducible cross-version benchmark driver:

`tools/perf/bench_versions.sh --threads 32 --particles 1024 --trace-time 1d-2 --npoiper 100 --npoiper2 128 --integmode 1`

Notes:
- The script creates temporary worktrees under `/tmp` for tags and applies print-suppression patches there (so older tags are not dominated by progress I/O).
- The script always uses `make clean && make -j32` for builds.
