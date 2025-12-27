# Non-SoA CPU performance (OpenMP) cross-version comparison

This note documents one focused performance check: **non-SoA** orbit tracing on CPU
with OpenMP, comparing the current `perf/restore-v131` branch against `main` and
older release tags where the same field mode exists.

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
| BOOZER | 14.298  | 9.492    | 1024 | 497 | 108  | 419 |
| CANFLUX| 26.134  | 20.872   | 1024 | 497 | 105  | 422 |
| MEISS  | 45.896  | 37.751   | 1024 | 497 | 108  | 419 |
| ALBERT | 34.056  | 23.623   | 1024 | 497 | 164  | 363 |

## Cross-version comparison (same benchmark case)

Numbers below are from a single run of:

`tools/perf/bench_versions.sh --refs HEAD,main,v1.5.1,v1.4.2,v1.3.1 --modes boozer,canflux,meiss,albert --threads 32 --particles 1024 --trace-time 1d-2 --npoiper 100 --npoiper2 128 --integmode 1`

Important:
- Older tags printed progress unconditionally. For apples-to-apples timing, the script **patches those worktrees in `/tmp`** to disable progress output (prints are replaced by constant `.false.` guards or `return` in progress helpers).
- `TRACE(s)` is `-` for `v1.3.1` because that version does not emit the modern internal timing line.

| REF    | MODE    | WALL(s) | TRACE(s) | N    | NEG | LOST | END |
|--------|---------|---------|----------|------|-----|------|-----|
| HEAD   | BOOZER  | 14.298  | 9.492    | 1024 | 497 | 108  | 419 |
| HEAD   | CANFLUX | 26.134  | 20.872   | 1024 | 497 | 105  | 422 |
| HEAD   | MEISS   | 45.896  | 37.751   | 1024 | 497 | 108  | 419 |
| HEAD   | ALBERT  | 34.056  | 23.623   | 1024 | 497 | 164  | 363 |
| main   | BOOZER  | 18.498  | 13.669   | 1024 | 497 | 108  | 419 |
| main   | CANFLUX | 35.255  | 29.892   | 1024 | 497 | 105  | 422 |
| main   | MEISS   | 46.145  | 37.998   | 1024 | 497 | 108  | 419 |
| main   | ALBERT  | 34.397  | 23.751   | 1024 | 497 | 164  | 363 |
| v1.5.1 | BOOZER  | 18.820  | 14.063   | 1024 | 474 | 126  | 424 |
| v1.5.1 | CANFLUX | 28.635  | 23.497   | 1024 | 474 | 128  | 422 |
| v1.5.1 | MEISS   | 48.514  | 40.204   | 1024 | 474 | 126  | 424 |
| v1.5.1 | ALBERT  | 34.868  | 24.208   | 1024 | 474 | 181  | 369 |
| v1.4.2 | BOOZER  | 18.781  | 14.094   | 1024 | 474 | 126  | 424 |
| v1.4.2 | CANFLUX | 28.993  | 23.897   | 1024 | 474 | 128  | 422 |
| v1.4.2 | MEISS   | 48.080  | 39.358   | 1024 | 474 | 126  | 424 |
| v1.4.2 | ALBERT  | 34.237  | 23.625   | 1024 | 474 | 180  | 370 |
| v1.3.1 | BOOZER  | 16.076  | -        | 1024 | 474 | 127  | 423 |
| v1.3.1 | CANFLUX | 25.575  | -        | 1024 | 474 | 127  | 423 |

## Automation

Run the reproducible cross-version benchmark driver:

`tools/perf/bench_versions.sh --threads 32 --particles 1024 --trace-time 1d-2 --npoiper 100 --npoiper2 128 --integmode 1`

Notes:
- The script creates temporary worktrees under `/tmp` for tags and applies print-suppression patches there (so older tags are not dominated by progress I/O).
- The script always uses `make clean && make -j32` for builds.
