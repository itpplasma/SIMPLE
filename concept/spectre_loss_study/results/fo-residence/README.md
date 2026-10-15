# Hybrid scheduling, GC step, FO residence and entry distance, 2026-10-15

NCSX fixture, 3.5 keV alphas, 256 matched volume starts, 2 ms, 16 threads,
one Apple M3 Ultra (load average 9–20 on 28 cores). The source is `a79b3cf`
plus the study-driver changes committed with this record. The FO-everywhere
reference loses 101/256 markers. Cost is core-seconds per simulated
particle-second (see `../boozer-cost/`). `summary_256.json` holds every case;
`residence_report.py` regenerates it from the run directories.

## Scheduling

The study driver used OpenMP's default static schedule, which caused the
42–45% thread efficiency previously reported. Production `trace_parallel`
already uses `schedule(dynamic, 1)`. With the same schedule in the driver,
hybrid wall time falls from 191 s to 94 s and marker outcomes are bitwise
identical. The hybrid's earlier wall-clock timings therefore overstate
production cost by about 2x; normalized costs are unaffected.

## GC step and integrator

| Hybrid (entry 3) | Losses | Status diff vs FO | Loss-curve max diff | Core-s/particle-s | Pure-GC core-s/particle-s | FO entries |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| midpoint, `npoiper2=1024` (study) | 98 | 9 | 0.016 | 3,862 | 390 | 979 |
| midpoint, 256 (example input) | 98 | 11 | 0.016 | 3,550 | 106 | 11,678 |
| Euler1, 1024 | 100 | 9 | 0.027 | 7,609 | 290 | 128,252 |
| Euler1, 256 | 100 | 11 | 0.023 | 8,661 | 116 | 1,296,041 |

Midpoint at 256, as in `examples/simple_spectre.in`, makes the GC part 3.7x
cheaper without a resolvable change in loss statistics. The hybrid total
improves only 8% because FO dominates. The coarser step makes the switch
chatter at interface 2 (4,118 episodes, median 0.06 us, plus 6,514 that enter
and leave within one step), which costs little FO time. Euler1 chatters far more
and doubles total cost; it is unsuitable for the SPECTRE hybrid.

## Where FO time goes

The driver now writes `study_episodes.csv`: one row per FO episode with
entry/exit time, accepted FO time, particle rho at entry and end kind. It
captures 893 of 979 counted entries; the remainder are exits and re-entries
within one microstep, merged into one episode. All logged FO time matches the
per-marker FO clock.

Midpoint, 1024, entry 3:

| Episode duration | Episodes | Share of FO time |
| --- | ---: | ---: |
| < 1 us | 104 | 0.0% |
| 1–10 us | 175 | 1.2% |
| 10–100 us | 500 | 16.3% |
| 0.1–1 ms | 79 | 19.1% |
| > 1 ms | 35 | 63.4% |

FO time is 70% at interface 1, 30% at interface 2, 0.1% at the outer wall
and none at the axis. Episodes still active at the 2 ms cutoff hold 65% of FO
time. 33 markers spend more than 90% of their life in FO; 26 markers hold half
of all FO time. **FO cost is dominated by particles whose orbits reside within
the entry/exit band of a sheet, not by crossings or wall hits.**

## Entry distance

Exit is at entry + 1 radial gyroradius. Current phase-discarding hybrid:

| Entry [radial gyroradii] | Losses | Status diff vs FO | Loss-curve max diff | FO time | Core-s/particle-s | Speedup vs FO |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.01 | 104 | 17 | 0.027 | 2.5% | 700 | 16.6 |
| 0.5 | 102 | 11 | 0.035 | 7.6% | 1,312 | 8.9 |
| 1 | 102 | 13 | 0.023 | 10.2% | 1,591 | 7.3 |
| 1.5 | 101 | 12 | 0.031 | 15.3% | 2,178 | 5.3 |
| 2 | 99 | 10 | 0.020 | 19.2% | 2,648 | 4.4 |
| 3 (default) | 98 | 9 | 0.016 | 30.0% | 3,862 | 3.0 |
| 4 | 104 | 9 | 0.016 | 38.2% | 4,797 | 2.4 |

FO time grows roughly linearly with the entry distance, as expected for a
resident population spread across the band. At 256 markers the binomial
uncertainty of the loss fraction is about 0.03, so the loss curves for entries
0.5–4 are statistically indistinguishable from FO. The earlier NCSX FO
default-versus-refined-grid pair differed in four statuses, a lower bound on
chaotic per-marker noise. Entry 0.01 differs in 17 statuses, the most of any
case. A 1,024-marker comparison of FO, entry 1 and entry 3 is the next
recorded result.

## Consequences

- A gyro-orbit within about one gyroradius of a sheet intersects it; FO is
  physically required there. Between about 1 and 3 gyroradii GC is valid in
  principle, and the cost difference is a factor of 2.4 (entry 1 versus 3).
- The phase-carrying switch (GC1+) is justified if it permits entry near
  1 gyroradius with loss agreement that the current switch cannot attain.
  If the current switch already agrees there, reduce the entry distance first.
- Residual FO of 2.5–10% of particle time remains intrinsic; it costs 4–8x
  pure GC.

## Reproduction

```sh
O=build/study/spectre-mac
python concept/spectre_loss_study/run.py --field ncsx --mode hybrid \
  --markers 256 --duration 2e-3 --threads 16 --output $O/hybrid
# then, each with --starts $O/hybrid/study_starts.csv:
#   --mode fo                         -> $O/fo
#   --entry E (0.01 0.5 1 1.5 2 4)    -> $O/entry_E
#   --npoiper2 256                    -> $O/ep_mid256
#   --integmode 1 [--npoiper2 256]    -> $O/hyb_eul1024, $O/hyb_eul256
python concept/spectre_loss_study/residence_report.py OUTPUT $O/fo CASE...
```

In `summary_256.json`, `gcprobe` is entry 0.01 and `ep_mid1024` is entry 3;
both also carry the episode breakdown.
