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

## Cause of the switch chatter at `npoiper2=256`

Short runs (64 markers, 200 us, entry 1.5) with per-episode exit diagnostics
and the new bulk-recovery counter (`sympl_bulk_recovery_stats`) show that the
chatter is not interface switching. Most short episodes start with
`iface = 0` deep inside a volume. They are FO recoveries after the GC
microstep exhausts its 200-iteration cap (255 of 277 recoveries; 22 persistent
teleports, no Newton failures). In those microsteps every rejection is the
`step_teleported` energy-shell guard, `|H - p^2| > 5e-4 p^2`: midpoint at 256
oscillates up to that band in some regions, so any substep raising the defect
is rejected, the substep shrinks to 1e-4–1e-7 of the microstep and the cap is
reached. FO recovery then resets the shell and exits almost at once through the
legacy rule. The guard's absolute band is deliberate (source comment), so it
is unchanged.

| `npoiper2` | Bulk recoveries | FO entries | Pure-GC core-s/particle-s | Median shell error |
| ---: | ---: | ---: | ---: | ---: |
| 256 | 277 | 348 | 133 | 1.7e-4 |
| 512 | 5 | 74 | 216 | 4.4e-5 |
| 1024 | 1 | 71 | 429 | 1.1e-5 |

`npoiper2=512` removes the chatter at half the GC cost of 1024. Because FO
dominates, the hybrid total changes by only about 10%.

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
case.

With 1,024 matched NCSX markers (`summary_1024.json`; FO loses 418):

| Entry | Losses | Final loss fraction minus FO, 95% bootstrap CI | Loss-curve max diff | FO time | Speedup vs FO |
| ---: | ---: | --- | ---: | ---: | ---: |
| 1 | 417 | -0.001 [-0.014, 0.013] | 0.007 | 10.6% | 7.1 |
| 3 | 423 | +0.005 [-0.010, 0.019] | 0.009 | 28.3% | 3.3 |

Both differ from FO in about 50 paired statuses (chaotic noise). G3V3L3Fi with
256 markers (`summary_g3v3_256.json`; FO loses 53): entries 1 and 1.5 match all
FO statuses at 5.3x and 4.0x; entry 3 differs only in marker 129 at 2.3x.
Marker 129 was not a selection criterion and remains phase-sensitive (entry 6
also loses it, entry 3 does not).

## Invariant gate and the new default

Loss statistics at this ensemble size do not resolve switch errors. The native
`test_spectre_fo_interfaces_rk45` does: in axisymmetric tok2vol it bounds the
relative spread of canonical p_zeta over all FO phases of a marker at 1e-5.

| Entry | FO phases | Max p_zeta spread | Test |
| ---: | ---: | ---: | --- |
| 3 | 24 | 4.5e-8 | pass |
| 2 | 27 | 7.6e-8 | pass |
| 1.5 | 30 | 8.3e-8 | pass |
| 1 | 38 | 1.2e-4 | fail |

Below about 1.5 radial gyroradii the gyro-orbit at entry reaches the sheet,
where the guiding-centre approximation no longer holds, and the switch error
rises by three orders of magnitude. The production default
`spectre_fo_entry_larmor` is therefore 1.5 (was 3): NCSX FO time 15% instead
of 30% and cost 1.8x lower; G3V3L3Fi 2.0x lower. With the default 1.5, the
affected native tests pass 24/25; `test_spectre_construction_grid` fails on
macOS before and after this change because it reads Linux `/proc` VmHWM.

## Consequences

- A gyro-orbit within about one gyroradius of a sheet intersects it; FO is
  physically required there. The invariant gate places the usable limit of the
  present switch at 1.5 gyroradii, now the default.
- A phase-carrying switch (GC1+) could at most move entry from 1.5 towards 1,
  worth about 1.3x on NCSX. That no longer justifies it on cost.
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
