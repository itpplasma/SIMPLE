# Decoded LCFS reporting baseline

These Albert coils records use the existing published NCSX fixture. The equilibrium
comes from the [simsopt test data](https://github.com/hiddenSymmetries/simsopt/blob/master/tests/test_files/wout_c09r00_fixedBoundary_0.5T_vacuum_ns201.nc)
and the coils from the [STELLOPT NCSX example](https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00).
The golden runner's in-tree converter retains every filament point, including
zero-current closing points.

The raw reference is SIMPLE c75e5295dcf8af1382ec2985ce748da4cd079919 built in
Release with deterministic floating-point settings and libneo
266f0a3c844315522d8b620fa4b96b69aedfeafc. An independent reference-only observer
reproduced start.dat, times_lost.dat and orbit_exit_code.dat byte for byte before
deriving observations. It classified all 17 reference LCFS losses: seven full
Euler-map events and ten already physical RK-active events. No candidate values
selected the event cohort, corrected cells or expected values.

For the seven Euler events, a field-value bisection of the canonical momentum
constraint recovers the physical radius at the advanced angles without using the
endpoint decoder's Newton derivatives. Physical momentum and pitch are derived at
both bracket endpoints, and an independent LCFS bisection supplies their reporting
fraction. Canonical step duration is dtaumin/sqrt(2), while the reported full
microstep duration is dtaumin/v0. Actual partial maps are outside this observation
contract. The observer changes no reference integration state or energy/recovery
decision, and RK-active records retain their raw physical observations.

Relative to those raw Release tables, 26 cells in times_lost.dat and three exit
times change. The corrected exit times are markers 18, 27 and 30. Every other
Release token is retained, including the already passing time for marker 11.
These files replace earlier Fast-context copies, so the whole-file change also
includes the raw Fast-to-Release context difference. The comparison retains
rtol=1.5e-7 and atol=1e-12. Meiss coils, historical physical-momentum baselines,
diagnostics and performance checks retain their existing contracts.

Selection requires exactly the audited SIMPLE reference above, unchanged raw
Release tables and published inputs, and actual Release caches with deterministic
settings and compiled libneo 266f0a3c844315522d8b620fa4b96b69aedfeafc in both builds.
The selector rejects other profiles, reference revisions and libneo consumer
contexts. Qualification of another context requires separate evidence.
