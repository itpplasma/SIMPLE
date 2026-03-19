# Classification Examples

Reproduces classification results in the style of the JPP 2020 paper:
Albert, Kasilov, Kernbichler, "Accelerated methods for direct computation
of fusion alpha particle losses within stellarator optimization",
J. Plasma Phys. 86, 815860201 (2020).

Uses a QH (quasi-helical) equilibrium from Landreman & Paul 2021.

## Input files

- `simple_s03.in` — Single surface at s=0.3, 10000 particles, trace_time=1s
- `simple_s06.in` — Single surface at s=0.6, 10000 particles, trace_time=1s
- `simple_volume.in` — Volume sampling, 100000 particles, trace_time=0.015s
- `simple.in` — Legacy volume config (5000 particles)

## Running

    make run        # runs all three cases under /tmp/simple_classification

Or manually:

    bash run_all.sh

## Plotting

    make plot       # generates paper-style plots from results

Scripts:
- `plot_paper_results.py` — Paper-style plots (Fig 6, Fig 8, Table 1)
- `plot_classification.py` — Legacy (s, J_perp) volume classification

## Output files

`class_parts.dat` columns:

1. particle index
2. starting radial coordinate s
3. proportional to perpendicular adiabatic invariant (magnetic moment)
4. result of j_parallel orbit classifier
5. result of topological (ideal) orbit classifier
6. fractal dimension classification (regular=1, chaotic=2)

Classification codes: 0=prompt loss, 1=regular/ideal, 2=chaotic/non-ideal.

`times_lost.dat` columns:

1. particle index
2. loss time (negative if confined, trace_time if regular)
3. trapping parameter theta_trap (Eq. 3.1 in paper)
4. starting s
5. perpendicular invariant J_perp
6-10. final state zend(1:5)
