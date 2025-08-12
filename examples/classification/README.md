class_parts.dat columns:

1: particle index
2: starting radial coordinate s
3: proportional to perpendicular adiabatic invariant (magnetic moment)
4: result of j_parallel orbitclassifier
5: result of topological (ideal) orbit classifier
6: fractal dimension classification (regular=1, chaotic=2)

Prompt losses: 0
Regular / ideal: 1
Chaotic / non-ideal: 2

according to classification.f90, check_orbit_type.f90 and simple_main.f90

tcut is used only for fractal dimension classification. Earlier it was recommended to use at t=0.1s. The new classifiers (topological, and j_parallel) work at t=0.015s.

Recipe: 

0) Minimize fraction of chaotic orbits. (10-100 times faster than usual run)
1) Look at s=0.25 and s=0.5. See at which J_perp particles are chaotic on s=0.25 and optimize such that at these J_perp values at s=0.5 there are regular regions. (10-100 times faster than usual run)
3) Optimize for loss fraction, tracing only chaotic orbits to the end. (only 1-10 times faster than usual run)

TODO:

Enable mode 3) also for new classifiers.