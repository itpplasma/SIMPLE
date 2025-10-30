Concept Integrator Prototypes
=============================

Created: 2019-03-07  
Author: Christopher Albert <albert@alumni.tugraz.at>

This directory collects the original Python proof-of-concept scripts for
guiding-center symplectic integrators described in [1].  The examples implement
axisymmetric guiding-center motion in canonicalized flux coordinates and were
the foundation for the production SIMPLE code.

Common functionality
--------------------

* ``common.py`` – derivative computation and Newton iteration helpers
* ``field_test.py`` – simple magnetic configuration for testing
* ``plotting.py`` – helper routines for plotting

Integrator implementations
--------------------------

* ``expl_impl_euler.py`` – explicit-implicit Euler method
* ``expl_impl_euler_optim.py`` – optimized explicit-implicit Euler method with
  user-supplied Jacobian
* ``expl_impl_euler_r_theta.py`` – explicit-implicit Euler in ``(r, theta)`` coordinates
* ``impl_expl_euler.py`` – implicit-explicit Euler method
* ``impl_midpoint.py`` – implicit midpoint rule
* ``impl_pendulum.py`` – nonlinear pendulum example
* ``verlet.py`` – Stoermer-Verlet / leap-frog method
* ``rk.py`` – adaptive Runge-Kutta 4/5 reference implementation (SciPy)
* ``pauli_particle.py`` – Pauli particle dynamics example
* ``discrete_variational.py`` – discrete variational integrator notebook helpers

Additional material
-------------------

* ``exportfig.py`` – export utilities for the published figures
* ``simple.in`` – sample SIMPLE configuration used by the prototypes
* ``iaea2019/`` – scripts and notebooks prepared for the 2019 IAEA meeting

Reference
---------

[1] C. G. Albert, S. V. Kasilov, and W. Kernbichler, *Symplectic integration with
non-canonical quadrature for guiding-center orbits in magnetic confinement
devices*, Mar. 2019, arXiv:1903.06885. Submitted to J. Comput. Phys.
