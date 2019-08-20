Created: 2019-03-07
Author:  Christopher Albert <albert@alumni.tugraz.at>

This dataset contains code examples for different symplectic integrators
with non-canonical quadrature points described in [1]. Here
guiding-center motion is implemented in its axisymmetric variant for 
tokamak magnetic fields in canonicalized flux coordinates.
All implementations except "expl_impl_euler_optim.py" use
SciPy library routines with approximate Jacobian for root-finding.

Common functionality
* common.py: computation of derivatives, and Newton iterations
* field_test.py: simple magnetic configuration for testing
* plotting.py: helper routines for plotting

Integrators:
* expl_impl_euler.py: Explicit-implicit Euler method
* expl_impl_euler_optim.py: Optimized explicit-implicit Euler method with user-supplied Jacobian
* impl_expl_euler.py: Implicit-explicit Euler method
* impl_midpoint.py: Implicit midpoint rule
* verlet.py: Stoermer-Verlet/leap-frog method
* rk.py: Adaptive Runge-Kutta 4/5 from SciPy library

[1] C. G. Albert, S. V. Kasilov, and W. Kernbichler, 
   Symplectic integration with non-canonical quadrature for guiding-center orbits in magnetic confinement devices, 
   Mar. 2019, arXiv:1903.06885. Submitted to J. Comp. Phys
