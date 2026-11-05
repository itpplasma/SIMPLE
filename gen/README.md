Axis map generator
==================

Use the same compiler as the installed FortSym. From the SIMPLE repository root:

```sh
cmake -S gen -B /scratch/axis-generator -G Ninja \
  -DCMAKE_PREFIX_PATH=/path/to/installed/fortsym
cmake --build /scratch/axis-generator --parallel 8
/scratch/axis-generator/generate_axis_maps src/generated
fprettify -c .fprettify src/generated/axis_map_generated.f90
```

For an existing FortSym build with compatible module files, the same small
generator can also be compiled against its public module/library interface:

```sh
gfortran -I/path/to/fortsym/build/include gen/generate_axis_maps.f90 \
  -L/path/to/fortsym/build/lib -Wl,-rpath,/path/to/fortsym/build/lib \
  -lfortsym -o /scratch/generate_axis_maps
/scratch/generate_axis_maps src/generated
fprettify -c .fprettify src/generated/axis_map_generated.f90
```

This emits `src/generated/axis_map_generated.f90`, module `axis_map_generated`,
with six scalar subroutines (each prefix is `axis_flux` or `axis_rho`):

- `_forward(q,theta,x,y)`;
- `_inverse(x,y,q)`;
- `_velocity(x,y,dq,dtheta,dx,dy)`.

Commit the emitted file after review; ordinary SIMPLE builds need no CAS.
Regenerate into an empty temporary directory and compare the module to check
reproducibility. No generated Fortran is handwritten.

Flux mode uses R=sqrt(q); rho mode uses R=q. The caller recovers theta with
atan2(y,x). Velocity derives its radial factor from the symbolic inverse-map
Jacobian, using R²=x*x+y*y. The caller passes the other RHS components unchanged.

Domains: both velocity routines require strictly positive radius. Both forward and
inverse maps extend to zero. Exact-axis evaluation must use the existing regular reconstruction,
never evaluate singular polar derivatives. Callers validate finite/nonnegative
inputs and select the mode from the active native coordinate contract. These
kernels do not infer chart semantics or set numerical thresholds. A native
radial floor q0 corresponds to Cartesian radius squared q0 (flux) or q0*q0
(rho). Use that converted radius squared throughout reconstruction and blending.

Independent verification belongs in the consuming tests: pull back the same
smooth translated/rotating Cartesian vector field through q=R² and q=R,
including nonzero translation, and compare to its analytic Cartesian RHS.
Symbolic differentiation and generation alone are not that behavioral oracle.

Event energy generator
----------------------

The same generator CMake project also builds `generate_event_energy`.
Run it with `src/generated` as the argument, then format
`src/generated/generated_event_energy.f90` with the repository formatter.
The emitted scalar conversion is E = E0*p^2, where E0 is the effective reference
kinetic energy in eV and p is speed divided by the reference speed. It does not
supply source normalization or change the orbit state.

For a compatible existing FortSym build:

```sh
gfortran -I/path/to/fortsym/build/include gen/generate_event_energy.f90 \
  -L/path/to/fortsym/build/lib -Wl,-rpath,/path/to/fortsym/build/lib \
  -lfortsym -o /scratch/generate_event_energy
/scratch/generate_event_energy src/generated
fo fmt src/generated/generated_event_energy.f90
```

For reproducibility, regenerate in an empty directory, format the result and
compare it with the checked-in file. Ordinary SIMPLE builds need no CAS.

Implicit-scheme residuals and Jacobians
---------------------------------------

The same generator CMake project builds `generate_orbit_systems`. Run it with
`src/generated` as the argument; it writes the `euler1` and midpoint kernels to
`src/generated/orbit_systems_euler1.inc` and
`src/generated/orbit_systems_midpoint.inc`, all others to
`src/generated/orbit_systems_kernels.inc`, and the module
`orbit_systems_generated` that includes all three to
`src/generated/orbit_systems_generated.f90`. `orbit_symplectic_euler1.f90` and
`simple_gpu.f90` include the `euler1` and midpoint kernels in their own scope:
called across modules, they made the OpenACC steps slower than the hand-written
systems with nvfortran 26.9, compiled in place they are as fast.

Each implicit symplectic scheme (`euler1`, `euler2`, `midpoint`, `gauss1` to
`gauss4`, `lobatto3`) is stated once as a residual in the canonical field
quantities H, p_theta, v_par, h_theta and h_phi, which are opaque functions of
the stage states. Their partial derivatives bind to the members of
`field_can_t` (the packed Hessian order of `get_derivatives2`), and every
`<scheme>_jacobian` is the symbolic derivative of `<scheme>_residual` with
respect to the Newton unknowns. `gc_rates_canonical` gives the guiding-centre
rates from the same field quantities; `gc_rates_canonical_jacobian` differentiates
those rates for the two-derivative integrators. All kernels carry `!$acc routine seq`, so
the CPU build and the OpenACC device code share them. The Runge-Kutta
coefficients stay arguments, supplied by `coeff_rk_gauss` and
`coeff_rk_lobatto`.

For a compatible existing FortSym build, which also provides FortGen:

```sh
gfortran -I/path/to/fortsym/build/include \
  -I/path/to/fortsym/build/_deps/fortgen-build/include \
  gen/generate_orbit_systems.f90 \
  -L/path/to/fortsym/build/lib -Wl,-rpath,/path/to/fortsym/build/lib \
  -lfortsym -lfortgen -o /scratch/generate_orbit_systems
/scratch/generate_orbit_systems src/generated <fortsym revision>
```

The same run writes two OpenACC kernels. `gpu_radial_mixed_second_derivatives`
(`src/generated/gpu_derivatives2_radial_mixed.inc`, included in
`simple_gpu.f90`) states v_par = (p_phi - A_phi/ro0)/h_phi,
H = v_par^2/2 + mu B + Phi_e and p_theta = h_theta v_par + A_theta/ro0 once and
gives the second derivatives with a radial or p_phi index, which the Euler
Newton system reads, from the values and gradients that `get_derivatives`
stored in `field_can_t`. v_par is a named intermediate of `fortsym_chain`, so
its derivatives reuse v_par and its lower derivatives. d2Phi_e is dropped
because the electric potential is linear in s. `gc_rates_boozer`
(`gc_rates_boozer.inc`, called by `gpu_rhs_canonical`) applies the
guiding-centre rates of `gc_rates_canonical` to a Boozer field; p_phi -
A_phi/ro0 stays a named difference, so no rewrite expands it where v_par is
close to zero. The CUDA rates in `simple_cuda_native_rk54.cu` stay
hand-written: the same kernel emitted for CUDA ran about 2 % slower there.
`emit_field_derivatives` can also emit `get_val`, `get_derivatives` and
`get_derivatives2`; the CPU chain stays hand-written because a different
rounding changes the chaotic orbits of the golden-record classifier cases.

The same run writes `src/generated/cuda_euler1_boozer.cuh`, the device code
of the direct-CUDA Euler step (`simple_cuda_euler.cu`). FortGen's Horner
emitter gives the quintic spline cell kernels in libneo's coefficient layout.
`euler1_boozer_system` takes the raw spline values of A_phi(s),
B_vartheta(rho), B_varphi(rho) and |B|(rho, vartheta, varphi) with
s = sigma r and rho = sqrt(s), where sigma = -1 is the axis continuation, and
returns the Euler residual, its Jacobian and the field quantities of the angle
advance. h = B/|B|, v_par, H and p_theta are chain intermediates, the residual
is the statement shared with `euler1_residual`, and its Jacobian is the chain
derivative. `euler1_boozer_advance` moves p_theta and the field quantities of
the last Newton evaluation to the converged point to first order and advances
the angles, as `sympl_euler1_extrapolate_field` and
`sympl_euler1_advance_angles` do on the CPU.

`gc_geometry_gradient_transport` transports the gradients of A_theta,
A_phi, h_theta, h_phi and B to the corrected Newton stage using the packed
Hessians already stored in `field_can_t`. Gauss endpoint Hermite extrapolation
uses those transported gradients, without another field evaluation.

The optional second argument records the FortSym revision in the banner of
every kernel. The emitted files are committed unformatted, so regenerating into
an empty directory with that FortSym revision reproduces it byte for byte.
`test_orbit_systems_generated` compares the kernels with the hand-derived
formulas they replaced on random field values.
