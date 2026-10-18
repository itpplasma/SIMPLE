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
