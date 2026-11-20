# cyl2vol test fixture

Two-volume axisymmetric torus in the cylinder limit (`Nfp=1`, `Ntor=0`,
`Mpol=4`, `Lrad=[8, 6]`, circular boundary `R0=100`, `a=1`, so `a/R0=0.01`)
generated with SPECTRE commit `0d8e45c4ed5ab807ec78b18222ce736b5084e2bd`:

    python concept/spectre_cylinder/make_input.py test/tests/data/tok2vol.toml cyl2vol.toml
    python scripts/run/calc_spectre_eq.py cyl2vol.toml   # SPECTRE; writes cyl2vol_opt.h5 -> cyl2vol.h5

`cyl2vol.toml` in this directory is the exact input.

## Design

In a straight cylinder each volume holds a Lundquist field, `curl B = mu B`:
`B_z = c_J J0(mu r) + c_Y Y0(mu r)`, `B_theta = c_J J1(mu r) + c_Y Y1(mu r)`,
with `c_Y = 0` in the core. The input was designed for

- core `r < 0.5`: `mu = 1.5 /m`, `B_z(0) = 5 T`; at the sheet
  `|B| = 4.660 T`, pitch `atan(B_theta/B_z) = 0.384 rad`;
- outer volume: `mu = 0.5 /m`; at the sheet `|B| = 1.02 x 4.660 T`, pitch
  `0.100 rad`, so the sheet current turns the field by 0.284 rad and raises
  its magnitude by 2 %;
- pressure step `p = [0.4387873309, 0]` (SPECTRE units, `mu0 p` in T^2)
  balancing `[[p + B^2/2]] = 0` at `r = 0.5`.

The fluxes follow from these constants (`make_input.py`, `lundquist.py` in
`concept/spectre_cylinder`): toroidal flux of volume k
`phiedge (tflux_k - tflux_{k-1})`, poloidal flux through a ribbon of length
`2 pi R0` in the `+theta` sense `phiedge (pflux_k - pflux_{k-1})`, where theta
is right-handed about `e_phi` (the boundary is `R = R0 + cos t`, `Z = -sin t`).
The sign of the poloidal flux was fixed by SPECTRE's interface force: only it
reproduces the designed pressure balance. `pflux_1` is unused in the core.

## Result

Force-minimized MRxMHD equilibrium, interface force residual `rms 1.2e-11`,
Beltrami field error `max 7.0e-4`. The interface is
`R = 100.0028 + 0.50004 cos t + 0.0038 cos 2t + ...`, to first order a
circle of radius 0.5 m whose centre is shifted outward by 6.6 mm. The
Lundquist field built from the fixture's own fluxes and m = 1 interface radius
reproduces SPECTRE's `|B|` within 0.96 % and its pitch within 5.2 mrad
everywhere (`test_spectre_cylinder_field`), i.e. to `O(a/R0)`; the same input
at `R0 = 50` gives 1.97 % and 10.4 mrad.

SHA-256:
- `cyl2vol.h5`: `07d36b08e64f35f1db10a3901e2c27b8a248ecd48b59e171786627040f286481`
