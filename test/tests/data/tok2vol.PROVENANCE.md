# tok2vol test fixtures

Two axisymmetric two-volume tokamaks (`Nfp=1`, `Ntor=0`, `Mpol=4`,
`Lrad=[8, 6]`, circular boundary `R0=10`, `a=1`, pressure `[1e-3, 0]`)
generated with SPECTRE commit `f51638520d2dbdee9b0b4e4174ace5b2349b7460`
through `scripts/run_case.sh` of `plasma/proj/spectre-orbits`:

    run_case.sh tok2vol.toml eq 2             # writes tok2vol_opt.h5 -> tok2vol.h5
    run_case.sh tok2vol_sheet.toml field 1    # writes tok2vol_sheet.h5

The `.toml` files in this directory are the exact inputs.

## tok2vol.h5

Force-minimized MRxMHD equilibrium with a current-carrying core
(`lconstraint=0`, Beltrami `mu=[0.07, 0]`, poloidal flux `0.225` in the
outer volume). The inner interface is Shafranov-shifted and mildly shaped
(`R = 10.299 + 0.449 cos(theta) + ...`). Interface force residual
`max 5.9e-9`, Beltrami field error `max 7.4e-4`. The interface carries a
sheet current, so the field direction jumps across it. Used by the orbit,
loss, crossing-map and validation tests.

## tok2vol_sheet.h5

Current-free field (`lconstraint=3` with zero currents) on a circular
interface at `r=0.5` with toroidal flux `0.2475`, giving a pure-magnitude
jump `[[B]]/B` up to `1.1e-2` outward and no direction jump. It is a solved
Beltrami field, not a force-balanced equilibrium (interface force residual
`max 2.1e-4`). Used by the tests of the pure-magnitude sheet models: the
CA15 relocation and the full-orbit interface entry.

## Chart validity

The SPEC interface blend must not fold: the Jacobian of
`(rho_g, theta, zeta) -> x` keeps one sign in every volume for both files
(checked on 4096 x 400 points per volume). The previous `tok2vol.h5`
(`lconstraint=3` force minimization of a current-free tokamak) had no
equilibrium with a pressure step: the minimizer squeezed the inner interface
until the volume-2 chart folded near `theta = pi/2, 3 pi/2`, where the field
direction then reversed across the interface.

SHA-256:
- `tok2vol.h5`: `bdd829513f9660e6ca27641a5be789fba869eb04d2a4fbcdb13ac9ab75d49072`
- `tok2vol_sheet.h5`: `e5f25b82d7d068ff9eb989293c994447a3b8e1a830907629ec60e72af16f47f5`
