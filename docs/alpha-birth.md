# Fusion-weighted alpha births

Set `startmode = 7` to draw equal-weight 3.5 MeV alpha markers from
`n_D(s) n_T(s) <sigma v>(T_rel(s)) dV`. Pitch `v_parallel/v` is uniform on
`[-1,1]`. Both angles span the full torus and are weighted by the local
absolute volume Jacobian. This mode uses the existing `&profiles` group,
also when collisions are disabled. Species 1 is D and species 2 is T;
densities are in m^-3 and temperatures are in eV. For Maxwellian species,
`T_rel = (3 Ti1 + 2 Ti2)/5`; equal ion temperatures give the usual `T_i`.

For example, add these settings to a VMEC or Boozer chartmap run:

```fortran
&config
startmode = 7
num_surf = 0
ntestpart = 100000
fusion_grid = 128, 32, 32
facE_al = 1
n_d = 4
n_e = 2
/
&profiles
profile_type = 'two_power'
ni1_scale = 5d19
ni2_scale = 5d19
ni1_p1 = 1
ni1_p2 = 1
ni2_p1 = 1
ni2_p2 = 1
Ti1_scale = 1d4
Ti2_scale = 1d4
/
```

`num_surf = 0` selects `[0,1]`, `1` selects `[0,sbeg(1)]`, and `2` selects
`[sbeg(1),sbeg(2)]` in normalized toroidal flux. Bounds must be ordered and
inside `[0,1]`. Unlike mode 5, this mode includes the core below `s=0.01`.
The axis itself is never a quadrature node. The reported power covers only
the selected radial interval and always covers the full torus.

VMEC runs sample in VMEC `(s,theta,phi)` with `magfie_vmec`'s Jacobian,
independently of the tracing coordinates. Boozer chartmap runs sample in
integration coordinates `(s,theta_B,zeta_B)` with the field-derived Boozer
Jacobian `sqrt(g) = torflux * (iota B_theta + B_zeta) / B^2`, then
convert to reference coordinates for `start.dat` (whose first coordinate is
`sqrt(s)` for these chartmaps). Jacobians in cm^3 are converted to m^3.
This volume element is physical only for Boozer-consistent chartmaps: the
field data and the chart geometry must describe the same equilibrium.
SIMPLE does not verify that consistency. In particular, the analytic test
chartmap's circular torus (`R0=1.5 m`, `a=0.5 m`) encloses about 2.8% less
volume than its field-derived Jacobian integrates to. Its sampling test
checks the field-derived measure, not the volume enclosed by that geometry.
SPECTRE's geometric radius and the normalized TEST field do not provide the
required physical flux coordinate and are rejected. Chartmap births require
`isw_field_type = 2`; compact GPU initialization is not supported.
Mode 5 and its output retain their existing behavior.

## Reactivity fit

The implementation uses H.-S. Bosch and G.M. Hale, *Improved formulas for
fusion cross-sections and thermal reactivities*, Nuclear Fusion **32** (1992)
611–631, [Eq. (12), Table VII](https://doi.org/10.1088/0029-5515/32/4/I07).
For D-T, `B_G = 34.3827`, `m_r c^2 = 1124656 keV`, and
`C_1 ... C_7 = (1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3,
1.35000e-2, -1.06750e-4, 1.36600e-5)`.
With `T` in keV:

```text
theta = T / [1 - T (C2 + T (C4 + T C6)) / (1 + T (C3 + T (C5 + T C7)))]
xi = [BG^2 / (4 theta)]^(1/3)
<sigma v> = C1 theta sqrt[xi / (mr c^2 T^3)] exp(-3 xi)
```

The last expression is in cm^3/s and is converted to m^3/s. The published
D-T fit range is 0.2–100 keV. Below 0.2 keV the formula is extrapolated to
allow cold profile edges; zero temperature gives zero source. No accuracy
claim is made for that extrapolation. Negative or nonfinite ion densities or
temperatures, an effective ion temperature above 100 keV, and zero total source
are errors.
`am1=2, am2=3, Z1=Z2=1` and the alpha parameters shown above are required.

## Resolution and marker normalization

`fusion_grid` gives the number of radial, poloidal and toroidal cells.
The source times the Jacobian is evaluated at each cell midpoint; cells are
drawn in proportion to their integrated source and positions are uniform
within the selected cell. Thus the sampled source is a piecewise constant
approximation. The power integral uses the same midpoint quadrature. Refine
all relevant grid dimensions until the radial distribution and total power
converge for the chosen equilibrium and profiles, especially for narrow peaks
or high angular harmonics. The default is `(128,32,32)`, with at most
8,388,608 cells and each dimension in `[1,1024]`.

Cell evaluations run in parallel with OpenMP after the VMEC or Boozer spline
tables are initialized. The CDF prefix sum remains ordered, so changing the
thread count preserves the cell weights and fixed-seed sampling.

`start.dat` retains five columns. The accompanying `alpha_source.dat` reports:

| Key | Meaning |
| --- | --- |
| `alpha_power_W` | `P_alpha = 3.5 MeV * integral(n_D n_T <sigma v> dV)` in W; uses SIMPLE's 3.5 MeV birth energy, not 3.52 MeV |
| `alpha_energy_MeV` | Energy per alpha used in the power normalization: 3.5 |
| `alpha_birth_rate_per_s` | Full birth rate in the selected interval |
| `marker_rate_per_s` | Birth rate divided by the number of sampled markers |
| `marker_power_W` | `P_alpha/N`, the common power represented by each marker |
| `marker_count` | Number of rows in `start.dat` |
| `s_bounds`, `fusion_grid` | Interval and quadrature resolution |

For heat-load normalization, sum `marker_power_W` over the relevant lost
markers (or divide by area for W/m^2). These are steady source rates, not
energies to divide by tracing time. For energy-resolved collisional loads,
multiply `marker_rate_per_s` by the marker's energy at impact. Every marker
has statistical weight `1/N`, so the ordinary loss fraction is already
source-weighted. Preserve this sidecar with the corresponding starts and
results; file-read or batch modes do not recompute its normalization.

## Verification and plot data

`make test TEST='test_fusion_'` runs fixed-seed chi-squared tests (20 bins,
19 degrees of freedom, threshold 65) for flat and peaked analytic profiles,
plus flat-profile tests on the downloaded VMEC equilibrium and the repository's
analytic Boozer chartmap. A separate VMEC case uses `n_D=5e19*(1-s) m^-3`,
`n_T=5e19 m^-3`, `T_D=T_T=10 keV` and three million markers. It asserts that
its expected chi-squared excess from omitting the Jacobian exceeds 65.
The first 100,000 draws exercise API output; the remaining draws stream into
the histogram without writing additional starts.
The peaked analytic case has densities proportional to `1-s`,
constant temperature and `dV/ds` proportional to `1+s`; its normalized source
has primitive `s - s^2/2 - s^3/3 + s^4/4` divided by `5/12`.
The tests also check angular moments, pitch, energy, fit values, independent
Simpson power quadrature, second-order grid convergence and invalid inputs.
VMEC bin volumes and total volume are checked against the independent wout
`vp` and `volume_p` data. In the test wout, `vp(1)` is padding and the
remaining entries are on the half mesh, with `dV/ds=4*pi^2*vp` in m^3.
VMEC histogram and power references use interpolated `vp`; chartmap
references use a separate quadrature of the same field-derived Jacobian.
OpenMP tests compare the complete CDF, power and fixed-seed histogram at
one and eight threads, and print three wall-clock CDF timings per thread count.

Small CSVs in `build/test/tests/fusion_flat/flat_radial.csv` and
`fusion_peaked/peaked_radial.csv` contain sampled counts and exact analytic
bin probabilities. `fusion_vmec/equilibrium_radial.csv` and
`fusion_chartmap/equilibrium_radial.csv` contain counts and reference volume
probabilities. `fusion_vmec_peaked/equilibrium_radial.csv` adds the peaked
VMEC case. Equilibrium CSVs also include probabilities without the Jacobian
and the reference bin-averaged `dV/ds` in m^3. Plot `count/sum(count)` against
the reference bin probabilities at the bin centers.
