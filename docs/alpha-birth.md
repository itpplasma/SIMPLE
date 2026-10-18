# Fusion-weighted alpha births

`startmode = 7` starts equal-weight 3.5 MeV alpha markers distributed like the
D-T fusion source, so that loss fractions and wall loads are weighted by where
alphas are born, and writes the power each marker represents.

## Setup

Add to a VMEC or Boozer-chartmap run:

```fortran
&config
startmode = 7
num_surf = 0               ! s in [0,1]; 1: [0,sbeg(1)]; 2: [sbeg(1),sbeg(2)]
ntestpart = 100000
fusion_grid = 128, 32, 32  ! midpoint cells in s, poloidal, toroidal
/
&profiles
profile_type = 'two_power'
ni1_scale = 5d19, ni1_p1 = 1, ni1_p2 = 1   ! deuterium, m^-3
ni2_scale = 5d19, ni2_p1 = 1, ni2_p2 = 1   ! tritium, m^-3
Ti1_scale = 1d4                            ! eV
Ti2_scale = 1d4
/
```

Profile syntax: {doc}`/user-guide/running`; `&profiles` is read even without
collisions. Required: species 1 = D, 2 = T (`am1 = 2`, `am2 = 3`,
`Z1 = Z2 = 1`, the defaults), the default alpha (`n_d = 4`, `n_e = 2`,
`facE_al = 1`), a VMEC file or a Boozer chartmap with `isw_field_type = 2`,
and `gpu_compact_init = .false.`.

## Source

The birth density per volume is

$$
S(s) = n_D(s)\, n_T(s)\, \langle\sigma v\rangle(T_\text{rel}(s)),
\qquad T_\text{rel} = \tfrac{3}{5} T_D + \tfrac{2}{5} T_T ,
$$

the Maxwellian relative temperature for D-T. Markers are drawn from
$S\,\sqrt{g}\;\mathrm{d}s\,\mathrm{d}\vartheta\,\mathrm{d}\varphi$ over the
selected $s$ interval and the full torus, with uniform pitch
$\lambda\in[-1,1]$ and $v/v_0 = 1$. VMEC runs use VMEC coordinates and the
VMEC Jacobian. Boozer-chartmap runs use Boozer coordinates with
$\sqrt{g} = \psi_a(\iota B_\vartheta + B_\varphi)/B^2$, $\psi_a$ the edge
toroidal flux over $2\pi$ (`torflux`), and convert to
reference coordinates for `start.dat` (radius $\rho=\sqrt{s}$). This
Jacobian is physical only if the chartmap geometry and field describe the same
equilibrium, which SIMPLE does not check. Unlike `startmode = 5`, the core
$s<0.01$ is included.

`fusion_grid` (default `128, 32, 32`; each entry in $[1,1024]$, at most
8,388,608 cells) sets the cells. $S\sqrt{g}$ is evaluated at each cell
midpoint, a cell is drawn with probability proportional to its integral and
the position is uniform within it. Refine until the radial distribution and
the power converge. Results do not depend on the thread count.

## Reactivity

$\langle\sigma v\rangle$ is the Bosch–Hale fit (Eq. 12, Table VII) with $T$ in
keV, $B_G = 34.3827\ \text{keV}^{1/2}$, $m_r c^2 = 1124656$ keV and
$C_{1..7} = (1.17302\times10^{-9}, 1.51361\times10^{-2}, 7.51886\times10^{-2},
4.60643\times10^{-3}, 1.35000\times10^{-2}, -1.06750\times10^{-4},
1.36600\times10^{-5})$:

$$
\theta = \frac{T}{1 - \dfrac{T(C_2 + T(C_4 + T C_6))}{1 + T(C_3 + T(C_5 + T C_7))}},
\qquad
\xi = \left(\frac{B_G^2}{4\theta}\right)^{1/3},
\qquad
\langle\sigma v\rangle = C_1\,\theta\,\sqrt{\frac{\xi}{m_r c^2\, T^3}}\,e^{-3\xi}
\ \text{cm}^3/\text{s}.
$$

Implemented in `src/fusion_source.f90`. The fit is valid for 0.2–100 keV;
below 0.2 keV it is extrapolated without accuracy claim ($T=0$ gives zero),
above 100 keV and for negative or non-finite profiles SIMPLE stops.

## Normalisation: `alpha_source.dat`

`start.dat` keeps its five columns. `alpha_source.dat` holds:

| Key | Meaning |
| --- | --- |
| `alpha_power_W` | $P_\alpha = E_\alpha \int S\,\mathrm{d}V$ in W, with $E_\alpha = 3.5$ MeV (not 3.52 MeV) |
| `alpha_energy_MeV` | $E_\alpha = 3.5$ |
| `alpha_birth_rate_per_s` | $\int S\,\mathrm{d}V$ over the selected interval |
| `marker_rate_per_s` | birth rate $/N$ |
| `marker_power_W` | $P_\alpha/N$, power carried by each marker |
| `marker_count` | $N$, rows of `start.dat` |
| `s_bounds`, `fusion_grid` | interval and cells |

Every marker has weight $1/N$, so the plain loss fraction is already
source-weighted. Lost alpha power is `marker_power_W` times the number of lost
markers (a steady rate; do not divide by `trace_time`). With collisions, use
`marker_rate_per_s` times each marker's energy at loss. Read-back modes
(`startmode` 2, 4) do not recompute `alpha_source.dat`; keep it with its
`start.dat`.

Reference: H.-S. Bosch, G. M. Hale, Nucl. Fusion 32, 611 (1992),
[doi:10.1088/0029-5515/32/4/I07](https://doi.org/10.1088/0029-5515/32/4/I07).
