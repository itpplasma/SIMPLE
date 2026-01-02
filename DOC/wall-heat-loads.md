# Wall Heat Load Physics

This document describes the physics behind the wall heat load calculation in `pysimple.engineering`.

## D-T Fusion Energy Balance

Each deuterium-tritium fusion reaction produces 17.6 MeV:

$$\text{D} + \text{T} \rightarrow \text{He}^{4} (3.5\,\text{MeV}) + n (14.1\,\text{MeV})$$

- **Alpha particle** ($\alpha$): 3.5 MeV kinetic energy
- **Neutron**: 14.1 MeV (escapes plasma, heats blanket)

The alpha power is therefore:

$$P_\alpha = \frac{P_\text{fusion}}{5}$$

For a 3 GW fusion reactor: $P_\alpha \approx 600$ MW.

## Alpha Particle Slowing Down

Alphas are born at $E_0 = 3.5$ MeV and slow down via Coulomb collisions with background electrons and ions. The characteristic slowing-down time is:

$$\tau_s = \frac{3(2\pi)^{3/2} \epsilon_0^2 m_\alpha}{e^4 Z_\alpha^2 n_e \ln\Lambda} \left(\frac{T_e}{m_e}\right)^{3/2}$$

For typical reactor parameters ($n_e = 10^{20}$ m$^{-3}$, $T_e = 10$ keV), $\tau_s \approx 0.1-1$ s.

The energy evolution follows:

$$E(t) = E_0 \left(1 - t/\tau_s\right)^{2/3} \quad \text{for } t < \tau_s$$

## Monte Carlo Heat Flux Calculation

SIMPLE traces $N$ test particles representing the alpha population. Each particle carries a fraction of total alpha power:

$$\text{Power per MC particle} = \frac{P_\alpha}{N}$$

### Energy Weighting

Particles that escape early (before significant slowing down) carry more energy than those that escape late. The normalized velocity at wall impact is:

$$p = \frac{v}{v_0} = \sqrt{\frac{E}{E_0}}$$

where $v_0$ is the birth velocity. The energy carried by particle $i$ at wall impact is:

$$E_i = p_i^2 \cdot E_0$$

Without collisions: $p = 1$ for all particles.
With collisions: $p < 1$ for particles that slowed down before escaping.

### Heat Flux Calculation

For a wall bin with area $A_\text{bin}$, the energy-weighted heat flux is:

$$q_\text{bin} = \frac{1}{A_\text{bin}} \sum_{i \in \text{bin}} p_i^2 \cdot \frac{P_\alpha}{N}$$

Or equivalently:

$$q_\text{bin} = \frac{\sum_{i \in \text{bin}} p_i^2}{N} \cdot \frac{P_\alpha}{A_\text{bin}}$$

### Energy Loss Fraction

The **particle loss fraction** is:

$$f_\text{particle} = \frac{n_\text{lost}}{N}$$

The **energy loss fraction** accounts for slowing down:

$$f_\text{energy} = \frac{1}{N} \sum_{i \in \text{lost}} p_i^2$$

Note: $f_\text{energy} \leq f_\text{particle}$ because $p_i^2 \leq 1$.

The lost power is:

$$P_\text{lost} = f_\text{energy} \cdot P_\alpha$$

## Typical Values

| Configuration | Peak Heat Flux | Reference |
|--------------|----------------|-----------|
| Infinity Two | ~2.5 MW/m$^2$ | Albert et al. (2024) |
| ARIES-CS | 5-18 MW/m$^2$ | Ku & Boozer (2008) |
| W7-X (experiments) | 0.1-1 MW/m$^2$ | Lazerson et al. (2021) |
| Material limits | 5-10 MW/m$^2$ | Engineering constraint |

## Code Implementation

The `WallHeatMap` class in `pysimple.engineering` implements this calculation:

```python
from pysimple.engineering import WallHeatMap

# For a 3 GW fusion reactor (600 MW alpha power)
heat_map = WallHeatMap.from_netcdf(
    "results.nc",
    total_alpha_power_MW=600.0,  # Required input
    major_radius=1000.0,  # cm
    minor_radius=100.0,   # cm
)

print(f"Particle loss fraction: {heat_map.loss_fraction:.1%}")
print(f"Energy loss fraction: {heat_map.energy_loss_fraction:.1%}")
print(f"Lost power: {heat_map.lost_power:.1f} MW")
print(f"Peak flux: {heat_map.peak_flux:.2f} MW/m^2")
```

## Data Format

The NetCDF results file from SIMPLE contains:

| Variable | Shape | Description |
|----------|-------|-------------|
| `zend` | (5, N) | Final phase space: s, theta, zeta, p=v/v0, lambda |
| `xend_cart` | (3, N) | Final Cartesian position (x, y, z) |
| `class_lost` | (N,) | 1 if particle hit wall, 0 otherwise |
| `times_lost` | (N,) | Time of wall impact (s) |

The velocity ratio `zend[3]` = $p = v/v_0$ is crucial for energy weighting.

## References

1. **Lazerson, S.A., et al.** (2021). "Modeling of energetic particle transport in optimized stellarators." *Plasma Physics and Controlled Fusion* 63, 125033.
   https://doi.org/10.1088/1361-6587/ac2b38

2. **Albert, C.G., et al.** (2024). "Alpha particle confinement in Infinity Two." *Journal of Plasma Physics*.
   https://doi.org/10.1017/S0022377824000163

3. **Ku, L.P. & Boozer, A.H.** (2008). "Stellarator alpha particle loss calculations for ARIES-CS." *Fusion Science and Technology* 54, 673-693.
   https://doi.org/10.13182/FST08-A1891

4. **Drevlak, M., et al.** (2014). "Fast particle confinement with optimized coil currents in the W7-X stellarator." *Nuclear Fusion* 54, 073002.
   https://doi.org/10.1088/0029-5515/54/7/073002

5. **BEAMS3D code documentation**: https://princetonuniversity.github.io/STELLOPT/BEAMS3D
