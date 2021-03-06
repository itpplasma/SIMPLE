# NEO-ORB

Physical time $t$ is represented by a normalized time parameter
$$
\tau = v_0 t\,,
$$
where
$$
v_0 = \sqrt{\frac{2T_\alpha}{m_\alpha}}\label{test}
$$
is the thermal velocity corresponding to the temperature $T_\alpha$ of the particles to trace (or energy in case of monoenergetic particles). Integration variables are
$$
\begin{align}
z^{1,2,3}&=x^{1,2,3}\,, \\
z^4 &= \frac{p}{m_\alpha v_0}\,, \\
z^5 &= \frac{p_\parallel}{p}\,. \\
\end{align}
$$
In the non-relativistic limit $z^4=1$ is constant over the whole orbit and $z^5=v_\parallel/v_0$.

To normalize length scales a reference Larmor radius $\rho_0$ is used with
$$
\rho_0 = v_0 \omega_c^{\mathrm{ref}}=\sqrt{\frac{2T_\alpha}{m_\alpha}}\frac{mc}{e B^{\mathrm{ref}}}.
$$
Time-stepping for RK4/5 uses the routine `velo_can` internally inside `sub_alpha_lifetime.f90`. For symplectic integrators quantities are rescaled in `init_sympl` in `neo_orb.f90`. Symplectic integrators are non-relativistic and work in four variables 
$$
\begin{align}
z_{\mathrm{sy}}^{1,2,3} &= x^{1,2,3}\,, \\
z_{\mathrm{sy}}^4 &= p_\varphi.
\end{align}
$$
The routine  `timestep_sympl_z` does the conversion back in every time-step such that the same variables $z$  as for the RK4/5 integrator can be used.

## Implementation

There is a trapping parameter
$$
\bar{\mu} = \left((1-v_\parallel^2/v^2)\frac{B_\max}{B}-1\right)\frac{B_\min}{B_\max-B_\min}
= \left(\frac{\mu}{\mu_{\mathrm{tp}}} - 1\right)B_\min/(B_\max-B_\min)
$$

For deeply trapped
$$
\mu_\mathrm{dt} = \frac{mv^2}{B_\mathrm{min}}
$$
So
$$
\bar{\mu} = \left(\frac{\mu}{\mu_{\mathrm{tp}}} - 1\right)/(mv^2\mu_\mathrm{dt}B_\max-1)
$$

$$
J_\perp =(1-\lambda^2)/B
$$

$$
(\eta B_\mathrm{max} - 1) B_\mathrm{min}/(B_\max-B_\min)
$$

$$
\eta = \frac{v_\perp^2}{v^2 B}
$$

