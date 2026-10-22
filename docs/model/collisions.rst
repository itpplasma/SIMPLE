Collisions
==========

This page derives the Monte Carlo Coulomb collision operator applied with
``swcoll = .true.``: coefficients, discrete kicks, substep bounds and random
streams. On this page only, :math:`p = v/v_0`.
Usage: :doc:`/user-guide/running`.

Scope and splitting
-------------------

A test ion scatters on a Maxwellian background of electrons and
two ion species with densities :math:`n_a` and temperatures :math:`T_a` from
the ``&profiles`` namelist. The test mass is ``n_d`` times the proton mass
and its charge number is ``n_e``. Defaults remain alpha particles (4, 2).
For 55 keV hydrogen use ``n_d = 1``, ``n_e = 1`` and
``particle_energy_eV = 55000``; for deuterium use ``n_d = 2``.
In the formulas below :math:`m_\alpha` and :math:`Z_\alpha` denote this
configured test ion, and :math:`E_\alpha` is its reference kinetic energy.
Collisions cannot be combined with classification.

Each orbit microstep :math:`\Delta\tau` is followed by a collision step of
the same length at fixed position (``collide`` in ``src/simple_main.f90``)
that changes :math:`p` and :math:`\lambda`; ``update_momentum`` re-encodes
:math:`J_\perp`, :math:`\bar v_\parallel`, :math:`\bar p_\varphi`
(:doc:`guiding-centre`). A failed collision step rejects the microstep as a
numerical failure (:doc:`integrators`).

Coefficients
------------

For species :math:`a \in \{1, 2, e\}` let :math:`v_{T,a} = \sqrt{2T_a/m_a}`,
:math:`x_a = p\,v_0/v_{T,a}`, :math:`n_e = Z_1n_1 + Z_2n_2`, :math:`Z_e = 1`,
:math:`E_\alpha = m_\alpha v_0^2/2 = 3.5\,\text{MeV}/\texttt{facE\_al}`, and

.. math::

   \nu_a = \frac{2\pi Z_\alpha^2 Z_a^2 e^4 n_a \ln\Lambda_a}{m_\alpha^2 v_0^3\, v_{T,a}}
   \quad [\text{cm}^{-1}].

With the Chandrasekhar function
:math:`G(x) = [\operatorname{erf}x - x\operatorname{erf}'x]/(2x^2)`,
``onseff`` evaluates (series for :math:`x < 0.01`, asymptotes for
:math:`x > 6`)

.. math::

   d_p = \frac{2G(x)}{x}, \qquad
   d_\lambda = \frac{\operatorname{erf}x - G(x)}{x}, \qquad
   d_F = \frac1x\frac{d(x^2 d_p)}{dx} = \frac{4e^{-x^2}}{\sqrt\pi} - d_p ,

and ``coleff`` sums the species:

.. math::

   D_p = \sum_a \nu_a d_p(x_a), \qquad
   D_\lambda = \frac{1}{p^2}\sum_a \nu_a d_\lambda(x_a), \qquad
   F_p = \sum_a \nu_a\Bigl[\frac{d_F(x_a)}{p} - 2p\frac{E_\alpha}{T_a}d_p(x_a)\Bigr]
       = \frac{1}{p^2}\frac{\partial(p^2 D_p)}{\partial p}
         - 2p\sum_a\frac{E_\alpha}{T_a}\nu_a d_p(x_a).

:math:`F_p` and :math:`D_p` are the Itô drift and diffusion of :math:`p`, and
:math:`D_\lambda` the pitch diffusion, of the Fokker–Planck equation for
:math:`f(p,\lambda)` per unit velocity volume

.. math::

   \frac{\partial f}{\partial\tau} = \sum_a \frac{1}{p^2}\frac{\partial}{\partial p}
     \Bigl[p^2\nu_a d_p(x_a)\Bigl(\frac{\partial f}{\partial p}
       + 2p\frac{E_\alpha}{T_a}f\Bigr)\Bigr]
   + D_\lambda\frac{\partial}{\partial\lambda}\Bigl[(1-\lambda^2)\frac{\partial f}{\partial\lambda}\Bigr],

whose speed flux vanishes for each species on the Maxwellian
:math:`f\propto\exp(-m_\alpha v^2/2T_a)`.

The Coulomb logarithms (``compute_collision_coeffs``; :math:`n` in
cm\ :sup:`-3`, :math:`T`, :math:`E_\alpha` in eV, :math:`A_i` = ``am1``,
``am2``) are the NRL electron–ion form and the NRL mixed ion–ion form with
:math:`m_\alpha = A_\alpha m_p`, :math:`T_\alpha \to E_\alpha` and only the
background-ion density under the root:

.. math::

   \ln\Lambda_e = 24 - \ln\frac{\sqrt{n_e}}{T_e}, \qquad
   \ln\Lambda_i = 23 - \ln\Bigl[\frac{Z_\alpha Z_i(A_\alpha + A_i)}{A_\alpha T_i + A_iE_\alpha}
                  \sqrt{\frac{n_iZ_i^2}{T_i}}\Bigr].

``init_collision_profiles`` tabulates :math:`\nu_a`, :math:`v_0/v_{T,a}` and
:math:`E_\alpha/T_a` on 101 equidistant points in :math:`[0,1]`, with floors
:math:`T_a \ge 1` eV and :math:`n_a \ge 10^6` m\ :sup:`-3` (both ion
densities zero is vacuum). Each step interpolates linearly at the clamped
normalised toroidal flux :math:`s` of the marker, which ``profile_s``
(``src/profile_radius.f90``) computes from the tracing state:
:math:`z_1` itself for ``integ_coords`` -1, 0, 1 and 2; for 3 and 4 the
radius of the source field after ``integ_to_ref``; for 5 the reference
radius; for 6 the enclosed toroidal flux of the SPECTRE chart
:math:`\rho_g`, normalised to ``tflux(Mvol)``. A ``rho_tor`` chartmap
radius gives :math:`s = \rho^2`, and ``vmec_extended``/``vmec_to_wall``
charts give :math:`s = (\rho/\rho_\mathrm{lcfs})^2`. Radii without
toroidal-flux meaning (geometric or poloidal chartmaps) are accepted only
for radially uniform ``&profiles``.

Discrete step
-------------

``stost`` (mode 1) advances :math:`(p,\lambda)` over :math:`\Delta\tau` in
substeps :math:`\delta`, with coefficients at the current :math:`p`:

.. math::

   \lambda' = \lambda\cos\alpha + \sqrt{1-\lambda^2}\,\sin\alpha\cos\psi,
   \quad \alpha = \sqrt{4D_\lambda\delta},\ \psi\sim\mathcal U[0,2\pi);
   \qquad
   p' = p + F_p\delta + \sqrt{2D_p\delta}\,\xi,\quad \xi\sim\mathcal U[-\sqrt3,\sqrt3].

The pitch rotation has :math:`\langle\Delta\lambda\rangle =
-2\lambda D_\lambda\delta` and :math:`\langle\Delta\lambda^2\rangle =
2(1-\lambda^2)D_\lambda\delta` to leading order, the generator of the
Lorentz operator, and keeps :math:`|\lambda| \le 1` exactly; :math:`\xi` has
zero mean and unit variance. Before the random draws,

.. math::

   \delta = \min\Bigl(\text{remaining},\ \frac{f^2}{D_\lambda},\
            \frac{f\,p}{|F_p|},\ \frac{f^2p^2}{6D_p}\Bigr), \qquad f = 0.1,

which bounds the rotation angle by :math:`2f` and the drift and the largest
noise amplitude :math:`\sqrt{6D_p\delta}` by :math:`fp`. ``stost`` accepts
:math:`0 < f \le 0.25` (``max_substep_fraction``); ``simple.x`` uses 0.1. A
non-finite or negative coefficient, :math:`p \le 0`, a substep that no longer
advances time, or more than 100,000 substeps restores the initial state and
returns an error. In vacuum neither the state nor the generator changes.

Monoenergetic Lorentz mode
--------------------------

``collision_model = 'lorentz_nustar'`` (with ``swcoll = .true.``) replaces the
Coulomb operator by pitch-angle scattering at fixed :math:`p`. The collision
frequency is :math:`\nu = \nu_* |\iota| v_0 / R_0`, with ``nu_star_standard``
:math:`= \nu_*`, ``lorentz_iota`` :math:`= \iota` and
``lorentz_major_radius_cm`` :math:`= R_0`; all three must be set (:math:`\nu_*
\ge 0`, :math:`R_0 > 0`, :math:`\iota \ne 0`). With path length
:math:`\tau = v_0 t`, one call to ``stost_lorentz`` applies the Euler kick
:math:`\Delta\lambda = \sqrt{2 D (1-\lambda^2)\Delta\tau}\,\xi - 2\lambda D \Delta\tau`,
:math:`D = \nu/(2 v_0)`, :math:`\xi = \pm 1`, reflected at :math:`|\lambda| = 1`.
It has no substep control: keep :math:`\nu \Delta t` small. ``p_abs`` is
unchanged. The particle energy can be set directly with ``particle_energy_eV``
(default ``-1``: use :math:`3.5\,\mathrm{MeV}/` ``facE_al``); a nondefault
``facE_al`` that disagrees with an explicit energy is rejected. Mode and energy
are written to ``orbits.nc`` and ``results.nc`` as global attributes.

Random streams
--------------

With ``deterministic = .true.``, ``reset_seed_if_deterministic``
(``src/params.f90``) reseeds the generator before each marker from
``ran_seed`` and :math:`k` = ``particle_index_offset`` + source row:

.. math::

   s_0 = (\texttt{ran\_seed} + 104729\,k) \bmod (2^{31}-2) + 1, \qquad
   s_i = (48271\,s_{i-1} + 7919\,i) \bmod (2^{31}-2) + 1 .

A marker thus draws the same kicks for any thread count or contiguous
shard; distinct seeds do not prove statistical independence.

References
----------

* J. D. Huba, NRL Plasma Formulary, NRL/PU/6790--19-640, Naval Research
  Laboratory (2019)

Tabulated files and exterior plasma
-----------------------------------

Set ``profile_file = 'plasma.dat'`` in ``&profiles`` to read a radial table.
The first line is the row count. Each following line contains
``s Te Ti1 Ti2 ni1 ni2``, with temperatures in eV and ion densities in
m\ :sup:`-3`. Radii must increase strictly from zero to a last radius at
least one; every value must be finite and temperatures/densities nonnegative.
The two ion species are selected by ``am1``, ``Z1``, ``am2``, ``Z2`` in
``&config``. A main ion plus an impurity therefore supplies both ion
backgrounds, with electron density fixed by quasineutrality.

Tables override the analytic profiles. Linear interpolation is used between
rows, with endpoint values held outside the table. Rows beyond ``s = 1``
define an explicit exterior plasma. To model an exterior vacuum, extend the
table to the desired exterior radius and set both ion densities to zero
there. The collision grid has 101 points across the table's entire radial
range, so resolve steep edge gradients in a radial convergence study.
Without a table, the existing grid over ``[0, 1]`` and held edge coefficients
are preserved. There is no automatic thermalization stop; ``trace_time``
sets the observation window.
