IMAS workflow
=============

``pysimple.imas`` reads IMAS IDS (equilibrium, fast-ion markers, bulk plasma,
wall), follows the markers as guiding centres in the frozen field and writes a
``distributions`` IDS with their initial and final states.

Method: :doc:`/model/guiding-centre`.

Installation and versions
-------------------------

``pip install -e '.[imas]'`` adds IMAS-Python (``>=2.3,<3``), Data Dictionary
(DD) 4.1.1, ``netCDF4``, ``scipy`` and libneo's equilibrium converter. Input IDS
of DD 3.x and 4.x are read; all IDS of one run must have the same DD version.
``run`` writes the DD version of its marker IDS; ``run_entry`` and
``simple-imas`` read and write ``dd_version`` (default 4.1.1) and let
IMAS-Python convert.

Running
-------

The three entry points do the same run. From a database entry (``.nc`` paths
use IMAS-Python's NetCDF backend, ``imas:`` URIs any other backend):

.. code-block:: bash

   simple-imas input.nc output.nc --duration 1e-3 --source-interval 0.01 --core-profiles

.. code-block:: python

   from pysimple.imas import run_entry

   out = run_entry("input.nc", "output.nc", duration=1e-3, source_interval=0.01,
                   core_profiles=True)

From IDS objects in memory, here with markers built in Python:

.. code-block:: python

   import imas
   import numpy as np
   from pysimple.imas import Markers, run, write_markers

   alphas = Markers(np.array([[10.5, 0.0, 0.0, 3.5e6, 0.3]]),  # R, phi, Z, E, pitch
                    np.array([1e18]), mass_number=4, nuclear_charge=2, charge=2,
                    rate=True)                                # weights in particles/s
   with imas.DBEntry("input.nc", "r") as entry:
       eq = entry.get("equilibrium")
   out = run(write_markers([alphas]), equilibrium=eq, duration=1e-3,
             source_interval=0.01)
   q = out.distribution[0].global_quantities
   print("confined fraction", q[1].particles_n / q[0].particles_n)

===========================  ==========================================================
``duration`` [s]             tracing time (``trace_time``); required
``source_interval`` [s]      required for ``distribution_sources``, forbidden otherwise
``time`` [s]                 initial time; default: first marker slice. Every IDS needs
                             a slice at exactly this time
``marker_ids``               ``distribution_sources`` (default) or ``distributions``;
                             ``occurrence`` (0) selects the marker occurrence
``core_profiles``            IDS object (``run``) or ``True`` to read it; collisions
                             are on whenever it is given
``conversion``               libneo resolution, integers ≥ 8: ``nrho`` (50), ``ntheta``
                             (48), ``nstep`` (3600), ``nlabel``, ``ntheta_int`` (500),
                             ``nsurfmax`` (10000), ``nsurf`` (2000)
other keywords               ``&config`` keys (:doc:`/reference/input-parameters`);
                             CLI: a JSON object in ``--settings``
===========================  ==========================================================

The actor owns ``n_d``, ``n_e``, ``trace_time``, ``ntestpart``, ``facE_al``,
``orbit_model``, the field, wall and profile inputs; setting them is an error.
An existing output entry is never overwritten. ``run(..., field="wout.nc")``
(or a chartmap) replaces the equilibrium IDS by a SIMPLE field file, e.g. a
stellarator, without ``core_profiles``. Repeat a run with doubled
``conversion`` resolution to check the converted equilibrium.

Inputs and conventions
----------------------

====================================================  ====================================================
IDS field (DD 4 names)                                SIMPLE quantity
====================================================  ====================================================
**equilibrium** ``time_slice[]``
``.profiles_2d[].grid.dim1``, ``.dim2``               R, Z [m]; uniform, ≥ 6 points; exactly one
                                                      grid with ``grid_type.index = 1``
``.profiles_2d[].psi``                                [Wb], shape (nR, nZ)
``.profiles_1d`` ``psi``, ``f``, ``pressure``, ``q``  [Wb], [T m], [Pa]; resampled to nR uniform
                                                      :math:`\psi` points from axis to boundary
``.global_quantities.psi_axis``, ``ip``, ...          G-EQDSK header [Wb], [A], [m]
``.boundary.outline`` (closed)                        LCFS [m]
``vacuum_toroidal_field.r0``                          :math:`R_0` [m]; :math:`B_0=F(\psi_\text{axis})/R_0`
**distribution_sources** ``source[]`` or              one population each
**distributions** ``distribution[]``
``.gyro_type``; ``.species.type.index``               2 (guiding centre); 2 or 3 (ion)
``.species.ion.element[0].a``                         ``n_d``, mass :math:`A\,m_p`
``.species.ion.z_ion`` (or ``.state.z_min``)          ``n_e``, charge
``.markers[]`` identifiers 4, 5, 3, 301, 403          :math:`R` [m], :math:`\varphi` [rad], :math:`Z` [m],
                                                      :math:`E` [eV], :math:`\lambda`; any column order
``.markers[].weights``                                :math:`w_i`: particles/s or particles
**core_profiles** ``profiles_1d[]``
``.grid.rho_tor_norm``, exactly 0 to 1                :math:`s=\rho^2`; linear interpolation in :math:`s`
``.electrons.temperature``                            :math:`T_e` [eV]
``.ion[i].temperature``, 1 or 2 ions                  :math:`T_i` [eV]
``.ion[i].density_thermal`` (or ``density``)          :math:`n_i` [m\ :sup:`-3`]
``.ion[i].element[0].a``, ``.z_ion``                  ``am1``/``am2``, ``Z1``/``Z2``
``.electrons.density``                                must equal :math:`\sum_i Z_i n_i`
====================================================  ====================================================

DD 3 stores equilibria in COCOS 11, DD 4 in COCOS 17. SIMPLE writes the COCOS 3
G-EQDSK that libneo reads (``python/pysimple/imas/equilibrium.py``):

.. math::

   \psi_\text{G} = \sigma\,\psi_\text{IDS}/2\pi,\quad
   q_\text{G} = -q_\text{IDS},\quad F_\text{G} = F_\text{IDS},\qquad
   \sigma = -1\ (\text{DD 3}),\ +1\ (\text{DD 4}),

with :math:`FF'` and :math:`p'` differentiated with respect to
:math:`\psi_\text{G}`. The directions of :math:`B_\varphi` and
:math:`B_\text{pol}` are kept, so identifier 403,
:math:`\lambda = \mathbf v\cdot\mathbf B/(vB)`, needs no sign change. libneo
turns the G-EQDSK into a Boozer chartmap with radius
:math:`\rho=\sqrt{\Phi/\Phi_\text{edge}}`, the IMAS ``rho_tor_norm``; orbits are
integrated in Boozer coordinates :math:`(s,\vartheta_B,\varphi_B)`. The inverse
chart map takes each marker from :math:`(R,\varphi,Z)` to
:math:`(\rho,\vartheta,\zeta)`, and :math:`v/v_0=\sqrt{E/E_0}` with
:math:`E_0 = 3.5` MeV, :math:`v_0=\sqrt{2E_0/(A\,m_p)}`.

Output and weights
------------------

``distributions.distribution[k]`` holds population :math:`k` with
``markers[0]`` at :math:`t_0` and ``markers[1]`` at :math:`t_0+T`
(:math:`T` = ``duration``), both with identifiers 4, 5, 3, 301, 403;
``markers[1]`` adds the private identifiers −1 (end time [s]) and −2 (exit code:
0 confined, 1 LCFS, 2 wall; :doc:`outputs`). A lost marker keeps its state at
the loss. The chartmap made from an equilibrium IDS ends at the LCFS
(``rho_lcfs = 1``) and SIMPLE queries the wall only beyond ``rho_lcfs``
(:doc:`wall-heat-loads`), so these markers end at the LCFS with code 1; code 2
needs a ``field=`` chartmap that extends to the wall.
Wall-intersection tracing is not yet supported for IDS equilibria: a wall IDS
(``wall=``, ``--wall``) raises an error before tracing, so omit it. The weights are

.. math::

   N_i = w_i\,\Delta t\ \text{(sources)},\qquad N_i = w_i\ \text{(distributions)},
   \qquad N_i' = N_i\,[\text{exit}_i = 0],

with :math:`\Delta t` = ``source_interval``; ``global_quantities[0]`` and
``[1]`` hold :math:`\sum_i N_i` and :math:`\sum_i N_i'` (``particles_n``) and
the matching kinetic energies :math:`\sum_i N_i E_i` [J] (``energy``).

The output is a finite cohort: all :math:`N_i` particles start together at
:math:`t_0` in the field and plasma of :math:`t_0`. It is not the steady-state
distribution of a continuous source. For a constant source with losses as the
only sink, :math:`\sum_i w_i\,(t_{\text{end},i}-t_0)` is a lower bound of the
steady-state particle number, reached when every marker is lost before
:math:`t_0+T`.

Rejected input
--------------

* Equilibria without exactly one uniform rectangular :math:`(R,Z)` grid (3D
  or GGD-only equilibria), open boundaries, non-monotonic ``profiles_1d.psi``.
  Use ``field=`` for 3D fields.
* Particle markers (``gyro_type = 1``), electrons, neutrals, multi-element
  species, non-integer charges, populations at different times.
* ``core_profiles`` with markers other than alphas (:math:`A=4`, :math:`Z=2`),
  because the collision operator uses alpha mass and charge
  (``src/collis_alphas.f90``); more than two bulk ions, rotation, nonzero
  ``density_fast``, non-quasineutral densities; ``core_profiles`` with ``field=``.
* Markers outside the chartmap, markers that SIMPLE skips or ends with a
  numerical exit code, and final states outside :math:`0\le\rho\le1` or
  :math:`|\lambda|\le1`: the run raises an error rather than count them as lost.

ASCOT5 compatibility
--------------------

SIMPLE reads the IDS fields that ASCOT5's ``a5py`` IMAS interface reads for
guiding-centre markers (identifiers 4, 5, 3, 301, 403; weights in particles/s),
axisymmetric equilibria (rectangular ``profiles_2d``, the same COCOS 11/17 to
COCOS 3 flux mapping). Both take
:math:`\lambda` relative to :math:`\mathbf B` and :math:`\varphi` in radians
from the IDS. Differences:

* ASCOT5 also reads particle markers and 3D equilibria; guiding-centre markers
  get gyrophase zero.
* ASCOT5 uses the nuclear mass and charge ``z_n``; SIMPLE uses :math:`A\,m_p`
  (0.7 % heavier for alphas) and ``z_ion``.
* ASCOT5 maps ``core_profiles`` to :math:`\rho_\text{pol}`, uses one ion
  temperature for all ions and accepts rotation.
* ASCOT5 reads the wall IDS; SIMPLE rejects it (see above).
* ASCOT5 writes histogram moments to ``distributions.profiles_2d``; SIMPLE
  writes marker snapshots.
