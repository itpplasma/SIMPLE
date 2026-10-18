Field sources and coordinates
=============================

This page tells you which files to give SIMPLE for each kind of magnetic
field and which integration coordinates work with it. Keys belong to the
``&config`` namelist of ``simple.in`` (:doc:`/reference/input-parameters`).

Method: :doc:`/model/coordinates`.

Input files
-----------

``field_input`` (default ``''``)
   The magnetic field: VMEC ``wout`` (``.nc``), Boozer chartmap (``.nc``),
   SPECTRE (``.h5``) or a coil file. The type is detected from the file.
``coord_input`` (default ``''``)
   The reference coordinates in which start points and results are given:
   a VMEC ``wout`` file or a chartmap. Empty means ``netcdffile``, then
   ``field_input``.
``netcdffile`` (default ``''``)
   VMEC ``wout`` file; fallback for both keys above. Where a VMEC file is
   required and all three are empty, ``wout.nc`` is used.

Integration coordinates
-----------------------

``integ_coords`` selects the coordinates in which orbits are integrated. Its
older alias ``isw_field_type`` (default 0) applies when ``integ_coords`` is
unset (default ``-1000``).

=====  ===========  ======================================  ==============================
Value  Name         Coordinates                             Use
=====  ===========  ======================================  ==============================
0      CANFLUX      :math:`(s, \vartheta_c, \varphi_c)`     canonical flux coordinates
                                                            from VMEC; symplectic
1      VMEC         :math:`(s, \vartheta, \varphi)`         VMEC; ``integmode = 0`` only
2      BOOZER       :math:`(s, \vartheta_B, \varphi_B)`     Boozer, from VMEC or a Boozer
                                                            chartmap; symplectic; standard
                                                            for nested flux surfaces
3      MEISS        :math:`(\rho, \vartheta_c, \varphi_c)`  canonical without flux
                                                            surfaces, e.g. coils; symplectic
4      ALBERT       :math:`(\psi, \vartheta_c, \varphi_c)`  as MEISS with
                                                            :math:`\psi \propto A_\vartheta`;
                                                            symplectic
5      REFCOORDS    reference coordinates                   ``field_input`` splined on
                                                            ``coord_input``; ``integmode = 0``
6      SPECTRE      :math:`(\rho_g, \vartheta, \zeta)`      per-volume SPECTRE/SPEC
-1     TEST         :math:`(r, \vartheta, \varphi)`         analytic tokamak for tests
=====  ===========  ======================================  ==============================

VMEC equilibria
---------------

.. code-block:: fortran

   &config
   netcdffile = 'wout.nc'
   isw_field_type = 2      ! Boozer
   /

SIMPLE splines VMEC's :math:`R`, :math:`Z` and :math:`\lambda` over :math:`s`
(order ``ns_s``, default 5) and the angles (order ``ns_tp``, default 5, on
:math:`m_\text{max}\,` ``multharm`` :math:`+1` poloidal and
:math:`n_\text{max}\,` ``multharm`` :math:`+1` toroidal points, ``multharm``
default 5), transforms to Boozer or canonical flux coordinates and splines the
canonical field components. ``axis_healing`` selects how the VMEC harmonics
are continued to the axis; use ``'polyfit'``, which replaces them below
``s_axis_heal`` (default 0.01) by a polynomial in :math:`s`.

``vmec_B_scale`` and ``vmec_RZ_scale`` (default 1) multiply the field and
the device size, for scans without a new equilibrium. The scaled field is an
MHD equilibrium only with consistently scaled pressure.

Start positions (``sbeg``, ``start.dat``) are given in VMEC coordinates
:math:`(s, \vartheta, \varphi)`.

Boozer chartmaps (GVEC, booz_xform)
-----------------------------------

A Boozer chartmap is a NetCDF file with geometry and Boozer field components
on a uniform grid; no VMEC file is needed:

.. literalinclude:: ../../examples/simple_chartmap.in
   :language: fortran
   :caption: examples/simple_chartmap.in

``sbeg`` stays the normalised toroidal flux :math:`s`; positions in
``start.dat`` and the outputs use the chartmap radius :math:`\rho=\sqrt{s}`.
Format, signs and scaling: :doc:`/boozer-chartmap-schema`. Converters:

.. code-block:: bash

   # GVEC (gvec Python package; installed by setup-venv.sh)
   python tools/gvec_to_boozer_chartmap.py parameter_final.ini State_final.dat \
       boozer_chartmap.nc
   # booz_xform
   python tools/booz_xform_to_boozer_chartmap.py boozmn_case.nc boozer_chartmap.nc

Both take ``--nrho``, ``--ntheta`` and ``--nphi`` (GVEC) or ``--nzeta``
(booz_xform). The GVEC converter reverses one angle (``--flip tor``, default,
or ``pol``) for SIMPLE's left-handed convention; its tested workflow is
``test/tests/run_gvec_qa_roundtrip.py``.

Coil fields
-----------

Coil fields follow particles outside the last closed flux surface (LCFS).
The coil file is text: the number of points, then :math:`x, y, z` (m) and
current (A) per line. Points :math:`i` and :math:`i+1` are joined by a
straight segment carrying the current of point :math:`i`; a zero current ends
a coil. Convert STELLOPT/MAKEGRID files with

.. code-block:: bash

   python test/tests/convert_coils_stellopt_to_simple.py coils.c09r00 coils.simple

``python -m libneo.convert_coils_to_simple`` drops the zero-current points and
so joins consecutive coils by spurious segments; do not use it for SIMPLE.

Coil fields need ``integ_coords = 3`` or ``4``, or ``integ_coords = 5`` with
``integmode = 0``. A VMEC file still provides start positions and reference
coordinates unless ``coord_input`` is a chartmap:

.. code-block:: fortran

   &config
   netcdffile = 'wout.nc'
   field_input = 'coils.simple'
   isw_field_type = 3      ! Meiss
   /

SIMPLE does not check that coils and equilibrium describe the same field.

Chartmaps beyond the LCFS
-------------------------

Wall losses need reference coordinates that continue beyond the LCFS. libneo
writes such a chartmap from a VMEC file:

.. code-block:: bash

   export PYTHONPATH=build/_deps/libneo-src/python
   python -m libneo.chartmap_cli from-vmec-extended wout.nc chartmap_wall.nc \
       --nrho 33 --ntheta 65 --nzeta 33 --boundary-offset 0.15

The global attribute ``rho_lcfs`` records the LCFS radius.
``chart_boundary_kind`` (default ``'auto'``) declares what the outer chartmap
boundary is: ``'lcfs'`` (loss, exit code 1), ``'wall'`` (exit code 2) or
``'domain'`` (numerical exit 101); ``'auto'`` means ``'lcfs'`` if
``rho_lcfs = 1``, else ``'domain'``. The subcommands ``from-stl`` (chartmap
inside an STL boundary) and ``from-vmec-to-wall`` (VMEC extended to an STL
wall) need the Python packages ``map2disc`` and ``shapely``.

SPECTRE equilibria
------------------

SPECTRE/SPEC multi-region relaxed MHD equilibria (HDF5) require
``integ_coords = 6``. SIMPLE integrates guiding-centre orbits inside each
volume and full orbits near the interfaces, where the field jumps:

.. literalinclude:: ../../examples/simple_spectre.in
   :language: fortran
   :caption: examples/simple_spectre.in

``sbeg`` is the SPECTRE radial label :math:`\rho_g \in [0, M_\text{vol}]`,
or :math:`s` with ``spectre_sbeg_is_toroidal_flux = .true.`` (``startmode = 1``). The interface
controls ``spectre_fo_entry_larmor``, ``spectre_fo_tolerance`` and
``spectre_fo_grid`` are explained in :doc:`/model/extensions`.

Each volume builds its own Meiss canonical coordinates on a grid of
``spectre_ncon_r`` (default 48) by ``spectre_ncon_th`` (48) by
``spectre_ncon_phi`` points. With ``spectre_ncon_phi = -1`` (default) SIMPLE
uses 32 points, or 8 for an axisymmetric field (all toroidal mode numbers
zero); a positive value is used as given. Memory scales with
:math:`M_\text{vol}\,n_r n_\vartheta n_\varphi`; the log line
``spectre_construction_grid`` reports the grid in use.
