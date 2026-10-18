Design notes
============

The directory ``DOC/`` of the repository contains the maintained technical
notes of the developers. They are included here as they are; the GPU note
is the page :doc:`/user-guide/gpu`. They describe the implementation in more
detail than the rest of this documentation and are updated together with the
code.

.. toctree::
   :maxdepth: 1

   design/coordinates-and-fields
   design/config
   design/neo-orb
   design/wall-heat-loads
   design/full-orbit-integration
   design/spectre-interface-crossing

The derivation of canonical and Boozer flux coordinates from VMEC data is
available as LaTeX source in ``DOC/canonical_and_boozer_flux_coords_via_VMEC.tex``
and the integrator notes in ``DOC/neo-orb.tex``.
