SIMPLE
======

**S**\ ymplectic **I**\ ntegration **M**\ ethods for **P**\ article **L**\ oss
**E**\ stimation.

SIMPLE traces guiding-centre orbits of energetic particles, by default
3.5 MeV fusion alphas, in 3D magnetic fields (VMEC, Boozer chartmaps from GVEC
or booz_xform, coils, SPECTRE) and records when, where and with which energy
they leave the plasma. It integrates in canonical coordinates with symplectic
schemes and classifies orbits as regular or chaotic. It is written in Fortran
with OpenMP, an optional CUDA backend and the Python package ``pysimple``.

.. figure:: /_static/figures/front.png
   :alt: Loss positions of alphas on the last closed flux surface of a
         quasi-axisymmetric stellarator and the confined fraction over time.
   :width: 90%

   Landreman–Paul QA equilibrium at reactor scale, 3000 collisionless 3.5 MeV
   alphas started at :math:`s = 0.5`: guiding-centre positions of the 194 lost
   markers on the LCFS (left) and confined fraction, split into passing and
   trapped, over 50 ms (right).

**Users start here:** :doc:`getting-started/installation`,
:doc:`getting-started/quickstart`, then the :doc:`tutorials/index` and the
task pages of the :doc:`user-guide/index`. Every input key is in
:doc:`reference/input-parameters`.

**Developers start here:** :doc:`model/index` for the equations as
implemented, :doc:`development/index` for build, layout and tests, and
:doc:`verification/index` for what the test suite checks.

Cite: :doc:`reference/publications`. Source: https://github.com/itpplasma/SIMPLE
(MIT License).

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   getting-started/installation
   getting-started/quickstart

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/index

.. toctree::
   :maxdepth: 2
   :caption: User guide

   user-guide/index

.. toctree::
   :maxdepth: 2
   :caption: Model and methods

   model/index

.. toctree::
   :maxdepth: 2
   :caption: Verification

   verification/index

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/index

.. toctree::
   :maxdepth: 2
   :caption: Development

   development/index
