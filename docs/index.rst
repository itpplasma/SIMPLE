.. SIMPLE documentation master file

SIMPLE
======

**S**\ ymplectic **I**\ ntegration **M**\ ethods for **P**\ article **L**\ oss
**E**\ stimation.

SIMPLE traces guiding-centre orbits of energetic particles in
three-dimensional magnetic fields and reports when and where they leave the
plasma. Its main application is the fast computation of fusion alpha-particle
losses in stellarators, both for the analysis of a given configuration and
inside a stellarator optimisation loop. The orbits are integrated with
symplectic methods in canonical coordinates, which improve long-time
behaviour at suitable time steps. Regular and chaotic
orbits can be classified early, which shortens the computation further.

.. figure:: /_static/figures/front.png
   :alt: Guiding-centre orbits of fusion alpha particles in a stellarator and
         the resulting loss pattern on the plasma boundary.
   :width: 90%

   Front figure of SIMPLE.

SIMPLE is written in Fortran with OpenMP parallelisation, an optional CUDA
backend for GPUs and a Python interface, ``pysimple``. It reads VMEC
equilibria, Boozer chartmaps (also from GVEC and booz_xform), coil files and
SPECTRE equilibria. It is distributed under the MIT License. The source code
is at https://github.com/itpplasma/SIMPLE.

How to read this documentation
------------------------------

* :doc:`getting-started/installation` and :doc:`getting-started/quickstart`
  get a first loss computation running in a few minutes.
* The :doc:`tutorials/index` walk through three complete applications:
  a loss fraction from VMEC, wall heat loads, and SIMPLE inside an
  optimisation loop.
* The :doc:`user-guide/index` explains inputs, runs and outputs task by task.
* :doc:`model/index` documents the physical model and the numerical methods
  with their equations.
* :doc:`verification/index` collects the published benchmarks and the test
  suite.
* :doc:`reference/index` lists every input parameter and the Python API.

If you use SIMPLE in a publication, cite the references listed in
:doc:`reference/publications`.

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

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
