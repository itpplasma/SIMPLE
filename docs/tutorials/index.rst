Tutorials
=========

The tutorials are written for physicists who know stellarator equilibria and
alpha-particle physics but have not used SIMPLE before. Each tutorial runs on
public input data in a few minutes on a workstation. Commands are run from
the repository root unless stated otherwise, after building SIMPLE as
described in :doc:`/getting-started/installation`.

.. toctree::
   :maxdepth: 1

   loss-fraction
   wall-heat-loads
   optimisation-loop

The numbers printed on these pages come from runs on a 32-thread Linux
workstation with a Release build and a fixed random seed
(``deterministic = .true.``). Another compiler, CPU or build type changes the
rounding, so individual orbits and loss times can differ. Loss fractions then
agree within their statistical error.
