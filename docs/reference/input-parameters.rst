Input parameters
================

A run is configured by the ``&config`` and ``&profiles`` namelists of
``simple.in``. The file ``examples/simple_full.in`` lists every key with its
default value and meaning. The test ``test_config_reference`` checks it
against the compiled defaults, so it is always complete and current. A run
needs only the keys that differ from the defaults.

:doc:`/input-reference` lists the parameter types, defaults and bounds
generated from the Fortran declarations. Guidance for choosing values is in
:doc:`/user-guide/running` and in ``DOC/config.md``
(:doc:`design-notes`).

.. literalinclude:: ../../examples/simple_full.in
   :language: fortran
   :caption: examples/simple_full.in
