Input parameters
================

A run reads the ``&config`` and ``&profiles`` namelists from ``simple.in``
and needs only the keys that differ from the defaults.
``examples/simple_full.in`` below lists every key with its default and
meaning; ``test_config_reference`` checks it against the compiled defaults.
:doc:`/input-reference` gives types, defaults and bounds generated from the
Fortran declarations, and :doc:`numerical-settings` explains the keys that set
step size, event location and random streams.

.. literalinclude:: ../../examples/simple_full.in
   :language: fortran
   :caption: examples/simple_full.in

Preflight validation
--------------------

Check keys, types and bounds before a run:

.. code-block:: bash

   python tools/check_input.py simple.in examples/*.in

The validator uses the generated schema ``docs/input-reference.json``. File
contents, fields and interacting options are checked by SIMPLE at run time.
