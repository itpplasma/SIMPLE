Python interface
================

.. include:: ../../python/README.rst
   :start-line: 3

Limitation
----------

With ``startmode = 1``, ``pysimple.sample()`` and ``sample_surface(n, s)``
return start points with invalid angles when the start surface (``s``, or
``sbeg`` for ``sample()``) equals that of the preceding ``init`` or sampling
call (``src/simple_api.f90``: the start-surface
arrays are reallocated, but the field line that fills them is not
re-integrated). Change ``s`` between calls, or use ``sample_volume`` or
``load_particles``.

Function reference: :doc:`/api/index`. Example scripts: :doc:`/examples`.
