Examples
========

The snippets below are pulled directly from the example scripts inside the
repository.  They rely exclusively on the public :mod:`simple` API, ensuring the
documentation never drifts from runnable code.

Simple session trace
--------------------

.. literalinclude:: ../examples/simple_api.py
   :language: python
   :caption: examples/simple_api.py
   :linenos:

Fast classification
-------------------

.. literalinclude:: ../examples/classify_fast.py
   :language: python
   :caption: examples/classify_fast.py
   :linenos:

Both examples expose helper functions (:func:`run_trace_example` and
:func:`classify_fast_example`) that are exercised directly in the pytest suite,
so the documentation mirrors the tested code paths one-to-one.
