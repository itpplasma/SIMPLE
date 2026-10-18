Plotting results
================

``simple-plot`` turns the files of a finished run into figures, including a
binomial confidence band on the loss fraction. It reads output only and never
runs SIMPLE. Install it with ``python -m pip install '.[plotting]'``.

.. code-block:: bash

   simple-plot losses run/ -o losses.png
   simple-plot losses run/ --confidence 0.99 -o losses-99.png
   simple-plot classes run/ -o classes.png
   simple-plot wall run/results.nc --color incidence -o wall.png
   simple-plot poincare run/ -o poincare.png

``-o`` saves the figure (format from the extension); without it the figure
opens interactively. Set ``MPLBACKEND=Agg`` on headless machines. A missing
file or variable gives an error and exit status 2.

Losses and their uncertainty
----------------------------

``losses`` reads ``confined_fraction.dat`` and plots, on a logarithmic time
axis, the confined total, its passing and trapped parts (same denominator,
so they add up to the total) and the lost fraction. Rows at :math:`t = 0` are
skipped; rows without resolved markers appear as gaps.

The shaded band is a pointwise Wilson interval for the lost fraction. With
:math:`n` the resolved markers (column 4), :math:`k = \text{round}(n(1 -
f_\text{pass} - f_\text{trap}))` the lost markers, :math:`p = k/n` and
:math:`z` the standard-normal quantile at :math:`(1 + c)/2` for confidence
:math:`c` (default 0.95, ``--confidence``),

.. math::

   p_\pm = \frac{p + \dfrac{z^2}{2n} \pm z\sqrt{\dfrac{p(1-p)}{n}
           + \dfrac{z^2}{4n^2}}}{1 + \dfrac{z^2}{n}} .

Implemented in ``python/pysimple/cli.py`` (``wilson_interval``). Unlike
:math:`p \pm z\sqrt{p(1-p)/n}`, it stays inside :math:`[0, 1]` and is
nonzero for :math:`k = 0`. It assumes independent, equally weighted markers
and covers sampling error only, not numerical error, and it is not a
simultaneous band for the whole curve. Unresolved markers are excluded, not
counted as lost; check ``unresolved_fraction.dat``.

Classes, wall hits, Poincaré sections
-------------------------------------

``classes``
   Counts the :math:`J_\parallel`, topological and fractal class codes,
   including 0 (unclassified), from ``class_parts.dat`` or, if absent,
   ``iclass`` in ``results.nc``.
``wall``
   R–Z footprint of markers with ``wall_hit = 1`` at ``wall_hit_cart``, all
   toroidal positions combined, in ``cart_units``. Colour is
   :math:`E/E_0 = (v/v_0)^2` from ``zend`` (default) or
   ``wall_hit_cos_incidence`` (``--color incidence``; 0 grazing, 1 normal).
   Needs ``results.nc`` from a run with ``wall_input`` and
   ``output_results_netcdf = .true.``; LCFS losses have no wall footprint.
``poincare``
   R–Z field-line sections from the ``poincare_vol<N>.dat`` files of
   ``spectre_poincare.x``, coloured by the volume-local :math:`s`.

File formats: :doc:`outputs`. Wall setup: :doc:`wall-heat-loads`.
