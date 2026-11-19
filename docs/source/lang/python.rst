dseams module reference
=======================

The API reference below covers the :mod:`dseams` module, generated from
its NumPy-style docstrings.

.. automodule:: dseams
   :members:
   :undoc-members:

Cubic ice verification
----------------------

On the eight-atom diamond cubic unit cell (nearest-neighbour spacing
2.75 angstrom), every atom is classified as cubic ice (label 0,
``counts[0] == 8``):

.. doctest::

   >>> import dseams
   >>> a = 2.75 * 4 / 3**0.5
   >>> sites = (
   ...     (0, 0, 0),
   ...     (0.5, 0.5, 0),
   ...     (0.5, 0, 0.5),
   ...     (0, 0.5, 0.5),
   ...     (0.25, 0.25, 0.25),
   ...     (0.75, 0.75, 0.25),
   ...     (0.75, 0.25, 0.75),
   ...     (0.25, 0.75, 0.75),
   ... )
   >>> xyz = [[x * a, y * a, z * a] for x, y, z in sites]
   >>> labels, counts = dseams.chill_plus(xyz, [a, a, a])
   >>> int(counts[0])
   8
