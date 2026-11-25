===================================
Bond-orientational order parameters
===================================


``Frame.steinhardt(order_l)`` computes the local Steinhardt parameter
:math:`q_l(i)` and its Lechner-Dellago neighbour-averaged counterpart
:math:`\bar{q}_l(i)` for every atom on the cutoff neighbour list. Supported
degrees are 3, 4, 6, 8, and 12.

Prerequisites
-------------

.. code:: bash

    pip install pydseamslib
    mkdir -p input/traj
    curl -sLo input/traj/exampleTraj.lammpstrj \
      https://raw.githubusercontent.com/d-SEAMS/seams-core/v2.14.3/input/traj/exampleTraj.lammpstrj

Compute Steinhardt parameters
-----------------------------

.. code:: python

    import pydseams as ds

    frame = ds.read("input/traj/exampleTraj.lammpstrj", atom_type=2)
    q3 = frame.steinhardt(3)
    mean_ql = sum(q3["ql"]) / len(q3["ql"])
    mean_ql_bar = sum(q3["ql_bar"]) / len(q3["ql_bar"])
    print(f"mean q3={mean_ql:.6f} mean q3_bar={mean_ql_bar:.6f}")
    assert len(q3["ql"]) == 250
    assert round(float(mean_ql), 6) == 0.374897
    assert round(float(mean_ql_bar), 6) == 0.187542

``q3["ql"]`` and ``q3["ql_bar"]`` are the per-atom arrays. On
``exampleTraj.lammpstrj``, their means are 0.374897 and 0.187542.

Definition
----------

The local parameter for atom :math:`i` is



.. math::

    q_l(i) = \sqrt{\frac{4\pi}{2l+1} \sum_{m=-l}^{l} |q_{lm}(i)|^2}

with complex bond-orientational moments



.. math::

    q_{lm}(i) = \frac{1}{N_b(i)} \sum_{j \in N(i)} Y_{lm}(\hat r_{ij}),

where :math:`N(i)` is the cutoff neighbour set of atom :math:`i` and :math:`Y_{lm}` are
spherical harmonics.
