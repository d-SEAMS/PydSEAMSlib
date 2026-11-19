======================
Tutorial: Classify ice
======================


Load a mixed ice-water LAMMPS dump, print CHILL+ counts, and compare
them to cage membership on the same frame.

Prerequisites
-------------

- Python 3.11+

- ``pydseams`` installed: ``pip install pydseamslib``

- The sample mixed TIP4P dump ``exampleTraj.lammpstrj``

Learning goals
--------------

In this tutorial you learn to:

1. Download the sample LAMMPS dump and load it with ``ds.read``

2. Classify every analysed oxygen with ``Frame.chill_plus``

3. Score hexagonal-cage (HC) and double-diamond-cage (DDC) membership
   with ``Frame.cages``

4. Read ``IceCounts`` and ``CageScore``

Step 1: Install pydseams
------------------------

.. code:: bash

    pip install pydseamslib

``pydseams`` is the primary import name. ``import pydseamslib`` is an alias
of ``pydseams``.

Step 2: Fetch the sample dump
-----------------------------

Download the sample TIP4P ice-water dump (``exampleTraj.lammpstrj``, 750
atoms: 250 oxygens of type 2 and 500 hydrogens of type 1):

.. code:: bash

    mkdir -p input/traj
    curl -sLo input/traj/exampleTraj.lammpstrj \
      https://raw.githubusercontent.com/d-SEAMS/seams-core/main/input/traj/exampleTraj.lammpstrj

Step 3: Load the frame
----------------------

.. code:: python

    import pydseams as ds

    frame = ds.read("input/traj/exampleTraj.lammpstrj")
    print(frame.n_atoms)
    print(frame.atom_type)

``ds.read`` selects the reader from the file suffix. For ``.lammpstrj``, it
calls ``Frame.from_file``, which automatically selects LAMMPS type 2
(oxygen) when present.

Expected output:

::

    250
    2

``read`` keeps the 250 oxygen atoms for analysis. Classification does not
write output files.

Step 4: Classify with CHILL+
----------------------------

.. code:: python

    counts = frame.chill_plus()
    print(counts)
    assert str(counts) == "IceCounts(interClathrate=9, water=241)"

Expected output:

::

    IceCounts(interClathrate=9, water=241)

CHILL+ classifies each oxygen from its four-neighbour bond geometry: 9
interfacial clathrate molecules (three eclipsed bonds) and 241 water
molecules.

``IceCounts`` behaves like a dictionary and an attribute record:
``counts.cubic`` and ``counts["cubic"]`` return the same integer, and
missing labels return ``0``.

Step 5: Score cages
-------------------

.. code:: python

    cages = frame.cages()
    print(cages)
    assert str(cages) == "CageScore(n_ih=0, n_ic=0, n_water=250)"

Expected output:

::

    CageScore(n_ih=0, n_ic=0, n_water=250)

``cages()`` scores complete hexagonal cages (HC, ice Ih) and
double-diamond cages (DDC, ice Ic). This mixed snapshot has no
complete HC or DDC, so every oxygen remains water under the cage
score. ``n_ih``, ``n_ic``, and ``n_water`` are atom counts, not cage counts.

Default ``cages()`` uses the seeded (hysteresis) construction. Pass
``seeded=False`` for cutoff-graph affiliation on this frame's six-rings.

``Frame.seeded_affiliation`` exposes that construction directly. Pass
``ring_adjacent=True`` on ``cages`` or ``seeded_affiliation`` to fill the last
vertex of a six-ring whose other five vertices carry a label.

Cage scores and CHILL+ answer different questions on the same cloud.
CHILL+ evaluates four-neighbour local coordination, whereas cages
require a closed HC or DDC polyhedron.

Step 6: Inspect the compiled module
-----------------------------------

``Frame`` delegates numerical work to the compiled module ``pydseams.yoda``:

.. code:: python

    from pydseams import yoda

    print(yoda.__doc__)
    print(type(frame.cloud))

Expected output includes ``d-SEAMS compiled surface (yoda)`` and
``<class 'pydseams.yoda.PointCloudDouble'>``.

Troubleshooting
---------------

``ModuleNotFoundError: No module named 'pydseams'``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install into the active Python interpreter:

.. code:: bash

    python -m pip install pydseamslib

``ValueError`` or missing dump file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Verify that ``input/traj/exampleTraj.lammpstrj`` exists in your working
directory, or pass the full path to the downloaded dump to ``ds.read``.

Counts differ on another trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The counts above apply to ``exampleTraj.lammpstrj``. A different
trajectory, 1-indexed ``frame``, ``cutoff``, or ``atom_type`` produces a
different histogram.

Next steps
----------

`Classify ASE Atoms <../howto/ase.rst>`_
    ``from_ase`` / ``to_ase``

`Features and ions <../howto/features.rst>`_
    ``IceFeaturizer`` / ``ion_environment``

`View a frame in solvis <../howto/solvis.rst>`_
    optional visualisation extra

`The yoda surface <../explanation/yoda-surface.rst>`_
    why helpers stay in Python

`Python surface <../reference/python.rst>`_
    live names

Summary
-------

In this tutorial you:

1. Installed ``pydseamslib`` and loaded a LAMMPS trajectory frame

2. Loaded oxygen atoms with ``ds.read``

3. Computed CHILL+ labels with ``frame.chill_plus()``

4. Scored HC / DDC cage membership with ``frame.cages()``
