======================
Tutorial: Classify ice
======================


Load a mixed ice-water LAMMPS dump, print CHILL+ counts, and compare
them to cage membership on the same frame.

Prerequisites
-------------

- Python 3.12+

- ``pydseams`` installed: ``pip install pydseamslib``

- The engine fixture dump, either:

  - from a PydSEAMSlib checkout: ``tests/data/exampleTraj.lammpstrj``

  - from `seams-core <https://github.com/d-SEAMS/seams-core>`_:
    ``input/traj/exampleTraj.lammpstrj`` (same file)

A wheel already links the C++ engine. The ``seams`` CLI lives in
seams-core; Lua is ``dseams`` in
`yodaStruct <https://github.com/d-SEAMS/yodaStruct>`_. This tutorial
does not run those front ends.

Learning Objectives
-------------------

By the end of this tutorial you will be able to:

1. Load one configuration with ``ds.read``

2. Classify every analysed oxygen with ``Frame.chill_plus``

3. Score hexagonal-cage and double-diamond-cage membership with
   ``Frame.cages``

4. Read ``IceCounts`` and ``CageScore``

5. Point at the ``seams`` CLI and the Lua ``dseams`` module when you
   want those front ends

Step 1: Install pydseams
------------------------

.. code:: bash

    pip install pydseamslib

``pydseams`` is the Python package. ``import pydseamslib`` is a
compatibility alias of ``pydseams``. Requires Python 3.12+.

Step 2: Locate the fixture dump
-------------------------------

From a PydSEAMSlib checkout, at the repository root, the vendored
copy is ``tests/data/exampleTraj.lammpstrj``.

That file is the seams-core engine fixture
``input/traj/exampleTraj.lammpstrj``. If you installed from PyPI only,
copy the seams-core path and pass that filename to ``ds.read``.

Step 3: Load the frame
----------------------

.. code:: python

    import pydseams as ds

    frame = ds.read("tests/data/exampleTraj.lammpstrj")
    print(frame.n_atoms)
    print(frame.atom_type)

``ds.read`` looks at the suffix and calls the matching constructor.
=.lammpstrj~ goes through ``Frame.from_file``, which keeps LAMMPS type
2 (oxygen) when that type is present.

Expected output:

::

    250
    2

``read`` keeps the 250 oxygen atoms. Classification does not write
files.

Step 4: Classify with CHILL+
----------------------------

.. code:: python

    print(frame.chill_plus())

Expected output:

::

    IceCounts(interClathrate=12, water=238)

These counts are the engine fixture result for seams-core
``input/traj/exampleTraj.lammpstrj``, vendored here as
``tests/data/exampleTraj.lammpstrj``. CHILL+ labels each oxygen from
its four neighbours: 12 interfacial clathrate (three eclipsed bonds)
and 238 water.

``IceCounts`` is a histogram. ``counts.cubic`` and ``counts["cubic"]``
are the same value. Missing labels read as ``0``.

Step 5: Score cages
-------------------

.. code:: python

    print(frame.cages())

Expected output:

::

    CageScore(n_ih=0, n_ic=0, n_water=250)

``cages()`` scores complete hexagonal cages (HC, ice Ih) and
double-diamond cages (DDC, ice Ic). This mixed snapshot has no
finished HC or DDC, so every oxygen stays water under the cage
score. ``n_ih``, ``n_ic``, and ``n_water`` are atom counts, not cage
counts.

Default ``cages()`` is the seeded (hysteresis) construction. That
path calls ``seeded_affiliation``. ``ring_adjacent`` fills the last
vertex of a six-ring whose other vertices carry a label. Pass
``seeded=0`` for cutoff-graph affiliation on this frame's
six-rings. For an ASE cloud that must keep ions, pass a sequence:
``Frame.from_ase(atoms, select=("O", "Na", "Cl"))``.

CHILL+ and cages answer different questions on the same cloud.
CHILL+ is a four-neighbour local label. Cages require a finished
HC or DDC.

Step 6: The same engine elsewhere
---------------------------------

``Frame`` calls the compiled module ``yoda``. Same session:

.. code:: python

    from pydseams import yoda

    print(yoda.__doc__)
    print(type(frame.cloud))

Expected output includes ``d-SEAMS compiled surface (yoda)`` and
``<class 'pydseams.yoda.PointCloudDouble'>``.

The ``seams`` CLI in seams-core and ``require("dseams")`` in yodaStruct
run ``libyodaLib``, the same engine ``yoda`` binds. Use those projects
for the command-line and Lua workflows. Stay on ``Frame`` for Python.

Troubleshooting
---------------

``ModuleNotFoundError: No module named 'pydseams'``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install into the active interpreter:

.. code:: bash

    python -m pip install pydseamslib

``FileNotFoundError`` on the dump
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From a checkout, run at the repository root or pass an absolute
path to ``tests/data/exampleTraj.lammpstrj``. From PyPI, copy
seams-core ``input/traj/exampleTraj.lammpstrj``.

Counts do not match the fixture
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The numbers above are the engine fixture result for that dump.
A different trajectory, frame index, cutoff, or ``atom_type``
produces a different histogram. ``ds.read(path, frame=1)`` is
1-indexed.

You want the CLI or Lua
~~~~~~~~~~~~~~~~~~~~~~~

Do not compile ``yoda`` from this tutorial. Install
`seams-core <https://github.com/d-SEAMS/seams-core>`_ for ``seams``,
and `yodaStruct <https://github.com/d-SEAMS/yodaStruct>`_ for
``require("dseams")``.

Next Steps
----------

`Classify ASE Atoms <../howto/ase.rst>`_
    ``from_ase`` / ``to_ase``

`Features and ions <../howto/features.rst>`_
    ``IceFeaturizer``, ``ring_adjacent``, sequence ``select``

`View a frame in solvis <../howto/solvis.rst>`_
    optional extra

`The yoda surface <../explanation/yoda-surface.rst>`_
    why helpers stay in Python

`Python surface <../reference/python.rst>`_
    live names

Summary
-------

You learned to:

1. Install ``pydseams`` and load the engine fixture dump

2. Call ``ds.read`` to get a ``Frame``

3. Print CHILL+ labels with ``frame.chill_plus()``

4. Print HC / DDC membership with ``frame.cages()``

5. Leave the ``seams`` CLI and Lua ``dseams`` module to their own
   repositories
