==========
Quickstart
==========



Install
-------

.. code:: bash

    pip install pydseamslib
    pip install 'pydseamslib[ase]'
    pip install 'pydseamslib[solvis]'
    pip install 'pydseamslib[mdanalysis]'

``pydseams`` is the package name to import (``import pydseams as ds``). The
``[ase]`` extra adds ASE support for ``from_ase`` and ``to_ase``, ``[solvis]``
adds ASE and solvis-tools for ``to_solvis``, and ``[mdanalysis]`` adds
MDAnalysis support for ``IceStates``. ``import pydseamslib`` and
``import dseams`` are included in the same installation.

With pixi, add the PyPI package:

.. code:: bash

    pixi add --pypi pydseamslib

Requires Python 3.11 or later. Prebuilt wheels bundle the d-SEAMS C++
engine so no compilation step is needed.

The full extra and source-build details are in the
`install how-to <howto/install.rst>`_.

Classify a frame
----------------

Build a 3x3x3 cubic diamond ice supercell from coordinates and classify
it with CHILL+ and topological cages:

.. code:: python

    import pydseams as ds

    a = 2.75 * 4 / 3**0.5
    basis = (
        (0, 0, 0),
        (0.5, 0.5, 0),
        (0.5, 0, 0.5),
        (0, 0.5, 0.5),
        (0.25, 0.25, 0.25),
        (0.75, 0.75, 0.25),
        (0.75, 0.25, 0.75),
        (0.25, 0.75, 0.75),
    )
    xyz = [
        [(i + x) * a, (j + y) * a, (k + z) * a]
        for i in range(3)
        for j in range(3)
        for k in range(3)
        for x, y, z in basis
    ]
    frame = ds.from_arrays(xyz, [3 * a, 3 * a, 3 * a])
    print(frame.chill_plus())
    print(frame.cages())

Expected output:

::

    IceCounts(cubic=216)
    CageScore(n_ih=0, n_ic=216, n_water=0)

To read a trajectory file instead, pass its path to ``ds.read(path)``
(supporting LAMMPS dumps, XYZ, and eOn ``.con``, or PDB/GRO/DCD when
chemfiles is linked). ``chill_plus`` and ``cages`` compute in memory and do
not write files. Pass ``all_atoms=True`` to ``ds.read`` when mixed-type
site analyses such as ``Frame.pairs`` and ``Frame.domain`` need every
species.

Seeded cages and ring-adjacent completion
-----------------------------------------

``Frame.seeded_affiliation`` runs the seeded hysteresis cage
construction that ``Frame.cages(seeded=True)`` wraps. Setting
``ring_adjacent=True`` fills the last vertex of a six-ring whose other
five vertices already carry a cage label:

.. code:: python

    print(frame.seeded_affiliation(ring_adjacent=True))
    print(frame.cages(seeded=True, ring_adjacent=True))

``IceFeaturizer`` turns ``ring_adjacent=True`` on by default.

Mixed ASE selections
--------------------

Pass a sequence to ``from_ase`` when an ASE ``Atoms`` object mixes water
and salt ions. The listed species remain in the cloud:

.. code:: python

    from ase import Atoms
    import pydseams as ds

    brine = Atoms(
        "O" * len(xyz) + "Na",
        positions=xyz + [[1.5 * a, 1.5 * a, 1.5 * a]],
        cell=[3 * a, 3 * a, 3 * a],
        pbc=True,
    )
    brine_frame = ds.from_ase(brine, select=("O", "Na"), bonded="cutoff")
    print(brine_frame.n_atoms)

The first entry is the analysed water species. The remaining species
keep their atomic numbers as ``c_type`` for ``ion_environment``. See the
`ASE how-to <howto/ase.rst>`_ and `features how-to <howto/features.rst>`_.
