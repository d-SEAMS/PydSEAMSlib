==================
Classify ASE Atoms
==================



Problem
-------

You have an ASE ``Atoms`` and want CHILL+ labels or cage flags on the
same configuration, then an ``Atoms`` back with those arrays.

Install the extra
-----------------

.. code:: bash

    pip install 'pydseamslib[ase]'

Without it, ``from_ase`` / ``to_ase`` raise ``ImportError`` and name that
command.

Classify ASE Atoms
------------------

.. code:: python

    import ase.io
    import pydseams as ds

    atoms = ase.io.read("water.lammpstrj", format="lammps-dump-text")
    frame = ds.from_ase(atoms)
    print(frame.chill_plus())
    labelled = frame.to_ase()

``from_ase`` is ``Frame.from_ase``. The cell must be nonsingular and
periodic in all three directions. Orthorhombic and general triclinic
cells are supported. ``from_ase`` keeps oxygen by default
(``select="O"``).

Select and bonding
------------------

Pass a symbol, an atomic number, a sequence of those, or ``None``
(every atom):

.. code:: python

    frame = ds.from_ase(atoms, select="O")
    frame = ds.from_ase(atoms, select=8)
    frame = ds.from_ase(atoms, select=("O", "Na", "Cl"))
    frame = ds.from_ase(atoms, select=None)

``select`` as a sequence keeps every listed species in the cloud
and analyses the first, so ions ride along for
``pydseams.features.ion_environment``. Bonding then follows the
same ``bonded`` rule as a single-species select.

``bonded="auto"`` uses hydrogen bonds only when the ``Atoms`` include ``H``
and the selected analysis cloud excludes ``H``. The default oxygen
selection on molecular water therefore uses hydrogen bonds, while
``select=None`` uses the cutoff neighbour list. Single-site models (mW)
also use the cutoff neighbour list:

.. code:: python

    frame = ds.from_ase(atoms, select="O", bonded="cutoff")

``bonded`` is ``"auto"``, ``"hbond"``, or ``"cutoff"``. Explicit
``bonded="hbond"`` requires a heavy-atom selection such as ``select="O"``.
For that selection, hydrogens stay in a side cloud so the analysed
species remain the CHILL / ring particles.

If ``atoms.arrays["mol-id"]`` exists, its molecule IDs associate each
hydrogen with its selected donor atom. Otherwise, ``from_ase`` assigns
each hydrogen to the nearest selected atom using ASE's periodic
minimum-image distance. Ordinary ``O H H`` water ordering does not need
extra metadata.

Labels on the way back
----------------------

``to_ase`` rebuilds an ``Atoms`` with the analysed positions and the
imported cell orientation and displacement (``pbc=True``). General cells
are represented internally in LAMMPS restricted-triclinic form. After
``chill_plus`` / ``chill``,
``atoms.arrays["ice_type"]`` holds the per-atom names. After
``cages()``, ``atoms.arrays["hc"]`` and ``atoms.arrays["ddc"]`` hold the
last ``CageScore``. ``atoms.info["dseams_n_atoms"]`` is the analysed
particle count.

A frame loaded from a LAMMPS dump (no ASE symbols) uses ``O`` as the
fallback species.

See also
--------

`View a frame in solvis <solvis.rst>`_
    wraps the same ``Atoms``

`Classify ice <../tutorials/classify-ice.rst>`_
    ``read`` / ``chill_plus`` / ``cages``

`Troubleshooting <troubleshooting.rst>`_
    periodicity, singular cells, missing H
