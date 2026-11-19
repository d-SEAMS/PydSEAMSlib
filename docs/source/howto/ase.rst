==================
Classify ASE Atoms
==================



Problem
-------

You have an ASE ``Atoms`` configuration and want CHILL+ labels or cage
flags on it, then return an ``Atoms`` object carrying those arrays.

Install the extra
-----------------

.. code:: bash

    pip install 'pydseamslib[ase]'

Without it, ``from_ase`` and ``to_ase`` raise ``ImportError``.

Classify ASE Atoms
------------------

Build or load an ASE ``Atoms`` object and pass it to ``ds.from_ase``:

.. code:: python

    from ase import Atoms
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
    atoms = Atoms("O" * len(xyz), positions=xyz, cell=[3 * a, 3 * a, 3 * a], pbc=True)

    frame = ds.from_ase(atoms, bonded="cutoff")
    print(frame.chill_plus())
    print(frame.cages())
    labelled = frame.to_ase()
    print("arrays:", sorted(labelled.arrays.keys()))

When reading a LAMMPS dump with ``ase.io.read(..., format="lammps-dump-text")``
where the file stores numeric types rather than element symbols, assign
element symbols first (for example ``atoms.set_chemical_symbols(["O" if int(t) == 2 else "H" for t in atoms.arrays["type"]])``).

``from_ase`` is ``Frame.from_ase``. The cell must be nonsingular and
periodic in all three directions. Orthorhombic and general triclinic
cells are supported. ``from_ase`` keeps oxygen by default
(``select="O"``).

Select and bonding
------------------

Pass a symbol, an atomic number, ``None`` (every atom), or a sequence:

.. code:: python

    frame_o = ds.from_ase(atoms, select="O", bonded="cutoff")
    frame_z = ds.from_ase(atoms, select=8, bonded="cutoff")
    frame_all = ds.from_ase(atoms, select=None)

    brine = Atoms(
        "O" * len(xyz) + "NaCl",
        positions=xyz + [[1.5 * a, 1.5 * a, 1.5 * a], [0.5 * a, 0.5 * a, 0.5 * a]],
        cell=[3 * a, 3 * a, 3 * a],
        pbc=True,
    )
    frame_mix = ds.from_ase(brine, select=("O", "Na", "Cl"), bonded="cutoff")
    print("mixed cloud n_atoms:", frame_mix.n_atoms)

A sequence keeps every listed species in the cloud and analyses the
first. Ions remain in the cloud for ``pydseams.features.ion_environment``.

``bonded="auto"`` uses hydrogen bonds only when the ``Atoms`` include ``H``
and the selected analysis cloud excludes ``H``. The default oxygen
selection on molecular water therefore uses hydrogen bonds, while
``select=None`` uses the cutoff neighbour list. Single-site models (such
as mW) also use the cutoff neighbour list:

.. code:: python

    frame_cutoff = ds.from_ase(atoms, select="O", bonded="cutoff")

``bonded`` accepts ``"auto"``, ``"hbond"``, or ``"cutoff"``. Explicit
``bonded="hbond"`` requires a heavy-atom selection such as ``select="O"``.
For that selection, hydrogens stay in a side cloud so the analysed
species remain the CHILL and ring particles.

If ``atoms.arrays["mol-id"]`` exists, its molecule IDs associate each
hydrogen with its selected donor atom. Otherwise, ``from_ase`` assigns
each hydrogen to the nearest selected atom using ASE's periodic
minimum-image distance.

Labels on the way back
----------------------

``to_ase`` rebuilds an ``Atoms`` object with the analysed positions and
the imported cell orientation and displacement (``pbc=True``). General
cells are represented internally in LAMMPS restricted-triclinic form.
After ``chill_plus`` or ``chill``, ``atoms.arrays["ice_type"]`` holds the
per-atom names. After ``cages()``, ``atoms.arrays["hc"]`` and
``atoms.arrays["ddc"]`` hold the last ``CageScore``.
``atoms.info["dseams_n_atoms"]`` stores the analysed particle count.

A frame loaded from a LAMMPS dump without element symbols uses ``O`` as
the fallback species.

See also
--------

`View a frame in solvis <solvis.rst>`_
    wraps the same ``Atoms``

`Classify ice <../tutorials/classify-ice.rst>`_
    ``read`` / ``chill_plus`` / ``cages``

`Troubleshooting <troubleshooting.rst>`_
    periodicity, singular cells, missing H
