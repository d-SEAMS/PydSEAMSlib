===============
Troubleshooting
===============


.. contents::


Installation
------------

``ModuleNotFoundError: No module named 'pydseams'``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The package is not installed in the active Python interpreter:

.. code:: bash

    python -c "import sys; print(sys.executable)"
    python -m pip install pydseamslib

``ImportError`` naming ``pip install 'pydseamslib[ase]'``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ASE is an optional extra:

.. code:: bash

    pip install 'pydseamslib[ase]'

``ImportError`` naming ``pip install 'pydseamslib[solvis]'``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``solvis-tools`` is an optional extra:

.. code:: bash

    pip install 'pydseamslib[solvis]'

Loading frames
--------------

``ValueError: ... has no atoms of type 1 or 2``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``atom_type`` is omitted, ``Frame.from_file`` looks for LAMMPS type 2
(oxygen) then type 1. Pass ``atom_type`` explicitly if your dump uses a
different type code, or use ``from_ase`` / ``from_arrays``.

``ValueError: from_ase requires a cell periodic in all three directions``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Enable periodic boundary conditions on all three ASE cell axes
(``atoms.pbc = True``). Partially periodic configurations do not define
three-dimensional periodic geometry.

``ValueError: from_ase cannot convert a singular cell``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Supply a cell with nonzero volume. Orthorhombic and general triclinic
cells are both accepted when nonsingular and fully periodic.

``ValueError: no atoms matched select=...``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``from_ase(..., select="O")`` found no oxygen atoms. Pass the chemical
symbol or atomic number present in your configuration, or ``select=None``.

``ValueError: bonded="hbond" requires a heavy-atom selection that excludes H``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ASE analysis selection includes hydrogen. Use ``select="O"`` (or
another heavy species) for hydrogen-bond topology, or pass
``bonded="cutoff"`` for an all-atom frame.

``ValueError: Hydrogen-bond analysis needs hydrogens``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``bonded="hbond"`` (or ``"auto"`` when a heavy-atom selection is used on a
water frame) requires hydrogens. Pass ASE ``Atoms`` with ``H``, a LAMMPS
dump containing hydrogen, or set ``bonded="cutoff"``.

``RuntimeError: chemfiles is not linked in this build of seams-core``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reading ``.pdb``, ``.gro``, or ``.dcd`` directly requires a build linked
against chemfiles (``available_readers()["chemfiles"]`` is ``True``).
Convert the trajectory to LAMMPS or XYZ (for example via ASE), or
build from source with chemfiles enabled.

``RuntimeError: readcon-core is not linked``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reading ``.con`` requires ``available_readers()["readcon"]`` to be ``True``.

``RuntimeError: this build has no readXYZ``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reading ``.xyz`` requires ``yoda.readXYZ``. Check
``available_readers()["xyz"]``.

Classification
--------------

Cage score is all water on a mixed snapshot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On mixed snapshots such as ``exampleTraj.lammpstrj``, CHILL+ can detect
local four-neighbour ice-like stars while ``cages()`` finds zero finished
HC or DDC polyhedra. See `Classify ice <../tutorials/classify-ice.rst>`_.

Counts differ from a documented dump
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Confirm the trajectory path, 1-indexed ``frame``, ``cutoff`` (default 3.5
Angstroms), and ``atom_type``.

Common messages
---------------

.. table::

    +--------------------------------------------+-------------------------+--------------------------------------------------------+
    | Error                                      | Cause                   | Action                                                 |
    +============================================+=========================+========================================================+
    | ``No module named 'pydseams'``             | Wrong interpreter       | ``python -m pip install pydseamslib``                  |
    +--------------------------------------------+-------------------------+--------------------------------------------------------+
    | ``pip install 'pydseamslib[ase]'``         | ASE extra missing       | install the extra                                      |
    +--------------------------------------------+-------------------------+--------------------------------------------------------+
    | ``periodic in all three directions``       | Partial ASE periodicity | set ``pbc=True`` on all three axes                     |
    +--------------------------------------------+-------------------------+--------------------------------------------------------+
    | ``cannot convert a singular cell``         | Zero-volume ASE cell    | supply a nonsingular cell                              |
    +--------------------------------------------+-------------------------+--------------------------------------------------------+
    | ``Hydrogen-bond analysis needs hydrogens`` | No H cloud              | ``bonded="cutoff"`` or pass H                          |
    +--------------------------------------------+-------------------------+--------------------------------------------------------+
    | ``chemfiles is not linked``                | Optional reader absent  | convert to LAMMPS/XYZ or check ``available_readers()`` |
    +--------------------------------------------+-------------------------+--------------------------------------------------------+

Getting help
------------

1. `FAQ <faq.rst>`_

2. `Python surface <../reference/python.rst>`_

3. `PydSEAMSlib <https://github.com/d-SEAMS/PydSEAMSlib>`_ issues

Include the ``pydseams`` version, Python version, file format, and exact
traceback when reporting an issue.
