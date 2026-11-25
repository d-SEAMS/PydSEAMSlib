=============
Paper figures
=============


The 2020 d-SEAMS paper figures use datasets hosted on figshare
(project 73545). This guide runs the same Python classification calls
on the bulk cubic ice, sI clathrate, and mixed TIP4P reference
structures.

Prerequisites
-------------

.. code:: bash

    pip install 'pydseamslib==2.14.2'
    mkdir -p input/traj
    curl -fsSL -o input/fetch.py \
      https://raw.githubusercontent.com/d-SEAMS/seams-core/v2.14.2/input/fetch.py
    curl -sLo input/traj/genice_sI.lammpstrj \
      https://raw.githubusercontent.com/d-SEAMS/seams-core/v2.14.2/input/traj/genice_sI.lammpstrj

Reference summary
-----------------

.. table::

    +-------------------------+----------------------------------------+-------------------------------+-------------------------------------------+
    | figure                  | dump                                   | call                          | result                                    |
    +=========================+========================================+===============================+===========================================+
    | cubic ice, CHILL+       | ``mW_cubic.lammpstrj`` frame 1, type 1 | ``chill_plus()``              | 4096 cubic                                |
    +-------------------------+----------------------------------------+-------------------------------+-------------------------------------------+
    | cubic ice, cages        | ``mW_cubic.lammpstrj`` frame 1, type 1 | ``cages()``                   | ``n_ic`` 4096                             |
    +-------------------------+----------------------------------------+-------------------------------+-------------------------------------------+
    | sI :math:`5^{12}` cages | ``genice_sI.lammpstrj``, type 1        | ``cages_by_signature("512")`` | 2 cages, 40 atoms                         |
    +-------------------------+----------------------------------------+-------------------------------+-------------------------------------------+
    | mixed frame             | ``exampleTraj.lammpstrj``, type 2      | cutoff CHILL+ vs cages        | 9 interfacial clathrate, no complete cage |
    +-------------------------+----------------------------------------+-------------------------------+-------------------------------------------+

Run the bulk and clathrate checks
---------------------------------

``Frame.chill_plus`` and ``Frame.cages`` run on the cubic mW cell, while
``dseams.userframe.Frame.cages_by_signature`` counts the sI :math:`5^{12}`
cages and their distinct vertex atoms. ``input/fetch.py`` downloads the
cubic mW cell that seams-core v2.14.2 pins and checks its SHA-256; a
verified local copy is reused:

.. code:: python

    import subprocess
    import sys

    import pydseams as ds
    from dseams import userframe

    subprocess.run(
        [sys.executable, "input/fetch.py", "mW_cubic.lammpstrj"],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    ice = ds.read("input/traj/mW_cubic.lammpstrj", frame=1, atom_type=1)
    print("cubic:", ice.chill_plus().cubic, "ddc:", ice.cages().n_ic)
    assert ice.chill_plus().cubic == 4096 and ice.cages().n_ic == 4096

    guest = userframe.read("input/traj/genice_sI.lammpstrj", atom_type=1)
    found = guest.cages_by_signature("512")
    print("5^12 cages:", len(found), "atoms:", found.n_atoms)
    assert len(found) == 2
    assert found.n_atoms == 40

Point the same calls (``chill_plus``, ``cages``, ``cages_by_signature``,
``find_prisms``, ``monolayer_rings``, and ``rdf_2d``) at the figshare
project 73545 trajectories to reproduce the nucleation, nanotube, and
monolayer figures.
