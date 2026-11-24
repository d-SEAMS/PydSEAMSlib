========================
Ice and clathrate labels
========================


You label a mixed frame with CHILL+ and cages, then count closed
clathrate polyhedra by ring-face signature on an sI clathrate crystal.

Prerequisites
-------------

.. code:: bash

    pip install 'pydseamslib==2.14.2'
    mkdir -p input/traj
    curl -sLo input/traj/exampleTraj.lammpstrj \
      https://raw.githubusercontent.com/d-SEAMS/seams-core/main/input/traj/exampleTraj.lammpstrj
    curl -sLo input/traj/genice_sI.lammpstrj \
      https://raw.githubusercontent.com/d-SEAMS/seams-core/main/input/traj/genice_sI.lammpstrj

The mixed TIP4P frame
---------------------

``input/traj/exampleTraj.lammpstrj`` has 750 atoms and 250 oxygens
(LAMMPS type 2). ``Frame.chill_plus`` uses the cutoff neighbour list,
while ``Frame.cages`` uses the seeded cage search:

.. code:: python

    import pydseams as ds

    frame = ds.read("input/traj/exampleTraj.lammpstrj", atom_type=2)
    print(frame.chill_plus())
    print(frame.cages())
    assert str(frame.chill_plus()) == "IceCounts(interClathrate=9, water=241)"
    assert str(frame.cages()) == "CageScore(n_ih=0, n_ic=0, n_water=250)"

To classify with the mutual four-nearest neighbour list instead, use
``dseams.userframe.Frame.chill_plus(graph="mutual")``:

.. code:: python

    from dseams import userframe

    uframe = userframe.read("input/traj/exampleTraj.lammpstrj", atom_type=2)
    print(uframe.chill_plus(graph="mutual"))
    assert str(uframe.chill_plus(graph="mutual")) == "IceCounts(clathrate=1, interClathrate=7, water=242)"

Clathrate cages by signature
----------------------------

A cage signature counts closed polyhedra of a given ring-face census,
such as ``"512"`` (12 pentagonal faces) or ``"51262"`` (12 pentagonal and
2 hexagonal faces). On ``input/traj/genice_sI.lammpstrj``,
``dseams.userframe.Frame.cages_by_signature("512")`` finds 2 cages
spanning 40 vertices.

The default face-edge length window is ``2.65,3.15`` in the coordinate
units of the file (angstroms for these dumps), covering standard
oxygen-oxygen contacts. Passing ``edge="none"`` disables the edge-length
filter so any ring accepted by the neighbour list can serve as a face.
On reference sI and sII clathrate cells, the default window and
``edge="none"`` return the same cages.

.. code:: python

    s1 = userframe.read("input/traj/genice_sI.lammpstrj", atom_type=1)
    found = s1.cages_by_signature("512")
    assert len(found) == 2 and found.n_atoms == 40

    found_unbounded = s1.cages_by_signature("512", edge="none")
    assert len(found_unbounded) == 2 and found_unbounded.n_atoms == 40

    found_6 = s1.cages_by_signature("512", max_ring_size=6)
    found_51262 = s1.cages_by_signature("51262", max_ring_size=8)
    print("5^12 cages:", len(found_6), "5^12 6^2 cages:", len(found_51262))
