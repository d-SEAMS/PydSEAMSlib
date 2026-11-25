===============
Primitive rings
===============


A primitive ring is a cycle that satisfies the Franzblau shortest-path
criterion: for every pair of vertices on the cycle, a shortest path
between them on the full graph lies along the cycle.

Prerequisites
-------------

.. code:: bash

    pip install pydseamslib
    mkdir -p input/traj
    curl -sLo input/traj/exampleTraj.lammpstrj \
      https://raw.githubusercontent.com/d-SEAMS/seams-core/v2.14.2/input/traj/exampleTraj.lammpstrj

Hydrogen-bond vs cutoff ring census
-----------------------------------

``pydseams.band.ring_census`` counts primitive rings on a ``Frame``'s
default bonded graph. When ``ds.read`` loads LAMMPS type 2 (oxygen) from
a dump containing hydrogens, the default graph is the hydrogen-bond
network; on ``exampleTraj.lammpstrj`` up to ring size 8, that gives 1
three-ring, 230 four-rings, 4 five-rings, and 192 six-rings.

``dseams.userframe.Frame.ring_census`` lets you choose the graph
explicitly. Passing ``graph="cutoff"`` builds the 3.5 A distance-cutoff
neighbour graph instead, giving 37 three-rings, 273 four-rings, 34
five-rings, and 167 six-rings on the same oxygens.

.. code:: python

    import pydseams as ds
    from pydseams.band import ring_census
    from dseams import userframe

    frame = ds.read("input/traj/exampleTraj.lammpstrj", atom_type=2)
    hbond_counts = ring_census(frame, max_size=8)
    print("hbond rings:", hbond_counts)
    assert hbond_counts == {3: 1, 4: 230, 5: 4, 6: 192}

    cutoff_counts = userframe.read(
        "input/traj/exampleTraj.lammpstrj", atom_type=2
    ).ring_census(8, graph="cutoff")
    print("cutoff rings:", cutoff_counts)
    assert cutoff_counts == {3: 37, 4: 273, 5: 34, 6: 167}
