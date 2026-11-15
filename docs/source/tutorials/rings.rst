===============
Primitive rings
===============


A primitive ring is a cycle that survives the Franzblau shortest-path test.
For every pair of vertices on the cycle, a shortest path between them lies on the cycle.
Distances are read on the whole graph, and only then are levels below the source cleared.
The hexagon with vertex 0 joined to vertices 1 and 4 is the worked example on the primitive-rings page.
``pydseams.band.ring_census`` counts those rings on the frame's default graph.
A type-2 read of this dump uses the hydrogen-bond list: one three-ring, 230 four-rings, four five-rings, and 192 six-rings.
``dseams.userframe.Frame.ring_census`` takes the graph by name; ``graph="cutoff"`` is the neighbour list inside 3.5 angstrom: 37 three-rings, 273 four-rings, 34 five-rings, and 167 six-rings.

``tests/dseams/docs/test_rings.py`` is the runnable copy of the hydrogen-bond census.

.. code:: python

    import pydseams as ds
    from pydseams.band import ring_census

    frame = ds.read("input/traj/exampleTraj.lammpstrj", atom_type=2)
    print(ring_census(frame, max_size=8))
    assert ring_census(frame, max_size=8) == {3: 1, 4: 230, 5: 4, 6: 192}
    from dseams import userframe

    cutoff = userframe.read("input/traj/exampleTraj.lammpstrj", atom_type=2).ring_census(8, graph="cutoff")
    print(cutoff)
    assert cutoff == {3: 37, 4: 273, 5: 34, 6: 167}

The CLI CHILL+ graph is the mutual four-nearest list, so a census of that graph is a different number.
