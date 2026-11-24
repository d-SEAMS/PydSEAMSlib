==============
Confined water
==============


d-SEAMS analyses quasi-one-dimensional ice nanotubes as stacked prism
blocks (``Frame.find_prisms``) and quasi-two-dimensional ice monolayers
using primitive rings (``Frame.monolayer_rings``) and an in-plane radial
distribution function (``Frame.rdf_2d``).

Prerequisites
-------------

.. code:: bash

    pip install 'pydseamslib==2.14.2'
    mkdir -p input/traj
    curl -sLo input/traj/exampleTraj.lammpstrj \
      https://raw.githubusercontent.com/d-SEAMS/seams-core/main/input/traj/exampleTraj.lammpstrj

Run confined-water analyses
---------------------------

Unlike bulk ``chill_plus`` and ``cages``, ``rdf_2d``, ``monolayer_rings``, and
``find_prisms`` write engine output files under the directory passed as
the first argument (``"output"`` below) in addition to returning results
to Python:

.. code:: python

    import pydseams as ds

    frame = ds.read("input/traj/exampleTraj.lammpstrj", atom_type=2)
    r, g = frame.rdf_2d("output", cutoff=6.0, binwidth=0.5)
    mono = frame.monolayer_rings("output", sheet_area=1.0, max_depth=6)
    prisms = frame.find_prisms(output_dir="output", max_depth=4)
    print("rdf_2d bins:", len(r))
    print("monolayer:", mono)
    print("prisms:", prisms)
    assert len(r) == len(g) and len(r) > 0

``sheet_area`` sets the area used to normalise the monolayer ring
coverage fraction. Pass the lateral area of the sheet you loaded; the
raw ring counts do not depend on ``sheet_area``.
