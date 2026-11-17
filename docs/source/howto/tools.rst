===============================
MDAnalysis, OVITO and notebooks
===============================


``pydseams.adapters`` connects the engine to the tools a trajectory
already lives in.

MDAnalysis
----------

``IceStates`` is an ``AnalysisBase`` subclass over an oxygen
``AtomGroup``; ions are a second group read against the assignment.
Boxes must be orthorhombic.

.. code:: python

    import MDAnalysis as mda
    from pydseams.adapters import IceStates

    u = mda.Universe("brine.gro", "brine.xtc")
    an = IceStates(u.select_atoms("name OW"), ions=u.select_atoms("name NA CL")).run()
    an.results.states        # frames x oxygens, codes in an.results.state_names
    an.results.features      # frames x len(an.results.names)
    an.results.ion_states    # frames x ions, codes in an.results.ion_state_names

OVITO
-----

``ice_states`` is a Python modifier function. Append it to a pipeline
and colour by the ``Ice state`` particle property (0 water, 1 cubic,
2 hexagonal, 3 mixed; ``-1`` outside the oxygen type).

.. code:: python

    from functools import partial
    from ovito.io import import_file
    from pydseams.adapters import ice_states

    pipeline = import_file("dump.lammpstrj")
    pipeline.modifiers.append(partial(ice_states, oxygen_type=1, cutoff=3.5))
    data = pipeline.compute()
    print(data.particles["Ice state"][...])

Triclinic OVITO cells pass through as LAMMPS bound spans and tilts.

Command line
------------

Without Python, ``seams cages dump.lammpstrj --per-atom labelled.lammpstrj``
appends a dump frame with a ``cage`` column (0 water, 1 hexagonal cage,
2 double-diamond cage, 3 both); ``fingerprint`` writes ``class`` or
``label`` and ``ions`` writes ``state``. Any visualiser that reads LAMMPS
dumps colours by that column.

Notebook
--------

``notebooks/classify_ice.ipynb`` installs the wheel, fetches one frame of
the cubic mW lattice and runs the assignment, CHILL+ and a fingerprint;
the README carries a Colab badge for it.
