====================================
Features for kinetic models and ions
====================================

Problem
-------

You have a nucleation trajectory and want one feature vector per frame
for a Markov state model or a committor fit, one integer state per
molecule per frame, and the first-shell class of each ion when ions
sit in the cloud.

Per-frame features
------------------

#skip_lint_start

.. code:: python

    from pydseams import Trajectory
    from pydseams.features import IceFeaturizer

    traj = Trajectory("nucleation.lammpstrj", frame=1, atom_type=2, cutoff=3.5)
    feat = IceFeaturizer(traj, ring_adjacent=True)
    X, S = feat.transform()
    print(feat.feature_names)

#skip_lint_end

``IceFeaturizer`` calls ``Frame.seeded_affiliation`` on each frame.
Ring-adjacent completion is the featurizer default: it fills the last
vertex of a six-ring whose other vertices already carry a label.
``Frame.cages`` and ``Frame.seeded_affiliation`` leave that flag off unless
you pass it.

``X`` lists cage counts (``n_ice``, ``n_ic``, ``n_ih``, ``n_mixed``), the
largest connected cage cluster ``n_max`` and the cluster count, the
cubicity, the ``chill_plus`` counts on the cutoff graph, the largest
``chill_plus`` bulk cluster, and the six-ring count. ``S`` lists
``STATE_WATER``, ``STATE_IC``, ``STATE_IH``, or ``STATE_MIXED`` per molecule.
Both arrays are plain NumPy and both are deterministic: a feature file
from one machine equals the same file from another to the last integer.

deeptime
--------

#skip_lint_start

.. code:: python

    from pydseams.features import discretize_nmax, to_deeptime

    dtrajs = discretize_nmax(X[:, 1], edges=[10, 50, 150, 400])
    msm = to_deeptime(X, lagtime=5)

#skip_lint_end

``discretize_nmax`` bins the largest-cluster series into integer states
for ``deeptime.markov.msm.MaximumLikelihoodMSM``. ``to_deeptime`` fits a
time-lagged independent component analysis (TICA) on the full vector
and returns the fitted model.

PyEMMA
------

#skip_lint_start

.. code:: python

    from pydseams.features import to_pyemma_featurizer

    featurizer = to_pyemma_featurizer(feat, topology="nucleation.pdb")

#skip_lint_end

The per-frame vector registers as a custom feature. PyEMMA is
unmaintained; deeptime succeeds it.

Ions
----

Ions sit outside the hydrogen-bond network. The cage assignment runs
on the water and the ions are read against it. Build the frame with
every species in the cloud and name the ion types.

Pass a sequence to ``Frame.from_ase`` when the ASE ``Atoms`` mix water and
salt. The listed species stay in the cloud. The first entry is the
analysed water:

#skip_lint_start

.. code:: python

    import pydseams as ds

    frame = ds.from_ase(atoms, select=("O", "Na", "Cl"), bonded="cutoff")

#skip_lint_end

``select=("O", "Na", "Cl")`` keeps oxygen, sodium, and chlorine. The
analysed species is oxygen (atomic number 8). Sodium and chlorine keep
their atomic numbers as ``c_type``. ``pydseams.features.ion_environment``
and ``Frame.ion_environment`` read those codes.

A LAMMPS dump stores integer type codes, not atomic numbers. Keep every
type with ``all_atoms`` on and pass those codes as ``ion_types``:

#skip_lint_start

.. code:: python

    from pydseams import Trajectory
    from pydseams.features import ION_ICE, IceFeaturizer, ion_environment

    traj = Trajectory("brine.lammpstrj", frame=1, atom_type=1, all_atoms=True)
    feat = IceFeaturizer(traj, ion_types=(3, 4))  # LAMMPS types
    x, states = feat.frame_features()
    ions, shell, fraction, ion_states = ion_environment(traj, states, (3, 4))
    trapped = ions[ion_states == ION_ICE]

#skip_lint_end

Those LAMMPS type codes classify each ion from its first water
shell. Water is type 1 on this dump; 3 and 4 are the ions, not atomic
numbers (Na is 11, Cl is 17). The return is the water count within
``cutoff``, the ice fraction of that shell, and a class. ``ION_ICE``
means every neighbour in the shell is ice. ``ION_LIQUID`` means none of
it is. ``ION_FRONT`` is the rest. ``ion_types`` adds ``n_ion_ice``,
``n_ion_front``, ``n_ion_liquid`` and the mean shell ice fraction to the
feature vector.

``Frame.ion_environment`` is the compiled first-shell class on one
frame. That method passes ``ring_adjacent`` through to
``seeded_affiliation``.

See also
--------

`Classify ice <../tutorials/classify-ice.rst>`_
    ``cages`` / ``seeded_affiliation``

`Classify ASE Atoms <ase.rst>`_
    ``from_ase`` sequence ``select``

`Python surface <../reference/python.rst>`_
    live names
