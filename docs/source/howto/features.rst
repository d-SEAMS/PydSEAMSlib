=====================================
Features for kinetic models and ions
=====================================

Problem
-------

A nucleation trajectory needs one feature vector per frame for a
Markov state model or a committor fit, one integer state per
molecule per frame, and, when ions sit in the cloud, a class for
what each ion sees.

Per-frame features
------------------

.. code:: python

    from pydseams import Trajectory
    from pydseams.features import IceFeaturizer

    traj = Trajectory("nucleation.lammpstrj", frame=1, atom_type=2, cutoff=3.5)
    feat = IceFeaturizer(traj)
    X, S = feat.transform()
    print(feat.feature_names)

``IceFeaturizer`` turns on ring-adjacent completion of the seeded
assignment by default (``ring_adjacent``). That is the same flag
``Frame.seeded_affiliation`` accepts. ``X`` has shape
``(n_frames, 15)``. ``S`` has shape ``(n_frames, n_molecules)``.

The columns of ``X`` name cage counts (``n_ice``, ``n_ic``,
``n_ih``, ``n_mixed``), the largest connected cage cluster
``n_max`` and the cluster count, the cubicity, the Moore and
Molinero ice labels (CHILL) plus counts on the cutoff graph, the
largest CHILL+ bulk cluster and the six-ring count. ``S`` stores
``STATE_WATER``, ``STATE_IC``, ``STATE_IH`` or ``STATE_MIXED`` per
molecule. Both arrays are plain NumPy and deterministic: a feature
file from one machine equals the same file from another to the last
integer.

deeptime
--------

.. code:: python

    from pydseams.features import discretize_nmax, to_deeptime

    dtrajs = discretize_nmax(X[:, 1], edges=[10, 50, 150, 400])
    msm = to_deeptime(X, lagtime=5)

``discretize_nmax`` bins the largest-cluster series into integer
states for ``deeptime.markov.msm.MaximumLikelihoodMSM``.
``to_deeptime`` fits a time lagged independent component analysis (TICA)
on the full vector and returns the estimator.

PyEMMA
------

.. code:: python

    from pydseams.features import to_pyemma_featurizer

    featurizer = to_pyemma_featurizer(feat, topology="nucleation.pdb")

The per-frame vector registers as a custom feature. deeptime is the
successor of PyEMMA; both adapters remain because both still appear
in the nucleation literature.

Ions
----

Ions sit outside the hydrogen-bond network. Cage assignment runs on
the water; the ions are then read against that assignment. Build the
frame with every species in the cloud and name the ion types.

For ASE input, pass a sequence to ``select`` so oxygen and the ions
stay in the cloud. Analysis uses the first listed species.

.. code:: python

    import pydseams as ds

    frame = ds.from_ase(atoms, select=("O", "Na", "Cl"), bonded="cutoff")

    from pydseams import Trajectory
    from pydseams.features import ION_ICE, IceFeaturizer, ion_environment

    traj = Trajectory("brine.lammpstrj", frame=1, atom_type=1, all_atoms=1)
    feat = IceFeaturizer(traj, ion_types=(3, 4))
    x, states = feat.frame_features()
    ions, shell, fraction, ion_states = ion_environment(traj, states, (3, 4))
    trapped = ions[ion_states == ION_ICE]

Those ion species then go through ``ion_environment``, which returns
four arrays per ion: the water count in the first shell (within
``cutoff``), the fraction of that shell carrying an ice label, and a
class. ``ION_ICE`` means every
shell water carries ice. ``ION_LIQUID`` means none of it does.
``ION_FRONT`` is the rest. With ``ion_types`` set the feature
vector gains ``n_ion_ice``, ``n_ion_front``, ``n_ion_liquid`` and
the mean shell ice fraction.

See also
--------

`Python surface <../reference/python.rst>`_
    ``IceFeaturizer``, ``ion_environment``, ``from_ase`` select

`API <../api.md>`_
    Napoleon autodoc of ``pydseams.features``

`Classify ASE Atoms <ase.rst>`_
    ``select`` as a symbol, number, sequence, or ``None``
