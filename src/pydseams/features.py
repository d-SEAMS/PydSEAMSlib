"""Per-frame and per-molecule ice features for kinetic modelling.

:class:`IceFeaturizer` walks a :class:`~pydseams.frame.Trajectory` and
returns what a Markov-state-model or committor pipeline consumes: one
feature vector per frame (cage counts, the largest connected cage
cluster, cubicity, CHILL+ counts, ring counts) and one integer state per
molecule per frame (water, ice Ic, ice Ih, mixed). The arrays are plain
NumPy, which is what deeptime takes directly; :func:`to_pyemma_featurizer`
registers the same per-frame vector as a custom feature with PyEMMA.
PyEMMA has been in maintenance mode since 2022 and deeptime is its
successor; both adapters are here because both still appear in the
nucleation literature.

Everything is deterministic and label-based, so a feature file from one
machine equals the same file from another to the last integer.
"""

from __future__ import annotations

import numpy as np

from . import yoda

STATE_WATER, STATE_IC, STATE_IH, STATE_MIXED = 0, 1, 2, 3
STATE_NAMES = ("water", "ic", "ih", "mixed")

FEATURE_NAMES = (
    "n_ice",
    "n_max",
    "n_clusters",
    "n_ic",
    "n_ih",
    "n_mixed",
    "cubicity",
    "chill_cubic",
    "chill_hexagonal",
    "chill_interfacial",
    "chill_clathrate",
    "chill_interclathrate",
    "chill_water",
    "chill_max",
    "six_rings",
)

ION_LIQUID, ION_FRONT, ION_ICE = 0, 1, 2
ION_STATE_NAMES = ("liquid", "front", "ice")
ION_FEATURE_NAMES = (
    "n_ion_ice",
    "n_ion_front",
    "n_ion_liquid",
    "ion_shell_ice",
)


def ion_environment(frame, states, ion_types, cutoff=None):
    """Ice content of the first water shell of every ion.

    Ions are not part of the hydrogen-bond network; the cage assignment
    runs on the water alone and the ions are read against it. A frame
    built with ``all_atoms=True`` or :meth:`Frame.from_arrays` carries
    every species in its cloud; ``frame.atom_type`` names the water. The
    engine's ``ionEnvironment`` does the work when the compiled module
    exposes it; the NumPy path below gives the same answer on older
    engines.

    Parameters
    ----------
    frame : Frame
        Frame whose cloud holds water and ions.
    states : array of int8
        Per-atom states from :meth:`IceFeaturizer.frame_features`;
        non-water entries are ignored.
    ion_types : iterable of int
        ``c_type`` codes of the ions.
    cutoff : float, optional
        First-shell radius in Angstrom. Default ``frame.cutoff``.

    Returns
    -------
    (indices, shell, fraction, ion_states)
        Cloud indices of the ions, the water count in each first shell,
        the fraction of that shell carrying an ice label, and one of
        ``ION_ICE`` (whole shell ice), ``ION_LIQUID`` (no ice) or
        ``ION_FRONT`` per ion. An ion with an empty shell is liquid.
    """
    cut = float(frame.cutoff if cutoff is None else cutoff)
    pts = frame.cloud.pts
    types = np.array([p.c_type for p in pts])
    ions = np.nonzero(np.isin(types, list(ion_types)))[0]
    states = np.asarray(states)
    if hasattr(yoda, "ionEnvironment"):
        env = yoda.ionEnvironment(
            frame.cloud,
            [bool(s != STATE_WATER) for s in states],
            [int(i) for i in ions],
            int(frame.atom_type),
            cut,
        )
        ion_states = np.array(
            [
                ION_ICE
                if s == yoda.IonState.ice
                else ION_FRONT
                if s == yoda.IonState.front
                else ION_LIQUID
                for s in env.state
            ],
            dtype=np.int8,
        )
        return (
            np.asarray(env.ion, dtype=int),
            np.asarray(env.shell, dtype=int),
            np.asarray(env.iceFraction, dtype=float),
            ion_states,
        )
    pos = np.array([[p.x, p.y, p.z] for p in pts], dtype=float)
    box = np.asarray(frame.box, dtype=float)[:3]
    water = np.nonzero(types == frame.atom_type)[0]
    ice_water = states[water] != STATE_WATER
    water_pos = pos[water]
    shell = np.zeros(len(ions), dtype=int)
    fraction = np.zeros(len(ions), dtype=float)
    ion_states = np.full(len(ions), ION_LIQUID, dtype=np.int8)
    for j, i in enumerate(ions):
        d = water_pos - pos[i]
        d -= box * np.round(d / box)
        near = np.sqrt((d**2).sum(axis=1)) < cut
        shell[j] = int(near.sum())
        if shell[j] == 0:
            continue
        fraction[j] = float(ice_water[near].mean())
        if fraction[j] >= 1.0:
            ion_states[j] = ION_ICE
        elif fraction[j] > 0.0:
            ion_states[j] = ION_FRONT
    return ions, shell, fraction, ion_states


def ion_features(frame, states, ion_types, cutoff=None):
    """Per-frame ion summary in the order of :data:`ION_FEATURE_NAMES`."""
    _, shell, fraction, ion_states = ion_environment(frame, states, ion_types, cutoff)
    return np.array(
        [
            int((ion_states == ION_ICE).sum()),
            int((ion_states == ION_FRONT).sum()),
            int((ion_states == ION_LIQUID).sum()),
            float(fraction[shell > 0].mean()) if (shell > 0).any() else 0.0,
        ],
        dtype=float,
    )


def count_frames(filename):
    """Stored frames in a LAMMPS dump (``ITEM: TIMESTEP`` headers)."""
    n = 0
    with open(filename, "rb") as fh:
        for line in fh:
            if line.startswith(b"ITEM: TIMESTEP"):
                n += 1
    return n


def _components(flags, rows):
    """Largest component and component count of flagged atoms over an
    index graph whose rows lead with the atom itself."""
    n = len(flags)
    parent = list(range(n))

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i in range(n):
        if not flags[i] or i >= len(rows):
            continue
        for j in rows[i][1:]:
            if 0 <= j < n and flags[j]:
                a, b = root(i), root(j)
                if a != b:
                    parent[b] = a
    sizes = {}
    for i in range(n):
        if flags[i]:
            r = root(i)
            sizes[r] = sizes.get(r, 0) + 1
    if not sizes:
        return 0, 0
    return max(sizes.values()), len(sizes)


class IceFeaturizer:
    """Turn a trajectory into per-frame features and per-molecule states.

    Parameters
    ----------
    trajectory : Trajectory
        An open :class:`~pydseams.frame.Trajectory`; frames are visited
        with :meth:`~pydseams.frame.Frame.load_frame`, so the
        incremental updaters stay warm.
    k : int, optional
        Neighbours in the four-nearest graphs. Default ``4``.
    ring_adjacent : bool, optional
        Ring-adjacent completion of the seeded assignment. Default
        ``True``.
    chill : bool, optional
        Also compute CHILL+ on the cutoff graph. Default ``True``.
    ion_types : iterable of int, optional
        ``c_type`` codes of ions present in the cloud. When given, the
        vector gains :data:`ION_FEATURE_NAMES` from
        :func:`ion_features`.
    """

    def __init__(self, trajectory, k=4, ring_adjacent=True, chill=True, ion_types=()):
        self.trajectory = trajectory
        self.k = k
        self.ring_adjacent = ring_adjacent
        self.chill = chill
        self.ion_types = tuple(int(t) for t in ion_types)

    @property
    def feature_names(self):
        names = list(FEATURE_NAMES)
        if self.ion_types:
            names += list(ION_FEATURE_NAMES)
        return names

    def frame_features(self):
        """Feature vector and per-molecule states of the loaded frame.

        Returns
        -------
        (features, states)
            ``features`` is a float array of length ``len(self.feature_names)``;
            ``states`` an ``int8`` array with one of ``STATE_*`` per
            molecule.
        """
        t = self.trajectory
        n = t.n_atoms
        score = t.seeded_affiliation(k=self.k, ring_adjacent=self.ring_adjacent)
        union = score.union
        six = int(score.n_six) if score.n_six is not None else 0
        hc = np.asarray(score.hc, dtype=bool)
        ddc = np.asarray(score.ddc, dtype=bool)
        states = np.zeros(n, dtype=np.int8)
        states[ddc & ~hc] = STATE_IC
        states[hc & ~ddc] = STATE_IH
        states[hc & ddc] = STATE_MIXED
        ice = hc | ddc
        n_max, n_clusters = _components(ice.tolist(), union)
        n_ic = int((states == STATE_IC).sum())
        n_ih = int((states == STATE_IH).sum())
        n_mixed = int((states == STATE_MIXED).sum())
        n_ice = int(ice.sum())
        cubicity = (n_ic + n_mixed) / n_ice if n_ice else 0.0

        chill = dict.fromkeys(
            (
                "cubic",
                "hexagonal",
                "interfacial",
                "clathrate",
                "interClathrate",
                "water",
            ),
            0,
        )
        chill_max = 0
        if self.chill:
            counts = t.chill_plus()
            for key in chill:
                chill[key] = int(counts.get(key, 0))
            bulk = [
                pt.iceType.name in ("cubic", "hexagonal", "reCubic", "reHex")
                for pt in cloud.pts
            ]
            chill_max, _ = _components(bulk, t.bonds_by_index)

        features = np.array(
            [
                n_ice,
                n_max,
                n_clusters,
                n_ic,
                n_ih,
                n_mixed,
                cubicity,
                chill["cubic"],
                chill["hexagonal"],
                chill["interfacial"],
                chill["clathrate"],
                chill["interClathrate"],
                chill["water"],
                chill_max,
                six,
            ],
            dtype=float,
        )
        if self.ion_types:
            features = np.concatenate(
                [features, ion_features(t, states, self.ion_types)]
            )
        return features, states

    def transform(self, frames=None):
        """Features and states for a range of frames.

        Parameters
        ----------
        frames : iterable of int or None
            1-based frame numbers; ``None`` walks the stored frames from
            the current one to the last.

        Returns
        -------
        (X, S)
            ``X`` has shape ``(n_frames, len(self.feature_names))``; ``S`` has
            shape ``(n_frames, n_molecules)`` with ``int8`` states.
        """
        t = self.trajectory
        if frames is None:
            frames = range(t.frame, count_frames(t.filename) + 1)
        X = []
        S = []
        for f in frames:
            if f != t.frame:
                t.load_frame(f)
            x, s = self.frame_features()
            X.append(x)
            S.append(s)
        return np.vstack(X), np.vstack(S)


def discretize_nmax(n_max, edges):
    """Map a largest-cluster series onto integer states by bin edges.

    Parameters
    ----------
    n_max : array_like
        Largest cage cluster per frame (``X[:, 1]``).
    edges : sequence of float
        Increasing bin edges; state ``i`` is ``edges[i] <= n_max < edges[i+1]``,
        below the first edge is state 0, at or above the last is the last.

    Returns
    -------
    numpy.ndarray of int
        A discrete trajectory for :class:`deeptime.markov.msm.MaximumLikelihoodMSM`
        or ``pyemma.msm.estimate_markov_model``.
    """
    return np.searchsorted(
        np.asarray(edges, dtype=float), np.asarray(n_max, dtype=float), side="right"
    ).astype(int)


def to_pyemma_featurizer(featurizer, topology):
    """Register the per-frame vector as a PyEMMA custom feature.

    PyEMMA's featurizer wants a function of an ``mdtraj`` trajectory; this
    ignores the coordinates it is given and reads the frames through the
    engine instead, so the feature order is the one of
    :attr:`FEATURE_NAMES`.
    """
    import pyemma  # noqa: F401  (optional, unmaintained since 2022)

    feat = pyemma.coordinates.featurizer(topology)

    def _fn(traj):
        X, _ = featurizer.transform(range(1, traj.n_frames + 1))
        return X.astype(np.float32)

    feat.add_custom_func(
        _fn, dim=len(FEATURE_NAMES), description=",".join(FEATURE_NAMES)
    )
    return feat


def to_deeptime(X, lagtime, dim=None):
    """TICA of a per-frame feature array with deeptime.

    Returns the fitted :class:`deeptime.decomposition.TICA` model; call
    ``model.transform(X)`` for the projection.
    """
    from deeptime.decomposition import TICA

    return TICA(lagtime=lagtime, dim=dim).fit(np.asarray(X, dtype=float)).fetch_model()
