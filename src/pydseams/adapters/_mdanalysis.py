"""d-SEAMS ice states as an MDAnalysis analysis.

MDAnalysis is imported when the class is built, so importing this module
without MDAnalysis installed is fine.
"""

from __future__ import annotations

import numpy as np

from pydseams.features import (
    ION_FEATURE_NAMES,
    ION_FRONT,
    ION_ICE,
    ION_LIQUID,
    ION_STATE_NAMES,
    STATE_NAMES,
    IceFeaturizer,
    ion_environment,
)
from pydseams.frame import Frame


def _box_from_dimensions(dimensions):
    a, b, c, alpha, beta, gamma = (float(x) for x in dimensions)
    if abs(alpha - 90.0) > 1e-6 or abs(beta - 90.0) > 1e-6 or abs(gamma - 90.0) > 1e-6:
        raise ValueError(
            "IceStates reads orthorhombic boxes; write the frame through "
            "MDAnalysis.transformations or pass a triclinic cell to Frame.from_arrays"
        )
    return [a, b, c]


def _analysis_base():
    from MDAnalysis.analysis.base import AnalysisBase

    return AnalysisBase


class IceStates:
    """Per-molecule ice state of an oxygen AtomGroup over a trajectory.

    Parameters
    ----------
    water : MDAnalysis.core.groups.AtomGroup
        The water oxygens (one atom per molecule).
    ions : MDAnalysis.core.groups.AtomGroup, optional
        Ions read against the assignment by their first water shell.
    cutoff : float, optional
        Neighbour cutoff in Angstrom. Default ``3.5``.
    k, ring_adjacent, chill
        As in :class:`pydseams.features.IceFeaturizer`.

    After ``run()``, ``results.states`` holds one ``int8`` row per frame
    with a :data:`pydseams.features.STATE_NAMES` code per oxygen,
    ``results.features`` the per-frame vector and ``results.names`` its
    labels; with ions, ``results.ion_states`` holds one row per frame with
    a :data:`pydseams.features.ION_STATE_NAMES` code per ion.
    """

    def __new__(
        cls, water, ions=None, cutoff=3.5, k=4, ring_adjacent=True, chill=True, **kwargs
    ):
        base = _analysis_base()

        class _IceStates(base):
            def __init__(self, water, ions, cutoff, k, ring_adjacent, chill, **kw):
                super().__init__(water.universe.trajectory, **kw)
                self.water = water
                self.ions = ions
                self.cutoff = float(cutoff)
                self.k = int(k)
                self.ring_adjacent = bool(ring_adjacent)
                self.chill = bool(chill)

            def _prepare(self):
                n_water = len(self.water)
                self._numbers = [1] * n_water
                if self.ions is not None and len(self.ions):
                    self._numbers = self._numbers + [2] * len(self.ions)
                self._states = []
                self._features = []
                self._ion_states = []
                self.results.names = None

            def _frame(self):
                ts = self.water.universe.trajectory.ts
                box = _box_from_dimensions(ts.dimensions)
                pos = np.asarray(self.water.positions, dtype=float)
                if self.ions is not None and len(self.ions):
                    pos = np.concatenate(
                        [pos, np.asarray(self.ions.positions, dtype=float)], axis=0
                    )
                box_arr = np.asarray(box, dtype=float)
                wrapped = np.remainder(pos, box_arr)
                return Frame.from_arrays(
                    wrapped, box, numbers=self._numbers, cutoff=self.cutoff
                )

            def _single_frame(self):
                frame = self._frame()
                feat = IceFeaturizer(
                    frame,
                    k=self.k,
                    ring_adjacent=self.ring_adjacent,
                    chill=self.chill,
                    ion_types=(),
                )
                x, states = feat.frame_features()
                names = list(feat.feature_names)
                x = np.asarray(x, dtype=float)
                if self.ions is not None and len(self.ions):
                    _, shell, fraction, ion_states = ion_environment(
                        frame, states, (2,), cutoff=self.cutoff
                    )
                    ion_x = np.array(
                        [
                            int((ion_states == ION_ICE).sum()),
                            int((ion_states == ION_FRONT).sum()),
                            int((ion_states == ION_LIQUID).sum()),
                            float(fraction[shell > 0].mean())
                            if (shell > 0).any()
                            else 0.0,
                        ],
                        dtype=float,
                    )
                    x = np.concatenate([x, ion_x])
                    names += list(ION_FEATURE_NAMES)
                    self._ion_states.append(np.asarray(ion_states, dtype=np.int8))
                if self.results.names is None:
                    self.results.names = names
                n_water = len(self.water)
                self._states.append(np.asarray(states[:n_water], dtype=np.int8))
                self._features.append(x)

            def _conclude(self):
                self.results.states = np.array(self._states, dtype=np.int8)
                self.results.features = np.array(self._features, dtype=float)
                self.results.state_names = STATE_NAMES
                if self._ion_states:
                    self.results.ion_states = np.array(self._ion_states, dtype=np.int8)
                    self.results.ion_state_names = ION_STATE_NAMES

        return _IceStates(water, ions, cutoff, k, ring_adjacent, chill, **kwargs)
