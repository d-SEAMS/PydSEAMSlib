"""Suffix dispatch onto the compiled I/O readers.

The compiled module stays thin. This helper picks LAMMPS, XYZ,
chemfiles, or readcon from the path suffix and returns a
:class:`~pydseams.frame.Frame`.
"""

from __future__ import annotations

from pathlib import Path

from .frame import Frame

_CHEMFILES = {".pdb", ".gro", ".dcd"}
_LAMMPS = {".lammpstrj", ".dump", ".lammps"}


def read(path, frame=1, **kwargs):
    """Load one configuration. Format follows the file suffix.

    Parameters
    ----------
    path : path-like
        Trajectory or structure file.
    frame : int, optional
        1-indexed frame for multi-frame formats (LAMMPS, chemfiles,
        ``.con``). Ignored for XYZ. Default ``1``.
    **kwargs
        Forwarded to the matching :class:`~pydseams.frame.Frame`
        constructor (``cutoff``, ``bonded``, ``atom_type``,
        ``region``, ``all_atoms``, ...). ``all_atoms=True`` retains every
        LAMMPS type for mixed-site analyses such as ionic pairs and domains.

    Returns
    -------
    Frame
        Analysable configuration.

    Notes
    -----
    Suffix dispatch:

    * ``.xyz`` -- :meth:`Frame.from_xyz`
    * ``.con`` -- :meth:`Frame.from_con`
    * ``.pdb``, ``.gro``, ``.dcd`` -- :meth:`Frame.from_chemfiles`
    * otherwise -- :meth:`Frame.from_file` (LAMMPS dump)

    ``.lammpstrj``, ``.dump``, and ``.lammps`` take the LAMMPS path.
    Builds without chemfiles or readcon raise ``RuntimeError`` when
    those suffixes are used.
    """
    suffix = Path(path).suffix.lower()
    if suffix == ".xyz":
        return Frame.from_xyz(path, **kwargs)
    if suffix == ".con":
        return Frame.from_con(path, frame=frame, **kwargs)
    if suffix in _CHEMFILES:
        return Frame.from_chemfiles(path, frame=frame, **kwargs)
    return Frame.from_file(path, frame=frame, **kwargs)


def has_chemfiles():
    """True when this build can read PDB, GRO, and DCD through chemfiles."""
    return bool(available_readers()["chemfiles"])


def has_readcon():
    """True when this build can read an eOn ``.con`` file through readcon-core."""
    return bool(available_readers()["readcon"])


def available_readers():
    """Report which compiled I/O readers this build linked.

    Returns
    -------
    dict of str to bool
        ``lammps`` is always ``True``. ``xyz``, ``chemfiles``, and
        ``readcon`` are ``True`` when the matching :mod:`pydseams.yoda`
        symbol exists (``readXYZ``, ``readChemfiles``, ``readCon``).
    """
    from . import yoda

    return {
        "lammps": True,
        "xyz": hasattr(yoda, "readXYZ"),
        "chemfiles": hasattr(yoda, "readChemfiles"),
        "readcon": hasattr(yoda, "readCon"),
    }
