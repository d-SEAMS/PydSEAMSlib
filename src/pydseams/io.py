"""Suffix dispatch onto the compiled I/O readers.

The compiled module stays thin. This helper picks LAMMPS, XYZ,
chemfiles, or readcon from the path suffix and returns a
:class:`~pydseams.frame.Frame`.
"""

from __future__ import annotations

from pathlib import Path

from .frame import Frame

_CHEMFILES = {".pdb", ".gro", ".dcd", ".xtc", ".trr"}
_LAMMPS = {".lammpstrj", ".dump", ".lammps"}


def iter_trajectory(
    path,
    *,
    atom_name,
    topology=None,
    residue_name=None,
    hydrogen_names=(),
    start=0,
    stop=None,
    stride=1,
    cutoff=None,
    bonded="auto",
):
    """Yield frames selected by atom and residue names through readcon.

    Coordinate-only formats such as XTC and TRR take a labelled topology
    with the same atom count and order. Indices follow Python slicing.
    Each selected residue contributes one centre and its named hydrogens.
    Coordinates are in angstroms. Frames retain the full periodic cell.
    Install the readcon extra to use this reader.
    """
    import numpy as np
    import readcon
    from ase.cell import Cell
    from ase.geometry import wrap_positions

    from .aseio import _dump_box

    if start < 0 or stride < 1 or (stop is not None and stop < 0):
        raise ValueError(
            "start and stop must be nonnegative and stride must be positive"
        )
    if isinstance(hydrogen_names, str):
        hydrogen_names = (hydrogen_names,)
    hydrogens = set(hydrogen_names)
    if atom_name in hydrogens:
        raise ValueError("centre and hydrogen names must be distinct")
    filename = str(path)
    total = readcon.read_chemfiles_nsteps(filename)
    limit = total if stop is None else min(stop, total)
    source = None if topology is None else str(topology)
    for first in range(start, limit, stride * 32):
        batch = readcon.read_chemfiles(
            filename,
            topology=source,
            skip=first,
            step=stride,
            stop=min(first + stride * 32, limit),
        )
        for offset, raw in enumerate(batch):
            residues = raw.metadata.get("chemfiles_residues", [])
            sites = raw.atoms
            centres, donors, centre_ids, donor_ids = [], [], [], []
            for ordinal, residue in enumerate(residues, 1):
                if residue_name is not None and residue["name"] != residue_name:
                    continue
                members = residue["atoms"]
                found = [i for i in members if sites[i].symbol == atom_name]
                if not found:
                    continue
                if len(found) != 1:
                    raise ValueError(
                        "selected atom name must identify one centre per residue"
                    )
                centres.append(found[0])
                centre_ids.append(ordinal)
                selected = [i for i in members if sites[i].symbol in hydrogens]
                missing = hydrogens.difference(sites[i].symbol for i in selected)
                if missing:
                    names = ", ".join(sorted(missing))
                    raise ValueError(
                        f"residue {ordinal} lacks requested hydrogen sites: {names}"
                    )
                donors.extend(selected)
                donor_ids.extend([ordinal] * len(selected))
            if not centres or (hydrogens and not donors):
                raise ValueError(
                    "no matching sites; check the atom names and the labelled topology"
                )
            if not all(raw.metadata.get("pbc", (False, False, False))):
                raise ValueError(
                    "trajectory analysis needs a three-dimensional periodic cell"
                )
            lattice = Cell(np.asarray(raw.metadata["lattice_vectors"]).T)
            restricted, rotation = lattice.standard_form("lower")
            if np.any(np.diag(restricted) < 0):
                restricted, rotation = -restricted, -rotation
            xyz = wrap_positions(raw.coords_array(), lattice) @ rotation.T
            box, low = _dump_box(restricted, (0.0, 0.0, 0.0))
            frame = Frame.from_arrays(
                xyz[centres],
                box,
                box_low=low,
                mol_ids=centre_ids,
                hydrogen_positions=xyz[donors] if donors else None,
                hydrogen_mol_ids=donor_ids if donors else None,
                cutoff=cutoff,
                bonded=bonded,
            )
            frame.frame = first + offset * stride + 1
            frame.cloud.currentFrame = frame.frame
            frame._cell_rotation = rotation
            yield frame


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
    if "atom_name" in kwargs:
        if frame < 1:
            raise ValueError("frame is one-indexed")
        result = next(
            iter_trajectory(path, start=frame - 1, stop=frame, **kwargs), None
        )
        if result is None:
            raise IndexError("frame index outside trajectory")
        return result
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
