"""Load and align MD trajectories (MDAnalysis).

Any topology/trajectory pair MDAnalysis can read (GROMACS, AMBER, CHARMM/NAMD,
OpenMM, ...) is supported. Several trajectory files (e.g. replicas) are read
back to back. NetCDF, DCD, XTC and TRR files are recognised by their contents
whatever their name (Amber replica-exchange files such as ``meld.nc.001``,
GROMACS backups such as ``#traj.xtc.1#``); other files by their extension, or
by the extension before a numeric suffix (``meld.mdcrd.001``). Coordinates are
returned as arrays of shape ``(n_frames, n_atoms, 3)`` in angstroms.
"""

from pathlib import Path

import numpy as np


def _sniff(path):
    """MDAnalysis format from the first bytes of a binary trajectory file, else None."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(8)
    except OSError:
        return None
    if head[:4] in (b"CDF\x01", b"CDF\x02"):         # NetCDF 3 (Amber, MELD)
        return "NC"
    if head[:4] in (b"T\0\0\0", b"\0\0\0T") and head[4:8] == b"CORD":   # record of 84 bytes
        return "DCD"
    if len(head) >= 4:                               # GROMACS magic numbers, big-endian
        return {1995: "XTC", 1993: "TRR"}.get(int.from_bytes(head[:4], "big"))
    return None


def _has_reader(path, format=None):
    from MDAnalysis.coordinates.core import get_reader_for

    try:
        get_reader_for(str(path), format=format)
    except ValueError:
        return False
    return True


def guess_format(path):
    """MDAnalysis format of a trajectory file, or None to leave it to MDAnalysis.

    The first bytes decide for NetCDF, DCD, XTC and TRR. Otherwise, for names
    with a numeric suffix (``meld.mdcrd.001``), the extension before it, if
    MDAnalysis has a reader for it.
    """
    fmt = _sniff(path)
    if fmt:
        return fmt
    suffixes = Path(path).suffixes
    if not suffixes or not suffixes[-1][1:].isdigit():
        return None
    while suffixes and suffixes[-1][1:].isdigit():
        suffixes.pop()
    if suffixes and _has_reader(path, suffixes[-1][1:].upper()):
        return suffixes[-1][1:].upper()
    return None


def load_universe(topology, trajectory, format=None):
    """Return an ``MDAnalysis.Universe`` for a topology and one or more trajectory files."""
    import MDAnalysis as mda

    files = ([str(t) for t in trajectory] if isinstance(trajectory, (list, tuple))
             else [str(trajectory)])
    formats = [format or guess_format(f) for f in files]
    for f, fmt in zip(files, formats, strict=True):
        if fmt is None and not _has_reader(f):
            raise ValueError(f"can't tell the trajectory format of {f}; give it with format= / "
                             "--format (e.g. NC, DCD, XTC, MDCRD)")
    if len(files) == 1:
        kwargs = {"format": formats[0]} if formats[0] else {}
        return mda.Universe(str(topology), files[0], **kwargs)
    if any(formats):
        return mda.Universe(str(topology), [(f, fmt) if fmt else f
                                            for f, fmt in zip(files, formats, strict=True)])
    return mda.Universe(str(topology), files)


def read_frames(universe, selection, start=0, stop=None, step=1):
    """Coordinates of ``selection`` for frames ``start:stop:step``.

    Returns ``(coords, frames, times)``: ``(n, n_atoms, 3)`` coordinates, the
    0-based index of each frame in the (chained) trajectory and its time in ps
    (NaN if the trajectory stores the same time for every frame, e.g. 0).
    """
    atoms = universe.select_atoms(selection)
    if len(atoms) == 0:
        raise ValueError(f"selection {selection!r} matched no atoms")
    frames = np.arange(universe.trajectory.n_frames)[start:stop:step]
    if len(frames) == 0:
        raise ValueError(f"no frames in range start={start} stop={stop} step={step}")
    out = np.empty((len(frames), len(atoms), 3), dtype=np.float64)
    times = np.empty(len(frames))
    for i, ts in enumerate(universe.trajectory[start:stop:step]):
        out[i] = atoms.positions
        times[i] = ts.time
    if len(times) > 1 and np.ptp(times) == 0:
        times[:] = np.nan
    return out, frames, times


def superpose(coords, reference):
    """Superpose every frame of ``coords`` ``(n_frames, n_atoms, 3)`` onto ``reference``.

    Returns ``(aligned, rotations, centres, reference_centre)`` with
    ``aligned[f] = (coords[f] - centres[f]) @ rotations[f] + reference_centre``.
    The same transform can be applied to other atoms of a frame with
    :func:`apply_transform`. Rotations from MDAnalysis (QCP, ``align.rotation_matrix``).
    """
    from MDAnalysis.analysis.align import rotation_matrix

    coords = np.asarray(coords, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    ref_c = reference.mean(axis=0)
    centres = coords.mean(axis=1)
    mobile = coords - centres[:, None, :]
    # rotation_matrix(a, b) gives R with b = R a (columns); our frames are rows: a @ R.T
    rot = np.array([rotation_matrix(m, reference - ref_c)[0].T for m in mobile]).reshape(-1, 3, 3)
    aligned = mobile @ rot + ref_c
    return aligned, rot, centres, ref_c


def apply_transform(coords, rotation, centre, reference_centre):
    """Apply one frame's transform from :func:`superpose` to other atoms of that frame."""
    return (coords - centre) @ rotation + reference_centre


def rmsd_to(coords, reference):
    """RMSD of each frame of ``coords`` to ``reference`` after optimal superposition."""
    from MDAnalysis.analysis.rms import rmsd

    reference = np.asarray(reference, dtype=np.float64)
    return np.array([rmsd(x, reference, center=True, superposition=True)
                     for x in np.asarray(coords, dtype=np.float64)])
