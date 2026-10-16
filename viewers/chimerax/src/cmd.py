"""``mdcif`` commands.

    open gb1.cif                                any mmCIF open (also File > Open, drag and drop):
                                                models #N.1, #N.2, ... get the file's mdcif data
    mdcif open gb1.cif                          the same, as mmCIF whatever the suffix
    mdcif modes #1                              list models, cluster populations and modes
    mdcif morph #1 frames 20                    morph through the models (ChimeraX ``morph``)
    mdcif animate #1.1 mode 1 amplitude 3       copy of model 1 oscillating along its mode 1
    mdcif animate #1.2 mode 1 scope global      ... along global PC1 (the clustering axis)
    mdcif porcupine #1.1 mode 1 scale 3         arrows on the mode atoms of model 1
    coordset stop #2                            stop an animation (ChimeraX ``coordset``)

``amplitude`` / ``scale`` are in standard deviations of the mode (k in
k * sigma * v, sigma in angstroms: sqrt(lambda) for PCA, with the elements'
standard masses taken out for quasi-harmonic modes). Text mmCIF only.

Structures get their data from ``attach`` (on open, on session restore, or the first
time a command needs it). An animated copy keeps its model's mode data, so ``mdcif
porcupine`` works on it too; both commands start from the deposited coordinates (the
copy's first coordinate set).
"""

import os

import numpy as np
from chimerax.atomic import AtomicStructureArg, AtomicStructuresArg
from chimerax.core.commands import (
    BoolArg,
    CmdDesc,
    ColorArg,
    EnumOf,
    FileNameArg,
    FloatArg,
    IntArg,
    OpenFileNameArg,
    run,
)
from chimerax.core.errors import UserError

from . import __version__, attach
from .mdcif_data import MdcifFile, arrow_mesh, arrows, mode_frames

ScopeArg = EnumOf(["cluster", "global"])
MIN_ANIMATE_FRAMES = 3                  # sin() is 0 at every frame of a 1- or 2-frame cycle


def _structures(models):
    from chimerax.atomic import AtomicStructure

    found = []
    for m in models:
        for sub in m.all_models():
            if isinstance(sub, AtomicStructure) and sub not in found:
                found.append(sub)
    return sorted(found, key=lambda s: s.id)


def _attached(session, structure):
    """The structure's mdcif data, attached now if the "add models" handler missed it."""
    data = getattr(structure, "mdcif_data", None)
    return data if data is not None else attach.lazy_attach(session, structure)


def _data(session, structure):
    data = _attached(session, structure)
    if data is None:
        path = attach.source_path(structure)
        if path is not None and not os.path.exists(path):
            raise UserError(f"#{structure.id_string}: its file {path} is gone (moved since the "
                            "session was saved?); open the mdcif file again")
        raise UserError(f"#{structure.id_string} has no mdcif data (open an mdcif .cif file "
                        "with 'open' or 'mdcif open')")
    return data, structure.mdcif_model


def _keys(structure):
    atoms = structure.atoms
    res = atoms.residues
    return [(str(c), str(n), str(a)) for c, n, a in zip(res.chain_ids, res.numbers, atoms.names,
                                                        strict=True)]


def _mode_set(session, structure, scope):
    data, model = _data(session, structure)
    try:
        return data.mode_set(model, scope == "global"), model
    except ValueError as exc:
        raise UserError(str(exc)) from None


def _rest_coords(structure):
    """Coordinates the modes are drawn from: the current ones, or an animated copy's first
    coordinate set (the deposited model; read by coord index, as ChimeraX keeps the slots of
    deleted atoms), so arrows and new animations don't start mid-swing."""
    atoms = structure.atoms
    if getattr(structure, "mdcif_copy", False) and 1 in structure.coordset_ids:
        return structure.coordset(1).xyzs[atoms.coord_indices]
    return atoms.coords


def _set_frames(structure, xyz):
    """Replace the coordinate sets with ``xyz`` (frames, atoms, 3; atom order).

    After atoms are deleted ChimeraX keeps their slots in the coordinate sets (coord_indices
    then has gaps), and ``add_coordsets`` would misplace the frames: fill by coord index.
    """
    ci = structure.atoms.coord_indices
    if np.array_equal(ci, np.arange(len(ci))):
        structure.add_coordsets(xyz, replace=True)
        return
    full = np.zeros((len(xyz), int(ci.max()) + 1, 3))
    full[:, ci] = xyz
    structure.remove_coordsets()
    for k, frame in enumerate(full, start=1):       # the first one becomes the active set
        structure.add_coordset(k, frame)


def mdcif(session):
    auto = "on" if attach.is_registered(session) else "off (restart ChimeraX after installing)"
    session.logger.info(
        f"ChimeraX-MDCIF {__version__}. Subcommands: mdcif open | modes | morph | animate | "
        f"porcupine. Attaching mdcif data to opened .cif files: {auto}.")


mdcif_desc = CmdDesc(synopsis="Report mdcif bundle status")


def mdcif_open(session, path):
    """Open a text mmCIF file with the mmCIF reader and attach its mdcif data to each model.

    Like plain 'open', but it checks the file first, reads it as mmCIF whatever its suffix,
    and attaches the data even to a file without _mdcif_ categories (no modes then).
    """
    try:
        data = MdcifFile(path)
    except (ValueError, OSError) as exc:            # OSError: a directory, unreadable file
        raise UserError(str(exc)) from None
    key = attach.path_key(path)
    attach.pending[key] = data                      # the "add models" handler attaches it
    try:
        models = run(session, f"open {FileNameArg.unparse(path)} format mmcif", log=False) or []
    finally:
        unused = attach.pending.pop(key, None)
    structures = _structures(models)
    if unused is not None and structures:           # handler not run (no .cif suffix, no init)
        attach.attach(session, structures, path, data=unused)
    return structures


mdcif_open_desc = CmdDesc(required=[("path", OpenFileNameArg)],
                          synopsis="Open an mdcif file and read its essential modes")


def mdcif_modes(session, structures=None):
    from chimerax.atomic import AtomicStructure

    if structures is None:
        structures = [m for m in session.models if isinstance(m, AtomicStructure)]
    seen = set()
    for s in structures:
        data = _attached(session, s)
        if data is not None and id(data) not in seen:
            seen.add(id(data))
            session.logger.info(data.summary())
    if not seen:
        session.logger.warning("no structures with mdcif data")


mdcif_modes_desc = CmdDesc(optional=[("structures", AtomicStructuresArg)],
                           synopsis="List the models and essential modes of mdcif structures")


def mdcif_morph(session, structures, frames=20, wrap=True):
    if frames < 1:
        raise UserError("frames must be at least 1")
    targets = sorted((s for s in structures if not getattr(s, "mdcif_copy", False)
                      and _attached(session, s) is not None), key=lambda s: s.id)
    if len(targets) < 2:
        raise UserError("need at least two models with mdcif data (animate copies don't count)")
    file_order = {}                                 # files in model-id order, then model number
    for s in targets:
        file_order.setdefault(id(s.mdcif_data), len(file_order))
    targets.sort(key=lambda s: (file_order[id(s.mdcif_data)], s.mdcif_model))
    spec = " ".join(f"#{s.id_string}" for s in targets)
    run(session, f"morph {spec} frames {frames} wrap {str(bool(wrap)).lower()} play true")


mdcif_morph_desc = CmdDesc(required=[("structures", AtomicStructuresArg)],
                           keyword=[("frames", IntArg), ("wrap", BoolArg)],
                           synopsis="Morph through the cluster representative models")


def mdcif_animate(session, structure, mode=1, amplitude=3.0, frames=30, scope="cluster"):
    if frames < MIN_ANIMATE_FRAMES:
        raise UserError(f"frames must be at least {MIN_ANIMATE_FRAMES}")
    ms, model = _mode_set(session, structure, scope)
    data = structure.mdcif_data
    try:
        xyz = mode_frames(_rest_coords(structure), _keys(structure), ms, mode, amplitude, frames)
    except ValueError as exc:
        raise UserError(str(exc)) from None
    base = os.path.basename(data.path)
    copy = structure.copy(f"{base} model {model} {scope} mode {mode}")
    _set_frames(copy, xyz)
    copy.position = structure.scene_position        # the copy is top level, the source may not be
    copy.mdcif_data, copy.mdcif_model, copy.mdcif_copy = data, model, True
    copy.mdcif_path = structure.mdcif_path          # data set before adding: the handler skips it
    session.models.add([copy])
    structure.display = False
    spec = f"#{copy.id_string}"
    run(session, f"coordset {spec} 1,{frames} loop 1000", log=False)
    session.logger.info(f"{spec}: {ms.scope} mode {mode} of model {model}, "
                        f"{100 * ms.variance_fraction[mode - 1]:.0f}% variance, "
                        f"+/-{abs(amplitude):g} sd ('coordset stop {spec}' to stop)")
    return copy


mdcif_animate_desc = CmdDesc(
    required=[("structure", AtomicStructureArg)],
    keyword=[("mode", IntArg), ("amplitude", FloatArg), ("frames", IntArg),
             ("scope", ScopeArg)],
    synopsis="Animate a model along one of its essential modes",
)


def mdcif_porcupine(session, structure, mode=1, scale=3.0, color=None, scope="cluster",
                    min_length=0.5, radius=0.15):
    from chimerax.core.models import Surface

    if radius <= 0:
        raise UserError("radius must be positive")
    if min_length < 0:
        raise UserError("minLength must not be negative")
    ms, model = _mode_set(session, structure, scope)
    try:
        starts, ends = arrows(_rest_coords(structure), _keys(structure), ms, mode, scale,
                              min_length)
    except ValueError as exc:
        raise UserError(str(exc)) from None
    va, na, ta = arrow_mesh(starts, ends, radius=radius, head_radius=2.7 * radius)
    if len(ta) == 0:
        raise UserError(f"no arrows to draw for {ms.scope} mode {mode} at {scale:g} sd "
                        f"with minLength {min_length:g} A")
    surf = Surface(f"{scope} mode {mode} porcupine", session)
    surf.set_geometry(va, na, ta)
    surf.color = (color.uint8x4() if color is not None
                  else np.array([255, 140, 0, 255], dtype=np.uint8))
    structure.add([surf])
    session.logger.info(f"#{surf.id_string}: {len(starts)} arrows, {ms.scope} mode {mode} of "
                        f"model {model}, {scale:g} sd")
    if not structure.visible:
        session.logger.warning(f"#{structure.id_string} is hidden, and so are its arrows "
                               f"('show #{structure.id_string} models' to show them)")
    return surf


mdcif_porcupine_desc = CmdDesc(
    required=[("structure", AtomicStructureArg)],
    keyword=[("mode", IntArg), ("scale", FloatArg), ("color", ColorArg), ("scope", ScopeArg),
             ("min_length", FloatArg), ("radius", FloatArg)],
    synopsis="Draw porcupine arrows for one mode of a model",
)
