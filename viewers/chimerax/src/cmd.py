"""``mdcif`` commands.

    mdcif open gb1.cif                          open the models (#N.1, #N.2, ...) and read the modes
    mdcif modes #1                              list models, cluster populations and modes
    mdcif morph #1 frames 20                    morph through the models (ChimeraX ``morph``)
    mdcif animate #1.1 mode 1 amplitude 3       copy of model 1 oscillating along its mode 1
    mdcif animate #1.2 mode 1 scope global      ... along global PC1 (the clustering axis)
    mdcif porcupine #1.1 mode 1 scale 3         arrows on the mode atoms of model 1

``amplitude`` / ``scale`` are in standard deviations of the mode (k in
k * sigma * v, sigma in angstroms: sqrt(lambda) for PCA, with the elements'
standard masses taken out for quasi-harmonic modes). Text mmCIF only.
"""

import numpy as np
from chimerax.atomic import AtomicStructureArg, AtomicStructuresArg
from chimerax.core.commands import (
    BoolArg,
    CmdDesc,
    ColorArg,
    EnumOf,
    FloatArg,
    IntArg,
    OpenFileNameArg,
    run,
)
from chimerax.core.errors import UserError

from . import __version__
from .mdcif_data import MdcifFile, arrow_mesh, arrows, mode_frames

ScopeArg = EnumOf(["cluster", "global"])


def _structures(models):
    from chimerax.atomic import AtomicStructure

    found = []
    for m in models:
        for sub in m.all_models():
            if isinstance(sub, AtomicStructure) and sub not in found:
                found.append(sub)
    return sorted(found, key=lambda s: s.id)


def _data(structure):
    data = getattr(structure, "mdcif_data", None)
    if data is None:
        raise UserError(f"#{structure.id_string} was not opened with 'mdcif open'")
    return data, structure.mdcif_model


def _keys(structure):
    atoms = structure.atoms
    res = atoms.residues
    return [(str(c), str(n), str(a)) for c, n, a in zip(res.chain_ids, res.numbers, atoms.names,
                                                        strict=True)]


def _mode_set(structure, scope):
    data, model = _data(structure)
    try:
        return data.mode_set(model, scope == "global"), model
    except ValueError as exc:
        raise UserError(str(exc)) from None


def mdcif(session):
    session.logger.info(
        f"ChimeraX-MDCIF {__version__}. Subcommands: mdcif open | modes | morph | animate | "
        "porcupine.")


mdcif_desc = CmdDesc(synopsis="Report mdcif bundle status")


def mdcif_open(session, path):
    """Open an mdcif file with the mmCIF reader and attach its modes to each model."""
    try:
        data = MdcifFile(path)
    except ValueError as exc:
        raise UserError(str(exc)) from None
    models = run(session, f'open "{path}" format mmcif', log=False) or []
    structures = _structures(models)
    if len(structures) != len(data.models):
        session.logger.warning(f"{len(structures)} structures opened for {len(data.models)} "
                               "models in the file; matching them in order")
    for s, m in zip(structures, data.models, strict=False):
        s.mdcif_data, s.mdcif_model = data, m
    session.logger.info(data.summary())
    return structures


mdcif_open_desc = CmdDesc(required=[("path", OpenFileNameArg)],
                          synopsis="Open an mdcif file and read its essential modes")


def mdcif_modes(session, structures=None):
    from chimerax.atomic import AtomicStructure

    if structures is None:
        structures = [m for m in session.models if isinstance(m, AtomicStructure)]
    seen = set()
    for s in structures:
        data = getattr(s, "mdcif_data", None)
        if data is not None and id(data) not in seen:
            seen.add(id(data))
            session.logger.info(data.summary())
    if not seen:
        session.logger.warning("no structures opened with 'mdcif open'")


mdcif_modes_desc = CmdDesc(optional=[("structures", AtomicStructuresArg)],
                           synopsis="List the models and essential modes of mdcif structures")


def mdcif_morph(session, structures, frames=20, wrap=True):
    targets = [s for s in structures if getattr(s, "mdcif_data", None) is not None]
    if len(targets) < 2:
        raise UserError("need at least two models opened with 'mdcif open'")
    targets.sort(key=lambda s: s.mdcif_model)
    spec = " ".join(f"#{s.id_string}" for s in targets)
    run(session, f"morph {spec} frames {frames} wrap {str(bool(wrap)).lower()} play true")


mdcif_morph_desc = CmdDesc(required=[("structures", AtomicStructuresArg)],
                           keyword=[("frames", IntArg), ("wrap", BoolArg)],
                           synopsis="Morph through the cluster representative models")


def mdcif_animate(session, structure, mode=1, amplitude=3.0, frames=30, scope="cluster"):
    ms, model = _mode_set(structure, scope)
    try:
        xyz = mode_frames(structure.atoms.coords, _keys(structure), ms, mode, amplitude, frames)
    except ValueError as exc:
        raise UserError(str(exc)) from None
    copy = structure.copy(f"{structure.name} {scope} mode {mode}")
    copy.add_coordsets(xyz, replace=True)
    session.models.add([copy])
    structure.display = False
    spec = f"#{copy.id_string}"
    try:
        run(session, f"coordset {spec} 1,{frames} loop 1000", log=False)
    except UserError:
        run(session, f"coordset {spec} 1,{frames}", log=False)
    session.logger.info(f"{spec}: {ms.scope} mode {mode} of model {model}, "
                        f"{100 * ms.variance_fraction[mode - 1]:.0f}% variance, "
                        f"+/-{amplitude:g} sd ('coordset {spec} stop' to stop)")
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

    ms, model = _mode_set(structure, scope)
    try:
        starts, ends = arrows(structure.atoms.coords, _keys(structure), ms, mode, scale,
                              min_length)
    except ValueError as exc:
        raise UserError(str(exc)) from None
    va, na, ta = arrow_mesh(starts, ends, radius=radius, head_radius=2.7 * radius)
    surf = Surface(f"{scope} mode {mode} porcupine", session)
    surf.set_geometry(va, na, ta)
    surf.color = (color.uint8x4() if color is not None
                  else np.array([255, 140, 0, 255], dtype=np.uint8))
    structure.add([surf])
    session.logger.info(f"#{surf.id_string}: {len(starts)} arrows, {ms.scope} mode {mode} of "
                        f"model {model}, {scale:g} sd")
    return surf


mdcif_porcupine_desc = CmdDesc(
    required=[("structure", AtomicStructureArg)],
    keyword=[("mode", IntArg), ("scale", FloatArg), ("color", ColorArg), ("scope", ScopeArg),
             ("min_length", FloatArg), ("radius", FloatArg)],
    synopsis="Draw porcupine arrows for one mode of a model",
)
