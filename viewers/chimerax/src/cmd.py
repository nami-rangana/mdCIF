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
    mdcif plot distance #2/A:10@CA #2/A:40@CA   the distance at every frame of animation #2
    mdcif plot angle #2/A:10@CA #2/A:11@CA #2/A:12@CA   ... and an angle, in a panel under it
    mdcif plot dihedral sel saveData d.csv      4 picked atoms (in picking order); a window
                                                with the GUI, saveData / saveImage headless too

``amplitude`` / ``scale`` are in standard deviations of the mode (k in
k * sigma * v, sigma in angstroms: sqrt(lambda) for PCA, with the elements'
standard masses taken out for quasi-harmonic modes). Text mmCIF only.

Structures get their data from ``attach`` (on open, on session restore, or the first
time a command needs it). An animated copy keeps its model's mode data, so ``mdcif
porcupine`` works on it too; both commands start from the deposited coordinates (the
copy's first coordinate set). ``mdcif plot`` measures on an animate copy (made for it from
a model, or the one given); one plot per copy holds every kind, a panel each (``new true``
starts another). Plots are in ``plot.py``, the Analyze window and its Geometry tab in
``analyze_tool.py`` / ``geometry_tab.py``.

The windows of the MDCIF menu (``menu.py``) work through the helpers here, which have no Qt:
``porcupine_geometry`` (the arrow mesh of 'mdcif porcupine') and ``porcupine_command`` (the
command to log), ``retune_animation`` (new frames for an animate copy, in place),
``retire_copy`` (close or put away a copy no window shows any more) and ``loop_command`` (the
'coordset' loop that plays a copy).
"""

import os

import numpy as np
from chimerax.atomic import AtomicStructureArg, AtomicStructuresArg, AtomsArg
from chimerax.core.commands import (
    BoolArg,
    CmdDesc,
    ColorArg,
    EnumOf,
    FileNameArg,
    FloatArg,
    IntArg,
    OpenFileNameArg,
    SaveFileNameArg,
    log_equivalent_command,
    run,
)
from chimerax.core.errors import UserError

from . import __version__, attach
from .mdcif_data import MEASURES, MdcifFile, arrow_mesh, arrows, measure, mode_frames
from .plot import A_KIND, IMAGE_SUFFIXES, KIND_OF, PLOT_MADE, UNITS, Plot, atom_label, plots
from .plot import save_image as _save_image
from .plot import write_csv as _write_csv

ScopeArg = EnumOf(["cluster", "global"])
KindArg = EnumOf(list(MEASURES))
MIN_ANIMATE_FRAMES = 3                  # sin() is 0 at every frame of a 1- or 2-frame cycle
ANIMATE_DEFAULTS = (1, "cluster", 3.0, 30)      # mode, scope, amplitude, frames
PORCUPINE_COLOR = (255, 140, 0, 255)            # 'mdcif porcupine' arrows without 'color'
ANIMATE_MADE = "_mdcif_animate_made"    # attribute of copies the Animate window made (not saved)


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


def _atom_keys(atoms):
    res = atoms.residues
    return [(str(c), str(n), str(a)) for c, n, a in zip(res.chain_ids, res.numbers, atoms.names,
                                                        strict=True)]


def _keys(structure):
    return _atom_keys(structure.atoms)


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
    where = "" if session.ui.is_gui else " (when ChimeraX has its window)"
    session.logger.info(
        f"ChimeraX-MDCIF {__version__}. Subcommands: mdcif open | modes | morph | animate | "
        f"porcupine | plot. Attaching mdcif data to opened .cif files: {auto}. Windows: menu "
        f"MDCIF > Porcupine, Animate, Analyze{where}.")


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
    copy = structure.copy(_copy_name(data, model, scope, mode))
    _set_frames(copy, xyz)
    copy.position = structure.scene_position        # the copy is top level, the source may not be
    copy.mdcif_data, copy.mdcif_model, copy.mdcif_copy = data, model, True
    copy.mdcif_path = structure.mdcif_path          # data set before adding: the handler skips it
    copy.mdcif_mode, copy.mdcif_scope, copy.mdcif_amplitude = mode, scope, float(amplitude)
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


def _copy_name(data, model, scope, mode):
    return f"{os.path.basename(data.path)} model {model} {scope} mode {mode}"


def animate_command(structure, mode, amplitude, frames, scope):
    """The 'mdcif animate' command for these settings, as text (for the Log)."""
    return (f"mdcif animate #{structure.id_string} mode {mode} amplitude {amplitude:g} "
            f"frames {frames} scope {scope}")


def loop_command(structure, pause_frames=1):
    """The 'coordset' command that plays an animate copy in a loop from frame 1, showing each
    frame for ``pause_frames`` graphics frames (1/N of full speed)."""
    pause = f" pauseFrames {pause_frames}" if pause_frames > 1 else ""
    return f"coordset #{structure.id_string} 1,{structure.num_coordsets} loop 1000{pause}"


def player(session, structure):
    """The 'coordset' player moving ``structure``, or None. ChimeraX 1.12 keeps the running
    ones in ``session._coord_set_players`` (std_commands/coordset.py); a stopped one leaves."""
    return next((p for p in getattr(session, "_coord_set_players", ())
                 if p.structure is structure), None)


def retune_animation(session, copy, mode=None, scope=None, amplitude=None, frames=None):
    """Give an animate copy the frames of other settings, in place (same model id); None keeps
    a setting. What the Animate window does when a setting changes.

    The frames start again from the copy's first coordinate set (the deposited model); the
    copy's ``mdcif_mode`` / ``mdcif_scope`` / ``mdcif_amplitude`` and its name follow. A playing
    copy keeps playing: with as many frames as before its loop just shows the new ones,
    otherwise the loop starts again at frame 1 at the same speed. A stopped copy stays on its
    frame. Returns the copy, or None (and changes nothing) when a plot measures on it, whose
    values would no longer match the frames: the caller makes a new copy then.
    """
    if not _alive(copy) or not getattr(copy, "mdcif_copy", False):
        raise UserError("not an 'mdcif animate' copy")
    if any(p.alive and p.structure is copy for p in plots(session)):
        return None
    wanted = [getattr(copy, name, None) if v is None else v for v, name in zip(
        (mode, scope, amplitude), ("mdcif_mode", "mdcif_scope", "mdcif_amplitude"), strict=True)]
    mode, scope, amplitude = (d if v is None else v for v, d in zip(
        wanted, ANIMATE_DEFAULTS[:3], strict=True))     # copies from before 0.3.0: defaults
    old_frames = copy.num_coordsets
    frames = old_frames if frames is None else int(frames)
    if frames < MIN_ANIMATE_FRAMES:
        raise UserError(f"frames must be at least {MIN_ANIMATE_FRAMES}")
    ms, model = _mode_set(session, copy, scope)
    try:
        xyz = mode_frames(_rest_coords(copy), _keys(copy), ms, mode, amplitude, frames)
    except ValueError as exc:
        raise UserError(str(exc)) from None
    playing = player(session, copy)
    shown = getattr(copy, "active_coordset_id", 1)
    restart = playing is not None and frames != old_frames
    if restart:
        pause = playing.pause_frames
        run(session, f"coordset stop #{copy.id_string}", log=False)
    _set_frames(copy, xyz)
    copy.mdcif_mode, copy.mdcif_scope, copy.mdcif_amplitude = mode, scope, float(amplitude)
    copy.name = _copy_name(copy.mdcif_data, model, scope, mode)
    if restart:
        run(session, loop_command(copy, pause), log=False)
    elif playing is None:                           # stay on the frame shown
        copy.active_coordset_id = min(max(shown, 1), frames)
    return copy


def retire_copy(session, copy, own):
    """Put away an animate copy that a plot or the Animate window has moved off: closed if
    ``own`` (they made it) and no plot measures on it, else stopped and hidden."""
    if not _alive(copy):
        return
    if own and not any(p.structure is copy for p in plots(session)):
        session.models.close([copy])
    else:
        run(session, f"coordset stop #{copy.id_string}", log=False)
        copy.display = False


def porcupine_geometry(session, structure, mode=1, scale=3.0, scope="cluster", min_length=0.5,
                       radius=0.15):
    """The arrows 'mdcif porcupine' draws, as a mesh, without drawing them:
    ``(va, na, ta, n_arrows, n_short, ms, model)`` = vertices, normals and triangles, the
    number of arrows and of those left out as shorter than ``min_length``, the mode set and the
    model number. UserError for what the command refuses (no arrows left, a bad radius, ...)."""
    if radius <= 0:
        raise UserError("radius must be positive")
    if min_length < 0:
        raise UserError("minLength must not be negative")
    ms, model = _mode_set(session, structure, scope)
    coords, keys = _rest_coords(structure), _keys(structure)
    try:
        starts, ends = arrows(coords, keys, ms, mode, scale, min_length)
        n_all = len(arrows(coords, keys, ms, mode, scale, 0.0)[0])
    except ValueError as exc:
        raise UserError(str(exc)) from None
    va, na, ta = arrow_mesh(starts, ends, radius=radius, head_radius=2.7 * radius)
    if len(ta) == 0:
        raise UserError(f"no arrows to draw for {ms.scope} mode {mode} at {scale:g} sd "
                        f"with minLength {min_length:g} A")
    return va, na, ta, len(starts), n_all - len(starts), ms, model


def porcupine_command(structure, mode=1, scale=3.0, color=None, scope="cluster", min_length=0.5,
                      radius=0.15):
    """The 'mdcif porcupine' command for these settings, as text (for the Log); ``color``: RGB
    or RGBA, 0-255 (None: the default orange)."""
    r, g, b = (int(c) for c in (PORCUPINE_COLOR if color is None else color)[:3])
    return (f"mdcif porcupine #{structure.id_string} mode {mode} scale {scale:g} color "
            f"#{r:02x}{g:02x}{b:02x} radius {radius:g} minLength {min_length:g} scope {scope}")


def arrow_surface(session, scope, mode):
    """An empty Surface for one mode's arrows that sessions save whole: a plain Surface saves
    only its place in the model tree (``SESSION_SAVE_DRAWING`` is False), and comes back from
    a session grey and empty."""
    from chimerax.core.models import Surface

    surf = Surface(f"{scope} mode {mode} porcupine", session)
    surf.SESSION_SAVE_DRAWING = True                # vertices, triangles, colour
    return surf


def mdcif_porcupine(session, structure, mode=1, scale=3.0, color=None, scope="cluster",
                    min_length=0.5, radius=0.15):
    va, na, ta, n_arrows, _, ms, model = porcupine_geometry(session, structure, mode, scale,
                                                            scope, min_length, radius)
    surf = arrow_surface(session, scope, mode)
    surf.set_geometry(va, na, ta)
    surf.color = (color.uint8x4() if color is not None
                  else np.array(PORCUPINE_COLOR, dtype=np.uint8))
    structure.add([surf])
    session.logger.info(f"#{surf.id_string}: {n_arrows} arrows, {ms.scope} mode {mode} of "
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


# --------------------------------------------------------------------------- plots

def _alive(structure):
    return structure is not None and not getattr(structure, "deleted", False)


def _matches(copy, mode=None, scope=None, amplitude=None, frames=None):
    """Whether an animate copy was made with these settings (None = any). Copies made before
    0.3.0 don't record theirs and never match."""
    have = (getattr(copy, "mdcif_mode", None), getattr(copy, "mdcif_scope", None),
            getattr(copy, "mdcif_amplitude", None), copy.num_coordsets)
    return have[0] is not None and all(w is None or w == h for w, h in zip(
        (mode, scope, amplitude, frames), have, strict=True))


def _animate(session, structure, mode=None, scope=None, amplitude=None, frames=None,
             mark=PLOT_MADE):
    """'mdcif animate' for a plot (``mark=ANIMATE_MADE``: for the Animate window), logged as
    that command; None takes its default. The copy is marked as made by them, so they may
    close it when they animate again."""
    mode, scope, amplitude, frames = (d if v is None else v for v, d in zip(
        (mode, scope, amplitude, frames), ANIMATE_DEFAULTS, strict=True))
    log_equivalent_command(session, animate_command(structure, mode, amplitude, frames, scope))
    copy = mdcif_animate(session, structure, mode, amplitude, frames, scope)
    setattr(copy, mark, True)
    return copy


def _values(structure, keys, kind):
    """``kind`` of the atoms ``keys`` at every coordinate set of ``structure``, in id order
    (read by coord index, as ChimeraX keeps the slots of deleted atoms)."""
    index = {}
    for i, k in enumerate(_keys(structure)):
        index.setdefault(k, i)
    missing = [k for k in keys if k not in index]
    if missing:
        raise UserError(f"#{structure.id_string} has no atom {atom_label(missing[0])}")
    ci = structure.atoms.coord_indices[[index[k] for k in keys]]
    xyz = np.array([structure.coordset(int(c)).xyzs[ci] for c in sorted(structure.coordset_ids)])
    return measure(xyz, kind)


def _of_file(structure, path):
    """Whether ``structure`` has the mdcif data of file ``path``. By path, not by data object:
    a restored session attaches the copies' data separately from the models'."""
    data = getattr(structure, "mdcif_data", None)
    return data is not None and attach.path_key(data.path) == attach.path_key(path)


def _model_structure(session, path, model):
    """The open structure of a file's model (not an animate copy), or None."""
    from chimerax.atomic import AtomicStructure

    for s in session.models.list(type=AtomicStructure):
        if (_of_file(s, path) and getattr(s, "mdcif_model", None) == model
                and not getattr(s, "mdcif_copy", False) and _alive(s)):
            return s
    return None


def _shown_copies(session, structure):
    """The shown animate copies of the model ``structure`` (same file and model), by id."""
    from chimerax.atomic import AtomicStructure

    path, model = structure.mdcif_data.path, structure.mdcif_model
    return sorted((s for s in session.models.list(type=AtomicStructure)
                   if getattr(s, "mdcif_copy", False) and _alive(s) and s.display
                   and _of_file(s, path) and getattr(s, "mdcif_model", None) == model),
                  key=lambda s: s.id)


def copy_for(session, structure, mode=None, scope=None, amplitude=None, frames=None,
             new=False):
    """The animation atoms of ``structure`` are measured on: ``(copy, source)``.

    Atoms on an animate copy made with the settings given (None = any) use that copy. Atoms on
    a model (unless ``new``) use the newest copy a plot animated from it with these settings,
    else a shown copy of that model with these settings (also one the user made with 'mdcif
    animate', or one restored from a session). Otherwise a new copy is animated ('mdcif
    animate', default mode 1, cluster, 3 sd, 30 frames); other settings on a copy animate the
    copy anew. ``source``: the model the copy animates, if known. Whether a plot may close the
    copy later depends on who made it (``Plot.own``), not on how it was found.
    """
    _data(session, structure)
    wanted = (mode, scope, amplitude, frames)
    live = [p for p in plots(session) if p.alive]
    if getattr(structure, "mdcif_copy", False):
        if _matches(structure, *wanted):
            return structure, next((p.source for p in reversed(live)
                                    if p.structure is structure and p.source is not None), None)
    elif not new:
        found = ([p.structure for p in live if p.source is structure
                  and _matches(p.structure, *wanted)]
                 or [s for s in _shown_copies(session, structure) if _matches(s, *wanted)])
        if found:
            return found[-1], structure
    return _animate(session, structure, *wanted), structure


def plot_for(session, structure, mode=None, scope=None, amplitude=None, frames=None,
             new=False):
    """The plot (found or new) that atoms of ``structure`` are measured in: the newest plot on
    the copy from ``copy_for``, whatever kinds it shows (a new kind gets a panel there),
    unless ``new`` (which on a model also means a new copy). A new plot joins
    ``plots(session)`` with its first series."""
    copy, source = copy_for(session, structure, mode, scope, amplitude, frames, new)
    if not new:
        found = [p for p in plots(session) if p.alive and p.structure is copy]
        if found:
            return found[-1]
    return Plot(copy, source)


def add_series(session, plot, keys):
    """Measure the atoms ``keys`` (2, 3 or 4: a distance, angle or dihedral) on the plot's
    copy and add them to the plot (registered in the session if new); logs the range. Atoms
    already in the plot, also in reverse order (the same measurement), are not added twice."""
    if not plot.alive:
        raise UserError(f"{plot.spec}, the animation of this plot, was closed")
    keys = [tuple(k) for k in keys]
    kind = KIND_OF.get(len(keys))
    if kind is None:
        raise UserError(f"a distance, angle or dihedral needs 2, 3 or 4 atoms, got {len(keys)}")
    twice = next((k for i, k in enumerate(keys) if k in keys[:i]), None)
    if twice is not None:                           # typed atoms can repeat, picked ones not
        raise UserError(f"{atom_label(twice)} is given twice: {A_KIND[kind]} needs "
                        f"{len(keys)} different atoms")
    same = next((s for s in plot.series if s.keys in (keys, keys[::-1])), None)
    if same is not None:
        session.logger.info(f"{plot.spec} {kind} {same.label} is already plotted")
        return same
    values = _values(plot.structure, keys, kind)
    try:
        series = plot.add(keys, values)
    except ValueError as exc:
        raise UserError(str(exc)) from None
    if plot not in plots(session):
        plots(session).append(plot)
    p = 2 if kind == "distance" else 1
    session.logger.info(f"{plot.spec} {kind} {series.label} along {plot.describe()}: "
                        f"min {values.min():.{p}f}, max {values.max():.{p}f}, "
                        f"mean {values.mean():.{p}f} {UNITS[kind]}")
    plot.changed()
    return series


def reanimate(session, plot, model=None, mode=None, scope=None, amplitude=None, frames=None):
    """Animate a plot again with other settings (None = keep) and measure its series anew.

    The new copy is made from the model's structure (for the same model: the plot's source,
    else the model, else the old copy). The old copy is then put away (``retire_copy``):
    closed if the plot made it and no other plot uses it, else stopped and hidden. Returns the
    new copy.
    """
    model = plot.model if model is None else model
    mode, scope, amplitude, frames = (old if v is None else v for v, old in zip(
        (mode, scope, amplitude, frames),
        (plot.mode, plot.scope, plot.amplitude, plot.n_frames), strict=True))
    source = plot.source if model == plot.model and _alive(plot.source) else None
    source = source or _model_structure(session, plot.path, model)
    if source is None and model == plot.model and plot.alive:
        source = plot.structure
    if source is None:
        raise UserError(f"no structure of model {model} of {os.path.basename(plot.path)} is "
                        "open")
    copy = _animate(session, source, mode, scope, amplitude, frames)
    try:
        values = [_values(copy, s.keys, s.kind) for s in plot.series]
    except UserError:
        session.models.close([copy])
        raise
    old, old_own = plot.structure, plot.own
    plot.set_animation(copy, source)
    for s, v in zip(plot.series, values, strict=True):
        s.values = v
    if old is not copy:
        retire_copy(session, old, old_own)
    plot.changed()
    return copy


def _check_target(option, path):
    """Refuse a file that can't be written before anything is animated or added."""
    full = os.path.abspath(os.path.expanduser(str(path)))
    folder = os.path.dirname(full)
    if os.path.isdir(full):
        raise UserError(f"{option} {path}: that is a folder, not a file")
    if not os.path.isdir(folder):
        raise UserError(f"{option} {path}: there is no folder {folder}")
    if not os.access(folder, os.W_OK) or (os.path.exists(full) and not os.access(full, os.W_OK)):
        raise UserError(f"{option} {path}: not allowed to write there")


def _save(session, plot, path, what):
    try:
        (_write_csv if what == "data" else _save_image)(plot, path)
    except (OSError, ValueError) as exc:
        raise UserError(f"could not save the plot {what}: {exc}") from None
    session.logger.info(f"{plot.spec} {' + '.join(plot.kinds)} plot: {what} saved to {path}")


def mdcif_plot(session, kind, atoms, mode=None, scope=None, amplitude=None, frames=None,
               new=False, save_data=None, save_image=None):
    """Plot a distance, angle or dihedral at every frame of a mode animation.

    ``atoms``: 2, 3 or 4 atoms of one structure, in the order given ('sel' keeps the picking
    order). The animation and plot are chosen as in ``plot_for``: one plot per animation, so
    an angle on the animation of plotted distances goes in a panel under them. With the GUI
    the plot's window opens (or is updated); ``save_data`` (CSV) / ``save_image`` (PNG, SVG,
    PDF) write the whole plot, every panel and series. Returns this series' values, one per
    frame.
    """
    structures = list(atoms.unique_structures)        # before the count: ':10@CA :40@CA'
    if len(structures) > 1:                         # with several models open
        ids = ", ".join(f"#{s.id_string}" for s in structures)
        raise UserError(f"the atoms must be in one structure, not {ids} (give the model, e.g. "
                        "#1.1/A:10@CA)")
    n = MEASURES[kind]
    if len(atoms) != n:
        raise UserError(f"{A_KIND[kind]} needs {n} atoms, got {len(atoms)} (used in the order "
                        "given; 'sel' keeps the picking order)")
    if save_image is not None and os.path.splitext(save_image)[1].lower() not in IMAGE_SUFFIXES:
        raise UserError(f"saveImage {save_image}: the image must be "
                        f"{', '.join(IMAGE_SUFFIXES)}")
    for option, path in (("saveData", save_data), ("saveImage", save_image)):
        if path is not None:                        # a bad file changes nothing
            _check_target(option, path)
    plot = plot_for(session, structures[0], mode, scope, amplitude, frames, new)
    series = add_series(session, plot, _atom_keys(atoms))
    if session.ui.is_gui:
        from .analyze_tool import show_plot

        show_plot(session, plot)
    if save_data is not None:
        _save(session, plot, save_data, "data")
    if save_image is not None:
        _save(session, plot, save_image, "image")
    return series.values


mdcif_plot_desc = CmdDesc(
    required=[("kind", KindArg), ("atoms", AtomsArg)],
    keyword=[("mode", IntArg), ("scope", ScopeArg), ("amplitude", FloatArg), ("frames", IntArg),
             ("new", BoolArg), ("save_data", SaveFileNameArg), ("save_image", SaveFileNameArg)],
    synopsis="Plot a distance, angle or dihedral along a mode animation",
)
