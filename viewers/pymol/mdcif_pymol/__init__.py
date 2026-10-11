"""PyMOL plugin for mdcif files.

Usage (models load as states of one object)::

    mdcif_load gb1.cif                          load + list models and modes
    mdcif_modes gb1                             list models, populations and modes again
    mdcif_morph gb1, steps=15                   new object gb1_morph cycling through the models
    mdcif_animate gb1, model=1, mode=1          new object gb1_m1_mode1 oscillating along mode 1
    mdcif_animate gb1, mode=1, scope=global     ... along global PC1 (the clustering axis)
    mdcif_porcupine gb1, model=1, mode=1        CGO arrows gb1_m1_mode1_porcupine

``amplitude`` / ``scale`` are in standard deviations of the mode (k in
k * sqrt(lambda) * v). Self-contained: needs only PyMOL and numpy.
"""

import os

import numpy as np

from .mdcif_data import MdcifFile, arrows, mode_frames, morph_frames

__version__ = "0.1.0"

_FILES = {}  # object name -> MdcifFile


def _data(obj):
    if obj not in _FILES:
        raise ValueError(f"{obj!r} was not loaded with mdcif_load")
    return _FILES[obj]


def _bool(v):
    return str(v).lower() in ("1", "true", "yes", "on")


def _is_global(scope):
    """``scope`` is "cluster" (the model's own modes) or "global" (the clustering PCs)."""
    scope = str(scope).lower()
    if scope not in ("cluster", "global"):
        raise ValueError(f"scope must be 'cluster' or 'global', not {scope!r}")
    return scope == "global"


def _state(cmd, obj, state):
    """Atom keys (chain, resi, name) and coordinates of one state, in object order."""
    keys = []
    cmd.iterate(obj, "keys.append((chain, resi, name))", space={"keys": keys})
    return keys, cmd.get_coords(obj, state=state)


def _new_states(cmd, src, model, name, frames):
    """Object ``name`` with one state per coordinate array in ``frames`` (src atom order)."""
    cmd.delete(name)
    for i, xyz in enumerate(frames, 1):
        cmd.create(name, src, source_state=model, target_state=i, quiet=1)
        cmd.load_coords(xyz, name, state=i)
    return name


def mdcif_load(path, name=None, quiet=0):
    """Load an mdcif file (models as states) and read its modes."""
    from pymol import cmd

    name = name or os.path.splitext(os.path.basename(path))[0]
    data = MdcifFile(path)
    cmd.load(path, name)
    cmd.dss(name, state=0)  # mdcif files carry no secondary structure records
    _FILES[name] = data
    if not int(quiet):
        print(data.summary())
    return name


def mdcif_modes(obj=None):
    """Print the models, cluster populations and modes of loaded mdcif objects."""
    for name in [obj] if obj else list(_FILES):
        print(f" {name}: " + _data(name).summary())


def mdcif_morph(obj, steps=10, wrap=1):
    """Object ``<obj>_morph``: linear morph through all models (cyclic by default)."""
    from pymol import cmd

    data = _data(obj)
    coords = [_state(cmd, obj, m)[1] for m in data.models]
    name = _new_states(cmd, obj, data.models[0], f"{obj}_morph",
                       morph_frames(coords, int(steps), wrap=_bool(wrap)))
    cmd.disable(obj)
    cmd.mplay()
    return name


def mdcif_animate(obj, model=1, mode=1, amplitude=3.0, frames=30, scope="cluster"):
    """Object oscillating along one mode of one model (``scope=global``: global PCs)."""
    from pymol import cmd

    data, model, mode = _data(obj), int(model), int(mode)
    use_global = _is_global(scope)
    ms = data.mode_set(model, use_global)
    keys, xyz = _state(cmd, obj, model)
    tag = f"global{mode}" if use_global else f"mode{mode}"
    name = _new_states(cmd, obj, model, f"{obj}_m{model}_{tag}",
                       mode_frames(xyz, keys, ms, mode, float(amplitude), int(frames)))
    cmd.disable(obj)
    cmd.mplay()
    print(f" {name}: {ms.scope} mode {mode} of model {model}, "
          f"{100 * ms.variance_fraction[mode - 1]:.0f}% variance, +/-{float(amplitude):g} sd")
    return name


def mdcif_porcupine(obj, model=1, mode=1, scale=3.0, color="orange", scope="cluster",
                    min_length=0.5, radius=0.15):
    """CGO arrows along one mode of one model (arrow = scale * sqrt(lambda) * v)."""
    from pymol import cmd
    from pymol.cgo import CONE, CYLINDER

    data, model, mode = _data(obj), int(model), int(mode)
    use_global = _is_global(scope)
    ms = data.mode_set(model, use_global)
    keys, xyz = _state(cmd, obj, model)
    starts, ends = arrows(xyz, keys, ms, mode, float(scale), float(min_length))
    rgb = list(cmd.get_color_tuple(color))
    r = float(radius)
    cgo = []
    for s, e in zip(starts, ends, strict=True):
        length = np.linalg.norm(e - s)
        head = min(1.0, 0.4 * length)
        head_r = min(2.7 * r, 0.5 * head + r)  # smaller heads on short arrows
        neck = e - (e - s) * head / length
        cgo += [CYLINDER, *s, *neck, r, *rgb, *rgb]
        cgo += [CONE, *neck, *e, head_r, 0.0, *rgb, *rgb, 1.0, 0.0]
    tag = f"global{mode}" if use_global else f"mode{mode}"
    name = f"{obj}_m{model}_{tag}_porcupine"
    cmd.delete(name)
    cmd.load_cgo(cgo, name)
    cmd.frame(model)
    print(f" {name}: {len(starts)} arrows, {ms.scope} mode {mode} of model {model}, "
          f"{float(scale):g} sd")
    return name


def __init_plugin__(app=None):
    from pymol import cmd

    for f in (mdcif_load, mdcif_modes, mdcif_morph, mdcif_animate, mdcif_porcupine):
        cmd.extend(f.__name__, f)
    print(f" mdcif_pymol {__version__} loaded: mdcif_load, mdcif_modes, mdcif_morph, "
          "mdcif_animate, mdcif_porcupine")
