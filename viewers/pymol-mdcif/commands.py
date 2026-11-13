"""The ``mdcif_*`` commands (ported from the ChimeraX bundle's cmd.py).

    load gb1.cif                                any mmCIF load (also File > Open, drag and
                                                drop): group gb1 with objects gb1_m1, gb1_m2, ...
    mdcif_load gb1.cif, gb1                     the same, as mmCIF whatever the suffix
    mdcif_modes gb1                             list models, cluster populations and modes
    mdcif_morph gb1, frames=20                  object gb1_morph through the models (mplay)
    mdcif_animate gb1_m1, mode=1, amplitude=3   copy gb1_m1_mode1 oscillating along mode 1
    mdcif_animate gb1_m2, mode=1, scope=global  ... along global PC1 (the clustering axis)
    mdcif_porcupine gb1_m1, mode=1, scale=3     CGO arrows gb1_m1_mode1_porcupine
    mstop / frame 8                             stop the animation / show one frame (PyMOL)
    mdcif_plot distance, gb1_m1///10/CA, gb1_m1///40/CA     the distance at every frame of a
                                                mode animation (a copy is made for it)
    mdcif_plot dihedral, pk1, pk2, pk3, pk4, save_data=d.csv    4 picked atoms, in order

``amplitude`` / ``scale`` are in standard deviations of the mode (k in k * sigma * v, sigma in
angstroms: sqrt(lambda) for PCA, with the elements' standard masses taken out for quasi-harmonic
modes). Text mmCIF only. The old calls on the group with ``model=`` (``mdcif_animate gb1,
model=2, mode=1``) still work.

Objects get their data from ``attach`` (on load, on session restore). An animate copy keeps its
model's mode data, so ``mdcif_porcupine`` works on it too; both start from the deposited
coordinates (the copy's state 1). Copies play on PyMOL's movie player (``mplay``, ``mstop``,
``frame``; speed ``movie_fps``): every object with several states follows one frame counter.
``mdcif_plot`` measures on an animate copy (made for it from a model, or the one given); one
plot per copy holds every kind, a panel each (``new=1`` starts another). Plots are in
``plot.py``, the Analyze window and its Geometry tab in ``analyze_tool.py`` /
``geometry_tab.py``.

The windows of the MDCIF menu (``menu.py``) work through the helpers here, which have no Qt:
``porcupine_cgo`` (the arrows of 'mdcif_porcupine') and ``porcupine_command`` (the command to
log), ``retune_animation`` (new frames for an animate copy, in place), ``retire_copy`` (delete or
disable a copy no window shows any more) and ``log_command`` (an equivalent command, shown as
PyMOL shows typed ones).
"""

import os

import numpy as np
from pymol import CmdException, cmd

from . import __version__, attach
from .mdcif_data import MEASURES, MdcifFile, arrows, measure, mode_frames, morph_frames
from .plot import A_KIND, IMAGE_SUFFIXES, KIND_OF, PLOT_MADE, UNITS, Plot, atom_label, plots
from .plot import save_image as _save_image
from .plot import write_csv as _write_csv

MIN_ANIMATE_FRAMES = 3                  # sin() is 0 at every frame of a 1- or 2-frame cycle
ANIMATE_DEFAULTS = (1, "cluster", 3.0, 30)      # mode, scope, amplitude, frames
PORCUPINE_COLOR = (255, 128, 0)         # 'mdcif_porcupine' arrows without 'color' (PyMOL orange)
ANIMATE_MADE = "animate"                # Entry.made of copies the Animate window made
SCOPES = ("cluster", "global")
COMMANDS = ("mdcif", "mdcif_load", "mdcif_open", "mdcif_modes", "mdcif_morph", "mdcif_animate",
            "mdcif_porcupine", "mdcif_plot", "mdcif_movie", "mdcif_window")


def _error(message):
    raise CmdException(message)


def message(exc):
    """The text of a CmdException (its str() starts with ' Error: ')."""
    return getattr(exc, "message", None) or str(exc)


def log_command(text):
    """Show ``text`` as PyMOL shows a typed command, and write it to the log file if one is open
    (what the windows do once a change has settled)."""
    print(f" PyMOL>{text}")
    try:
        cmd.log(text + "\n")
    except Exception:                               # no log file, or an old PyMOL
        pass


def gui_window():
    """PyMOL's main window (a QMainWindow), or None without the GUI."""
    try:
        from pymol import gui

        return gui.get_qtwindow()
    except Exception:
        return None


# --------------------------------------------------------------------------- arguments

def _bool(v):
    return str(v).strip().lower() in ("1", "true", "yes", "on")


def _int(v, what):
    try:
        f = float(v)
    except (TypeError, ValueError):
        _error(f"{what} must be a whole number, not {v!r}")
    if f != int(f):
        _error(f"{what} must be a whole number, not {v!r}")
    return int(f)


def _float(v, what):
    try:
        return float(v)
    except (TypeError, ValueError):
        _error(f"{what} must be a number, not {v!r}")


def _given(v):
    return v is not None and str(v).strip() not in ("", "None", "none")


def _opt(v, convert, what):
    return convert(v, what) if _given(v) else None


def _scope(scope):
    """``scope`` is "cluster" (the model's own modes) or "global" (the clustering PCs)."""
    s = str(scope).strip().lower()
    if s not in SCOPES:
        _error(f"scope must be 'cluster' or 'global', not {scope!r}")
    return s


def _rgb(color):
    """RGB 0-1 from a PyMOL colour name, 0xRRGGBB, #RRGGBB or (r, g, b) (0-1 or 0-255)."""
    if isinstance(color, (tuple, list, np.ndarray)):
        rgb = [float(c) for c in color][:3]
        return tuple(c / 255 for c in rgb) if max(rgb) > 1 else tuple(rgb)
    text = str(color).strip()
    hexa = text[1:] if text.startswith("#") else text[2:] if text.lower().startswith("0x") else ""
    if hexa:
        try:
            return tuple(int(hexa[k:k + 2], 16) / 255 for k in (0, 2, 4))
        except ValueError:
            _error(f"unknown colour {color!r}")
    try:
        rgb = cmd.get_color_tuple(text)
    except Exception:
        rgb = None
    if not rgb or len(rgb) < 3:
        _error(f"unknown colour {color!r}")
    return tuple(float(c) for c in rgb[:3])


def hex_color(rgb):
    """``0xrrggbb`` (a colour PyMOL commands take) from RGB 0-255."""
    r, g, b = (int(round(float(c))) for c in list(rgb)[:3])
    return f"0x{r:02x}{g:02x}{b:02x}"


# --------------------------------------------------------------------------- objects

def _keys(name):
    """Atom keys (chain, resi, name) of an object, in its atom order."""
    keys = []
    cmd.iterate(f"%{name}", "keys.append((chain, resi, name))", space={"keys": keys})
    return keys


def _rest_coords(name):
    """Coordinates the modes are drawn from: state 1 (an animate copy's first state is its
    deposited model), so arrows and new animations don't start mid-swing."""
    return cmd.get_coords(f"%{name}", 1)


def _check_data(e):
    if e.data is None:
        _error(f"{e.name}: its file {e.path} is gone (moved since the session was saved?); "
               "load the mdcif file again")
    return e


def model_entry(path, model):
    """The entry of a file's model (not an animate copy), or None."""
    for e in attach.entries():
        if (not e.copy and e.model == model
                and attach.path_key(e.path) == attach.path_key(path)):
            return e
    return None


def resolve(obj, model=None):
    """The entry of ``obj``: a model or an animate copy, or the group of a file with ``model``
    (default: its first model). For a model object ``model`` may name another model of its
    file."""
    obj = str(obj).strip()
    model = _opt(model, _int, "model")
    e = attach.entry(obj)
    if e is not None:
        if model is None or model == e.model:
            return _check_data(e)
        if e.copy:
            _error(f"{obj} is an animation of model {e.model}, not of model {model}")
        other = model_entry(e.path, model)
        if other is None:
            _error(f"{os.path.basename(e.path)} has no model {model} open")
        return _check_data(other)
    members = [x for x in attach.entries() if x.group == obj and not x.copy]
    if members:
        if model is None:
            return _check_data(members[0])
        found = next((x for x in members if x.model == model), None)
        if found is None:
            models = ", ".join(str(x.model) for x in members)
            _error(f"{obj} has no model {model} (models {models})")
        return _check_data(found)
    if obj in cmd.get_names("all"):
        _error(f"{obj} has no mdcif data (load an mdcif .cif file with 'load' or "
               "'mdcif_load')")
    _error(f"no object {obj!r}")


def mode_set(e, scope):
    """``(mode set, model)`` of an entry for ``scope``."""
    _check_data(e)
    try:
        return e.data.mode_set(e.model, scope == "global"), e.model
    except ValueError as exc:
        _error(str(exc))


def _model_name(e):
    """The model object a copy animates (its name when the copy was made)."""
    if not e.copy:
        return e.name
    return e.source or (f"{e.group}_m{e.model}" if e.group else e.name)


def _tag(scope, mode):
    return f"global{mode}" if scope == "global" else f"mode{mode}"


def _in_plot(e):
    return any(p.alive and p.entry is e for p in plots())


def _copy_name(base):
    """A name for a new animate copy: ``base``, unless an object of that name must be kept (it
    is not an animate copy, or a plot measures on it); then ``base_2``, ``base_3``, ..."""
    taken = set(cmd.get_names("all"))
    for k in range(1, 10000):
        name = base if k == 1 else f"{base}_{k}"
        if name not in taken:
            return name
        e = attach.entry(name)
        if e is not None and e.copy and not _in_plot(e):
            return name                             # replaced, as PyMOL replaces objects
    _error(f"no free name for {base}")


def _new_states(name, src, frames):
    """Object ``name``: the atoms of ``src`` (its state 1), one state per array of ``frames``
    (``src`` atom order)."""
    cmd.delete(name)
    for k, xyz in enumerate(frames, start=1):
        cmd.create(name, f"%{src}", source_state=1, target_state=k, quiet=1, zoom=0)
        cmd.load_coords(xyz, f"%{name}", state=k)
    return name


def _set_frames(name, frames):
    """Replace the states of object ``name`` by ``frames`` (its atom order), in place."""
    old = cmd.count_states(f"%{name}")
    for k in range(old + 1, len(frames) + 1):
        cmd.create(name, f"%{name}", source_state=1, target_state=k, quiet=1, zoom=0)
    if old > len(frames):
        cmd.delete_states(name, f"{len(frames) + 1}-{old}")
    for k, xyz in enumerate(frames, start=1):
        cmd.load_coords(xyz, f"%{name}", state=k)


def _play():
    """Play the animations from frame 1 (PyMOL's movie player, every multi-state object)."""
    cmd.frame(1)
    cmd.mplay()


def playing():
    try:
        return bool(cmd.get_movie_playing())
    except Exception:
        return False


# --------------------------------------------------------------------------- commands

def mdcif():
    """Report the plugin's status."""
    hook = "on" if attach.is_hooked() else "off"
    where = "" if gui_window() is not None else " (when PyMOL has its window)"
    print(f" mdcif_pymol {__version__}. Commands: mdcif_load | mdcif_modes | mdcif_morph | "
          f"mdcif_animate | mdcif_porcupine | mdcif_plot | mdcif_movie | mdcif_window. "
          f"Attaching mdcif data to loaded .cif files: {hook}. Windows: menu MDCIF > "
          f"Porcupine, Animate, Analyze{where}.")


def mdcif_load(path, name=None, quiet=0):
    """Load an mdcif file: one object per model, in a group (also for a file without modes, or
    a suffix other than .cif). Returns the model objects."""
    path = str(path).strip()
    try:
        data = MdcifFile(path)
    except (ValueError, OSError) as exc:            # OSError: a directory, unreadable file
        _error(str(exc))
    objects = attach.load_mdcif(path, name if _given(name) else None, data=data,
                                quiet=_int(quiet, "quiet"))
    if not objects:
        _error(f"{os.path.basename(path)}: its models don't match the states PyMOL read")
    return objects


def mdcif_modes(obj=None):
    """List the models, cluster populations and modes of mdcif objects (all without ``obj``)."""
    found = attach.entries()
    if _given(obj):
        try:
            names = set(cmd.get_object_list(f"({obj})") or [])
        except Exception:
            names = {str(obj)}
        found = [e for e in found if e.name in names or e.group == str(obj)]
    seen = []
    for e in found:
        if e.data is not None and all(e.data is not d for d in seen):
            seen.append(e.data)
            print(e.data.summary())
            objects = [f"{x.name} = model {x.model}" for x in attach.entries()
                       if x.data is e.data and not x.copy]
            copies = [f"{x.name} ({x.scope} mode {x.mode} of model {x.model})"
                      for x in attach.entries() if x.data is e.data and x.copy]
            print(" objects: " + ", ".join(objects) + (
                f"; animations: {', '.join(copies)}" if copies else ""))
    if not seen:
        print(" mdcif: no objects with mdcif data")


def mdcif_morph(obj, frames=20, wrap=1, steps=None, quiet=0):
    """Object ``<group>_morph``: a linear morph through the models in ``obj`` (a group or
    model objects; animate copies don't count), ``frames`` per step, back to the first model
    with ``wrap``. Plays with PyMOL's movie player. Returns its name."""
    if _given(steps):
        frames = steps
    frames = _int(frames, "frames")
    if frames < 1:
        _error("frames must be at least 1")
    try:
        names = cmd.get_object_list(f"({obj})") or []
    except Exception:
        names = []
    targets = [e for e in (attach.entry(n) for n in names) if e is not None and not e.copy]
    if len(targets) < 2:
        _error("need at least two models with mdcif data (animate copies don't count)")
    order = {}                                      # files in load order, then model number
    for e in attach.entries():
        order.setdefault(attach.path_key(e.path), len(order))
    targets.sort(key=lambda e: (order[attach.path_key(e.path)], e.model))
    keys = _keys(targets[0].name)
    if any(_keys(e.name) != keys for e in targets[1:]):
        _error("the models must have the same atoms (models of one file)")
    xyz = morph_frames([_rest_coords(e.name) for e in targets], frames, wrap=_bool(wrap))
    name = f"{targets[0].group or targets[0].name}_morph"
    _new_states(name, targets[0].name, xyz)
    for e in targets:
        cmd.disable(e.name)
    _play()
    if not _int(quiet, "quiet"):
        models = ", ".join(str(e.model) for e in targets)
        print(f" {name}: {len(xyz)} frames, models {models}"
              + (" and back" if _bool(wrap) else "") + f", {frames} frames per step "
              "('mstop' to stop)")
    return name


def mdcif_animate(obj, model=None, mode=1, amplitude=3.0, frames=30, scope="cluster", quiet=0):
    """Object ``<model>_mode<N>`` (``_global<N>``) oscillating along one mode of a model, from
    its deposited coordinates; the model is disabled and the animation plays. Returns its
    name."""
    mode, frames = _int(mode, "mode"), _int(frames, "frames")
    amplitude, scope = _float(amplitude, "amplitude"), _scope(scope)
    if frames < MIN_ANIMATE_FRAMES:
        _error(f"frames must be at least {MIN_ANIMATE_FRAMES}")
    e = resolve(obj, model)
    ms, model = mode_set(e, scope)
    try:
        xyz = mode_frames(_rest_coords(e.name), _keys(e.name), ms, mode, amplitude, frames)
    except ValueError as exc:
        _error(str(exc))
    source = _model_name(e)
    name = _copy_name(f"{source}_{_tag(scope, mode)}")
    _new_states(name, e.name, xyz)
    copy = attach.register_entry(attach.Entry(
        name, e.data, model, e.path, group=e.group, copy=True, mode=mode, scope=scope,
        amplitude=amplitude, source=source))
    cmd.disable(e.name)
    _play()
    if not _int(quiet, "quiet"):
        print(f" {name}: {ms.scope} mode {mode} of model {model}, "
              f"{100 * ms.variance_fraction[mode - 1]:.0f}% variance, +/-{abs(amplitude):g} sd "
              "('mstop' to stop)")
    return copy.name


def animate_command(name, mode, amplitude, frames, scope):
    """The 'mdcif_animate' command for these settings, as text (for the log)."""
    return (f"mdcif_animate {name}, mode={mode}, amplitude={amplitude:g}, frames={frames}, "
            f"scope={scope}")


def retune_animation(copy, mode=None, scope=None, amplitude=None, frames=None):
    """Give an animate copy the frames of other settings, in place; None keeps a setting. What
    the Animate window does when a setting changes.

    The frames start again from the copy's state 1 (the deposited model); the entry's mode,
    scope and amplitude follow, and so does the name (``gb1_m1_mode2``) when that name is free.
    The movie player goes on as it was. Returns the copy, or None (and changes nothing) when a
    plot measures on it, whose values would no longer match the frames: the caller makes a new
    copy then.
    """
    if copy is None or not copy.alive or not copy.copy:
        _error("not an 'mdcif_animate' copy")
    if _in_plot(copy):
        return None
    wanted = [have if v is None else v for v, have in zip(
        (mode, scope, amplitude), (copy.mode, copy.scope, copy.amplitude), strict=True)]
    mode, scope, amplitude = (d if v is None else v for v, d in zip(
        wanted, ANIMATE_DEFAULTS[:3], strict=True))     # copies without settings: defaults
    frames = copy.n_frames if frames is None else int(frames)
    if frames < MIN_ANIMATE_FRAMES:
        _error(f"frames must be at least {MIN_ANIMATE_FRAMES}")
    ms, _ = mode_set(copy, scope)
    try:
        xyz = mode_frames(_rest_coords(copy.name), _keys(copy.name), ms, mode, amplitude,
                          frames)
    except ValueError as exc:
        _error(str(exc))
    shown = cmd.get_state()
    _set_frames(copy.name, xyz)
    copy.mode, copy.scope, copy.amplitude = int(mode), scope, float(amplitude)
    base = f"{_model_name(copy)}_{_tag(scope, mode)}"
    if not copy.name.startswith(base) and base not in cmd.get_names("all"):
        attach.rename(copy, base)
    if not playing() and shown > frames:
        cmd.frame(frames)
    return copy


def retire_copy(e, own):
    """Put away an animate copy that a plot or the Animate window has moved off: deleted if
    ``own`` (they made it) and no plot measures on it, else disabled."""
    if e is None or not e.alive:
        return
    if own and not _in_plot(e):
        cmd.delete(e.name)
    else:
        cmd.disable(e.name)


def _arrow_cgo(starts, ends, radius, rgb):
    """CGO arrows: a cylinder and a cone each, smaller heads on short arrows."""
    from pymol.cgo import CONE, CYLINDER

    cgo = []
    r = float(radius)
    for s, e in zip(starts, ends, strict=True):
        length = float(np.linalg.norm(e - s))
        if length < 1e-6:
            continue
        head = min(1.0, 0.4 * length)
        head_r = min(2.7 * r, 0.5 * head + r)
        neck = e - (e - s) * head / length
        cgo += [CYLINDER, *s, *neck, r, *rgb, *rgb]
        cgo += [CONE, *neck, *e, head_r, 0.0, *rgb, *rgb, 1.0, 0.0]
    return cgo


def porcupine_cgo(e, mode=1, scale=3.0, scope="cluster", min_length=0.5, radius=0.15,
                  rgb=(1.0, 0.5, 0.0)):
    """The arrows 'mdcif_porcupine' draws, without drawing them: ``(cgo, n_arrows, n_short,
    mode_set)``, ``n_short`` the arrows left out as shorter than ``min_length``.
    CmdException for what the command refuses (no arrows left, a bad radius, ...)."""
    if radius <= 0:
        _error("radius must be positive")
    if min_length < 0:
        _error("min_length must not be negative")
    ms, _ = mode_set(e, scope)
    coords, keys = _rest_coords(e.name), _keys(e.name)
    try:
        starts, ends = arrows(coords, keys, ms, mode, scale, min_length)
        n_all = len(arrows(coords, keys, ms, mode, scale, 0.0)[0])
    except ValueError as exc:
        _error(str(exc))
    cgo = _arrow_cgo(starts, ends, radius, rgb)
    if not cgo:
        _error(f"no arrows to draw for {ms.scope} mode {mode} at {scale:g} sd with min_length "
               f"{min_length:g} A")
    return cgo, len(starts), n_all - len(starts), ms


def porcupine_name(e, scope, mode):
    return f"{e.name}_{_tag(scope, mode)}_porcupine"


def porcupine_command(name, mode=1, scale=3.0, color=None, scope="cluster", min_length=0.5,
                      radius=0.15):
    """The 'mdcif_porcupine' command for these settings, as text (for the log); ``color``:
    RGB 0-255 (None: the default orange)."""
    color = hex_color(PORCUPINE_COLOR if color is None else color)
    return (f"mdcif_porcupine {name}, mode={mode}, scale={scale:g}, color={color}, "
            f"scope={scope}, min_length={min_length:g}, radius={radius:g}")


def mdcif_porcupine(obj, model=None, mode=1, scale=3.0, color="orange", scope="cluster",
                    min_length=0.5, radius=0.15, name=None, quiet=0):
    """CGO arrows along one mode of a model (arrow = scale * sigma * v), object
    ``<obj>_mode<N>_porcupine`` (replaced if there); on an animate copy from its state 1.
    Returns the arrows' name."""
    mode, scale = _int(mode, "mode"), _float(scale, "scale")
    min_length, radius = _float(min_length, "min_length"), _float(radius, "radius")
    scope = _scope(scope)
    e = resolve(obj, model)
    cgo, n, _, ms = porcupine_cgo(e, mode, scale, scope, min_length, radius, _rgb(color))
    name = str(name).strip() if _given(name) else porcupine_name(e, scope, mode)
    cmd.delete(name)
    cmd.load_cgo(cgo, name, 1, zoom=0)
    if not _int(quiet, "quiet"):
        print(f" {name}: {n} arrows, {ms.scope} mode {mode} of model {e.model}, {scale:g} sd")
    return name


# --------------------------------------------------------------------------- movies

MOVIE_SUFFIXES = (".gif", ".mp4", ".webm", ".mov")
VIDEO_SUFFIXES = (".mp4", ".webm", ".mov")      # need ffmpeg


def ffmpeg_exe():
    """The ffmpeg PyMOL would use (its own folders and PATH), or None."""
    try:
        from pymol.movie import find_exe

        return find_exe("ffmpeg")
    except Exception:
        import shutil

        return shutil.which("ffmpeg")


def _gif(images, path, fps):
    """An animated GIF (looping) from PNG images, one shared palette."""
    import io

    from PIL import Image

    frames = [Image.open(io.BytesIO(b)).convert("RGB") for b in images]
    sample = frames[::max(1, len(frames) // 8)][:8]     # the palette sees every part of the cycle
    w, h = frames[0].size
    sheet = Image.new("RGB", (w * len(sample), h))
    for k, im in enumerate(sample):
        sheet.paste(im, (k * w, 0))
    palette = sheet.quantize(colors=255)
    dither = getattr(getattr(Image, "Dither", Image), "NONE", 0)
    out = [im.quantize(palette=palette, dither=dither) for im in frames]
    out[0].save(path, save_all=True, append_images=out[1:], duration=round(1000 / fps), loop=0,
                optimize=False, disposal=1)


def _video(images, path, fps, exe):
    """MP4 / WebM / MOV from PNG images with ffmpeg."""
    import subprocess
    import tempfile

    with tempfile.TemporaryDirectory(prefix="mdcif_movie") as folder:
        for k, b in enumerate(images, start=1):
            with open(os.path.join(folder, f"frame{k:04d}.png"), "wb") as fh:
                fh.write(b)
        args = [exe, "-y", "-loglevel", "error", "-framerate", f"{fps:g}", "-i",
                os.path.join(folder, "frame%04d.png"), "-vf",
                "pad=ceil(iw/2)*2:ceil(ih/2)*2"]
        if path.lower().endswith((".mp4", ".mov")):
            args += ["-pix_fmt", "yuv420p"]
        done = subprocess.run(args + [path], capture_output=True, text=True)
        if done.returncode != 0:
            _error(f"ffmpeg could not write {path}: {done.stderr.strip()[-300:]}")


def mdcif_movie(filename, obj=None, width=0, height=0, ray=None, quiet=0):
    """Record one cycle of an animation as a movie: frames 1 to the number of states of
    ``obj`` (any object: an animate copy, a morph, ...; default: the newest shown animate copy,
    else every frame of PyMOL's movie), at ``movie_fps``. GIF always (Pillow); MP4, WebM and MOV
    with ffmpeg. ``ray``: ray-trace each frame (default: only without the GUI). The animation
    goes on as it was. Returns the path."""
    path = os.path.abspath(os.path.expanduser(str(filename).strip()))
    suffix = os.path.splitext(path)[1].lower()
    if suffix not in MOVIE_SUFFIXES:
        _error(f"{filename}: a movie must be {', '.join(MOVIE_SUFFIXES)}")
    exe = ffmpeg_exe() if suffix in VIDEO_SUFFIXES else None
    if suffix in VIDEO_SUFFIXES and exe is None:
        _error(f"{suffix} needs ffmpeg, which PyMOL can't find (File > Export Movie offers to "
               "install it); a .gif works without it")
    _check_target("filename", path)
    if _given(obj):
        name = str(obj).strip()
        if name not in cmd.get_names("objects"):
            _error(f"no object {name!r}")
        n = cmd.count_states(f"%{name}")
    else:
        enabled = set(cmd.get_names("objects", enabled_only=1))
        copies = [c for c in attach.entries() if c.copy and c.name in enabled]
        n = copies[-1].n_frames if copies else cmd.count_frames()
    if n < 1:
        _error("nothing to record")
    ray = (0 if gui_window() is not None else 1) if not _given(ray) else _int(ray, "ray")
    width, height = _int(width, "width"), _int(height, "height")
    fps = max(cmd.get_setting_float("movie_fps"), 0.1)
    was, shown = playing(), cmd.get_frame()
    cmd.mstop()
    images = []
    try:
        for k in range(1, n + 1):
            cmd.frame(k)
            try:
                images.append(cmd.png(None, width, height, ray=ray))
            except CmdException:
                if ray or not (width or height):
                    raise
                ray = 1                                 # incentive PyMOL 3.1 draws only at
                print(" mdcif: PyMOL can't draw pictures of this size without ray tracing; "
                      "ray-tracing the frames")         # the viewport's size
                images.append(cmd.png(None, width, height, ray=ray))
    finally:
        cmd.frame(shown)
        if was:
            cmd.mplay()
    if any(not b for b in images):
        _error("PyMOL gave no image (is the viewer shown?)")
    if suffix == ".gif":
        try:
            _gif(images, path, fps)
        except ImportError:
            _error("a GIF needs Pillow in PyMOL's Python (python -m pip install pillow)")
    else:
        _video(images, path, fps, exe)
    if not _int(quiet, "quiet"):
        print(f" mdcif: movie saved to {path} ({n} frames at {fps:g} frames per second)")
    return path


# --------------------------------------------------------------------------- plots

def _matches(copy, mode=None, scope=None, amplitude=None, frames=None):
    """Whether an animate copy was made with these settings (None = any)."""
    have = (copy.mode, copy.scope, copy.amplitude, copy.n_frames)
    return have[0] is not None and all(w is None or w == h for w, h in zip(
        (mode, scope, amplitude, frames), have, strict=True))


def _of_file(e, path):
    """Whether entry ``e`` has the data of file ``path`` (by path: a restored session reads
    the files again)."""
    return attach.path_key(e.path) == attach.path_key(path)


def shown_copies(e):
    """The enabled animate copies of the model ``e`` (same file and model), oldest first."""
    enabled = set(cmd.get_names("objects", enabled_only=1))
    return [c for c in attach.entries()
            if c.copy and c.name in enabled and _of_file(c, e.path) and c.model == e.model]


def _animate(e, mode=None, scope=None, amplitude=None, frames=None, made=PLOT_MADE):
    """'mdcif_animate' for a plot (``made=ANIMATE_MADE``: for the Animate window), logged as
    that command; None takes its default. The copy is marked as made by them, so they may
    delete it when they animate again. Returns the copy's entry."""
    mode, scope, amplitude, frames = (d if v is None else v for v, d in zip(
        (mode, scope, amplitude, frames), ANIMATE_DEFAULTS, strict=True))
    log_command(animate_command(e.name, mode, amplitude, frames, scope))
    copy = attach.entry(mdcif_animate(e.name, None, mode, amplitude, frames, scope))
    copy.made = made
    return copy


def _values(e, keys, kind):
    """``kind`` of the atoms ``keys`` at every state of object ``e``."""
    index = {}
    for i, k in enumerate(_keys(e.name)):
        index.setdefault(k, i)
    missing = [k for k in keys if k not in index]
    if missing:
        _error(f"{e.name} has no atom {atom_label(missing[0])}")
    rows = [index[k] for k in keys]
    xyz = np.array([cmd.get_coords(f"%{e.name}", k)[rows] for k in range(1, e.n_frames + 1)])
    return measure(xyz, kind)


def copy_for(e, mode=None, scope=None, amplitude=None, frames=None, new=False):
    """The animation atoms of entry ``e`` are measured on: ``(copy, source)`` entries.

    Atoms on an animate copy made with the settings given (None = any) use that copy. Atoms on
    a model (unless ``new``) use the newest copy a plot animated from it with these settings,
    else an enabled copy of that model with these settings (also one the user made with
    'mdcif_animate', or one restored from a session). Otherwise a new copy is animated
    ('mdcif_animate', default mode 1, cluster, 3 sd, 30 frames); other settings on a copy
    animate anew. ``source``: the model the copy animates, if known.
    """
    _check_data(e)
    wanted = (mode, scope, amplitude, frames)
    live = [p for p in plots() if p.alive]
    if e.copy:
        if _matches(e, *wanted):
            source = next((p.source for p in reversed(live)
                           if p.entry is e and p.source is not None), None)
            return e, source or attach.entry(e.source)
    elif not new:
        found = ([p.entry for p in live if p.source is e and _matches(p.entry, *wanted)]
                 or [c for c in shown_copies(e) if _matches(c, *wanted)])
        if found:
            return found[-1], e
    copy = _animate(e, *wanted)
    return copy, (e if not e.copy else attach.entry(copy.source))


def plot_for(e, mode=None, scope=None, amplitude=None, frames=None, new=False):
    """The plot (found or new) that atoms of entry ``e`` are measured in: the newest plot on
    the copy from ``copy_for``, whatever kinds it shows (a new kind gets a panel there),
    unless ``new``. A new plot joins ``plots()`` with its first series."""
    copy, source = copy_for(e, mode, scope, amplitude, frames, new)
    if not new:
        found = [p for p in plots() if p.alive and p.entry is copy]
        if found:
            return found[-1]
    return Plot(copy, source)


def add_series(plot, keys, quiet=0):
    """Measure the atoms ``keys`` (2, 3 or 4: a distance, angle or dihedral) on the plot's
    copy and add them to the plot (registered if new); prints the range. Atoms already in the
    plot, also in reverse order (the same measurement), are not added twice."""
    if not plot.alive:
        _error(f"{plot.spec}, the animation of this plot, was deleted")
    keys = [tuple(str(v) for v in k) for k in keys]
    kind = KIND_OF.get(len(keys))
    if kind is None:
        _error(f"a distance, angle or dihedral needs 2, 3 or 4 atoms, got {len(keys)}")
    twice = next((k for i, k in enumerate(keys) if k in keys[:i]), None)
    if twice is not None:
        _error(f"{atom_label(twice)} is given twice: {A_KIND[kind]} needs {len(keys)} "
               "different atoms")
    same = next((s for s in plot.series if s.keys in (keys, keys[::-1])), None)
    if same is not None:
        print(f" {plot.spec} {kind} {same.label} is already plotted")
        return same
    values = _values(plot.entry, keys, kind)
    try:
        series = plot.add(keys, values)
    except ValueError as exc:
        _error(str(exc))
    if plot not in plots():
        plots().append(plot)
    if not quiet:
        p = 2 if kind == "distance" else 1
        print(f" {plot.spec} {kind} {series.label} along {plot.describe()}: min "
              f"{values.min():.{p}f}, max {values.max():.{p}f}, mean {values.mean():.{p}f} "
              f"{UNITS[kind]}")
    plot.changed()
    return series


def reanimate(plot, model=None, mode=None, scope=None, amplitude=None, frames=None):
    """Animate a plot again with other settings (None = keep) and measure its series anew.

    The new copy is made from the model's object (for the same model: the plot's source, else
    the model, else the old copy). The old copy is then put away (``retire_copy``): deleted if
    the plot made it and no other plot uses it, else disabled. Returns the new copy.
    """
    model = plot.model if model is None else model
    mode, scope, amplitude, frames = (old if v is None else v for v, old in zip(
        (mode, scope, amplitude, frames),
        (plot.mode, plot.scope, plot.amplitude, plot.n_frames), strict=True))
    source = plot.source if (model == plot.model and plot.source is not None
                             and plot.source.alive) else None
    source = source or model_entry(plot.path, model)
    if source is None and model == plot.model and plot.alive:
        source = plot.entry
    if source is None:
        _error(f"no object of model {model} of {os.path.basename(plot.path)} is loaded")
    copy = _animate(source, mode, scope, amplitude, frames)
    try:
        values = [_values(copy, s.keys, s.kind) for s in plot.series]
    except CmdException:
        cmd.delete(copy.name)
        raise
    old, old_own = plot.entry, plot.own
    plot.set_animation(copy, source if not source.copy else attach.entry(copy.source))
    for s, v in zip(plot.series, values, strict=True):
        s.values = v
    if old is not copy:
        retire_copy(old, old_own)
    plot.changed()
    return copy


def _check_target(option, path):
    """Refuse a file that can't be written before anything is animated or added."""
    full = os.path.abspath(os.path.expanduser(str(path)))
    folder = os.path.dirname(full)
    if os.path.isdir(full):
        _error(f"{option} {path}: that is a folder, not a file")
    if not os.path.isdir(folder):
        _error(f"{option} {path}: there is no folder {folder}")
    if not os.access(folder, os.W_OK) or (os.path.exists(full) and not os.access(full, os.W_OK)):
        _error(f"{option} {path}: not allowed to write there")


def save_plot_file(plot, path, what):
    """Write the plot's ``what`` ("data": CSV, "image": PNG/SVG/PDF) to ``path``."""
    path = os.path.expanduser(str(path))
    try:
        (_write_csv if what == "data" else _save_image)(plot, path)
    except ImportError:
        _error("saving a plot image needs matplotlib in PyMOL's Python (see the tutorial: "
               "python -m pip install matplotlib)")
    except (OSError, ValueError) as exc:
        _error(f"could not save the plot {what}: {exc}")
    print(f" {plot.spec} {' + '.join(plot.kinds)} plot: {what} saved to {path}")


def _atoms_of(selection):
    """``[(object, (chain, resi, name)), ...]`` of a selection, in atom order."""
    out = []
    try:
        cmd.iterate(f"({selection})", "out.append((model, chain, resi, name))",
                    space={"out": out})
    except Exception:
        _error(f"{selection!r} is not a valid selection")
    return [(m, (c, r, n)) for m, c, r, n in out]


def plot_atoms(kind, atoms):
    """The object and atom keys, in order, of the atoms for a ``kind`` measurement: one
    selection per atom (each one atom: ``pk1``, ``/gb1_m1//A/10/CA``), or one selection of
    exactly as many atoms as needed (taken in atom order)."""
    n = MEASURES[kind]
    atoms = [str(a).strip() for a in atoms if _given(a)]
    if len(atoms) == 1 and n > 1:
        picked = _atoms_of(atoms[0])
        if len(picked) != n:
            _error(f"{A_KIND[kind]} needs {n} atoms, {atoms[0]} has {len(picked)} (give one "
                   "selection per atom, in order: e.g. pk1, pk2)")
    else:
        if len(atoms) != n:
            _error(f"{A_KIND[kind]} needs {n} atoms, got {len(atoms)} (one selection per atom, "
                   "in order)")
        picked = []
        for k, sel in enumerate(atoms, start=1):
            found = _atoms_of(sel)
            if len(found) != 1:
                _error(f"atom {k} ({sel}) selects {len(found)} atoms, not 1")
            picked.append(found[0])
    objects = list(dict.fromkeys(o for o, _ in picked))
    if len(objects) > 1:
        _error(f"the atoms must be in one object, not {', '.join(objects)} (give the object, "
               "e.g. /gb1_m1//A/10/CA)")
    return objects[0], [k for _, k in picked]


def mdcif_plot(kind, atom1=None, atom2=None, atom3=None, atom4=None, mode=None, scope=None,
               amplitude=None, frames=None, new=0, save_data=None, save_image=None, quiet=0):
    """Plot a distance, angle or dihedral at every frame of a mode animation.

    The atoms (2, 3 or 4, in order) are on a model or an animate copy; the animation and plot
    are chosen as in ``plot_for``: one plot per animation, so an angle on the animation of
    plotted distances goes in a panel under them. With the GUI the plot's window opens (or is
    updated); ``save_data`` (CSV) / ``save_image`` (PNG, SVG, PDF) write the whole plot, every
    panel and series. Returns this series' values, one per frame.
    """
    kind = str(kind).strip().lower()
    if kind not in MEASURES:
        _error(f"unknown measurement {kind!r} (distance, angle or dihedral)")
    obj, keys = plot_atoms(kind, (atom1, atom2, atom3, atom4))
    e = attach.entry(obj)
    if e is None:
        _error(f"{obj} has no mdcif data (load an mdcif .cif file with 'load')")
    mode, frames = _opt(mode, _int, "mode"), _opt(frames, _int, "frames")
    amplitude = _opt(amplitude, _float, "amplitude")
    scope = _scope(scope) if _given(scope) else None
    if frames is not None and frames < MIN_ANIMATE_FRAMES:
        _error(f"frames must be at least {MIN_ANIMATE_FRAMES}")
    save_data = str(save_data).strip() if _given(save_data) else None
    save_image = str(save_image).strip() if _given(save_image) else None
    if save_image is not None and os.path.splitext(save_image)[1].lower() not in IMAGE_SUFFIXES:
        _error(f"save_image {save_image}: the image must be {', '.join(IMAGE_SUFFIXES)}")
    for option, path in (("save_data", save_data), ("save_image", save_image)):
        if path is not None:                        # a bad file changes nothing
            _check_target(option, path)
    plot = plot_for(e, mode, scope, amplitude, frames, _bool(new))
    series = add_series(plot, keys, quiet=_int(quiet, "quiet"))
    if gui_window() is not None:
        from .analyze_tool import show_plot

        show_plot(plot)
    if save_data is not None:
        save_plot_file(plot, save_data, "data")
    if save_image is not None:
        save_plot_file(plot, save_image, "image")
    return [float(v) for v in series.values]


def mdcif_window(name="analyze"):
    """Open (or bring to the front) a window of the MDCIF menu: porcupine, animate or
    analyze."""
    from .menu import WINDOWS, show_window

    key = str(name).strip().lower()
    if key not in WINDOWS:
        _error(f"no window {name!r} (porcupine, animate or analyze)")
    if gui_window() is None:
        _error("the MDCIF windows need PyMOL's window (not available with -c)")
    return show_window(key)
