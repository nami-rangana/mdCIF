"""Attach mdcif data to PyMOL objects, however the file was loaded.

``load file.cif`` (also File > Open, drag and drop, ``fetch`` of a saved file) goes through
``pymol.importing.loadfunctions["cif"]``, where ``register`` puts ``load_hook``. Ordinary mmCIF
files cost a byte scan for "_mdcif_" (about 0.15 s per 300 MB from the page cache) and go on to
PyMOL's own reader unchanged. An mdcif file is read by PyMOL's reader into a temporary object,
whose states (one per model) become one object each, in a group named after the file:

    gb1                 group (the name 'load' would give the object)
      gb1_m1            model 1 (state 1 of the file), 'dss' run, its title the cluster label
      gb1_m2            model 2, ...

PyMOL objects can't carry Python attributes, so the mdcif data lives in a registry here,
``name -> Entry``: the ``MdcifFile``, the model number, the file's path and group, and for
'mdcif_animate' copies the mode, scope and amplitude they play and the model object they animate.
An entry is alive while its object exists (and was not replaced by another object of the same
name). Sessions keep the registry (``session_save`` / ``session_restore``, PyMOL's session tasks);
a restored session reads the files again.
"""

import os

REGISTRY = {}                           # object name -> Entry, in the order they were made
CIF_SUFFIXES = (".cif", ".mmcif")
MARK = b"_mdcif_"
SESSION_KEY = "mdcif_pymol"


def path_key(path):
    return os.path.realpath(os.path.expanduser(str(path)))


def _cmd():
    from pymol import cmd

    return cmd


class Entry:
    """The mdcif data of one PyMOL object: a model of a file or an 'mdcif_animate' copy.

    ``data`` is None when the file was gone at session restore (``path`` says which file).
    ``source``: for a copy, the name of the model object it animates; ``made``: who made a copy
    ("plot", "animate" for the Animate window, None for the command), so they know they may
    delete it when they animate again.
    """

    def __init__(self, name, data, model, path, group=None, copy=False, mode=None, scope=None,
                 amplitude=None, source=None, made=None):
        self.name, self.data, self.model, self.path = name, data, int(model), path
        self.group = group
        self.copy, self.mode, self.scope = bool(copy), mode, scope
        self.amplitude = None if amplitude is None else float(amplitude)
        self.source, self.made = source, made
        self.n_atoms = None

    def __repr__(self):
        what = f"copy of model {self.model}" if self.copy else f"model {self.model}"
        return f"<mdcif Entry {self.name}: {what}>"

    @property
    def alive(self):
        """Whether the object is still there, and still this entry's."""
        if REGISTRY.get(self.name) is not self:
            return False
        cmd = _cmd()
        if self.name not in cmd.get_names("objects"):
            return False
        return self.n_atoms is None or cmd.count_atoms(f"%{self.name}") == self.n_atoms

    @property
    def deleted(self):
        return not self.alive

    @property
    def n_frames(self):
        return _cmd().count_states(f"%{self.name}") if self.alive else 0

    def fields(self):
        """What a session keeps (plain values)."""
        return {"path": self.path, "model": self.model, "group": self.group, "copy": self.copy,
                "mode": self.mode, "scope": self.scope, "amplitude": self.amplitude,
                "source": self.source}


def register_entry(entry):
    """Put ``entry`` in the registry under its name (replacing whatever was there)."""
    REGISTRY.pop(entry.name, None)                  # a new entry goes to the end
    REGISTRY[entry.name] = entry
    entry.n_atoms = _cmd().count_atoms(f"%{entry.name}")
    return entry


def entry(name):
    """The alive entry of object ``name``, or None."""
    found = REGISTRY.get(str(name))
    return found if found is not None and found.alive else None


def entries():
    """The alive entries, oldest first; dead ones are dropped."""
    for name in [n for n, e in REGISTRY.items() if not e.alive]:
        del REGISTRY[name]
    return list(REGISTRY.values())


def rename(entry, new):
    """Rename the entry's object to ``new`` (``cmd.set_name``) and follow in the registry."""
    cmd = _cmd()
    old = entry.name
    cmd.set_name(old, new)
    del REGISTRY[old]
    entry.name = new
    REGISTRY[new] = entry
    for e in REGISTRY.values():                     # copies that animate a renamed model
        if e.source == old:
            e.source = new


def has_mdcif_categories(path, chunk=1 << 22):
    """Whether the file contains "_mdcif_" (a byte scan in 4 MB blocks, no parsing)."""
    keep = len(MARK) - 1
    tail = b""
    try:
        with open(path, "rb") as fh:
            while block := fh.read(chunk):
                if MARK in block or MARK in tail + block[:keep]:
                    return True
                tail = block[-keep:]
    except OSError:
        return False
    return False


# --------------------------------------------------------------------------- loading

_previous = {}                          # what loadfunctions["cif"] was before the hook


def _raw_load(path, name, quiet=1):
    """PyMOL's own mmCIF reader (the hook, or whoever had "cif" before it, put aside)."""
    import pymol.importing as importing

    cmd = _cmd()
    mine = importing.loadfunctions.get("cif")
    hooked = mine is load_hook
    if hooked:
        if "cif" in _previous:
            importing.loadfunctions["cif"] = _previous["cif"]
        else:
            del importing.loadfunctions["cif"]
    try:
        cmd.load(path, name, format="cif", quiet=quiet, zoom=0)
    finally:
        if hooked:
            importing.loadfunctions["cif"] = mine


def load_mdcif(path, name=None, data=None, quiet=0, zoom=-1):
    """Load an mdcif file as one object per model in a group; returns the model objects.

    ``data`` (an ``MdcifFile``) skips reading the file again. A file whose states don't match
    its models (or without models) is left as PyMOL's reader gives it, with a warning.
    """
    from .mdcif_data import MdcifFile

    cmd = _cmd()
    path = os.path.abspath(os.path.expanduser(str(path)))
    if data is None:
        data = MdcifFile(path)
    name = str(name or os.path.splitext(os.path.basename(path))[0]).strip()
    tmp = cmd.get_unused_name("_mdcif_load")
    _raw_load(path, tmp)
    n_states = cmd.count_states(f"%{tmp}")
    if n_states != len(data.models):
        cmd.delete(name)
        cmd.set_name(tmp, name)
        print(f" mdcif: {os.path.basename(path)} has {len(data.models)} models but PyMOL read "
              f"{n_states} states; loaded as one object {name} without mdcif data")
        return []
    for old in [e.name for e in entries() if e.group == name and not e.copy]:
        cmd.delete(old)                             # the same file loaded again
    cmd.delete(name)
    objects = []
    for state, model in enumerate(data.models, start=1):
        obj = f"{name}_m{model}"
        cmd.delete(obj)
        cmd.create(obj, f"%{tmp}", source_state=state, target_state=1, quiet=1, zoom=0)
        cmd.set_title(obj, 1, data.model_label(model))
        cmd.dss(obj, state=1)                       # mdcif files carry no secondary structure
        objects.append(obj)
    cmd.delete(tmp)
    cmd.group(name, " ".join(objects))
    for obj, model in zip(objects, data.models, strict=True):
        register_entry(Entry(obj, data, model, path, group=name))
    if int(zoom) != 0 and (int(zoom) > 0 or cmd.get_setting_int("auto_zoom")):
        cmd.zoom(name)
    if not int(quiet):
        print(data.summary())
        print(" mdcif: " + ", ".join(f"{o} = model {m}" for o, m in zip(
            objects, data.models, strict=True)) + f" (group {name})")
    return objects


def load_hook(**kw):
    """``loadfunctions["cif"]``: an mdcif file loads as one object per model (``load_mdcif``),
    anything else goes to PyMOL's reader (or the hook that was there before) unchanged."""
    import pymol.internal

    filename, state = kw.get("filename"), int(kw.get("state") or 0)
    if (state == 0 and not int(kw.get("partial") or 0) and isinstance(filename, str)
            and has_mdcif_categories(filename)):
        from .mdcif_data import MdcifFile

        try:
            data = MdcifFile(filename)
        except (ValueError, OSError) as exc:        # not readable here (BinaryCIF, ...)
            print(f" mdcif: {exc}; loading it as plain mmCIF")
        else:
            load_mdcif(filename, kw.get("object"), data=data, quiet=kw.get("quiet", 1),
                       zoom=kw.get("zoom", -1))
            return None
    previous = _previous.get("cif")
    if previous is not None:
        return _call_previous(previous, kw)
    return pymol.internal._load(kw["oname"], kw["finfo"], kw["state"], kw["ftype"],
                                kw["finish"], kw["discrete"], kw["quiet"], kw["multiplex"],
                                kw["zoom"], kw["mimic"], kw["plugin"], kw["object_props"],
                                kw["atom_props"], _self=kw["_self"])


def _call_previous(func, kw):
    """Call another load function with the arguments its signature names (as PyMOL does)."""
    import inspect

    sig = inspect.signature(func, follow_wrapped=False)
    if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()):
        return func(**kw)
    return func(**{n: kw[n] for n in sig.parameters if n in kw})


def register_hook():
    """Put ``load_hook`` in ``loadfunctions["cif"]`` (once)."""
    from pymol.importing import loadfunctions

    current = loadfunctions.get("cif")
    if current is not None and not _is_hook(current):    # an older copy of the plugin's hook
        _previous["cif"] = current                       # is replaced, anything else kept
    loadfunctions["cif"] = load_hook


def _is_hook(func):
    return func is load_hook or (getattr(func, "__name__", "") == "load_hook"
                                 and getattr(func, "__module__", "").endswith(".attach"))


def unregister_hook():
    from pymol.importing import loadfunctions

    if loadfunctions.get("cif") is load_hook:
        if "cif" in _previous:
            loadfunctions["cif"] = _previous.pop("cif")
        else:
            del loadfunctions["cif"]


def is_hooked():
    try:
        from pymol.importing import loadfunctions
    except ImportError:
        return False
    return loadfunctions.get("cif") is load_hook


# --------------------------------------------------------------------------- sessions

def session_save(session, _self=None):
    """Session task: the registry, as plain values."""
    session[SESSION_KEY] = {e.name: e.fields() for e in entries()}
    return None


def session_restore(session, _self=None):
    """Session task: entries again for the saved objects, each file read once; a file that is
    gone leaves entries without data (the commands then say so). Must return True."""
    from .mdcif_data import MdcifFile

    saved = session.get(SESSION_KEY) if isinstance(session, dict) else None
    if not saved:
        return True
    cmd = _cmd()
    names = set(cmd.get_names("objects"))
    files, gone = {}, set()
    attached = []
    for name, f in saved.items():
        if name not in names:
            continue
        path = f.get("path")
        key = path_key(path)
        if key not in files and key not in gone:
            try:
                files[key] = MdcifFile(path)
            except (ValueError, OSError):
                gone.add(key)
        data = files.get(key)
        register_entry(Entry(name, data, f.get("model", 1), path, group=f.get("group"),
                             copy=f.get("copy", False), mode=f.get("mode"),
                             scope=f.get("scope"), amplitude=f.get("amplitude"),
                             source=f.get("source")))
        attached.append(name)
    for key in sorted(gone):
        print(f" mdcif: {key} is gone (moved since the session was saved?): its objects have "
              "no modes; load the mdcif file again")
    if attached:
        print(f" mdcif: data attached to {' '.join(attached)}")
    return True


def _register_task(tasks, func):
    """Add ``func`` to a session task list, replacing an older copy (a reloaded plugin)."""
    tasks[:] = [t for t in tasks if not (getattr(t, "__name__", None) == func.__name__ and
                                         getattr(t, "__module__", "").endswith(".attach"))]
    tasks.append(func)


def register_sessions():
    import pymol

    _register_task(pymol._session_save_tasks, session_save)
    _register_task(pymol._session_restore_tasks, session_restore)


def register():
    """The load hook and the session tasks (the plugin's init)."""
    register_hook()
    register_sessions()
