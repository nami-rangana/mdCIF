"""Attach mdcif data to structures, however the file was opened.

``open file.cif`` (also File > Open and drag and drop) reads an mdcif file with ChimeraX's own
mmCIF reader. The bundle's custom init puts a handler on the "add models" trigger that finds
the mdcif content of new structures and attaches it the way ``mdcif open`` does:

    s.mdcif_data    MdcifFile of the file (not saved in sessions)
    s.mdcif_model   model number in the file          | registered attributes: saved in
    s.mdcif_path    absolute path of the file         | sessions, so a restored session
    s.mdcif_copy    True for 'mdcif animate' copies   | gets its data back from the file

The structures of one file are matched to its models in id order (#1.1 = model 1), or by their
saved ``mdcif_model`` (restored sessions, animate copies). Ordinary mmCIF files only cost a
byte scan for "_mdcif_" (about 0.15 s per 300 MB from the page cache) and are left alone.
"""

import os

REGISTRANT = "MDCIF"
ATTRS = (("mdcif_path", str), ("mdcif_model", int), ("mdcif_copy", bool))
CIF_SUFFIXES = (".cif", ".mmcif")
MARK = b"_mdcif_"

pending = {}                            # path_key -> MdcifFile parsed by 'mdcif open'
_handlers = {}                          # id(session) -> "add models" trigger handler


def path_key(path):
    return os.path.realpath(os.path.expanduser(str(path)))


def register(session):
    """Register the saved attributes and the "add models" handler (bundle custom init)."""
    from chimerax.atomic import AtomicStructure
    from chimerax.core.models import ADD_MODELS

    for name, kind in ATTRS:
        AtomicStructure.register_attr(session, name, REGISTRANT, attr_type=kind)
    unregister(session)
    _handlers[id(session)] = session.triggers.add_handler(
        ADD_MODELS, lambda _trigger, models: on_add_models(session, models))


def unregister(session):
    handler = _handlers.pop(id(session), None)
    if handler is not None:
        session.triggers.remove_handler(handler)


def is_registered(session):
    return id(session) in _handlers


def source_path(structure):
    """The text mmCIF file a structure came from: its saved ``mdcif_path``, else ChimeraX's
    ``filename`` (set by the mmCIF reader, saved in sessions); None for anything else."""
    path = getattr(structure, "mdcif_path", None) or getattr(structure, "filename", None)
    if not isinstance(path, str) or not path.lower().endswith(CIF_SUFFIXES):
        return None
    return path


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


def _in_ihm_model(structure):
    p = structure.parent
    while p is not None:
        if type(p).__name__ == "IHMModel":     # 'open ... format ihm' (chimerax.ihm)
            return True
        p = p.parent
    return False


def attach(session, structures, path, data=None):
    """Attach the mdcif data of ``path`` to ``structures`` (all opened from it).

    Returns the ``MdcifFile``, or None if the file has no mdcif data or the structures can't
    be matched to its models (with a warning). ``data`` (from 'mdcif open') skips the byte
    scan, so any text mmCIF file is accepted.
    """
    from .mdcif_data import MdcifFile

    log = session.logger
    name = os.path.basename(path)
    if data is None:
        if not has_mdcif_categories(path):
            return None
        try:
            data = MdcifFile(path)
        except (ValueError, OSError) as exc:
            log.warning(f"mdcif: {exc}")
            return None
    structures = sorted(structures, key=lambda s: s.id)
    ids = " ".join(f"#{s.id_string}" for s in structures)
    if any(_in_ihm_model(s) for s in structures):
        log.warning(f"mdcif: {ids} opened as IHM; for the mdcif commands open {name} as "
                    "mmCIF ('open FILE' or 'mdcif open FILE')")
        return None
    saved = [getattr(s, "mdcif_model", None) for s in structures]
    if all(m in data.models for m in saved):     # restored session, animate copies
        models, fresh = saved, False
    elif len(structures) == len(data.models):
        models, fresh = data.models, True
    else:
        s0 = structures[0]
        how = (f"one structure with {s0.num_coordsets} coordinate sets ('coordsets true')"
               if len(structures) == 1 and s0.num_coordsets > 1
               else f"{len(structures)} structure(s) ({ids})")
        log.warning(f"mdcif: {name} has {len(data.models)} models but was opened as {how}; "
                    "the mdcif commands need one structure per model ('open FILE')")
        return None
    for s, m in zip(structures, models, strict=True):
        s.mdcif_data, s.mdcif_model, s.mdcif_path = data, m, os.path.abspath(path)
        if fresh and len(structures) > 1:       # submodels: name them after their cluster
            s.name = data.model_label(m)
    if fresh:
        log.info(data.summary())
    else:
        log.info(f"mdcif: data of {name} attached to {ids}")
    return data


def on_add_models(session, models):
    """"add models" handler: attach mdcif data to new structures from mdcif files."""
    from chimerax.atomic import AtomicStructure

    by_path = {}
    for m in models:
        if isinstance(m, AtomicStructure) and getattr(m, "mdcif_data", None) is None:
            path = source_path(m)
            if path is not None:
                by_path.setdefault(path, []).append(m)
    for path, structures in by_path.items():
        attach(session, structures, path, data=pending.pop(path_key(path), None))


def lazy_attach(session, structure):
    """Attach data to a structure the handler missed (opened before the bundle's init ran,
    a session saved without the attributes, ...); returns the data or None."""
    from chimerax.atomic import AtomicStructure

    path = source_path(structure)
    if path is None:
        return None
    group = [structure]
    if getattr(structure, "mdcif_model", None) is None and len(structure.id) > 1:
        group = [m for m in structure.parent.child_models()     # its file's other models
                 if isinstance(m, AtomicStructure) and source_path(m) == path
                 and getattr(m, "mdcif_data", None) is None]
    return attach(session, group, path)
