"""3D views of an mdcif file in a notebook (py3Dmol): models, porcupine plots, mode movies.

:class:`MdcifViewer` reads a file once and builds ``py3Dmol.view`` objects, so everything
py3Dmol does still works on them (``rotate``, ``zoom``, ``show``, ``write_html``)::

    viewer = MdcifViewer("gb1_meld_300K.cif")
    viewer.show()                                # all models, coloured as in the landscape plot
    viewer.porcupine(1).show()                   # model 1, one panel per mode
    viewer.porcupine(1, mode=1, scope="global")  # one panel: global mode 1 drawn on model 1
    viewer.movie(1).show()                       # model 1 swinging along each of its modes
    save_png(viewer.porcupine(1, mode=1), "figures/porcupine.png", width_in=3.5)

A mode is drawn as the spec displays it: atom i moves by k * sigma * v_i, where sigma is the
standard deviation of the motion along the mode (the square root of
``ModeSet.cartesian_variance``) and k is a number of standard deviations (``scale`` for arrows,
``amplitude`` for movies). Atoms without a vector component of their own move with their
residue's mode atom. Mode signs are arbitrary, so every arrow could point the other way too.
py3Dmol is an optional dependency (``viz`` extra); ``save_png`` also needs Pillow and a Chrome.
"""

import math
import os
import re
import shutil
import subprocess
import tempfile
from collections import namedtuple
from pathlib import Path

import numpy as np

from .plot import NOISE, OTHER, SERIES, TEXT, TEXT_2
from .reader import MdcifData, read_mdcif

Atom = namedtuple("Atom", "name resname chain resnum element")
BACKBONE = ("N", "CA", "C", "O")       # enough for a cartoon; keeps the movies small
SCOPES = ("cluster", "global")

# Chrome (or Chromium, or Edge) for save_png, when neither `chrome` nor $CHROME names one
CHROME_NAMES = ("google-chrome", "google-chrome-stable", "chromium", "chromium-browser",
                "chrome", "msedge")
CHROME_PATHS = ("/apps/ubuntu/24.04/bin/google-chrome-stable",          # HiPerGator (UF)
                "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
                "C:/Program Files/Google/Chrome/Application/chrome.exe",
                "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe")


def _py3dmol():
    try:
        import py3Dmol
    except ImportError as e:
        raise ImportError("mdcif.view needs py3Dmol (PyPI: py3Dmol, conda-forge: py3dmol)") from e
    return py3Dmol


def _sigma(ms):
    # standard deviation (angstrom) of the Cartesian motion along each mode
    return np.sqrt(ms.cartesian_variance)


def _point(v):
    return {"x": float(v[0]), "y": float(v[1]), "z": float(v[2])}


def _label(view, text, panel=None, y=8, color=TEXT):
    # caption in the top left corner of a view (or of one panel of a grid)
    where = {} if panel is None else {"viewer": panel}
    view.addLabel(text, {"useScreen": True, "position": {"x": 8, "y": y, "z": 0},
                         "alignment": "topLeft", "inFront": True, "fontSize": 13,
                         "fontColor": color, "backgroundColor": "white",
                         "backgroundOpacity": 0.8}, **where)


def pdb_text(atoms, xyz, model=None):
    """PDB ATOM records for ``Atom``s and their ``(n, 3)`` coordinates (for py3Dmol).

    With ``model``, the records are wrapped in MODEL/ENDMDL, so several can be joined into
    the frames of a movie.
    """
    lines = [f"MODEL {model:8d}\n"] if model is not None else []
    for i, (atom, (x, y, z)) in enumerate(zip(atoms, xyz, strict=True), 1):
        name, resname, chain, resnum, element = atom
        name = name if len(name) == 4 else " " + name
        lines.append(f"ATOM  {i % 100000:5d} {name:<4s} {resname:3s} {chain[:1]:1s}"
                     f"{int(resnum):4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          "
                     f"{element:>2s}\n")
    return "".join(lines) + ("ENDMDL\n" if model is not None else "END\n")


class MdcifViewer:
    """3D views (py3Dmol) of the models and modes of an mdcif file.

    ``source`` is a path (text mmCIF or BinaryCIF) or the ``MdcifData`` that
    :func:`~mdcif.reader.read_mdcif` returned. The drawing methods (``models``,
    ``porcupine``, ``movie``) each return a new ``py3Dmol.view``: call ``.show()`` on it.
    ``scope="cluster"`` uses a model's own modes, ``scope="global"`` the global ones (the
    PCs the frames were clustered on), drawn on that model. Modes are numbered from 1,
    largest motion first.
    """

    def __init__(self, source):
        self.data = source if isinstance(source, MdcifData) else read_mdcif(source)
        self._atoms, self._keys, self._xyz = {}, {}, {}
        for _, model in self.data.system._all_models():
            m = int(model._id)
            atoms, keys, xyz = [], [], []
            for a in model.get_atoms():
                res = a.asym_unit.residue(a.seq_id)
                atoms.append(Atom(a.atom_id, res.comp.id, a.asym_unit.strand_id,
                                  res.auth_seq_id, a.type_symbol))
                keys.append((a.asym_unit._id, a.seq_id, a.atom_id))   # matches across models
                xyz.append((a.x, a.y, a.z))
            self._atoms[m], self._keys[m], self._xyz[m] = atoms, keys, np.array(xyz, float)
        self.model_ids = sorted(self._atoms)

    def __repr__(self):
        own = sorted(self.data.mode_sets)
        glob = self.data.global_modes
        return (f"<MdcifViewer: models {self.model_ids}, own modes on models {own}, "
                f"{glob.n_modes if glob else 'no'} global modes>")

    # ------------------------------------------------------------------ data

    def atoms(self, model):
        """The atoms of ``model``, as ``Atom(name, resname, chain, resnum, element)``."""
        return list(self._atoms[self._model(model)])

    def coordinates(self, model):
        """``(n_atoms, 3)`` coordinates of ``model`` (Å)."""
        return self._xyz[self._model(model)].copy()

    def mode_set(self, model, scope="cluster"):
        """``model``'s own ``ModeSet`` (``scope="cluster"``) or the global one."""
        model = self._model(model)
        if scope not in SCOPES:
            raise ValueError(f"scope must be 'cluster' or 'global', not {scope!r}")
        if scope == "global":
            if self.data.global_modes is None:
                raise ValueError("the file has no global modes")
            return self.data.global_modes
        if model not in self.data.mode_sets:
            raise ValueError(f"model {model} has no modes of its own")
        return self.data.mode_sets[model]

    def color(self, model, scope="cluster"):
        """``model``'s colour: as in the landscape plot (``mdcif.plot``); grey for global."""
        if scope == "global":
            return TEXT_2
        i = self.model_ids.index(self._model(model))
        return SERIES[i] if i < len(SERIES) else OTHER

    def caption(self, model, mode, scope="cluster"):
        """Panel caption: the mode's σ (Å) and its share of the motion."""
        ms = self.mode_set(model, scope)
        k = self._modes(ms, mode)[0]
        which = "global mode" if scope == "global" else "mode"
        text = f"{which} {k}: σ = {_sigma(ms)[k - 1]:.1f} Å"
        share = ms.variance_fraction[k - 1]
        return text if np.isnan(share) else f"{text}, {100 * share:.0f}%"

    def directions(self, model, mode, scope="cluster"):
        """Mode ``mode`` as a unit-vector component for every atom of ``model``.

        Returns ``(d, own)``: ``(n_atoms, 3)`` components (atoms without their own take their
        residue's mode atom's) and a boolean mask of the atoms with their own component.
        """
        model, ms = self._model(model), self.mode_set(model, scope)
        k = self._modes(ms, mode)[0]
        stored = self._keys[ms.model_id]     # the modes are stored on atoms of ms.model_id
        exact, by_residue = {}, {}
        for i, v in zip(ms.atom_indices, ms.eigenvectors[k - 1], strict=True):
            exact[stored[i]] = v
            by_residue.setdefault(stored[i][:2], v)
        keys = self._keys[model]
        d = np.zeros((len(keys), 3))
        own = np.zeros(len(keys), dtype=bool)
        for i, key in enumerate(keys):
            if key in exact:
                d[i], own[i] = exact[key], True
            elif key[:2] in by_residue:
                d[i] = by_residue[key[:2]]
        if not own.any():
            raise ValueError(f"no atom of model {model} matches the atoms of the modes")
        return d, own

    def arrows(self, model, mode, scope="cluster", scale=3.0, min_length=0.5):
        """Porcupine arrows ``(starts, ends)``, ``(n, 3)`` each, on the atoms with a component.

        Each arrow is ``scale`` σ times the atom's share of the unit vector; arrows shorter
        than ``min_length`` Å are left out.
        """
        model, ms = self._model(model), self.mode_set(model, scope)
        k = self._modes(ms, mode)[0]
        d, own = self.directions(model, k, scope)
        start = self._xyz[model][own]
        vec = scale * _sigma(ms)[k - 1] * d[own]
        keep = np.linalg.norm(vec, axis=1) >= min_length
        return start[keep], start[keep] + vec[keep]

    def mode_frames(self, model, mode, scope="cluster", amplitude=3.0, n_frames=30):
        """``(n_frames, n_atoms, 3)``: ``model`` swinging along ``mode`` to ±``amplitude`` σ.

        Frame t is x + amplitude σ sin(2π t / n_frames) d: out one way, back, out the other.
        """
        model, ms = self._model(model), self.mode_set(model, scope)
        k = self._modes(ms, mode)[0]
        d, _ = self.directions(model, k, scope)
        swing = amplitude * _sigma(ms)[k - 1] * np.sin(2 * np.pi * np.arange(n_frames) / n_frames)
        return self._xyz[model] + swing[:, None, None] * d

    # ------------------------------------------------------------------ views

    def models(self, *models, width=640, height=360, legend=True):
        """All models (or the given ones) together as cartoons, one colour each.

        The models of a file share one frame (that of the reference entry, if it was built
        with one), so they can be compared as they are.
        """
        models = [self._model(m) for m in models] or self.model_ids
        view = _py3dmol().view(width=width, height=height)
        for i, m in enumerate(models):
            view.addModel(pdb_text(self._atoms[m], self._xyz[m]), "pdb")
            view.setStyle({"model": i}, {"cartoon": {"color": self.color(m)}})
            if legend:
                _label(view, f"model {m}", y=8 + 24 * i, color=self.color(m))
        view.zoomTo()
        return view

    def show(self, *models, **kwargs):
        """Show all models together, or only the given ones (``models(...).show()``)."""
        return self.models(*models, **kwargs).show()

    def porcupine(self, model=1, mode=None, scope="cluster", scale=3.0, min_length=0.5,
                  color=None, width=None, height=None):
        """Porcupine plot: ``model`` as a grey cartoon, one arrow per mode atom.

        ``mode``: a mode number (one panel), a list of them, or None for all stored modes;
        panels rotate together. Arrows as in :meth:`arrows`, in ``color`` (default: the
        model's colour).
        """
        model, ms = self._model(model), self.mode_set(model, scope)
        modes = self._modes(ms, mode)
        color = color or self.color(model, scope)
        pdb = pdb_text(self._atoms[model], self._xyz[model])
        view = self._grid(len(modes), width, height)
        for p, k in enumerate(modes):
            panel = (0, p)
            view.addModel(pdb, "pdb", viewer=panel)
            for start, end in zip(*self.arrows(model, k, scope, scale, min_length),
                                  strict=True):
                view.addArrow({"start": _point(start), "end": _point(end), "radius": 0.2,
                               "radiusRatio": 2.0, "mid": 0.75, "color": color}, viewer=panel)
            _label(view, self.caption(model, k, scope), panel)
        view.setStyle({"cartoon": {"color": NOISE}})
        view.zoomTo()
        view.zoom(0.85)                    # room for arrows that reach past the protein
        return view

    def movie(self, model=1, mode=None, scope="cluster", amplitude=3.0, n_frames=16,
              interval=100, backbone_only=True, color=None, width=None, height=None):
        """Animation: ``model`` swinging along each mode out to ±``amplitude`` σ and back.

        ``mode`` as in :meth:`porcupine`; ``interval`` ms per frame. With ``backbone_only``
        only N, CA, C and O are drawn (other atoms make the page larger, not the cartoon
        better). The modes are straight lines, so a large swing can stretch bonds.
        """
        model, ms = self._model(model), self.mode_set(model, scope)
        modes = self._modes(ms, mode)
        atoms = self._atoms[model]
        keep = np.array([a.name in BACKBONE for a in atoms])
        if not backbone_only or not keep.any():
            keep[:] = True
        drawn = [a for a, k in zip(atoms, keep, strict=True) if k]
        view = self._grid(len(modes), width, height)
        for p, k in enumerate(modes):
            frames = self.mode_frames(model, k, scope, amplitude, n_frames)[:, keep]
            pdb = "".join(pdb_text(drawn, xyz, model=i + 1) for i, xyz in enumerate(frames))
            # assignBonds off: residues are joined by number, not bond length, so the ribbon
            # stays whole where the straight-line motion stretches the chain
            view.addModelsAsFrames(pdb, "pdb", {"assignBonds": False}, viewer=(0, p))
            _label(view, self.caption(model, k, scope), (0, p))
        view.setStyle({"cartoon": {"color": color or self.color(model, scope)}})
        view.animate({"loop": "forward", "interval": interval})
        view.zoomTo()
        view.zoom(0.8)                     # room for the far ends of the swing
        return view

    # ------------------------------------------------------------------ helpers

    def _model(self, model):
        if int(model) not in self._atoms:
            raise ValueError(f"no model {model}; the file has models {self.model_ids}")
        return int(model)

    @staticmethod
    def _modes(ms, mode):
        # mode numbers to draw: one, several, or None for all
        if mode is None:
            return list(range(1, ms.n_modes + 1))
        modes = [int(mode)] if np.ndim(mode) == 0 else [int(k) for k in mode]
        for k in modes:
            if not 1 <= k <= ms.n_modes:
                raise ValueError(f"mode {k} not in 1..{ms.n_modes}")
        return modes

    @staticmethod
    def _grid(n, width, height):
        # n panels side by side that rotate together
        width = width or (450 if n == 1 else 300 * n)
        height = height or (400 if n == 1 else 330)
        return _py3dmol().view(width=width, height=height, viewergrid=(1, n), linked=True)


def find_chrome(chrome=None):
    """Path of a Chrome, Chromium or Edge: ``chrome``, $CHROME, PATH, usual places; or None."""
    for c in (chrome, os.environ.get("CHROME"), *map(shutil.which, CHROME_NAMES),
              *CHROME_PATHS):
        if c and Path(c).exists():
            return str(c)
    return None


def save_png(view, path, width_in=7.0, dpi=600, chrome=None):
    """Save a py3Dmol view as a PNG ``width_in`` inches wide at ``dpi``; returns the path.

    py3Dmol draws in the browser, so Python never sees the picture: the view is opened as a
    web page in a headless Chrome (see :func:`find_chrome`), which saves a screenshot.
    3Dmol.js comes from the internet (jsdelivr). The page is rendered at a whole-number zoom
    and then downsized, since a fractional zoom leaves seams between the panels of a grid.
    """
    from PIL import Image

    program = find_chrome(chrome)
    if program is None:
        raise RuntimeError("no Chrome, Chromium or Edge found: pass chrome= or set CHROME")
    path = Path(path)
    with tempfile.TemporaryDirectory() as tmp:
        page, shot = Path(tmp) / "view.html", Path(tmp) / "view.png"
        view.write_html(str(page))
        w, h = (int(v) for v in re.search(r"width: (\d+)px; height: (\d+)px",
                                          page.read_text()).groups())
        zoom = math.ceil(width_in * dpi / w)
        subprocess.run([program, "--headless=new", "--no-sandbox", "--hide-scrollbars",
                        "--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader",
                        f"--user-data-dir={tmp}", f"--window-size={w},{h}",
                        f"--force-device-scale-factor={zoom}", "--virtual-time-budget=8000",
                        f"--screenshot={shot}", page.as_uri()],
                       capture_output=True, timeout=300)
        if not shot.exists():
            raise RuntimeError(f"{program} did not save a screenshot")
        image = Image.open(shot).convert("RGB")
    image = image.resize((round(width_in * dpi), round(h * width_in * dpi / w)),
                         Image.Resampling.LANCZOS)
    if image.getextrema() == ((255, 255),) * 3:
        raise RuntimeError("the saved view is blank: is the internet reachable (3Dmol.js)?")
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path, dpi=(dpi, dpi))
    return path
