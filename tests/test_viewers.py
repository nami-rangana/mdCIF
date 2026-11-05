"""Viewer-side code: the shared mdcif_data module, the PyMOL plugin, the ChimeraX commands.

The viewers must not import the mdcif package, so they are loaded from their own
directories. The ChimeraX commands run here against minimal stand-ins for the
chimerax modules (ChimeraX itself is not installable with pip): this checks the
command logic, not the ChimeraX API.
"""

import csv
import importlib.util
import os
import sys
import types
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYMOL_DIR = ROOT / "viewers" / "pymol"
CHIMERAX_SRC = ROOT / "viewers" / "chimerax" / "src"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def vdata():
    return _load("viewer_mdcif_data", PYMOL_DIR / "mdcif_pymol" / "mdcif_data.py")


def test_viewer_copies_are_identical():
    a = (PYMOL_DIR / "mdcif_pymol" / "mdcif_data.py").read_bytes()
    b = (CHIMERAX_SRC / "mdcif_data.py").read_bytes()
    assert a == b, "viewers/*/mdcif_data.py differ: copy one over the other"


def test_tokenizer_handles_quotes_and_text_fields(vdata, tmp_path):
    path = tmp_path / "t.cif"
    path.write_text(
        "data_x\n_mdcif_mode_set.id 1\n_mdcif_mode_set.details\n;two\nlines\n;\n"
        "loop_\n_mdcif_mode.id\n_mdcif_mode.method\n_mdcif_mode.note\n"
        "1 'it''s fine' .\n2 \"a b\" ?   # comment\n3 x 'q'\n"
        "_other.item value\n")
    t = vdata.read_tables(path, ["_mdcif_mode_set", "_mdcif_mode"])
    assert t["_mdcif_mode_set"] == {"id": ["1"], "details": ["two\nlines"]}
    assert t["_mdcif_mode"]["method"] == ["it''s fine", "a b", "x"]
    assert t["_mdcif_mode"]["note"] == [None, None, "q"]
    assert "_other" not in t


def test_viewer_reader_matches_core_reader(vdata, adk_mdcif):
    from mdcif.reader import read_mdcif

    path, summary = adk_mdcif
    core = read_mdcif(path)
    f = vdata.MdcifFile(path)
    assert f.models == list(range(1, summary.n_models + 1))
    g, cg = f.global_set, core.global_modes
    assert np.allclose(g.eigenvalues, cg.eigenvalues, rtol=1e-5)
    assert np.allclose(g.variance_fraction, cg.variance_fraction, rtol=1e-5)
    assert g.scope == "global" and all(f.cluster_set(m).scope == "cluster" for m in core.mode_sets)
    assert np.allclose(g.vectors, cg.eigenvectors, atol=1e-5)
    assert np.allclose(g.origin, cg.origin, atol=1e-3)
    for m, ms in core.mode_sets.items():
        vs = f.cluster_set(m)
        assert np.allclose(vs.vectors, ms.eigenvectors, atol=1e-5)
        assert len(vs.keys) == len(ms.atom_indices) and {k[2] for k in vs.keys} == {"CA"}
    assert np.array_equal(f.projection["model"], core.projection["model_id"])
    assert np.allclose(f.projection["coords"], core.projection["coords"], atol=1e-3)
    assert "global set" in f.summary()


def test_displacements_mode_frames_and_arrows(vdata, adk_mdcif):
    f = vdata.MdcifFile(adk_mdcif[0])
    ms = f.cluster_set(1)
    keys = f.model_keys[1]
    d, own = vdata.displacements(keys, ms, 1)
    assert own.sum() == len(ms.keys)
    # side-chain atoms move with their residue's CA
    res = {}
    for k, v in zip(keys, d, strict=True):
        res.setdefault(k[:2], []).append(v)
    assert all(np.allclose(vs, vs[0]) for vs in res.values())
    xyz = np.random.default_rng(0).normal(size=(len(keys), 3))
    frames = vdata.mode_frames(xyz, keys, ms, 1, amplitude=2.0, n_frames=20)
    assert frames.shape == (20, len(keys), 3) and np.allclose(frames[0], xyz)
    peak = np.abs(frames - xyz).max()
    assert ms.sigma(1) == pytest.approx(np.sqrt(ms.eigenvalues[0]))   # PCA: sqrt(lambda)
    assert peak == pytest.approx(2.0 * ms.sigma(1) * np.abs(d).max(), rel=0.05)
    starts, ends = vdata.arrows(xyz, keys, ms, 1, scale=3.0, min_length=0.0)
    assert len(starts) == len(ms.keys)
    lengths = np.linalg.norm(ends - starts, axis=1)
    assert np.sqrt((lengths ** 2).sum()) == pytest.approx(3.0 * ms.sigma(1))
    va, na, ta = vdata.arrow_mesh(starts[:5], ends[:5], sides=8)
    assert va.shape == na.shape and ta.max() < len(va)
    assert np.allclose(np.linalg.norm(na, axis=1), 1.0, atol=1e-5)
    with pytest.raises(ValueError):
        vdata.displacements(keys, ms, ms.n_modes + 1)



def test_mass_weighted_modes_are_drawn_in_angstrom(vdata, adk_mdcif_mw):
    """Quasi-harmonic sets: arrows in A from the eigenvalue and the elements' masses."""
    from mdcif.reader import read_mdcif

    f = vdata.MdcifFile(adk_mdcif_mw[0])
    ms = f.cluster_set(1)
    assert ms.mass_weighted and not f.global_set.mass_weighted
    core = read_mdcif(adk_mdcif_mw[0]).mode_sets[1]
    assert np.allclose(ms.cartesian_variance, core.cartesian_variance, rtol=1e-6)
    assert "amu A^2" in ms.summary()
    keys = f.model_keys[1]
    xyz = np.zeros((len(keys), 3))
    starts, ends = vdata.arrows(xyz, keys, ms, 1, scale=3.0, min_length=0.0)
    rms = np.sqrt((np.linalg.norm(ends - starts, axis=1) ** 2).sum())
    assert rms == pytest.approx(3.0 * np.sqrt(ms.cartesian_variance[0]), rel=1e-5)
    assert rms == pytest.approx(3.0 * np.sqrt(ms.eigenvalues[0] / 12.011), rel=1e-5)  # CA only
    ms.cartesian_variance[:] = np.nan          # as for an element without a standard mass
    with pytest.raises(ValueError, match="standard mass"):
        vdata.arrows(xyz, keys, ms, 1)


def test_morph_frames(vdata):
    a, b = np.zeros((4, 3)), np.ones((4, 3))
    frames = vdata.morph_frames([a, b], steps=4)
    assert frames.shape == (8, 4, 3) and np.allclose(frames[2], 0.5)
    assert vdata.morph_frames([a, b], steps=4, wrap=False).shape == (5, 4, 3)


def test_measure(vdata):
    distances = pytest.importorskip("MDAnalysis.lib.distances")
    assert vdata.measure([[[0, 0, 0], [3, 4, 0]]], "distance") == pytest.approx([5.0])
    assert vdata.measure([[[1, 0, 0], [0, 0, 0], [0, 2, 0]]], "angle") == pytest.approx([90.0])
    assert vdata.measure([[[1, 0, 0], [0, 0, 0], [-1, 0, 0]]], "angle") == pytest.approx([180.0])
    p = np.random.default_rng(2).normal(size=(50, 4, 3))
    ref = np.degrees(distances.calc_dihedrals(*(p[:, i] for i in range(4))))   # IUPAC sign
    one_by_one = np.array([vdata.measure(f[None], "dihedral")[0] for f in p])
    np.testing.assert_allclose(one_by_one, ref, atol=1e-3)        # MDAnalysis: float32
    ref_angle = np.degrees(distances.calc_angles(*(p[:, i] for i in range(3))))
    np.testing.assert_allclose(vdata.measure(p[:, :3], "angle"), ref_angle, atol=1e-3)
    # a dihedral swinging across 180: unwrapped along the frames, frame 1 in (-180, 180]
    t = np.radians(np.arange(150, 260, 10))
    sweep = np.zeros((len(t), 4, 3))
    sweep[:, 0], sweep[:, 2] = (1, 0, 0), (0, 0, 1)
    sweep[:, 3] = np.column_stack([np.cos(t), np.sin(t), np.ones_like(t)])
    d = vdata.measure(sweep, "dihedral")
    assert d[0] == pytest.approx(150) and np.allclose(np.diff(d), 10)    # 150 ... 250
    wrapped = np.degrees(distances.calc_dihedrals(*sweep.transpose(1, 0, 2)))   # -180 ... 180
    assert np.abs(np.diff(wrapped)).max() > 300
    np.testing.assert_allclose(180 - (180 - d) % 360, wrapped, atol=1e-3)
    flat = [[[1, 0, 0], [0, 0, 0], [0, 0, 1], [-1, 0, 1]]]       # trans: +180, never -180
    assert vdata.measure(flat, "dihedral")[0] == 180.0
    with pytest.raises(ValueError, match="unknown"):
        vdata.measure(p, "torsion")
    with pytest.raises(ValueError, match="needs points"):
        vdata.measure(p[:, :3], "dihedral")


# --------------------------------------------------------------------------- PyMOL

@pytest.fixture(scope="module")
def pymol_plugin():
    pytest.importorskip("pymol")
    sys.path.insert(0, str(PYMOL_DIR))
    try:
        import mdcif_pymol
    finally:
        sys.path.remove(str(PYMOL_DIR))
    from pymol import cmd

    mdcif_pymol.__init_plugin__()
    yield mdcif_pymol, cmd
    cmd.reinitialize()


def test_pymol_plugin_commands(pymol_plugin, adk_mdcif):
    mp, cmd = pymol_plugin
    path, summary = adk_mdcif
    name = mp.mdcif_load(str(path), "adk", quiet=1)
    assert cmd.count_states(name) == summary.n_models
    morph = mp.mdcif_morph(name, steps=3)
    assert cmd.count_states(morph) == 3 * summary.n_models
    anim = mp.mdcif_animate(name, model=2, mode=1, amplitude=2.0, frames=8)
    assert cmd.count_states(anim) == 8
    moved = cmd.get_coords(anim, state=3) - cmd.get_coords(name, state=2)
    assert np.abs(moved).max() > 0.1
    g = mp.mdcif_animate(name, model=1, mode=2, frames=6, scope="global")
    assert g.endswith("_global2") and cmd.count_states(g) == 6
    porc = mp.mdcif_porcupine(name, model=1, mode=1)
    assert porc in cmd.get_names("objects")
    with pytest.raises(ValueError):
        mp.mdcif_animate(name, scope="sideways")


# --------------------------------------------------------------------------- ChimeraX

class _Residues:
    def __init__(self, keys):
        self.chain_ids = [k[0] for k in keys]
        self.numbers = [int(k[1]) for k in keys]


class _Atoms:
    def __init__(self, keys, coords):
        self.names = [k[2] for k in keys]
        self.residues = _Residues(keys)
        self.coords = coords
        self.coord_indices = np.arange(len(keys))

    def __len__(self):
        return len(self.names)


class _Picked:
    """What AtomsArg gives a command: atoms in the order typed and their structures."""

    def __init__(self, *pairs):                 # (structure, key) pairs
        keys = [k for _, k in pairs]
        self.names = [k[2] for k in keys]
        self.residues = _Residues(keys)
        self.unique_structures = list(dict.fromkeys(s for s, _ in pairs))

    def __len__(self):
        return len(self.names)


def _pick(structure, *indices):
    return _Picked(*[(structure, structure.keys[i]) for i in indices])


class _Model:
    _next = [1]

    def __init__(self, name=""):
        self.name, self.children, self.display, self.parent = name, [], True, None
        self.deleted = False
        self.id = (_Model._next[0],)
        self.scene_position = object()          # stand-in for a Place
        _Model._next[0] += 1

    @property
    def id_string(self):
        return ".".join(map(str, self.id))

    @property
    def visible(self):
        return self.display

    def all_models(self):
        out = [self]
        for c in self.children:
            out.extend(c.all_models())
        return out

    def child_models(self):
        return list(self.children)

    def add(self, models):
        for k, m in enumerate(models, start=len(self.children) + 1):
            m.parent, m.id = self, self.id + (k,)
        self.children.extend(models)


class _Structure(_Model):
    registered = []                             # AtomicStructure.register_attr calls

    def __init__(self, keys, coords, name="s", filename=None):
        super().__init__(name)
        self.keys, self.atoms, self.coordsets = keys, _Atoms(keys, coords), None
        self.filename, self.active_coordset_id = filename, 1

    @classmethod
    def register_attr(cls, session, name, registrant, attr_type=None):
        cls.registered.append((name, registrant, attr_type))

    def copy(self, name):
        c = _Structure(self.keys, self.atoms.coords.copy(), name)
        c.atoms.coord_indices = self.atoms.coord_indices.copy()
        return c

    def add_coordsets(self, xyz, replace=False):
        assert xyz.shape[1] == len(self.atoms)              # as ChimeraX checks
        self.coordsets = xyz

    def remove_coordsets(self):
        self.coordsets = None

    def add_coordset(self, cs_id, xyz):
        assert cs_id == self.num_coordsets + (self.coordsets is not None)
        self.coordsets = xyz[None] if self.coordsets is None else np.concatenate(
            [self.coordsets, xyz[None]])

    @property
    def num_coordsets(self):
        return 1 if self.coordsets is None else len(self.coordsets)

    @property
    def coordset_ids(self):
        return list(range(1, self.num_coordsets + 1))

    def coordset(self, cs_id):
        xyz = self.atoms.coords if self.coordsets is None else self.coordsets[cs_id - 1]
        return types.SimpleNamespace(xyzs=xyz)


class _Surface(_Model):
    SESSION_SAVE_DRAWING = False                # as ChimeraX's Model: geometry not in sessions

    def __init__(self, name, session):
        super().__init__(name)

    def set_geometry(self, va, na, ta):
        self.geometry = (va, na, ta)


@pytest.fixture()
def chimerax_cmd(monkeypatch):
    """Import viewers/chimerax/src/cmd.py against stand-in chimerax modules."""
    ran = []

    class UserError(Exception):
        pass

    commands = types.ModuleType("chimerax.core.commands")
    for n in ("BoolArg", "ColorArg", "FloatArg", "IntArg", "OpenFileNameArg", "SaveFileNameArg"):
        setattr(commands, n, n)
    commands.log_equivalent_command = lambda session, text: session.logger.commands.append(text)
    commands.CmdDesc = lambda **kw: kw
    commands.EnumOf = lambda values: ("EnumOf", tuple(values))
    commands.FileNameArg = types.SimpleNamespace(unparse=lambda path: f"<{path}>")
    commands.StringArg = types.SimpleNamespace(unparse=lambda text: f'"{text}"')
    commands.run = lambda session, text, log=True: (ran.append(text), session.run(text))[1]
    modules = {
        "chimerax": types.ModuleType("chimerax"),
        "chimerax.core": types.ModuleType("chimerax.core"),
        "chimerax.core.commands": commands,
        "chimerax.core.errors": types.ModuleType("chimerax.core.errors"),
        "chimerax.core.models": types.ModuleType("chimerax.core.models"),
        "chimerax.core.toolshed": types.ModuleType("chimerax.core.toolshed"),
        "chimerax.core.triggerset": types.ModuleType("chimerax.core.triggerset"),
        "chimerax.atomic": types.ModuleType("chimerax.atomic"),
    }
    modules["chimerax.core.errors"].UserError = UserError
    modules["chimerax.core.models"].Surface = _Surface
    modules["chimerax.core.models"].ADD_MODELS = "add models"
    toolshed = types.SimpleNamespace(triggers=_Triggers())
    modules["chimerax.core.toolshed"].BundleAPI = object
    modules["chimerax.core.toolshed"].TOOLSHED_BUNDLE_UNINSTALLED = "bundle uninstalled"
    modules["chimerax.core.toolshed"].get_toolshed = lambda: toolshed
    modules["chimerax.core.triggerset"].DEREGISTER = _DEREGISTER
    modules["chimerax.atomic"].AtomicStructure = _Structure
    modules["chimerax.atomic"].AtomicStructureArg = "AtomicStructureArg"
    modules["chimerax.atomic"].AtomicStructuresArg = "AtomicStructuresArg"
    modules["chimerax.atomic"].AtomsArg = "AtomsArg"
    for k, v in modules.items():
        monkeypatch.setitem(sys.modules, k, v)
    pkg = types.ModuleType("mdcif_chimerax")
    pkg.__path__ = [str(CHIMERAX_SRC)]
    monkeypatch.setitem(sys.modules, "mdcif_chimerax", pkg)
    for sub in ("__init__", "mdcif_data", "attach", "plot", "cmd", "menu"):
        name = "mdcif_chimerax" if sub == "__init__" else f"mdcif_chimerax.{sub}"
        spec = importlib.util.spec_from_file_location(name, CHIMERAX_SRC / f"{sub}.py")
        module = importlib.util.module_from_spec(spec)
        monkeypatch.setitem(sys.modules, name, module)
        spec.loader.exec_module(module)
        if sub == "__init__":
            pkg = module
    return sys.modules["mdcif_chimerax.cmd"], ran, UserError


class _Logger:
    def __init__(self):
        self.lines, self.commands = [], []      # commands: log_equivalent_command

    def info(self, text, **kw):
        self.lines.append(text)

    warning = info


_DEREGISTER = "delete handler"                  # chimerax.core.triggerset.DEREGISTER


class _Triggers:
    def __init__(self):
        self.handlers = {}

    def add_handler(self, name, func):
        handler = (name, func)
        self.handlers.setdefault(name, []).append(handler)
        return handler

    def remove_handler(self, handler):
        self.handlers[handler[0]].remove(handler)

    def activate_trigger(self, name, data):
        for handler in list(self.handlers.get(name, [])):
            if handler[1](name, data) == _DEREGISTER:
                self.remove_handler(handler)


class _Session:
    """Stand-in session: ``open`` makes one structure per model of the file (random
    coordinates) under a group model and fires "add models", as ChimeraX's mmCIF reader does."""

    def __init__(self, vdata, seed=1, bundle_init=True):
        self.logger, self.models_added, self.triggers = _Logger(), [], _Triggers()
        self.models = types.SimpleNamespace(add=self.add_models, list=self.list_models,
                                            close=self.close_models)
        self.ui = types.SimpleNamespace(is_gui=False)
        self._vdata, self._rng = vdata, np.random.default_rng(seed)
        if bundle_init:
            sys.modules["mdcif_chimerax"].bundle_api.initialize(self, None)

    def add_models(self, models):
        self.models_added.extend(models)
        self.triggers.activate_trigger("add models", [a for m in models for a in m.all_models()])

    def list_models(self, type=None):
        found = [a for m in self.models_added for a in m.all_models() if not a.deleted]
        return [m for m in found if type is None or isinstance(m, type)]

    def close_models(self, models):
        for m in models:
            for a in m.all_models():
                a.deleted = True

    def structures_for(self, path, n=None):
        f = self._vdata.MdcifFile(path)
        models = f.models if n is None else f.models[:n]
        return [_Structure(f.model_keys[m], self._rng.normal(size=(len(f.model_keys[m]), 3)),
                           name=Path(path).name, filename=str(path)) for m in models]

    def open_file(self, path, n=None):
        """Plain 'open FILE'; ``n`` structures instead of one per model."""
        structures = self.structures_for(path, n)
        group = _Model(f"{Path(path).name} group")
        group.add(structures)
        self.add_models([group])
        return structures

    def run(self, text):
        if text.startswith("open <"):
            return self.open_file(text[len("open <"):text.index(">")])
        return None


def _open_standins(cmd, vdata, path, seed=1, bundle_init=True):
    """``mdcif open`` on stand-in structures."""
    session = _Session(vdata, seed, bundle_init)
    return vdata.MdcifFile(path), session, cmd.mdcif_open(session, str(path))


def _summaries(session, path):
    return [line for line in session.logger.lines if line.startswith(f"{path}: ")]


def test_chimerax_commands(chimerax_cmd, vdata, adk_mdcif):
    cmd, ran, UserError = chimerax_cmd
    path = adk_mdcif[0]
    f, session, opened = _open_standins(cmd, vdata, path)
    assert [s.mdcif_model for s in opened] == f.models
    assert [s.name for s in opened] == [f.model_label(m) for m in f.models]
    assert {s.mdcif_path for s in opened} == {str(path)}
    assert ran[0] == f"open <{path}> format mmcif"          # path quoted by FileNameArg
    assert len(_summaries(session, path)) == 1               # attached and logged once
    assert ("mdcif_model", "MDCIF", int) in _Structure.registered
    cmd.mdcif_modes(session, opened)
    assert len(_summaries(session, path)) == 2
    cmd.mdcif(session)                                       # bare 'mdcif': status, the menu
    assert session.logger.lines[-1].endswith(
        "Windows: menu MDCIF > Porcupine, Animate, Analyze (when ChimeraX has its window).")

    n_added = len(session.models_added)
    copy = cmd.mdcif_animate(session, opened[0], mode=1, amplitude=-2.0, frames=12)
    assert copy.coordsets.shape == (12, len(f.model_keys[1]), 3)
    assert session.models_added[n_added:] == [copy] and opened[0].display is False
    assert ran[-1] == f"coordset #{copy.id_string} 1,12 loop 1000"
    assert session.logger.lines[-1].endswith(
        f"+/-2 sd ('coordset stop #{copy.id_string}' to stop)")
    assert copy.position is opened[0].scene_position
    assert (copy.mdcif_data, copy.mdcif_model) == (opened[0].mdcif_data, 1)
    assert copy.mdcif_copy is True and copy.mdcif_path == str(path)
    assert (copy.mdcif_mode, copy.mdcif_scope, copy.mdcif_amplitude) == (1, "cluster", -2.0)
    assert ("mdcif_amplitude", "MDCIF", float) in _Structure.registered
    assert len(_summaries(session, path)) == 2               # the handler left the copy alone
    copy2 = cmd.mdcif_animate(session, opened[1], mode=1, frames=6, scope="global")
    assert copy2.coordsets.shape[0] == 6 and copy2.name != copy.name
    assert copy2.name == f"{path.name} model 2 global mode 1"

    surf = cmd.mdcif_porcupine(session, opened[0], mode=1, scale=3.0, min_length=0.0)
    va, na, ta = surf.geometry
    assert len(va) > 0 and surf in opened[0].children
    assert "is hidden" in session.logger.lines[-1]           # animate hid opened[0]
    copy.atoms.coords = copy.coordsets[3]                    # the copy mid-swing
    on_copy = cmd.mdcif_porcupine(session, copy, mode=1, scale=3.0, min_length=0.0)
    np.testing.assert_array_equal(on_copy.geometry[0], va)   # drawn from coordset 1

    cmd.mdcif_morph(session, [copy, *opened[::-1]], frames=5)
    assert ran[-1].startswith("morph #") and "frames 5" in ran[-1]
    assert ran[-1].split(" frames")[0].split()[1:] == [f"#{s.id_string}" for s in opened]
    with pytest.raises(UserError, match="no mdcif data"):
        cmd.mdcif_animate(session, _Structure([], np.zeros((0, 3))), mode=1)


def test_chimerax_plain_open_attaches(chimerax_cmd, vdata, adk_mdcif, tmp_path):
    cmd, ran, UserError = chimerax_cmd
    path = adk_mdcif[0]
    f = vdata.MdcifFile(path)
    session = _Session(vdata)
    opened = session.open_file(str(path))                    # 'open FILE', File > Open, ...
    assert [s.mdcif_model for s in opened] == f.models
    assert [s.name for s in opened] == [f.model_label(m) for m in f.models]
    assert len(_summaries(session, path)) == 1
    copy = cmd.mdcif_animate(session, opened[0], mode=1, frames=4)
    assert copy.mdcif_data is opened[0].mdcif_data and ran == [ran[-1]]

    plain = tmp_path / "plain.cif"                           # an ordinary mmCIF file
    plain.write_text(Path(path).read_text().split("_mdcif_")[0].rsplit("loop_", 1)[0])
    n_lines = len(session.logger.lines)
    others = session.open_file(str(plain))
    assert all(getattr(s, "mdcif_data", None) is None for s in others)
    assert len(session.logger.lines) == n_lines              # silent
    with pytest.raises(UserError, match="no mdcif data"):
        cmd.mdcif_porcupine(session, others[0])
    forced = cmd.mdcif_open(session, str(plain))             # 'mdcif open' still attaches
    assert [s.mdcif_model for s in forced] == f.models
    assert len(_summaries(session, plain)) == 1

    one = session.open_file(str(path), n=1)                  # e.g. 'coordsets true'
    assert getattr(one[0], "mdcif_data", None) is None
    assert "the mdcif commands need one structure per model" in session.logger.lines[-1]


def test_chimerax_lazy_and_restored(chimerax_cmd, vdata, adk_mdcif):
    cmd, _, UserError = chimerax_cmd
    path = adk_mdcif[0]
    f = vdata.MdcifFile(path)
    session = _Session(vdata, bundle_init=False)             # no handler: opened before init
    opened = session.open_file(str(path))
    assert all(getattr(s, "mdcif_data", None) is None for s in opened)
    copy = cmd.mdcif_animate(session, opened[1], mode=1, frames=4)   # attached on first use
    assert [s.mdcif_model for s in opened] == f.models and copy.mdcif_model == 2
    assert len(_summaries(session, path)) == 1
    cmd.mdcif_open(session, str(path))                       # attached without the handler
    assert len(_summaries(session, path)) == 2

    restored = _Session(vdata)                               # a session restore: the saved
    structures = restored.structures_for(path)               # attributes come back, the data
    for s, m in zip(structures, f.models[::-1], strict=True):    # are read again
        s.mdcif_model, s.mdcif_path, s.filename = m, str(path), None
    group = _Model("group")
    group.add(structures)
    restored.add_models([group])
    assert [s.mdcif_model for s in structures] == f.models[::-1]
    assert all(s.mdcif_data is structures[0].mdcif_data for s in structures)
    assert restored.logger.lines[-1].startswith(f"mdcif: data of {path.name} attached to #")
    assert structures[0].name == path.name                   # names are not reset

    gone = restored.structures_for(path, n=1)[0]             # the file moved since the save
    gone.mdcif_model, gone.mdcif_path = 1, str(path.parent / "moved.cif")
    with pytest.raises(UserError, match="is gone"):
        cmd.mdcif_porcupine(restored, gone)


def test_chimerax_animate_after_deleting_atoms(chimerax_cmd, vdata, adk_mdcif):
    cmd, _, _ = chimerax_cmd
    _, session, opened = _open_standins(cmd, vdata, adk_mdcif[0])
    source = opened[0]
    source.atoms.coord_indices = 2 * np.arange(len(source.atoms))   # deleted atoms' slots
    copy = cmd.mdcif_animate(session, source, mode=1, frames=4)
    ci = copy.atoms.coord_indices
    assert copy.coordsets.shape == (4, 2 * len(ci) - 1, 3)
    np.testing.assert_array_equal(copy.coordset(1).xyzs[ci], source.atoms.coords)
    np.testing.assert_array_equal(cmd._rest_coords(copy), source.atoms.coords)
    again = cmd.mdcif_animate(session, copy, mode=2, frames=3)    # a copy of the copy
    np.testing.assert_array_equal(again.coordset(1).xyzs[ci], source.atoms.coords)


def test_chimerax_morph_orders_by_file(chimerax_cmd, vdata, adk_mdcif):
    cmd, ran, _ = chimerax_cmd
    _, session, first = _open_standins(cmd, vdata, adk_mdcif[0], seed=1)
    second = cmd.mdcif_open(session, str(adk_mdcif[0]))
    assert first[0].mdcif_data is not second[0].mdcif_data   # one MdcifFile per open
    mixed = [s for pair in zip(second, first, strict=False) for s in pair]
    cmd.mdcif_morph(session, mixed)
    order = ran[-1].split(" frames")[0].split()[1:]
    assert order == [f"#{s.id_string}" for s in first + second]


def test_chimerax_input_checks(chimerax_cmd, vdata, adk_mdcif, tmp_path):
    cmd, ran, UserError = chimerax_cmd
    _, session, opened = _open_standins(cmd, vdata, adk_mdcif[0])
    n_ran, n_added = len(ran), len(session.models_added)
    with pytest.raises(UserError):
        cmd.mdcif_open(session, str(tmp_path))               # a directory: OSError
    for frames in (0, -5, 2):
        with pytest.raises(UserError, match="frames"):
            cmd.mdcif_animate(session, opened[0], frames=frames)
    with pytest.raises(UserError, match="frames"):
        cmd.mdcif_morph(session, opened, frames=0)
    for kw in ({"radius": 0}, {"radius": -0.1}, {"min_length": -1},
               {"scale": 0, "min_length": 0}, {"min_length": 1e6}):
        with pytest.raises(UserError):
            cmd.mdcif_porcupine(session, opened[0], **kw)
    assert len(ran) == n_ran and len(session.models_added) == n_added   # nothing ran or added
    assert all(not s.children for s in opened)


def test_mdcif_byte_scan(chimerax_cmd, adk_mdcif, tmp_path):
    attach = chimerax_cmd[0].attach
    assert attach.has_mdcif_categories(adk_mdcif[0])
    big = tmp_path / "split.cif"                             # the mark across a block edge
    big.write_bytes(b"x" * 13 + b"_mdcif_mode_set.id 1\n")
    assert attach.has_mdcif_categories(big, chunk=16)
    big.write_bytes(b"x" * 13 + b"_mdcif" + b"\n" * 40)
    assert not attach.has_mdcif_categories(big, chunk=16)
    assert not attach.has_mdcif_categories(tmp_path / "missing.cif")


def _ca(structure):
    return [i for i, k in enumerate(structure.keys) if k[2] == "CA"]


def _distance(structure, i, j):
    xyz = structure.coordsets
    return np.linalg.norm(xyz[:, i] - xyz[:, j], axis=1)


def test_chimerax_plot(chimerax_cmd, vdata, adk_mdcif, tmp_path):
    cmd, ran, UserError = chimerax_cmd
    plot_module = sys.modules["mdcif_chimerax.plot"]
    path = adk_mdcif[0]
    _, session, opened = _open_standins(cmd, vdata, path)
    model, ca = opened[0], _ca(opened[0])
    n_added = len(session.models_added)
    d1 = cmd.mdcif_plot(session, "distance", _pick(model, ca[9], ca[39]))
    copy = session.models_added[-1]                          # animated for the plot (defaults)
    assert len(session.models_added) == n_added + 1 and copy.mdcif_copy and model.display is False
    assert session.logger.commands == [f"mdcif animate #{model.id_string} mode 1 amplitude 3 "
                                       "frames 30 scope cluster"]
    assert ran[-1] == f"coordset #{copy.id_string} 1,30 loop 1000"
    assert d1.shape == (30,) and np.ptp(d1) > 0
    np.testing.assert_allclose(d1, _distance(copy, ca[9], ca[39]))
    [plot] = plot_module.plots(session)
    assert plot.structure is copy and plot.source is model and plot.kinds == ["distance"]
    label = plot.series[0].label
    assert label == f"/4AKE:{model.keys[ca[9]][1]}@CA - /4AKE:{model.keys[ca[39]][1]}@CA"
    assert session.logger.lines[-1].startswith(f"#{copy.id_string} distance {label} along "
                                               "cluster mode 1 of model 1, +/-3 sd, 30 frames")
    assert session.logger.lines[-1].endswith(f"mean {d1.mean():.2f} A")

    # the model again (no settings: any of its plot copies) and the copy itself: same plot
    cmd.mdcif_plot(session, "distance", _pick(model, ca[11], ca[37]))
    cmd.mdcif_plot(session, "distance", _pick(copy, ca[0], ca[50]))
    cmd.mdcif_plot(session, "distance", _pick(copy, ca[0], ca[50]))      # not added twice,
    cmd.mdcif_plot(session, "distance", _pick(copy, ca[50], ca[0]))      # also not reversed
    assert session.logger.lines[-1].endswith("is already plotted")
    assert len(session.models_added) == n_added + 1 and len(plot.series) == 3
    assert [s.color for s in plot.series] == plot_module.COLORS[:3]
    np.testing.assert_allclose(plot.series[2].values, _distance(copy, ca[0], ca[50]))
    # the whole plot (every series) as CSV and as an image
    data, image = tmp_path / "d.csv", tmp_path / "d.png"
    again = cmd.mdcif_plot(session, "distance", _pick(model, ca[9], ca[39]),
                           save_data=str(data), save_image=str(image))
    np.testing.assert_array_equal(again, d1)                 # the first series, not added again
    lines = data.read_text().splitlines()
    version = sys.modules["mdcif_chimerax"].__version__
    assert version == "0.5.0" and lines[0] == f"# mdcif plot (ChimeraX-MDCIF {version})"
    for item in (f"# file: {path}", "# model: 1 (", "# scope: cluster", "# mode: 1",
                 "# amplitude: 3 sd", "# frames: 30", "# kind: distance (A)"):
        assert any(line.startswith(item) for line in lines if line.startswith("#")), item
    rows = list(csv.reader(line for line in lines if not line.startswith("#")))
    assert rows[0] == ["frame"] + [s.label for s in plot.series]
    table = np.array(rows[1:], float)
    assert table.shape == (30, 4) and np.array_equal(table[:, 0], np.arange(1, 31))
    np.testing.assert_allclose(table[:, 1:].T, [s.values for s in plot.series], rtol=1e-5)
    png = image.read_bytes()
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert (int.from_bytes(png[16:20], "big"), int.from_bytes(png[20:24], "big")) == (1800, 1050)
    for suffix, magic in ((".svg", b"<?xml"), (".pdf", b"%PDF")):
        plot_module.save_image(plot, tmp_path / f"d{suffix}")
        assert (tmp_path / f"d{suffix}").read_bytes().startswith(magic)
    assert plot_module.frame_ticks(30) == [1, 5, 10, 15, 20, 25, 30]     # frame 1 labelled
    assert plot_module.frame_ticks(12) == [1, 2, 4, 6, 8, 10, 12]
    assert plot_module.frame_ticks(3) == [1, 2, 3]
    assert [plot_module.frame_minor_step(n) for n in (3, 12, 30, 100, 1000)] == [None, 1, 1, 2,
                                                                                20]
    from matplotlib.figure import Figure

    figure = Figure()
    marker = plot_module.draw(figure, plot, current=4)
    assert len(figure.axes) == 1 and len(marker.artists) == 4    # the line + a dot each
    assert all(a.get_animated() for a in marker.artists)     # left out of full draws: blitted
    from matplotlib.lines import Line2D

    assert all(isinstance(a, Line2D) for a in marker.artists)   # no scatter: 1.8 ms per draw
    assert marker.line is marker.panels["distance"][0] and len(marker.dots) == 3
    assert plot_module.draw(Figure(), plot_module.Plot(copy), current=4) is None

    # other kinds: panels of the same plot (one plot per animation), also picked on the model
    angle = cmd.mdcif_plot(session, "angle", _pick(model, ca[5], ca[6], ca[7]))
    dihedral = cmd.mdcif_plot(session, "dihedral", _pick(copy, *ca[5:9]))
    assert len(session.models_added) == n_added + 1 and plot_module.plots(session) == [plot]
    assert plot.kinds == ["distance", "angle", "dihedral"]
    assert session.logger.lines[-1].startswith(f"#{copy.id_string} dihedral /4AKE:")
    assert session.logger.lines[-1].endswith(" deg")
    assert np.all((angle > 0) & (angle < 180))
    np.testing.assert_allclose(dihedral, vdata.measure(copy.coordsets[:, ca[5:9]], "dihedral"))
    cmd.mdcif_plot(session, "angle", _pick(copy, ca[7], ca[6], ca[5]))   # reversed: the same
    cmd.mdcif_plot(session, "dihedral", _pick(copy, *ca[8:4:-1]))
    assert [len(plot.series_of(k)) for k in plot.kinds] == [3, 1, 1]
    assert [s.color for s in plot.series_of("angle")] == plot_module.COLORS[:1]   # per panel
    both, image3 = tmp_path / "t.csv", tmp_path / "t.png"
    cmd.mdcif_plot(session, "angle", _pick(copy, *ca[5:8]), save_data=str(both),
                   save_image=str(image3))
    assert session.logger.lines[-1] == (f"#{copy.id_string} distance + angle + dihedral plot: "
                                        f"image saved to {image3}")
    text = both.read_text()
    assert "# kind: distance (A), angle (deg), dihedral (deg)\n" in text
    assert "# dihedrals unwrapped from frame 1" in text and "unwrapped" not in data.read_text()
    rows = list(csv.reader(line for line in text.splitlines() if not line.startswith("#")))
    assert rows[0] == ["frame"] + [s.label for s in plot.series]       # 2, 3, 4 atoms
    png = image3.read_bytes()                                # 3 panels: 3.5 + 2 * 2.2 in high
    assert (int.from_bytes(png[16:20], "big"), int.from_bytes(png[20:24], "big")) == (1800, 2370)
    figure = Figure()
    marker = plot_module.draw(figure, plot, current=4)
    assert len(figure.axes) == 3 and len(marker.artists) == 3 + 5
    assert [ax.get_ylabel() for ax in figure.axes] == ["Distance (Å)", "Angle (°)",
                                                       "Dihedral (°)"]
    marker.set(7)
    for kind, (line, _) in marker.panels.items():
        assert list(line.get_xdata()) == [7, 7]
        np.testing.assert_array_equal(marker.points(kind),
                                      [[7, s.values[6]] for s in plot.series_of(kind)])
    marker.set(31)                                           # past the end: no dots
    assert np.isnan(marker.points("angle")[:, 1]).all()
    # new true: another plot of the same copy; other settings: a new animation (of the copy)
    cmd.mdcif_plot(session, "distance", _pick(copy, ca[1], ca[2]), new=True)
    assert len(plot_module.plots(session)) == 2 and len(session.models_added) == n_added + 1
    cmd.mdcif_plot(session, "distance", _pick(copy, ca[1], ca[2]), mode=2, frames=10)
    copy2 = session.models_added[-1]
    assert copy2 is not copy and (copy2.mdcif_mode, copy2.num_coordsets) == (2, 10)
    assert plot_module.plots(session)[-1].structure is copy2

    # bad input: nothing animated, nothing added
    n_added, n_plots = len(session.models_added), len(plot_module.plots(session))
    with pytest.raises(UserError, match="a dihedral needs 4 atoms, got 3"):
        cmd.mdcif_plot(session, "dihedral", _pick(copy, *ca[:3]))
    with pytest.raises(UserError, match="must be in one structure"):
        cmd.mdcif_plot(session, "distance", _Picked((opened[0], opened[0].keys[ca[0]]),
                                                    (opened[1], opened[1].keys[ca[1]])))
    with pytest.raises(UserError, match="must be in one structure"):    # before the count
        cmd.mdcif_plot(session, "distance", _Picked(*[(s, s.keys[ca[0]]) for s in opened[:3]]))
    with pytest.raises(UserError, match="image must be"):
        cmd.mdcif_plot(session, "distance", _pick(model, ca[3], ca[4]), save_image="d.jpg")
    with pytest.raises(UserError, match="no mdcif data"):
        cmd.mdcif_plot(session, "distance",
                       _pick(_Structure([("A", "1", "CA"), ("A", "2", "CA")], np.zeros((2, 3))),
                             0, 1))
    locked = tmp_path / "locked"
    locked.mkdir(mode=0o500)
    targets = [({"save_data": str(tmp_path / "missing" / "d.csv")}, "there is no folder"),
               ({"save_image": str(tmp_path / "missing" / "d.png")}, "there is no folder"),
               ({"save_data": str(tmp_path)}, "that is a folder")]
    if not os.access(locked, os.W_OK):                      # (root can write anywhere)
        targets.append(({"save_image": str(locked / "d.png")}, "not allowed to write there"))
    for kw, message in targets:                              # checked before animating
        with pytest.raises(UserError, match=message):
            cmd.mdcif_plot(session, "distance", _pick(model, ca[3], ca[4]), mode=3, **kw)
    with pytest.raises(UserError, match="is given twice"):
        cmd.add_series(session, plot, [model.keys[ca[3]], model.keys[ca[3]]])
    assert len(session.models_added) == n_added
    assert len(plot_module.plots(session)) == n_plots and len(plot.series) == 5


def test_chimerax_plot_typed_atoms(chimerax_cmd):
    """The window's atom boxes: residue number + atom name, matched to a structure's keys."""
    plot_module = sys.modules["mdcif_chimerax.plot"]
    parse, match = plot_module.parse_atom, plot_module.match_key
    for text in ("38CA", "38 CA", "38@CA", ":38@CA", " 38ca ", "038CA"):
        assert parse(text) == (None, 38, "CA"), text
    for text in ("A:38CA", "/A:38@CA", "A:38 CA", "/A:38ca"):
        assert parse(text) == ("A", 38, "CA"), text
    assert parse("38@1HB") == (None, 38, "1HB") and parse("-2@CA") == (None, -2, "CA")
    assert parse("4AKE:5@OD1") == ("4AKE", 5, "OD1") and parse("12C4'") == (None, 12, "C4'")
    for text in ("abc", "CA38", "38", "38@", "381HB CA", "#2/A:38@CA", "A38CA", "", "  "):
        with pytest.raises(ValueError, match="e.g. 38CA"):
            parse(text)

    keys = [("A", "38", "N"), ("A", "38", "CA"), ("A", "38", "CB"), ("A", "39", "CA"),
            ("B", "39", "CA"), ("b", "40", "CA"), ("A", "41", "Zn")]
    assert match(parse("38CA"), keys) == ("A", "38", "CA")              # one chain: no chain
    assert match(parse("A:38@CB"), keys) == ("A", "38", "CB")
    assert match(parse("B:39CA"), keys) == ("B", "39", "CA")
    assert match(parse("B:40CA"), keys) == ("b", "40", "CA")            # case only: taken
    assert match(parse("41ZN"), keys) == ("A", "41", "Zn")
    with pytest.raises(ValueError, match="39CA is in chains A, B: write A:39CA for chain A"):
        match(parse("39CA"), keys)
    with pytest.raises(ValueError, match=r"no atom /A:38@XX \(residue 38 has N, CA, CB\)"):
        match(parse("38XX"), keys)
    with pytest.raises(ValueError, match="no residue :99"):
        match(parse("99CA"), keys)
    with pytest.raises(ValueError, match="no residue /C:38"):
        match(parse("C:38CA"), keys)

    # '381HB' is residue 381, atom HB: the window then hints at the '@' rule
    hint, rule = plot_module.digit_hint, "; a name that starts with a digit needs @, e.g. 38@1HB"
    keys += [("A", "38", "1HB"), ("A", "9", "CA")]
    assert hint("381HB", parse("381HB"), keys) == ("; did you mean 38@1HB? (a name that starts "
                                                   "with a digit needs @)")
    assert hint("A:381HB", parse("A:381HB"), keys).startswith("; did you mean A:38@1HB?")
    for text in ("382HB", "99CA"):                           # no such residue: the rule
        assert hint(text, parse(text), keys) == rule, text
    for text in ("38@XX", "38XX", "5XX"):                    # '@' typed; residue 38 is there;
        assert hint(text, parse(text), keys) == "", text     # one digit
    no_digits = [k for k in keys if not k[2][:1].isdigit()]  # no such names: no reminder
    assert hint("99CA", parse("99CA"), no_digits) == ""


def test_chimerax_plot_panels(chimerax_cmd, vdata, adk_mdcif):
    """Kinds by atom count, eight colours per panel, figure sizes."""
    cmd, _, _ = chimerax_cmd
    plot_module = sys.modules["mdcif_chimerax.plot"]
    _, session, opened = _open_standins(cmd, vdata, adk_mdcif[0])
    ca = _ca(opened[0])
    copy = cmd.mdcif_animate(session, opened[0], frames=12)
    plot = plot_module.Plot(copy)
    keys = [copy.keys[i] for i in ca]
    for k in range(8):
        plot.add(keys[k:k + 2], np.zeros(12))
    with pytest.raises(ValueError, match="at most 8 distances"):
        plot.add(keys[20:22], np.zeros(12))
    assert plot.add(keys[:4], np.zeros(12)).color == plot_module.COLORS[0]   # own panel
    assert plot.kinds == ["distance", "dihedral"] and plot.series[-1].kind == "dihedral"
    with pytest.raises(ValueError, match="needs 2, 3 or 4 atoms, got 5"):
        plot.add(keys[:5], np.zeros(12))
    assert plot_module.figure_size(1) == plot_module.IMAGE_SIZE == plot_module.figure_size(0)
    assert plot_module.figure_size(3) == pytest.approx((6.0, 7.9))
    assert plot_module.header(plot)[7] == "kind: distance (A), dihedral (deg)"
    from matplotlib.figure import Figure

    # legends: as many columns as fit across the panel
    assert plot_module.legend_layout(Figure(figsize=(6, 3.5)), plot) == {"distance": 2}
    assert plot_module.legend_layout(Figure(figsize=(16, 3.5)), plot) == {"distance": 4}
    assert plot_module.legend_layout(Figure(figsize=(3.2, 2)), plot) == {"distance": 1}


def test_chimerax_plot_reanimate(chimerax_cmd, vdata, adk_mdcif):
    cmd, ran, UserError = chimerax_cmd
    plot_module = sys.modules["mdcif_chimerax.plot"]
    _, session, opened = _open_standins(cmd, vdata, adk_mdcif[0])
    ca = _ca(opened[0])
    cmd.mdcif_plot(session, "distance", _pick(opened[0], ca[9], ca[39]))
    [plot] = plot_module.plots(session)
    first, heard = plot.structure, []
    plot.listeners.append(heard.append)                      # as a window does
    new = cmd.reanimate(session, plot, mode=2, frames=12)
    assert first.deleted and plot.structure is new and heard == [plot]   # its own copy: closed
    assert (plot.mode, plot.n_frames, new.num_coordsets, plot.source) == (2, 12, 12, opened[0])
    np.testing.assert_allclose(plot.series[0].values, _distance(new, ca[9], ca[39]))
    other = cmd.reanimate(session, plot, model=2)            # from model 2's structure
    assert other.mdcif_model == 2 and plot.source is opened[1] and new.deleted
    assert (plot.model, plot.mode, plot.n_frames) == (2, 2, 12)
    with pytest.raises(UserError, match="no structure of model 3"):
        opened[2].deleted = True
        cmd.reanimate(session, plot, model=3)
    assert plot.structure is other and not other.deleted

    user = cmd.mdcif_animate(session, opened[0], mode=1, frames=8)     # a copy the user made
    cmd.mdcif_plot(session, "angle", _pick(user, *ca[5:8]))
    uplot = plot_module.plots(session)[-1]
    assert uplot.structure is user and not uplot.own
    cmd.reanimate(session, uplot, amplitude=1.5)
    assert not user.deleted and user.display is False                  # stopped and hidden
    assert f"coordset stop #{user.id_string}" in ran and uplot.amplitude == 1.5

    session.models.close([uplot.structure])                  # the animation is closed
    assert not uplot.alive and uplot not in plot_module.plots(session)  # no window: dropped
    with pytest.raises(UserError, match="was closed"):
        cmd.add_series(session, uplot, uplot.series[0].keys)
    plot.listeners.clear()
    session.models.close([plot.structure])
    plot.listeners.append(heard.append)                      # a window keeps a closed plot
    assert plot in plot_module.plots(session) and len(plot.series[0].values) == 12
    revived = cmd.reanimate(session, plot)                   # animated again from model 2
    assert plot.alive and revived.mdcif_model == 2 and plot.source is opened[1]


def test_chimerax_plot_reuses_shown_copies(chimerax_cmd, vdata, adk_mdcif):
    """Atoms on a model are measured on a copy already on screen (one the user made, one a
    closed window's plot made, one restored from a session) instead of a second copy."""
    cmd, ran, _ = chimerax_cmd
    plot_module = sys.modules["mdcif_chimerax.plot"]
    _, session, opened = _open_standins(cmd, vdata, adk_mdcif[0])
    model, ca = opened[0], _ca(opened[0])
    user = cmd.mdcif_animate(session, model, mode=2, frames=10)       # 'mdcif animate' by hand
    n_added = len(session.models_added)
    cmd.mdcif_plot(session, "distance", _pick(model, ca[9], ca[39]))  # no settings: any copy
    [plot] = plot_module.plots(session)
    assert len(session.models_added) == n_added and plot.structure is user
    assert plot.source is model and not plot.own
    np.testing.assert_allclose(plot.series[0].values, _distance(user, ca[9], ca[39]))
    cmd.reanimate(session, plot, frames=12)                  # the user's copy: stopped, hidden
    assert not user.deleted and user.display is False and f"coordset stop #{user.id_string}" in ran
    assert plot.own and len(session.models_added) == n_added + 1

    cmd.mdcif_plot(session, "angle", _pick(model, *ca[5:8]), mode=1)  # no copy with mode 1
    made = session.models_added[-1]
    assert len(session.models_added) == n_added + 2 and made.mdcif_mode == 1
    plot_module.plots(session).clear()                       # windows closed / session restored
    cmd.mdcif_plot(session, "dihedral", _pick(model, *ca[5:9]), mode=1)
    found = plot_module.plots(session)[-1]
    assert len(session.models_added) == n_added + 2          # found on screen again
    assert found.structure is made and found.own and found.source is model
    plot_module.plots(session).clear()
    made.display = False                                     # hidden copies are not used
    cmd.mdcif_plot(session, "dihedral", _pick(model, *ca[5:9]), mode=1)
    assert len(session.models_added) == n_added + 3
    cmd.mdcif_plot(session, "dihedral", _pick(model, *ca[5:9]), mode=1, new=True)
    assert len(session.models_added) == n_added + 4          # new true: always a new copy
    other = cmd.mdcif_animate(session, opened[1], mode=1, frames=30)  # another model's copy
    shown = cmd._shown_copies(session, model)
    assert other not in shown and made not in shown and user not in shown and len(shown) == 3

    # a restored session attaches the copies' data apart from the models': matched by file
    other.mdcif_data = vdata.MdcifFile(adk_mdcif[0])
    assert other.mdcif_data is not opened[1].mdcif_data
    assert cmd._shown_copies(session, opened[1]) == [other]
    plot_module.plots(session).clear()
    n_added = len(session.models_added)
    cmd.mdcif_plot(session, "distance", _pick(opened[1], ca[9], ca[39]))
    restored = plot_module.plots(session)[-1]
    assert len(session.models_added) == n_added and restored.structure is other
    cmd.reanimate(session, restored, model=1)                # the model found by file, too
    assert restored.source is model and len(session.models_added) == n_added + 1


# --------------------------------------------------------------------------- ChimeraX: windows

def test_chimerax_porcupine_geometry(chimerax_cmd, vdata, adk_mdcif):
    """The Porcupine window draws through ``porcupine_geometry``: the command's mesh."""
    cmd, _, UserError = chimerax_cmd
    _, session, opened = _open_standins(cmd, vdata, adk_mdcif[0])
    s = opened[0]
    n_short_seen = 0
    for mode, scale, min_length, radius in ((1, 3.0, 0.5, 0.15), (2, 2.0, 0.0, 0.3),
                                            (1, -4.0, 1.5, 0.1)):
        va, na, ta, n, n_short, ms, model = cmd.porcupine_geometry(session, s, mode, scale,
                                                                   "cluster", min_length, radius)
        surf = cmd.mdcif_porcupine(session, s, mode=mode, scale=scale, min_length=min_length,
                                   radius=radius)
        for mine, drawn in zip((va, na, ta), surf.geometry, strict=True):
            np.testing.assert_array_equal(mine, drawn)
        assert surf.SESSION_SAVE_DRAWING is True             # arrows come back from sessions
        kept = vdata.arrows(s.atoms.coords, s.keys, ms, mode, scale, min_length)[0]
        every = vdata.arrows(s.atoms.coords, s.keys, ms, mode, scale, 0.0)[0]
        assert (n, n_short) == (len(kept), len(every) - len(kept))
        assert (ms.scope, model) == ("cluster", 1)
        assert session.logger.lines[-1].startswith(f"#{surf.id_string}: {n} arrows, cluster mode")
        n_short_seen += n_short
    assert n_short_seen > 0
    glob = cmd.porcupine_geometry(session, opened[1], 1, 3.0, "global", 0.5, 0.15)
    assert (glob[5].scope, glob[6]) == ("global", 2)
    with pytest.raises(UserError, match="radius must be positive"):
        cmd.porcupine_geometry(session, s, radius=0)
    with pytest.raises(UserError, match="no arrows to draw"):
        cmd.porcupine_geometry(session, s, min_length=1e6)

    # the command the window logs for each mode
    assert cmd.porcupine_command(s) == (f"mdcif porcupine #{s.id_string} mode 1 scale 3 color "
                                        "#ff8c00 radius 0.15 minLength 0.5 scope cluster")
    rgba = np.array([42, 120, 214, 255], dtype=np.uint8)
    assert cmd.porcupine_command(s, 2, 2.5, rgba, "global", 0.25, 0.2) == (
        f"mdcif porcupine #{s.id_string} mode 2 scale 2.5 color #2a78d6 radius 0.2 "
        "minLength 0.25 scope global")


def test_chimerax_retune_animation(chimerax_cmd, vdata, adk_mdcif):
    """The Animate window changes its copy's frames in place, unless a plot measures on it."""
    cmd, ran, UserError = chimerax_cmd
    path = adk_mdcif[0]
    _, session, opened = _open_standins(cmd, vdata, path)
    model, ca = opened[0], _ca(opened[0])
    rest, data = model.atoms.coords.copy(), model.mdcif_data
    copy = cmd.mdcif_animate(session, model, mode=1, amplitude=3.0, frames=12)
    n_added = len(session.models_added)
    copy.active_coordset_id = 5                              # stopped on frame 5
    assert cmd.retune_animation(session, copy, mode=2, amplitude=-1.5, frames=8) is copy
    np.testing.assert_allclose(copy.coordsets, vdata.mode_frames(
        rest, model.keys, data.mode_set(1), 2, -1.5, 8))
    assert (copy.mdcif_mode, copy.mdcif_scope, copy.mdcif_amplitude) == (2, "cluster", -1.5)
    assert copy.num_coordsets == 8 and copy.active_coordset_id == 5
    assert copy.name == f"{path.name} model 1 cluster mode 2"
    assert len(session.models_added) == n_added               # in place: no new model
    cmd.retune_animation(session, copy, scope="global")      # None keeps the others
    np.testing.assert_allclose(copy.coordsets, vdata.mode_frames(
        rest, model.keys, data.mode_set(1, True), 2, -1.5, 8))
    assert copy.name == f"{path.name} model 1 global mode 2" and copy.mdcif_amplitude == -1.5

    # a retuned copy is found by its new settings, from the model and as itself
    assert cmd.copy_for(session, model, 2, "global", -1.5, 8) == (copy, model)
    assert cmd.copy_for(session, copy, 2, "global", -1.5, 8)[0] is copy
    assert len(session.models_added) == n_added

    # playing: as many frames, the loop goes on; other frames, it starts again at its speed
    session._coord_set_players = [types.SimpleNamespace(structure=copy, pause_frames=3)]
    n_ran = len(ran)
    cmd.retune_animation(session, copy, amplitude=2.0)
    assert len(ran) == n_ran and copy.mdcif_amplitude == 2.0
    cmd.retune_animation(session, copy, frames=10)
    assert ran[n_ran:] == [f"coordset stop #{copy.id_string}",
                           f"coordset #{copy.id_string} 1,10 loop 1000 pauseFrames 3"]
    session._coord_set_players = []
    copy.active_coordset_id = 10
    cmd.retune_animation(session, copy, frames=4)            # stopped past the new end
    assert copy.active_coordset_id == 4

    # bad settings change nothing
    before = copy.coordsets.copy()
    for kw, message in (({"frames": 2}, "frames must be at least 3"),
                        ({"mode": 99}, "mode 99 not in")):
        with pytest.raises(UserError, match=message):
            cmd.retune_animation(session, copy, **kw)
    np.testing.assert_array_equal(copy.coordsets, before)
    assert (copy.mdcif_mode, copy.num_coordsets) == (2, 4)
    with pytest.raises(UserError, match="not an 'mdcif animate' copy"):
        cmd.retune_animation(session, model, mode=1)

    # a plot measures on the copy: refused (None), the copy unchanged
    cmd.mdcif_plot(session, "distance", _pick(copy, ca[9], ca[39]))
    assert sys.modules["mdcif_chimerax.plot"].plots(session)[-1].structure is copy
    assert cmd.retune_animation(session, copy, mode=1) is None
    np.testing.assert_array_equal(copy.coordsets, before)
    assert copy.mdcif_mode == 2

    assert cmd.loop_command(copy) == f"coordset #{copy.id_string} 1,4 loop 1000"
    assert cmd.loop_command(copy, 2) == f"coordset #{copy.id_string} 1,4 loop 1000 pauseFrames 2"
    assert cmd.animate_command(model, 2, 1.5, 20, "global") == (
        f"mdcif animate #{model.id_string} mode 2 amplitude 1.5 frames 20 scope global")


def test_chimerax_retire_copy(chimerax_cmd, vdata, adk_mdcif):
    """A window closes the copies it made and no plot uses; others are stopped and hidden."""
    cmd, ran, _ = chimerax_cmd
    plot_module = sys.modules["mdcif_chimerax.plot"]
    _, session, opened = _open_standins(cmd, vdata, adk_mdcif[0])
    ca = _ca(opened[0])
    own = cmd._animate(session, opened[0], frames=6, mark=cmd.ANIMATE_MADE)
    assert session.logger.commands[-1] == (f"mdcif animate #{opened[0].id_string} mode 1 "
                                           "amplitude 3 frames 6 scope cluster")
    assert getattr(own, cmd.ANIMATE_MADE) and not getattr(own, plot_module.PLOT_MADE, False)
    user = cmd.mdcif_animate(session, opened[0], mode=2, frames=6)     # 'mdcif animate' by hand
    cmd.retire_copy(session, own, own=True)
    assert own.deleted
    cmd.retire_copy(session, user, own=False)
    assert not user.deleted and user.display is False
    assert ran[-1] == f"coordset stop #{user.id_string}"
    plotted = cmd._animate(session, opened[1], frames=6, mark=cmd.ANIMATE_MADE)
    cmd.mdcif_plot(session, "distance", _pick(plotted, ca[0], ca[5]))
    cmd.retire_copy(session, plotted, own=True)               # a plot measures on it: kept
    assert not plotted.deleted and plotted.display is False
    n_ran = len(ran)
    cmd.retire_copy(session, own, own=True)                   # closed already: nothing
    assert len(ran) == n_ran


class _QObject:
    """Stand-in for the Qt objects the MDCIF menu touches (names, deleteLater)."""

    def __init__(self, name=""):
        self._name, self.deleted_later = name, False

    def objectName(self):
        return self._name

    def setObjectName(self, name):
        self._name = name

    def deleteLater(self):
        self.deleted_later = True


class _Action(_QObject):
    def __init__(self, text, callback, tip):
        super().__init__()
        self.text, self.callback, self.tip = text, callback, tip


class _Menu(_QObject):
    def __init__(self, name):
        super().__init__(name)
        self._actions, self._menu_action = [], _QObject(name)

    def actions(self):
        return list(self._actions)

    def addAction(self, action):
        self._actions.append(action)

    def removeAction(self, action):
        self._actions.remove(action)

    def menuAction(self):
        return self._menu_action


class _MenuBar:
    def __init__(self):
        self.kids, self.shown = [], [_QObject("Help")]      # shown: the bar's entries

    def children(self):
        return list(self.kids)

    def removeAction(self, action):
        self.shown.remove(action)


class _MainWindow:
    """ChimeraX's ``MainWindow.add_menu_entry`` (ui/gui.py, 1.12): a top-level menu is found by
    object name among the menu bar's children, else made just before Help. A menu Qt has not
    deleted yet (``deleteLater``) is still a child."""

    def __init__(self):
        self.bar = _MenuBar()

    def menuBar(self):
        return self.bar

    def add_menu_entry(self, menu_names, entry_name, callback, *, tool_tip=None,
                       insertion_point=None, shortcut=None):
        [name] = menu_names
        menu = next((c for c in self.bar.children() if c.objectName() == name), None)
        if menu is None:
            menu = _Menu(name)
            self.bar.kids.append(menu)
            self.bar.shown.insert(len(self.bar.shown) - 1, menu.menuAction())
        action = _Action(entry_name, callback, tool_tip)
        menu.addAction(action)
        return action


def _mdcif_menus(main_window):
    """The MDCIF menus a user sees in the menu bar."""
    return [m for m in main_window.bar.kids
            if m.objectName() == "MDCIF" and m.menuAction() in main_window.bar.shown]


def test_chimerax_menu(chimerax_cmd):
    _, ran, _ = chimerax_cmd
    menu = sys.modules["mdcif_chimerax.menu"]
    labels = ["Porcupine...", "Animate...", "Analyze..."]

    def gui(main_window):
        ui = types.SimpleNamespace(is_gui=True, main_window=main_window, triggers=_Triggers())
        return types.SimpleNamespace(ui=ui, run=lambda text: None)

    nogui = types.SimpleNamespace(ui=types.SimpleNamespace(is_gui=False))
    menu.register(nogui)                                     # no window: nothing
    assert not menu.is_registered(nogui)
    menu.unregister(nogui)

    # installed during a session: the window exists, the entries come at once, before Help
    window = _MainWindow()
    session = gui(window)
    menu.register(session)
    [first] = _mdcif_menus(window)
    assert [a.text for a in first.actions()] == labels and all(a.tip for a in first.actions())
    assert [e.objectName() for e in window.bar.shown] == ["MDCIF", "Help"]
    assert session.ui.triggers.handlers.get("ready", []) == []
    first.actions()[0].callback(False)                       # QAction.triggered(checked)
    assert ran[-1] == 'ui tool show "MDCIF Porcupine"'
    first.actions()[1].callback()
    assert ran[-1] == 'ui tool show "MDCIF Animate"'

    menu.register(session)                                   # again: still one menu
    [second] = _mdcif_menus(window)
    assert second is not first and len(second.actions()) == 3
    assert first.deleted_later and first.actions() == [] and first.objectName() != "MDCIF"

    # 'toolshed uninstall': gone at once, not just scheduled; 'devel install' right after
    # (the old menu not deleted yet) makes one new menu with three entries
    old_actions = second.actions()
    menu.unregister(session)
    assert _mdcif_menus(window) == [] and not menu.is_registered(session)
    assert all(m.objectName() != "MDCIF" for m in window.bar.kids)
    assert all(a.deleted_later for a in old_actions) and second.actions() == []
    assert [e.objectName() for e in window.bar.shown] == ["Help"]
    menu.register(session)
    [third] = _mdcif_menus(window)
    assert third not in (first, second) and [a.text for a in third.actions()] == labels
    menu.unregister(session)

    # at startup the bundle is initialized before the main window: entries on 'ready'
    session = gui(None)
    menu.register(session)
    assert len(session.ui.triggers.handlers["ready"]) == 1 and menu.is_registered(session)
    session.ui.main_window = window = _MainWindow()
    session.ui.triggers.activate_trigger("ready", None)
    [m] = _mdcif_menus(window)
    assert [a.text for a in m.actions()] == labels
    menu.register(session)                                   # again after 'ready'
    assert len(_mdcif_menus(window)) == 1 and len(_mdcif_menus(window)[0].actions()) == 3
    assert session.ui.triggers.handlers["ready"] == []
    menu.unregister(session)

    session = gui(None)                                      # unloaded before the GUI is ready
    menu.register(session)
    menu.unregister(session)
    assert session.ui.triggers.handlers["ready"] == []
    session.ui.main_window = window = _MainWindow()
    session.ui.triggers.activate_trigger("ready", None)
    assert _mdcif_menus(window) == []


def test_chimerax_uninstall(chimerax_cmd, monkeypatch):
    """'toolshed uninstall' (ChimeraX 1.12 calls no ``finish`` then): the toolshed's "bundle
    uninstalled" trigger removes the MDCIF menu and the "add models" handler and closes the
    MDCIF windows without logging what they had pending."""
    init = sys.modules["mdcif_chimerax"]
    attach, menu = sys.modules["mdcif_chimerax.attach"], sys.modules["mdcif_chimerax.menu"]
    toolshed = sys.modules["chimerax.core.toolshed"].get_toolshed()
    uninstalled = toolshed.triggers.handlers.setdefault("bundle uninstalled", [])

    class Tool:
        def __init__(self, name):
            self.tool_name, self.deleted, self.discarded = name, False, False

        def discard_pending(self):
            self.discarded = True

        def delete(self):
            self.deleted = True

    tools = [Tool("MDCIF Porcupine"), Tool("MDCIF Animate"), Tool("MDCIF Analyze"), Tool("Log")]
    window = _MainWindow()
    session = types.SimpleNamespace(
        ui=types.SimpleNamespace(is_gui=True, main_window=window, triggers=_Triggers()),
        triggers=_Triggers(), run=lambda text: None,
        tools=types.SimpleNamespace(list=lambda: list(tools)))
    bundle = types.SimpleNamespace(name="ChimeraX-MDCIF")
    init.bundle_api.initialize(session, bundle)
    init.bundle_api.initialize(session, bundle)              # again (reinstall): one handler
    assert len(uninstalled) == 1 and attach.is_registered(session)
    assert len(_mdcif_menus(window)) == 1

    toolshed.triggers.activate_trigger("bundle uninstalled", types.SimpleNamespace(name="X"))
    assert attach.is_registered(session) and len(_mdcif_menus(window)) == 1   # not this one
    toolshed.triggers.activate_trigger("bundle uninstalled", bundle)
    assert not attach.is_registered(session) and not menu.is_registered(session)
    assert _mdcif_menus(window) == [] and session.triggers.handlers["add models"] == []
    assert [t.deleted for t in tools] == [True, True, True, False]
    assert all(t.discarded for t in tools[:3])
    assert uninstalled == []                                 # the handler removed itself

    init.bundle_api.initialize(session, bundle)              # finish (unload) removes it too
    init.bundle_api.finish(session, bundle)
    assert uninstalled == [] and not attach.is_registered(session)

    # uninstalled without the trigger (a restart pending): the "add models" handler stops
    # itself once the bundle's own modules can't be imported, and only then
    attach.register(session)

    def missing(name):
        def on_add_models(session, models):
            raise ModuleNotFoundError(f"No module named '{name}'", name=name)
        return on_add_models

    monkeypatch.setattr(attach, "on_add_models", missing("numpy"))
    with pytest.raises(ModuleNotFoundError):
        session.triggers.activate_trigger("add models", [])
    assert attach.is_registered(session)
    monkeypatch.setattr(attach, "on_add_models", missing(f"{attach.__package__}.mdcif_data"))
    session.triggers.activate_trigger("add models", [])
    assert not attach.is_registered(session)
