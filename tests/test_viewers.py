"""Viewer-side code: the shared mdcif_data module, the PyMOL plugin, the ChimeraX commands.

The viewers must not import the mdcif package, so they are loaded from their own
directories. The ChimeraX commands run here against minimal stand-ins for the
chimerax modules (ChimeraX itself is not installable with pip): this checks the
command logic, not the ChimeraX API.
"""

import importlib.util
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


class _Model:
    _next = [1]

    def __init__(self, name=""):
        self.name, self.children, self.display, self.parent = name, [], True, None
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
        self.filename = filename

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
    for n in ("BoolArg", "ColorArg", "FloatArg", "IntArg", "OpenFileNameArg"):
        setattr(commands, n, n)
    commands.CmdDesc = lambda **kw: kw
    commands.EnumOf = lambda values: ("EnumOf", tuple(values))
    commands.FileNameArg = types.SimpleNamespace(unparse=lambda path: f"<{path}>")
    commands.run = lambda session, text, log=True: (ran.append(text), session.run(text))[1]
    modules = {
        "chimerax": types.ModuleType("chimerax"),
        "chimerax.core": types.ModuleType("chimerax.core"),
        "chimerax.core.commands": commands,
        "chimerax.core.errors": types.ModuleType("chimerax.core.errors"),
        "chimerax.core.models": types.ModuleType("chimerax.core.models"),
        "chimerax.core.toolshed": types.ModuleType("chimerax.core.toolshed"),
        "chimerax.atomic": types.ModuleType("chimerax.atomic"),
    }
    modules["chimerax.core.errors"].UserError = UserError
    modules["chimerax.core.models"].Surface = _Surface
    modules["chimerax.core.models"].ADD_MODELS = "add models"
    modules["chimerax.core.toolshed"].BundleAPI = object
    modules["chimerax.atomic"].AtomicStructure = _Structure
    modules["chimerax.atomic"].AtomicStructureArg = "AtomicStructureArg"
    modules["chimerax.atomic"].AtomicStructuresArg = "AtomicStructuresArg"
    for k, v in modules.items():
        monkeypatch.setitem(sys.modules, k, v)
    pkg = types.ModuleType("mdcif_chimerax")
    pkg.__path__ = [str(CHIMERAX_SRC)]
    monkeypatch.setitem(sys.modules, "mdcif_chimerax", pkg)
    for sub in ("__init__", "mdcif_data", "attach", "cmd"):
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
        self.lines = []

    def info(self, text, **kw):
        self.lines.append(text)

    warning = info


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
        for _, func in list(self.handlers.get(name, [])):
            func(name, data)


class _Session:
    """Stand-in session: ``open`` makes one structure per model of the file (random
    coordinates) under a group model and fires "add models", as ChimeraX's mmCIF reader does."""

    def __init__(self, vdata, seed=1, bundle_init=True):
        self.logger, self.models_added, self.triggers = _Logger(), [], _Triggers()
        self.models = types.SimpleNamespace(add=self.add_models)
        self._vdata, self._rng = vdata, np.random.default_rng(seed)
        if bundle_init:
            sys.modules["mdcif_chimerax"].bundle_api.initialize(self, None)

    def add_models(self, models):
        self.models_added.extend(models)
        self.triggers.activate_trigger("add models", [a for m in models for a in m.all_models()])

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
