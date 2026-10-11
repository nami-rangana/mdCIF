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
    assert peak == pytest.approx(2.0 * np.sqrt(ms.eigenvalues[0]) * np.abs(d).max(), rel=0.05)
    starts, ends = vdata.arrows(xyz, keys, ms, 1, scale=3.0, min_length=0.0)
    assert len(starts) == len(ms.keys)
    lengths = np.linalg.norm(ends - starts, axis=1)
    assert np.sqrt((lengths ** 2).sum()) == pytest.approx(3.0 * np.sqrt(ms.eigenvalues[0]))
    va, na, ta = vdata.arrow_mesh(starts[:5], ends[:5], sides=8)
    assert va.shape == na.shape and ta.max() < len(va)
    assert np.allclose(np.linalg.norm(na, axis=1), 1.0, atol=1e-5)
    with pytest.raises(ValueError):
        vdata.displacements(keys, ms, ms.n_modes + 1)


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


class _Model:
    _next = [1]

    def __init__(self, name=""):
        self.name, self.children, self.display = name, [], True
        self.id = (_Model._next[0],)
        _Model._next[0] += 1

    @property
    def id_string(self):
        return ".".join(map(str, self.id))

    def all_models(self):
        out = [self]
        for c in self.children:
            out.extend(c.all_models())
        return out

    def add(self, models):
        self.children.extend(models)


class _Structure(_Model):
    def __init__(self, keys, coords, name="s"):
        super().__init__(name)
        self.keys, self.atoms, self.coordsets = keys, _Atoms(keys, coords), None

    def copy(self, name):
        return _Structure(self.keys, self.atoms.coords.copy(), name)

    def add_coordsets(self, xyz, replace=False):
        self.coordsets = xyz


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
    modules["chimerax.core.toolshed"].BundleAPI = object
    modules["chimerax.atomic"].AtomicStructure = _Structure
    modules["chimerax.atomic"].AtomicStructureArg = "AtomicStructureArg"
    modules["chimerax.atomic"].AtomicStructuresArg = "AtomicStructuresArg"
    for k, v in modules.items():
        monkeypatch.setitem(sys.modules, k, v)
    pkg = types.ModuleType("mdcif_chimerax")
    pkg.__path__ = [str(CHIMERAX_SRC)]
    monkeypatch.setitem(sys.modules, "mdcif_chimerax", pkg)
    for sub in ("__init__", "mdcif_data", "cmd"):
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


class _Session:
    def __init__(self, structures):
        self.logger, self.models_added = _Logger(), []
        self._structures = structures
        self.models = types.SimpleNamespace(add=self.models_added.extend)

    def run(self, text):
        if text.startswith("open"):
            parent = _Model("group")
            parent.children = self._structures
            return [parent]
        return None


def test_chimerax_commands(chimerax_cmd, vdata, adk_mdcif):
    cmd, ran, UserError = chimerax_cmd
    path = adk_mdcif[0]
    f = vdata.MdcifFile(path)
    rng = np.random.default_rng(1)
    structures = [_Structure(f.model_keys[m], rng.normal(size=(len(f.model_keys[m]), 3)),
                             name=f"model {m}") for m in f.models]
    session = _Session(structures)
    opened = cmd.mdcif_open(session, str(path))
    assert [s.mdcif_model for s in opened] == f.models
    assert ran[0].startswith(f'open "{path}" format mmcif')
    cmd.mdcif_modes(session, opened)
    assert any("global set" in line for line in session.logger.lines)

    copy = cmd.mdcif_animate(session, opened[0], mode=1, amplitude=2.0, frames=12)
    assert copy.coordsets.shape == (12, len(f.model_keys[1]), 3)
    assert session.models_added == [copy] and opened[0].display is False
    assert ran[-1].startswith(f"coordset #{copy.id_string} 1,12")
    copy = cmd.mdcif_animate(session, opened[1], mode=1, frames=6, scope="global")
    assert copy.coordsets.shape[0] == 6

    surf = cmd.mdcif_porcupine(session, opened[0], mode=1, scale=3.0, min_length=0.0)
    va, na, ta = surf.geometry
    assert len(va) > 0 and surf in opened[0].children

    cmd.mdcif_morph(session, opened, frames=5)
    assert ran[-1].startswith("morph #") and "frames 5" in ran[-1]
    with pytest.raises(UserError):
        cmd.mdcif_animate(session, _Structure([], np.zeros((0, 3))), mode=1)
