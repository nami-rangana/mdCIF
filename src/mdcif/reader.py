"""Read mdcif files back with python-ihm.

``ihm.reader.read`` with extra handlers for the ``_mdcif_`` categories,
returning the ``ihm.System`` plus the simulation conditions, the per-model
and global ``ModeSet``s and the frame projections. What the file leaves out
because it follows from other items is filled in here: the scope of a mode set
(global = no ensemble), mass weighting (method quasi-harmonic), the mode order
(largest eigenvalue first), variance fractions (eigenvalue / total variance) and,
for quasi-harmonic sets, the Cartesian variance (from the standard masses of the
mode atoms' elements, as in :func:`mdcif.trajectory.element_masses`).

``_mdcif_mode_vector.atom_id`` and ``_mdcif_mode_origin.atom_id`` refer to
``_atom_site.id``. python-ihm does not keep atom_site ids, so they are mapped
back to (model, atom index) assuming consecutive numbering in model order,
which is how python-ihm (and therefore ``mdcif.writer``) writes them.
"""

from dataclasses import dataclass, field
from pathlib import Path

import ihm
import ihm.reader
import numpy as np

from .modes import ModeSet
from .trajectory import standard_masses


class _SimulationHandler(ihm.reader.Handler):
    category = "_mdcif_simulation"

    def __call__(self, id, protocol_id, software_id, force_field, water_model, ensemble,
                 temperature: float, pressure: float, timestep: float, length: float,
                 num_replicas: int, frame_interval: float, enhanced_sampling, details):
        self.sysr.system._mdcif_simulation = {
            k: v for k, v in dict(
                force_field=force_field, water_model=water_model, ensemble=ensemble,
                temperature_K=temperature, pressure_bar=pressure, timestep_fs=timestep,
                length_ns=length, n_replicas=num_replicas, save_interval_ps=frame_interval,
                enhanced_sampling=enhanced_sampling, details=details).items()
            if v is not None}


def _store(handler, name):
    return handler.sysr.system.__dict__.setdefault(name, {})


class _ModeSetHandler(ihm.reader.Handler):
    category = "_mdcif_mode_set"

    def __call__(self, id: int, model_id: int, ensemble_id: int, method, atom_selection,
                 num_frames: int, total_variance: float):
        _store(self, "_mdcif_mode_sets")[id] = dict(
            scope="global" if ensemble_id is None else "cluster", model_id=model_id,
            ensemble_id=ensemble_id, method=method, mass_weighted=method == "quasi-harmonic",
            selection=atom_selection, n_frames=num_frames, total_variance=total_variance)


class _ModeHandler(ihm.reader.Handler):
    category = "_mdcif_mode"

    def __call__(self, id: int, mode_set_id: int, eigenvalue: float, model_id: int):
        if mode_set_id is None and model_id is not None:
            raise ValueError("mdcif v0.1 file (_mdcif_mode.model_id, no _mdcif_mode_set); "
                             "rebuild it with this version of mdcif")
        _store(self, "_mdcif_modes")[id] = dict(mode_set_id=mode_set_id, eigenvalue=eigenvalue)


class _ModeVectorHandler(ihm.reader.Handler):
    category = "_mdcif_mode_vector"

    def __call__(self, mode_id: int, atom_id: int, dx: float, dy: float, dz: float):
        _store(self, "_mdcif_vectors").setdefault(mode_id, []).append((atom_id, dx, dy, dz))


class _ModeOriginHandler(ihm.reader.Handler):
    category = "_mdcif_mode_origin"

    def __call__(self, mode_set_id: int, atom_id: int, x: float, y: float, z: float):
        _store(self, "_mdcif_origins").setdefault(mode_set_id, []).append((atom_id, x, y, z))


class _ProjectionHandler(ihm.reader.Handler):
    category = "_mdcif_projection"

    def __call__(self, frame: int, time: float, model_id: int, proj_1: float, proj_2: float,
                 proj_3: float):
        rows = self.sysr.system.__dict__.setdefault("_mdcif_projection_rows", [])
        rows.append((frame, time, model_id, proj_1, proj_2, proj_3))


HANDLERS = [_SimulationHandler, _ModeSetHandler, _ModeHandler, _ModeVectorHandler,
            _ModeOriginHandler, _ProjectionHandler]


@dataclass
class MdcifData:
    system: object                 # ihm.System
    simulation: dict               # _mdcif_simulation items
    mode_sets: dict                # model_id -> ModeSet of that model's cluster
    global_modes: ModeSet | None   # PCA of all frames (clustering axes), or None
    projection: dict | None        # frame, time, model_id (0 = none), coords (n, 3)
    mode_set_info: dict = field(default_factory=dict)  # mode_set id -> _mdcif_mode_set row


def _nan(v):
    return np.nan if v is None else v


def read_mdcif(path) -> MdcifData:
    """Read an mdcif file.

    Each ``ModeSet`` gets extra attributes ``atom_indices`` (index, within the
    atom list of model ``model_id``, of every atom the vector components refer
    to), ``model_id`` and ``scope``; the global set also has ``origin``.
    """
    path = Path(path)
    if path.suffix.lower() == ".bcif":
        with open(path, "rb") as fh:
            system = ihm.reader.read(fh, format="BCIF", handlers=HANDLERS)[0]
    else:
        with open(path, encoding="utf-8") as fh:
            system = ihm.reader.read(fh, handlers=HANDLERS)[0]

    offsets, model_atoms = {}, {}
    n = 0
    for _, model in system._all_models():
        atoms = list(model.get_atoms())
        model_atoms[int(model._id)] = atoms
        offsets[int(model._id)] = (n, len(atoms))
        n += len(atoms)

    def atom_indices(model_id, atom_ids, what):
        start, count = offsets[model_id]
        idx = np.asarray(atom_ids) - 1 - start
        if np.any((idx < 0) | (idx >= count)):
            raise ValueError(f"{what}: atom ids outside model {model_id}")
        return idx

    infos = getattr(system, "_mdcif_mode_sets", {})
    modes = getattr(system, "_mdcif_modes", {})
    vectors = getattr(system, "_mdcif_vectors", {})
    origins = getattr(system, "_mdcif_origins", {})
    by_set = {}
    # mode 1 of a set = its largest eigenvalue
    for mode_id, m in sorted(modes.items(), key=lambda kv: (kv[1]["mode_set_id"],
                                                             -kv[1]["eigenvalue"], kv[0])):
        by_set.setdefault(m["mode_set_id"], []).append((mode_id, m))
    mode_sets, global_modes = {}, None
    for set_id, rows in by_set.items():
        if set_id not in infos:
            raise ValueError(f"modes refer to unknown mode set {set_id}")
        info = infos[set_id]
        model_id = info["model_id"]
        vecs, atom_idx = [], None
        for mode_id, _ in rows:
            comp = sorted(vectors.get(mode_id, []))
            idx = atom_indices(model_id, [a for a, *_ in comp], f"mode {mode_id}")
            if atom_idx is None:
                atom_idx = idx
            elif not np.array_equal(atom_idx, idx):
                raise ValueError(f"mode set {set_id}: modes use different atoms")
            vecs.append([c[1:] for c in comp])
        origin = None
        if set_id in origins:
            comp = sorted(origins[set_id])
            idx = atom_indices(model_id, [a for a, *_ in comp], f"mode set {set_id} origin")
            if not np.array_equal(idx, atom_idx):
                raise ValueError(f"mode set {set_id}: origin and vectors use different atoms")
            origin = np.array([c[1:] for c in comp], dtype=float)
        evals = np.array([m["eigenvalue"] for _, m in rows], dtype=float)
        vecs = np.array(vecs, dtype=float)
        cart = None
        if info["mass_weighted"]:  # NaN if a mode atom's element has no standard mass
            masses = standard_masses([model_atoms[model_id][i].type_symbol for i in atom_idx])
            cart = evals / ((vecs ** 2).sum(axis=2) @ masses)
        total = info["total_variance"]
        ms = ModeSet(
            eigenvalues=evals, eigenvectors=vecs,
            variance_fraction=evals / total if total else np.full(len(evals), np.nan),
            method=info["method"], mass_weighted=info["mass_weighted"],
            selection=info["selection"], n_frames=info["n_frames"], origin=origin,
            total_variance=total, cartesian_variance=cart)
        ms.atom_indices = atom_idx
        ms.model_id = model_id
        ms.scope = info["scope"]
        if info["scope"] == "global":
            global_modes = ms
        else:
            mode_sets[model_id] = ms

    projection = None
    rows = getattr(system, "_mdcif_projection_rows", None)
    if rows:
        rows.sort(key=lambda r: r[0])
        projection = {  # always on the global set
            "mode_set_id": next((i for i, s in infos.items() if s["scope"] == "global"), None),
            "frame": np.array([r[0] for r in rows], dtype=int),
            "time": np.array([_nan(r[1]) for r in rows], dtype=float),
            "model_id": np.array([r[2] or 0 for r in rows], dtype=int),
            "coords": np.array([[_nan(v) for v in r[3:]] for r in rows], dtype=float),
        }
    return MdcifData(system=system, simulation=getattr(system, "_mdcif_simulation", {}),
                     mode_sets=mode_sets, global_modes=global_modes, projection=projection,
                     mode_set_info=infos)
