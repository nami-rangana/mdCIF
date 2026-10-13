"""Read mdcif files back with python-ihm.

``ihm.reader.read`` with extra handlers for the ``_mdcif_`` categories,
returning the ``ihm.System`` plus the simulation conditions, the per-model
and global ``ModeSet``s and the frame projections.

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

    def __call__(self, id: int, scope, model_id: int, ensemble_id: int, method, mass_weighted,
                 atom_selection, num_frames: int, total_variance: float, details):
        _store(self, "_mdcif_mode_sets")[id] = dict(
            scope=scope, model_id=model_id, ensemble_id=ensemble_id, method=method,
            mass_weighted=self.get_bool(mass_weighted), selection=atom_selection,
            n_frames=num_frames, total_variance=total_variance, details=details)


class _ModeHandler(ihm.reader.Handler):
    category = "_mdcif_mode"

    def __call__(self, id: int, mode_set_id: int, mode_index: int, eigenvalue: float,
                 variance_fraction: float, cartesian_variance: float, model_id: int):
        if mode_set_id is None and model_id is not None:
            raise ValueError("mdcif v0.1 file (_mdcif_mode.model_id, no _mdcif_mode_set); "
                             "rebuild it with this version of mdcif")
        _store(self, "_mdcif_modes")[id] = dict(mode_set_id=mode_set_id, mode_index=mode_index,
                                                eigenvalue=eigenvalue,
                                                variance_fraction=variance_fraction,
                                                cartesian_variance=cartesian_variance)


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

    def __call__(self, id: int, mode_set_id: int, frame: int, time: float, model_id: int,
                 proj_1: float, proj_2: float, proj_3: float):
        rows = self.sysr.system.__dict__.setdefault("_mdcif_projection_rows", [])
        rows.append((mode_set_id, frame, time, model_id, proj_1, proj_2, proj_3))


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

    offsets = {}
    n = 0
    for _, model in system._all_models():
        count = sum(1 for _ in model.get_atoms())
        offsets[int(model._id)] = (n, count)
        n += count

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
    for mode_id, m in sorted(modes.items(), key=lambda kv: (kv[1]["mode_set_id"],
                                                             kv[1]["mode_index"])):
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
        cart = [m["cartesian_variance"] for _, m in rows]
        ms = ModeSet(
            eigenvalues=np.array([m["eigenvalue"] for _, m in rows]),
            eigenvectors=np.array(vecs, dtype=float),
            variance_fraction=np.array([_nan(m["variance_fraction"]) for _, m in rows]),
            method=info["method"], mass_weighted=bool(info["mass_weighted"]),
            selection=info["selection"], n_frames=info["n_frames"], origin=origin,
            total_variance=info["total_variance"],
            cartesian_variance=(None if all(v is None for v in cart)
                                else np.array([_nan(v) for v in cart], dtype=float)))
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
        rows.sort(key=lambda r: r[1])
        projection = {
            "mode_set_id": rows[0][0],
            "frame": np.array([r[1] for r in rows], dtype=int),
            "time": np.array([_nan(r[2]) for r in rows], dtype=float),
            "model_id": np.array([r[3] or 0 for r in rows], dtype=int),
            "coords": np.array([[_nan(v) for v in r[4:]] for r in rows], dtype=float),
        }
    return MdcifData(system=system, simulation=getattr(system, "_mdcif_simulation", {}),
                     mode_sets=mode_sets, global_modes=global_modes, projection=projection,
                     mode_set_info=infos)
