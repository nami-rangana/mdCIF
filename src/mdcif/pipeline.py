"""Trajectory -> mdcif file, in three steps (``build_mdcif`` / ``mdcif build`` runs all three).

1. :func:`load_trajectory` reads the analysis atoms (default CA) of frames
   ``start:stop:step``, superposes them on their average and runs the global
   PCA; ``Trajectory.projections`` are the frames' coordinates on the PCs.
2. :func:`cluster` gives one label per frame: built-in HDBSCAN on PCs 1-3, any
   scikit-learn-style estimator, or labels computed with any other tool. The
   medoid of each cluster is its representative.
3. :func:`save_mdcif` puts each medoid's deposited atoms in the common frame,
   computes each cluster's own modes (frames superposed on the medoid),
   optionally superposes everything on a reference PDB entry, and writes the
   file with python-ihm.

Steps 2 and 3 work on the in-memory ``Trajectory``, so a notebook can try
several clusterings and write each without re-reading the trajectory.
"""

import dataclasses
import time
from dataclasses import dataclass, field
from importlib import metadata
from pathlib import Path

import numpy as np

from .clustering import ClusterResult, clusters_from_labels, hdbscan_labels
from .modes import ModeSet, anm_modes, cluster_modes, cumulative_overlap, project
from .provenance import Provenance
from .reference import load_reference, match_chains, reference_coordinates
from .trajectory import apply_transform, load_universe, read_frames, superpose
from .writer import chain_layout, chain_sequences, write_mdcif

DEFAULT_DEPOSIT_SELECTION = "protein and not name H* and not resname ACE NME NHE NH2"


@dataclass
class Trajectory:
    """Analysed frames of a trajectory, superposed on their average, with the global PCA."""

    universe: object              # MDAnalysis Universe (re-read for the medoids' atoms)
    analysis: object              # AtomGroup of the analysis atoms
    deposit: object               # AtomGroup of the atoms written for each model
    mode_atoms: np.ndarray        # index in ``deposit`` of every analysis atom
    selection: str
    frames: np.ndarray            # (n,) 0-based index of each analysed frame in the trajectory
    times: np.ndarray             # (n,) time in ps (NaN if unknown)
    raw: np.ndarray               # (n, n_atoms, 3) analysis-atom coordinates as read
    aligned: np.ndarray           # (n, n_atoms, 3) superposed on the average structure
    rotations: np.ndarray         # per-frame transform raw -> aligned (see trajectory.superpose)
    centres: np.ndarray
    ref_centre: np.ndarray
    global_modes: ModeSet         # PCA of all analysed frames (with origin)
    projections: np.ndarray       # (n, n_global_modes) coordinates on the global modes, angstrom

    @property
    def n_frames(self):
        return len(self.frames)

    def deposit_coordinates(self, i):
        """Coordinates of the deposited atoms in analysed frame ``i``, in the common frame."""
        self.universe.trajectory[int(self.frames[i])]
        return apply_transform(self.deposit.positions.astype(np.float64), self.rotations[i],
                               self.centres[i], self.ref_centre)

    def frame_label(self, i):
        """'frame N' or 'replica R frame N' (1-based R, 0-based N), plus the time if known."""
        frame, time_ps = int(self.frames[i]), self.times[i]
        label = f"frame {frame}"
        readers = getattr(self.universe.trajectory, "readers", None)
        if readers and len(readers) > 1:
            start = 0
            for r, reader in enumerate(readers):
                if frame < start + reader.n_frames:
                    label = f"replica {r + 1} frame {frame - start}"
                    break
                start += reader.n_frames
        if np.isfinite(time_ps):
            label += f", t = {time_ps / 1000:.3f} ns"
        return label


def load_trajectory(topology, trajectories, *, selection="name CA",
                    deposit_selection=DEFAULT_DEPOSIT_SELECTION, start=0, stop=None, step=1,
                    n_global_modes=3, traj_format=None) -> Trajectory:
    """Read, superpose and run the global PCA (``n_global_modes`` kept, at least 3)."""
    u = load_universe(topology, trajectories, format=traj_format)
    analysis = u.select_atoms(selection)
    deposit = u.select_atoms(deposit_selection)
    if len(deposit) == 0:
        raise ValueError(f"deposit selection {deposit_selection!r} matched no atoms")
    position = {ix: k for k, ix in enumerate(deposit.ix)}
    missing = [ix for ix in analysis.ix if ix not in position]
    if missing:
        raise ValueError(f"{len(missing)} analysis atoms are not in the deposit selection")
    raw, frames, times = read_frames(u, selection, start, stop, step)
    # Global alignment on the average structure (frame 0 first, then the mean).
    aligned = superpose(raw, raw[0])[0]
    aligned, rot, centres, ref_centre = superpose(raw, aligned.mean(axis=0))
    # Plain Cartesian PCA, never mass-weighted, so that projections are well defined.
    global_modes = cluster_modes(aligned, n_modes=max(3, n_global_modes), selection=selection)
    return Trajectory(universe=u, analysis=analysis, deposit=deposit,
                      mode_atoms=np.array([position[ix] for ix in analysis.ix]),
                      selection=selection, frames=frames, times=times, raw=raw,
                      aligned=aligned, rotations=rot, centres=centres, ref_centre=ref_centre,
                      global_modes=global_modes, projections=project(aligned, global_modes))


def _software_of(obj):
    root = type(obj).__module__.split(".")[0]
    try:
        version = metadata.version("scikit-learn" if root == "sklearn" else root)
    except metadata.PackageNotFoundError:
        version = None
    return (root, version, None)


def _changed_params(estimator):
    """The simple-valued parameters of a scikit-learn-style estimator that differ from defaults."""
    if not hasattr(estimator, "get_params"):
        return {}
    try:
        defaults = type(estimator)().get_params()
    except TypeError:  # constructor needs arguments
        defaults = {}
    return {k: v for k, v in estimator.get_params().items()
            if isinstance(v, (int, float, str, bool)) and v != defaults.get(k, object())
            and k not in ("copy", "n_jobs", "verbose")}  # how it runs, not what it does


def cluster(traj: Trajectory, clusterer=None, *, labels=None, n_pcs=3, features="pca",
            method=None, software=None, details="", min_cluster_size=None, min_samples=None,
            cluster_selection="eom") -> ClusterResult:
    """Cluster the frames of ``traj``; use any clustering you prefer.

    - default: HDBSCAN (scikit-learn) on PCs 1..``n_pcs`` (``min_cluster_size``,
      ``min_samples``, ``cluster_selection``);
    - ``clusterer``: any object with ``fit_predict(X)`` (scikit-learn API), run on
      the ``features``: ``"pca"`` (PCs 1..``n_pcs``) or ``"cartesian"`` (aligned
      coordinates);
    - ``labels``: one label per analysed frame (-1 = noise) from any other tool.

    ``method`` names the algorithm in the file (default: the estimator's class
    name); well-known names get the matching IHM clustering method. ``software``
    (a name, or ``(name, version, url)``) names the tool that made ``labels``.
    """
    if features not in ("pca", "cartesian"):
        raise ValueError(f"features must be 'pca' or 'cartesian', not {features!r}")
    if features == "pca":
        if not 1 <= n_pcs <= traj.global_modes.n_modes:
            raise ValueError(f"n_pcs must be 1..{traj.global_modes.n_modes}")
        X, feature = traj.projections[:, :n_pcs], f"PCA projection (global PCs 1-{n_pcs})"
        ihm_feature = "other"
    else:
        X, feature = traj.aligned.reshape(traj.n_frames, -1), "aligned Cartesian coordinates"
        ihm_feature = "RMSD"
    params = {}
    if isinstance(software, str):
        software = (software, None, None)
    if labels is not None:
        if clusterer is not None:
            raise ValueError("give either clusterer or labels, not both")
        method = method or "user-supplied labels"
    elif clusterer is None:
        labels, params = hdbscan_labels(X, min_cluster_size, min_samples, cluster_selection)
        method = "HDBSCAN"
        software = ("scikit-learn", metadata.version("scikit-learn"), "https://scikit-learn.org")
    else:
        labels = clusterer.fit_predict(X)
        method = method or type(clusterer).__name__
        software = _software_of(clusterer)
        params = _changed_params(clusterer)
    text = ", ".join(f"{k} {v}" for k, v in params.items())
    details = "; ".join(x for x in (f"{method} on {feature}" + (f", {text}" if text else ""),
                                    details) if x)
    return clusters_from_labels(traj.aligned, labels, method=method, feature=feature,
                                ihm_feature=ihm_feature, software=software, details=details,
                                params={**params, **({"n_pcs": n_pcs} if features == "pca"
                                                     else {})})


@dataclass
class BuildSummary:
    output: Path
    n_frames: int
    n_analysis_atoms: int
    n_deposited_atoms: int
    clusters: ClusterResult
    mode_sets: list
    global_modes: ModeSet | None = None
    n_models: int = 0
    anm_overlap: list = field(default_factory=list)  # per model: cumulative ANM overlap/mode
    reference_label: str = ""
    reference_rmsd: list = field(default_factory=list)  # per model, CA RMSD to the reference
    frame_labels: list = field(default_factory=list)  # per model, where its medoid came from
    timings: dict = field(default_factory=dict)

    def report(self):
        c = self.clusters
        lines = [f"{self.output}: {self.n_frames} frames, {self.n_analysis_atoms} analysis atoms, "
                 f"{self.n_deposited_atoms} deposited atoms"]
        if self.global_modes is not None:
            g = self.global_modes
            lines.append("  global PCA var% " + "/".join(f"{100 * v:.0f}"
                                                        for v in g.variance_fraction)
                         + f" (cum {100 * g.cumulative_variance_fraction[-1]:.0f}%)")
        lines.append(f"  clustering: {c.details}, {c.n_clusters} clusters"
                     + (f", noise {100 * c.noise_fraction:.1f}%" if c.noise_fraction else ""))
        if self.n_models < c.n_clusters:
            lines.append(f"  deposited the {self.n_models} most populated clusters")
        total = len(c.labels)
        for i in range(self.n_models):
            ms = self.mode_sets[i]
            modes = ("no modes" if ms is None else
                     "modes var% " + "/".join(f"{100 * v:.0f}" for v in ms.variance_fraction)
                     + f" (cum {100 * ms.cumulative_variance_fraction[-1]:.0f}%)")
            anm = (", ANM cumulative overlap " + "/".join(f"{v:.2f}" for v in self.anm_overlap[i])
                   if i < len(self.anm_overlap) and self.anm_overlap[i] is not None else "")
            ref = (f", RMSD to {self.reference_label} {self.reference_rmsd[i]:.2f} A"
                   if self.reference_rmsd else "")
            medoid = (self.frame_labels[i] if self.frame_labels
                      else f"frame {c.representatives[i]}")
            lines.append(f"  model {i + 1}: {c.populations[i]:5d} frames "
                         f"({100 * c.populations[i] / total:4.1f}%), medoid {medoid}, "
                         f"precision {c.precision[i]:.2f} A{ref}, {modes}{anm}")
        if self.timings:
            lines.append("  timings: " + ", ".join(f"{k} {v:.1f}s"
                                                   for k, v in self.timings.items()))
        return "\n".join(lines)


def _subsample(n, max_points):
    if max_points is None or n <= max_points:
        return np.arange(n)
    return np.unique(np.linspace(0, n - 1, max_points).round().astype(int))


def _rigid(ms, r, centre, ref_centre):
    """Copy of a ModeSet moved by a rigid transform (the input is left untouched)."""
    if ms is None:
        return None
    origin = None if ms.origin is None else apply_transform(ms.origin, r, centre, ref_centre)
    return dataclasses.replace(ms, eigenvectors=ms.eigenvectors @ r, origin=origin)


def save_mdcif(output, traj: Trajectory, clusters: ClusterResult, provenance: Provenance, *,
               n_modes=None, mass_weighted=None, max_models=None, min_cluster_frames=None,
               reference=None, max_projection_points=5000, anm=False) -> BuildSummary:
    """Write ``traj`` + ``clusters`` to ``output`` (.cif or .bcif).

    ``reference`` is a PDB ID or PDBx/mmCIF file of an experimental structure of
    the same molecule (see :mod:`mdcif.reference`); ``max_projection_points``
    caps the frames stored in ``_mdcif_projection`` (``0`` stores none); ``anm``
    reports the cumulative overlap of each cluster's modes with the 10 lowest
    ANM modes of its medoid.
    """
    t = {}
    t0 = time.perf_counter()
    n_modes = n_modes or provenance.modes.n_modes or 3
    if mass_weighted is None:
        mass_weighted = provenance.modes.mass_weighted
    n_models = clusters.n_clusters if not max_models else min(max_models, clusters.n_clusters)
    reps = clusters.representatives[:n_models]
    model_coords = np.array([traj.deposit_coordinates(f) for f in reps])
    masses = traj.analysis.masses if mass_weighted else None
    min_frames = min_cluster_frames or max(10, n_modes + 1)
    mode_sets, overlaps = [], []
    for c, f in enumerate(reps):
        members = clusters.members(c)
        if len(members) < min_frames:
            mode_sets.append(None)
            overlaps.append(None)
            continue
        ms = cluster_modes(traj.raw[members], n_modes=n_modes, masses=masses,
                           reference=traj.aligned[f], selection=traj.selection)
        mode_sets.append(ms)
        overlaps.append(cumulative_overlap(ms, anm_modes(traj.aligned[f], max(n_modes, 10)))
                        if anm else None)
    global_modes = traj.global_modes
    t["modes"] = time.perf_counter() - t0

    # Reference entry: same chains/numbering, superpose everything on its first model.
    ref, ref_rmsd = None, None
    if reference is not None:
        t0 = time.perf_counter()
        ref = load_reference(reference)
        chains, atom_res = chain_layout(traj.deposit)
        cmap = match_chains(ref, chain_sequences(chains))
        is_ca = np.array([a.name == "CA" for a in traj.deposit])
        keys = [(cmap[cid][0], cmap[cid][1] + idx)
                for (cid, idx), ca in zip(atom_res, is_ca, strict=True) if ca]
        found, ref_xyz = reference_coordinates(ref, keys)
        ca_idx = np.flatnonzero(is_ca)[found]
        if len(ca_idx) < 3:
            raise ValueError(f"fewer than 3 CA atoms in common with {ref.label}")
        _, r, cen, rc = superpose(model_coords[:1, ca_idx], ref_xyz)
        model_coords = apply_transform(model_coords, r[0], cen[0], rc)
        global_modes = _rigid(global_modes, r[0], cen[0], rc)
        mode_sets = [_rigid(ms, r[0], cen[0], rc) for ms in mode_sets]
        ref_rmsd = [float(np.sqrt(((superpose(m[None, ca_idx], ref_xyz)[0][0] - ref_xyz) ** 2)
                                  .sum(axis=1).mean())) for m in model_coords]
        t["reference"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    provenance.clustering.n_clusters = clusters.n_clusters
    provenance.clustering.n_models_deposited = n_models
    labels = [traj.frame_label(f) for f in reps]
    projection = None
    if max_projection_points != 0:
        keep = _subsample(traj.n_frames, max_projection_points)
        projection = {"frames": traj.frames[keep], "times": traj.times[keep],
                      "labels": clusters.labels[keep], "coords": traj.projections[keep, :3]}
    write_mdcif(output, traj.deposit, model_coords, clusters, mode_sets, provenance,
                traj.mode_atoms, n_frames_analysed=traj.n_frames, frame_labels=labels,
                global_modes=global_modes, projection=projection, reference=ref,
                reference_rmsd=ref_rmsd)
    t["write"] = time.perf_counter() - t0
    return BuildSummary(output=Path(output), n_frames=traj.n_frames,
                        n_analysis_atoms=len(traj.analysis), n_deposited_atoms=len(traj.deposit),
                        clusters=clusters, mode_sets=mode_sets, global_modes=global_modes,
                        n_models=n_models, anm_overlap=overlaps,
                        reference_label=f"{ref.label} model 1" if ref else "",
                        reference_rmsd=ref_rmsd or [], frame_labels=labels, timings=t)


def build_mdcif(topology, trajectories, provenance: Provenance, output, *, clusterer=None,
                labels=None, selection=None, deposit_selection=DEFAULT_DEPOSIT_SELECTION,
                start=0, stop=None, step=1, traj_format=None, n_pcs=3, features="pca",
                method=None, min_cluster_size=None, min_samples=None, cluster_selection="eom",
                **save_kwargs) -> BuildSummary:
    """All three steps: :func:`load_trajectory`, :func:`cluster`, :func:`save_mdcif`.

    ``labels`` (one per analysed frame, after ``start:stop:step``) or
    ``clusterer`` replace the built-in HDBSCAN; ``save_kwargs`` go to
    :func:`save_mdcif`.
    """
    t = {}
    t0 = time.perf_counter()
    traj = load_trajectory(topology, trajectories,
                           selection=selection or provenance.modes.selection or "name CA",
                           deposit_selection=deposit_selection, start=start, stop=stop,
                           step=step, n_global_modes=max(3, n_pcs), traj_format=traj_format)
    t["read+pca"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    software = None
    if labels is not None:  # name the tool that made them, from the provenance YAML
        method = method or provenance.clustering.method or None
        software = provenance.clustering.software or None
    clusters = cluster(traj, clusterer, labels=labels, n_pcs=n_pcs, features=features,
                       method=method, software=software, min_cluster_size=min_cluster_size,
                       min_samples=min_samples, cluster_selection=cluster_selection)
    t["cluster"] = time.perf_counter() - t0
    summary = save_mdcif(output, traj, clusters, provenance, **save_kwargs)
    summary.timings = {**t, **summary.timings}
    return summary
