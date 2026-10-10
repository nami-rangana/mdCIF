"""Write mdcif files with python-ihm.

Builds an ``ihm.System`` (entities, assembly, models, model groups,
ensembles, protocol, software, datasets) and passes extra ``ihm.dumper.Dumper``
subclasses for the ``_mdcif_`` categories. With a reference PDB entry
(:mod:`mdcif.reference`) its entities, chains, source, sequence references and
citations are reused. Output format follows the file suffix: ``.cif`` ->
mmCIF, ``.bcif`` -> BinaryCIF.
"""

import re
from importlib import metadata
from pathlib import Path

import ihm
import ihm.analysis
import ihm.dataset
import ihm.dumper
import ihm.format_bcif
import ihm.location
import ihm.model
import ihm.protocol
import ihm.representation
import numpy as np

from . import __version__
from .reference import match_chains

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
# Force-field protonation/disulfide variants -> standard residue
RESNAME_ALIASES = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "HISD": "HIS", "HISE": "HIS", "HISH": "HIS", "HIS1": "HIS", "HIS2": "HIS",
    "CYX": "CYS", "CYM": "CYS", "CYS2": "CYS", "ASH": "ASP", "GLH": "GLU",
    "LYN": "LYS", "LSN": "LYS", "ARN": "ARG", "TYM": "TYR",
}


def residue_code(resname):
    """One-letter code for an MD residue name (handles variants and Amber N/C termini)."""
    name = resname.upper()
    if name not in THREE_TO_ONE and name not in RESNAME_ALIASES and len(name) == 4 \
            and name[0] in "NC":
        name = name[1:]  # Amber terminal residues, e.g. NALA / CLYS
    name = RESNAME_ALIASES.get(name, name)
    return THREE_TO_ONE.get(name)


def _elements(atoms):
    try:
        elements = list(atoms.elements)
    except Exception:  # MDAnalysis NoDataError when the topology has no elements
        elements = [""] * len(atoms)
    out = []
    for el, name in zip(elements, atoms.names, strict=True):
        el = (el or "").strip()
        if not el:
            el = re.sub(r"[^A-Za-z]", "", name)[:1]
        out.append(el.upper())
    return out


def _chain_ids(atoms):
    for attr in ("chainIDs", "segids"):
        try:
            ids = [str(c).strip() for c in getattr(atoms, attr)]
        except Exception:
            continue
        if any(ids):
            return [c or "A" for c in ids]
    return ["A"] * len(atoms)


def chain_layout(atoms):
    """Residues of each chain of an AtomGroup, in order of appearance.

    Returns ``(chains, atom_res)``: ``{chain_id: {"resnames": [...], "resids": [...]}}``
    and, for every atom, ``(chain_id, 1-based residue index within its chain)``.
    """
    chains, atom_res, last = {}, [], {}
    for atom, cid in zip(atoms, _chain_ids(atoms), strict=True):
        res = chains.setdefault(cid, {"resnames": [], "resids": []})
        if last.get(cid) != atom.resindex:
            last[cid] = atom.resindex
            res["resnames"].append(atom.resname)
            res["resids"].append(int(atom.resid))
        atom_res.append((cid, len(res["resnames"])))
    return chains, atom_res


def chain_sequences(chains):
    """``{chain_id: one-letter sequence}``; raises for residues without a standard code."""
    out = {}
    for cid, res in chains.items():
        seq = []
        for rn in res["resnames"]:
            code = residue_code(rn)
            if code is None:
                raise ValueError(f"unsupported residue {rn!r} in chain {cid}; "
                                 "exclude it with the deposit selection")
            seq.append(code)
        out[cid] = "".join(seq)
    return out


class MDTrajectoryDataset(ihm.dataset.Dataset):
    """The full MD trajectory the deposited models were extracted from."""

    data_type = "Other"


def _software(name, classification, description, location=None, version=None):
    return ihm.Software(name=name, classification=classification, description=description,
                        location=location, version=version)


def _pkg_version(dist):
    try:
        return metadata.version(dist)
    except metadata.PackageNotFoundError:
        return None


def _entities_and_asyms(chains, provenance, reference):
    """``(entities, asyms, seq_offsets)`` keyed by MD chain ID."""
    seqs = chain_sequences(chains)
    if reference is not None:
        cmap = reference.chain_map or match_chains(reference, seqs)
        asyms = {cid: cmap[cid][0] for cid in chains}
        offsets = {cid: cmap[cid][1] for cid in chains}
        entities = list(dict.fromkeys(a.entity for a in asyms.values()))
        return entities, asyms, offsets
    entities, asyms = {}, {}
    for cid, seq in seqs.items():
        if seq not in entities:
            entities[seq] = ihm.Entity(list(seq), description=provenance.name or None)
        resids = chains[cid]["resids"]
        asyms[cid] = ihm.AsymUnit(
            entities[seq], details=f"Chain {cid}", strand_id=cid,
            auth_seq_id_map={i + 1: r for i, r in enumerate(resids)})
    return list(entities.values()), asyms, dict.fromkeys(chains, 0)


def build_system(atoms, model_coords, clusters, mode_sets, provenance, mode_atoms,
                 n_frames_analysed=None, frame_labels=None, global_modes=None,
                 projection=None, reference=None, reference_rmsd=None):
    """Build the ``ihm.System`` for one mdcif file.

    :param atoms: MDAnalysis AtomGroup of the deposited atoms (topology only).
    :param model_coords: ``(n_models, n_atoms, 3)`` coordinates of each representative;
        model ``c`` is cluster ``c`` (``n_models`` may be below the number of clusters).
    :param clusters: :class:`~mdcif.clustering.ClusterResult`.
    :param mode_sets: one :class:`~mdcif.modes.ModeSet` (or None) per model.
    :param provenance: :class:`~mdcif.provenance.Provenance`.
    :param mode_atoms: indices into ``atoms`` of the atoms the modes are defined on,
        in the same order as the mode eigenvectors.
    :param frame_labels: optional text per model saying where its medoid came from
        (default ``"frame <index>"``).
    :param global_modes: optional ModeSet of all frames (the clustering axes), with origin.
    :param projection: optional dict ``frames``, ``times``, ``labels``, ``coords``: frame
        coordinates along the global modes and the cluster of each frame (-1 = noise).
    :param reference: optional :class:`~mdcif.reference.Reference` PDB entry.
    :param reference_rmsd: CA RMSD of each model to the reference's first model.
    """
    prov, sim = provenance, provenance.simulation
    chains, atom_res = chain_layout(atoms)
    elements = _elements(atoms)
    n_models = len(model_coords)
    entities, asyms, offsets = _entities_and_asyms(chains, prov, reference)

    title = prov.title or prov.name
    if not title and reference is not None:
        title = f"MD ensemble of {reference.system.title or reference.label}"
    system = ihm.System(title=title or "mdcif",
                        id=re.sub(r"[^0-9A-Za-z_-]", "_", prov.name or "mdcif"))
    system.authors.extend(prov.authors)
    system.entities.extend(entities)
    system.asym_units.extend(dict.fromkeys(asyms.values()))
    if reference is not None:
        # Fresh copies: they describe the experimental entry, not this file, so none is
        # primary (python-ihm would keep the 'primary' ID read from the entry).
        for c in reference.system.citations:
            system.citations.append(ihm.Citation(
                pmid=c.pmid, title=c.title, journal=c.journal, volume=c.volume,
                page_range=c.page_range, year=c.year, authors=c.authors, doi=c.doi,
                is_primary=False, journal_astm=c.journal_astm, journal_csd=c.journal_csd,
                journal_issn=c.journal_issn))

    engine_sw = _software(sim.engine or "unknown", "molecular dynamics",
                          "MD engine used for the simulation", version=sim.engine_version or None)
    mdcif_sw = _software("mdcif", "data processing",
                         "PCA, clustering, essential-mode analysis and mmCIF output",
                         version=__version__)
    system.software.extend([engine_sw, mdcif_sw])
    cluster_sw = None
    if clusters.software:  # (name, version, url); None for labels from an unnamed tool
        name, version, url = clusters.software
        cluster_sw = _software(name, "data processing", f"{clusters.method} clustering",
                               url, version)
        system.software.append(cluster_sw)
    system.software.append(_software("MDAnalysis", "data processing", "Trajectory reading",
                                     "https://www.mdanalysis.org", _pkg_version("MDAnalysis")))

    assembly = ihm.Assembly(list(dict.fromkeys(asyms.values())), name="Simulated system")
    rep = ihm.representation.Representation(
        [ihm.representation.AtomicSegment(a, rigid=False) for a in dict.fromkeys(asyms.values())],
        name="All-atom")

    # --- datasets: starting structure + full trajectory
    datasets = []
    if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", prov.source_structure or ""):
        datasets.append(ihm.dataset.PDBDataset(
            ihm.location.PDBLocation(prov.source_structure.upper()),
            details="Starting structure of the simulation"))
    if prov.trajectory_doi:
        repo = ihm.location.Repository(doi=prov.trajectory_doi, url=prov.trajectory_url or None)
        datasets.append(MDTrajectoryDataset(
            ihm.location.InputFileLocation(None, repo=repo, details="Full MD trajectory"),
            details="Full MD trajectory"))
    elif prov.trajectory_accession:
        loc = ihm.location.DatabaseLocation(
            prov.trajectory_accession,
            details=" ".join(x for x in (prov.trajectory_database, prov.trajectory_url) if x))
        datasets.append(MDTrajectoryDataset(loc, details="Full MD trajectory"))
    dgroup = ihm.dataset.DatasetGroup(datasets, name="Simulation inputs")

    # --- protocol: MD step + clustering post-process (+ comparison with the reference)
    n_frames = int(n_frames_analysed or len(clusters.labels))
    protocol = ihm.protocol.Protocol(name="Molecular dynamics")
    protocol.steps.append(ihm.protocol.Step(
        assembly=assembly, dataset_group=dgroup, method="Molecular dynamics",
        name="Production MD", num_models_begin=0, num_models_end=n_frames,
        software=engine_sw, ensemble=True, description=sim.details or None))
    analysis = ihm.analysis.Analysis()
    noise = (f"{100 * clusters.noise_fraction:.1f}% of frames in no cluster (noise)"
             if clusters.noise_fraction > 0 else "")
    kept = (f"{n_models} most populated of {clusters.n_clusters} clusters deposited"
            if n_models < clusters.n_clusters else "")
    feature = clusters.ihm_feature
    cstep = ihm.analysis.ClusterStep(
        feature=feature, num_models_begin=n_frames, num_models_end=clusters.n_clusters,
        assembly=assembly, software=cluster_sw,
        details="; ".join(x for x in (clusters.details, prov.frames_used, noise, kept) if x))
    analysis.steps.append(cstep)
    ref_label = ""
    if reference is not None and reference_rmsd is not None:
        ref_label = f"{reference.label} model 1"
        method = f" ({reference.method})" if reference.method else ""
        ref_group = None
        if reference.pdb_id:
            ref_group = ihm.dataset.DatasetGroup([ihm.dataset.PDBDataset(
                ihm.location.PDBLocation(reference.pdb_id),
                details=f"Experimental reference structure{method}")],
                name="Reference structure")
        analysis.steps.append(ihm.analysis.ValidationStep(
            feature="RMSD", num_models_begin=n_models, num_models_end=n_models,
            assembly=assembly, dataset_group=ref_group, software=mdcif_sw,
            details=f"CA RMSD of each deposited model to {ref_label}{method}, after "
                    "superposition: " + ", ".join(f"model {i + 1} {r:.2f} A"
                                                  for i, r in enumerate(reference_rmsd))))
    protocol.analyses.append(analysis)

    # --- one model / model group / ensemble per deposited cluster, ordered by population
    state = ihm.model.State(name="MD conformational ensemble")
    models, ensembles = [], []
    total = len(clusters.labels)
    for c in range(n_models):
        frame = (frame_labels[c] if frame_labels
                 else f"frame {clusters.representatives[c]}")
        model = ihm.model.Model(assembly, protocol, rep, name=f"Cluster {c + 1} medoid ({frame})")
        xyz = model_coords[c]
        for i, (cid, res_idx) in enumerate(atom_res):
            model.add_atom(ihm.model.Atom(
                asym_unit=asyms[cid], seq_id=offsets[cid] + res_idx, atom_id=atoms[i].name,
                type_symbol=elements[i], x=float(xyz[i, 0]), y=float(xyz[i, 1]),
                z=float(xyz[i, 2])))
        pop = int(clusters.populations[c])
        details = f"{pop} of {total} frames ({pop / total:.1%})"
        if ref_label:
            details += f"; medoid CA RMSD to {ref_label} {reference_rmsd[c]:.2f} A"
        group = ihm.model.ModelGroup([model], name=f"Cluster {c + 1}", details=details)
        group.representatives.append(ihm.model.ModelRepresentative(model, "medoid"))
        state.append(group)
        ensemble = ihm.model.Ensemble(
            group, num_models=pop, post_process=cstep,
            clustering_method=clusters.ihm_method,
            clustering_feature=feature, name=f"Cluster {c + 1}",
            precision=float(clusters.precision[c]), superimposed=True,
            details="precision = mean RMSD of cluster frames to the medoid (angstrom)")
        system.ensembles.append(ensemble)
        models.append(model)
        ensembles.append(ensemble)
    system.state_groups.append(ihm.model.StateGroup([state]))

    # (scope, model whose atoms the vectors refer to, ensemble, ModeSet, write origin, details)
    sets = []
    if global_modes is not None:
        sets.append(("global", models[0], None, global_modes, True,
                     "PCA of all analysed frames superposed on their average"
                     + (f"; the frames were clustered on modes 1-{clusters.params['n_pcs']}"
                        if "n_pcs" in clusters.params else "")))
    for c in range(n_models):
        if c < len(mode_sets) and mode_sets[c] is not None:
            sets.append(("cluster", models[c], ensembles[c], mode_sets[c], False,
                         f"frames of cluster {c + 1} superposed on its medoid"))
    system._mdcif = {"simulation": sim, "protocol": protocol, "engine": engine_sw,
                     "mode_sets": sets, "mode_atoms": np.asarray(mode_atoms, dtype=int),
                     "projection": projection, "models": models}
    return system


def _none_if_empty(v):
    return None if v in ("", None) else v


class _SimulationDumper(ihm.dumper.Dumper):
    def dump(self, system, writer):
        info = getattr(system, "_mdcif", None)
        if not info:
            return
        s = info["simulation"]
        ens = (s.ensemble or "").upper()
        ens = ens if ens in ("NVE", "NVT", "NPT", "NPAT") else ("other" if ens else None)
        with writer.loop("_mdcif_simulation", [
                "id", "protocol_id", "software_id", "force_field", "water_model", "ensemble",
                "temperature", "pressure", "timestep", "length", "num_replicas",
                "frame_interval", "enhanced_sampling", "details"]) as lp:
            lp.write(id=1, protocol_id=info["protocol"]._id,
                     software_id=info["engine"]._id,
                     force_field=_none_if_empty(s.force_field),
                     water_model=_none_if_empty(s.water_model), ensemble=ens,
                     temperature=s.temperature_K, pressure=s.pressure_bar,
                     timestep=s.timestep_fs, length=s.length_ns, num_replicas=s.n_replicas,
                     frame_interval=s.save_interval_ps,
                     enhanced_sampling=_none_if_empty(s.enhanced_sampling),
                     details=_none_if_empty(s.details))


class _ModeDumper(ihm.dumper.Dumper):
    def dump(self, system, writer):
        info = getattr(system, "_mdcif", None)
        if not info or not info["mode_sets"]:
            return
        # _atom_site.id is a running count over all models, in output order.
        offsets, n = {}, 0
        for _, model in system._all_models():
            offsets[id(model)] = n
            n += sum(1 for _ in model.get_atoms())
        if isinstance(writer, ihm.format_bcif.BinaryCifWriter):
            def fmt(v):
                return float(v)
        else:
            # The mmCIF writer prints floats with 3 decimals: too coarse for unit vectors.
            def fmt(v):
                return f"{v:.6g}" if abs(v) >= 1e-4 or v == 0 else f"{v:.4e}"
        atom_idx = info["mode_atoms"]
        with writer.loop("_mdcif_mode_set", [
                "id", "scope", "model_id", "ensemble_id", "method", "mass_weighted",
                "atom_selection", "num_frames", "total_variance", "details"]) as lp:
            for set_id, (scope, model, ensemble, ms, _, details) in enumerate(
                    info["mode_sets"], 1):
                lp.write(id=set_id, scope=scope, model_id=model._id,
                         ensemble_id=ensemble._id if ensemble is not None else None,
                         method=ms.method, mass_weighted="YES" if ms.mass_weighted else "NO",
                         atom_selection=ms.selection, num_frames=ms.n_frames,
                         total_variance=(fmt(ms.total_variance)
                                         if ms.total_variance is not None else None),
                         details=details)
        rows = []
        with writer.loop("_mdcif_mode", [
                "id", "mode_set_id", "mode_index", "eigenvalue", "variance_fraction",
                "cumulative_variance_fraction"]) as lp:
            mode_id = 0
            for set_id, (_, model, _, ms, _, _) in enumerate(info["mode_sets"], 1):
                cum = ms.cumulative_variance_fraction
                for i in range(ms.n_modes):
                    mode_id += 1
                    lp.write(id=mode_id, mode_set_id=set_id, mode_index=i + 1,
                             eigenvalue=fmt(ms.eigenvalues[i]),
                             variance_fraction=fmt(ms.variance_fraction[i]),
                             cumulative_variance_fraction=fmt(cum[i]))
                    rows.append((mode_id, offsets[id(model)], ms.eigenvectors[i]))
        with writer.loop("_mdcif_mode_vector", ["mode_id", "atom_id", "dx", "dy", "dz"]) as lp:
            for mode_id, offset, vec in rows:
                for k, a in enumerate(atom_idx):
                    lp.write(mode_id=mode_id, atom_id=offset + int(a) + 1,
                             dx=fmt(vec[k, 0]), dy=fmt(vec[k, 1]), dz=fmt(vec[k, 2]))
        with writer.loop("_mdcif_mode_origin", ["mode_set_id", "atom_id", "x", "y", "z"]) as lp:
            for set_id, (_, model, _, ms, write_origin, _) in enumerate(info["mode_sets"], 1):
                if not write_origin or ms.origin is None:
                    continue
                for k, a in enumerate(atom_idx):
                    lp.write(mode_set_id=set_id, atom_id=offsets[id(model)] + int(a) + 1,
                             x=float(ms.origin[k, 0]), y=float(ms.origin[k, 1]),
                             z=float(ms.origin[k, 2]))
        proj = info.get("projection")
        if proj is None or info["mode_sets"][0][0] != "global":
            return
        models = info["models"]
        n_cols = min(3, proj["coords"].shape[1])
        with writer.loop("_mdcif_projection", [
                "id", "mode_set_id", "frame", "time", "model_id", "proj_1", "proj_2",
                "proj_3"]) as lp:
            for i in range(len(proj["frames"])):
                label = int(proj["labels"][i])
                p = [float(proj["coords"][i, j]) if j < n_cols else None for j in range(3)]
                time = proj.get("times")
                lp.write(id=i + 1, mode_set_id=1, frame=int(proj["frames"][i]),
                         time=(float(time[i]) if time is not None and np.isfinite(time[i])
                               else None),
                         model_id=models[label]._id if 0 <= label < len(models) else None,
                         proj_1=p[0], proj_2=p[1], proj_3=p[2])


MDCIF_DUMPERS = [_SimulationDumper, _ModeDumper]


def write_system(path, system):
    """Write an already built system; ``.bcif`` gives BinaryCIF, anything else mmCIF."""
    path = Path(path)
    if path.suffix.lower() == ".bcif":
        with open(path, "wb") as fh:
            ihm.dumper.write(fh, [system], format="BCIF", dumpers=MDCIF_DUMPERS)
    else:
        with open(path, "w", encoding="utf-8") as fh:
            ihm.dumper.write(fh, [system], dumpers=MDCIF_DUMPERS)
    return path


def write_mdcif(path, atoms, model_coords, clusters, mode_sets, provenance, mode_atoms,
                **kwargs):
    """Write cluster representatives, provenance and modes to ``path``.

    ``kwargs`` are passed to :func:`build_system`.
    """
    system = build_system(atoms, model_coords, clusters, mode_sets, provenance, mode_atoms,
                          **kwargs)
    write_system(path, system)
    return system
