"""Provenance metadata read from a YAML file (see examples/provenance_template.yaml).

These dataclasses are what ``writer.py`` turns into IHM categories and the
``_mdcif_simulation`` extension category.
"""

import dataclasses
from dataclasses import dataclass, field


@dataclass
class SimulationInfo:
    engine: str = ""
    engine_version: str = ""
    force_field: str = ""
    water_model: str = ""
    ensemble: str = ""
    temperature_K: float | None = None
    pressure_bar: float | None = None
    timestep_fs: float | None = None
    length_ns: float | None = None
    n_replicas: int = 1
    enhanced_sampling: str = "none"
    save_interval_ps: float | None = None
    details: str = ""


@dataclass
class ClusteringInfo:
    software: str = ""
    method: str = ""
    feature: str = "RMSD"
    selection: str = "name CA"
    n_clusters: int | None = None
    cutoff_angstrom: float | None = None
    n_models_deposited: int | None = None


@dataclass
class ModeAnalysisInfo:
    method: str = "PCA"
    selection: str = "name CA"
    mass_weighted: bool = False
    n_modes: int = 3
    per_cluster: bool = True


@dataclass
class Provenance:
    title: str = ""
    authors: list[str] = field(default_factory=list)
    trajectory_doi: str = ""
    trajectory_url: str = ""
    trajectory_database: str = ""   # e.g. "MDposit" when there is no DOI
    trajectory_accession: str = ""  # accession in that database
    name: str = ""
    source_structure: str = ""      # PDB ID the simulation started from
    reference_structure: str = ""   # experimental PDB entry used as reference (see reference.py)
    composition: str = ""
    frames_used: str = ""
    alignment_selection: str = "name CA"
    simulation: SimulationInfo = field(default_factory=SimulationInfo)
    clustering: ClusteringInfo = field(default_factory=ClusteringInfo)
    modes: ModeAnalysisInfo = field(default_factory=ModeAnalysisInfo)


def _fill(cls, data):
    """Instantiate dataclass ``cls`` from ``data``, ignoring unknown keys and nulls."""
    names = {f.name for f in dataclasses.fields(cls)}
    return cls(**{k: v for k, v in (data or {}).items() if k in names and v is not None})


def provenance_from_dict(data) -> Provenance:
    """Build a ``Provenance`` from the nested dict layout of the YAML template."""
    data = data or {}
    dataset = data.get("dataset") or {}
    system = data.get("system") or {}
    analysis = data.get("analysis") or {}
    top = {**dataset, **system, **analysis}
    prov = _fill(Provenance, top)
    prov.authors = list(dataset.get("authors") or [])
    prov.simulation = _fill(SimulationInfo, data.get("simulation"))
    prov.clustering = _fill(ClusteringInfo, data.get("clustering"))
    prov.modes = _fill(ModeAnalysisInfo, data.get("modes"))
    return prov


def provenance_to_dict(prov: Provenance) -> dict:
    """Inverse of :func:`provenance_from_dict` (the YAML template layout)."""
    return {
        "dataset": {"title": prov.title, "authors": list(prov.authors),
                    "trajectory_doi": prov.trajectory_doi, "trajectory_url": prov.trajectory_url,
                    "trajectory_database": prov.trajectory_database,
                    "trajectory_accession": prov.trajectory_accession},
        "system": {"name": prov.name, "source_structure": prov.source_structure,
                   "reference_structure": prov.reference_structure,
                   "composition": prov.composition},
        "simulation": dataclasses.asdict(prov.simulation),
        "analysis": {"frames_used": prov.frames_used,
                     "alignment_selection": prov.alignment_selection},
        "clustering": dataclasses.asdict(prov.clustering),
        "modes": dataclasses.asdict(prov.modes),
    }


def load_provenance(path) -> Provenance:
    """Read a provenance YAML file into a ``Provenance`` object."""
    import yaml

    with open(path, encoding="utf-8") as fh:
        return provenance_from_dict(yaml.safe_load(fh))
