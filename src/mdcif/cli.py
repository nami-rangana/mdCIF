"""Command-line interface: ``mdcif build | validate | info | landscape``."""

from pathlib import Path

import typer

app = typer.Typer(help="Build, validate and inspect mdcif files.", no_args_is_help=True)


def _read_labels(path):
    """One integer label per analysed frame (-1 = noise): .npy, or text with one per line."""
    import numpy as np

    labels = np.load(path) if path.suffix == ".npy" else np.loadtxt(path, dtype=float)
    labels = np.asarray(labels).ravel()
    if not np.all(labels == np.round(labels)):
        raise typer.BadParameter(f"{path}: labels must be integers (-1 = noise)")
    return labels.astype(int)


@app.command()
def build(
    topology: Path,
    trajectory: list[Path] = typer.Argument(..., help="One or more trajectory files (replicas)."),
    provenance: Path = typer.Option(..., help="Provenance YAML (see examples/)."),
    output: Path = typer.Option(Path("output.cif"), help="Output .cif or .bcif file."),
    labels: Path = typer.Option(
        None, help="Your own clustering instead of HDBSCAN: one integer label per analysed "
                   "frame (after --start/--stop/--step), -1 = noise; .npy or text."),
    labels_method: str = typer.Option(
        None, help="Name of the algorithm behind --labels (default: provenance "
                   "clustering.method)."),
    n_pcs: int = typer.Option(3, help="Number of global PCs HDBSCAN clusters on."),
    features: str = typer.Option("pca", help="What HDBSCAN clusters: pca | cartesian."),
    min_cluster_size: int = typer.Option(
        None, help="HDBSCAN: smallest cluster, in frames (default: 1% of the frames, "
                   "at least 25)."),
    min_samples: int = typer.Option(
        None, help="HDBSCAN: density smoothing (default: min_cluster_size)."),
    cluster_selection: str = typer.Option(
        "eom", help="HDBSCAN: eom (fewer, larger clusters) | leaf (finer substates)."),
    max_models: int = typer.Option(None, help="Deposit only the N most populated clusters."),
    n_modes: int = typer.Option(None, help="Essential modes per model (default: provenance, 3)."),
    selection: str = typer.Option(None, help="Analysis atoms (default: provenance, 'name CA')."),
    deposit_selection: str = typer.Option(
        "protein and not name H* and not resname ACE NME NHE NH2",
        help="Atoms written for each model."),
    start: int = typer.Option(0, help="First frame used (0-based, e.g. to skip equilibration)."),
    stop: int = typer.Option(None, help="Stop before this frame."),
    step: int = typer.Option(1, help="Use every step-th frame."),
    traj_format: str = typer.Option(
        None, "--format", help="MDAnalysis trajectory format, if it can't be told from the "
                               "file (NetCDF, DCD, XTC and TRR are recognised by their "
                               "contents, e.g. meld.nc.001)."),
    reference: str = typer.Option(
        None, help="Experimental structure of the same molecule: PDB ID (downloaded from "
                   "RCSB) or PDBx/mmCIF file. Its entity, chains, source and citations are "
                   "reused, the models are superposed on it and their RMSD is recorded."),
    max_projection_points: int = typer.Option(
        5000, help="Frames stored in _mdcif_projection (regular subsample; 0 = none)."),
    mass_weighted: bool = typer.Option(None, help="Quasi-harmonic (mass-weighted) modes."),
    anm: bool = typer.Option(False, help="Report the cumulative overlap of each cluster's "
                                         "modes with the 10 lowest ANM modes of its medoid."),
) -> None:
    """Cluster a trajectory, compute global and per-cluster modes and write an mdcif file."""
    from .pipeline import build_mdcif
    from .provenance import load_provenance

    prov = load_provenance(provenance)
    summary = build_mdcif(topology, trajectory, prov, output,
                          labels=_read_labels(labels) if labels else None, method=labels_method,
                          n_pcs=n_pcs, features=features, min_cluster_size=min_cluster_size,
                          min_samples=min_samples, cluster_selection=cluster_selection,
                          selection=selection, deposit_selection=deposit_selection, start=start,
                          stop=stop, step=step, traj_format=traj_format, n_modes=n_modes,
                          mass_weighted=mass_weighted, max_models=max_models,
                          reference=reference or prov.reference_structure or None,
                          max_projection_points=max_projection_points, anm=anm)
    typer.echo(summary.report())


@app.command()
def validate(path: Path) -> None:
    """Validate a file against the PDBx, IHM and mdcif dictionaries."""
    from .validate import MdcifValidationError, validate_mdcif

    try:
        validate_mdcif(path)
    except MdcifValidationError as exc:
        typer.echo(f"{path}: {len(exc.problems)} problem(s)")
        for p in exc.problems:
            typer.echo(f"  {p}")
        raise typer.Exit(1) from None
    typer.echo(f"{path}: valid (PDBx + IHM + mdcif_ext, mode checks)")


def _var(ms):
    return "/".join(f"{100 * v:.0f}" for v in ms.variance_fraction)


@app.command()
def info(path: Path) -> None:
    """Summarize models, provenance and modes stored in an mdcif file."""
    from .reader import read_mdcif

    data = read_mdcif(path)
    system = data.system
    typer.echo(f"{path}: {system.title}")
    if data.simulation:
        typer.echo("  simulation: " + ", ".join(f"{k}={v}" for k, v in data.simulation.items()
                                                 if k != "details"))
    g = data.global_modes
    if g is not None:
        typer.echo(f"  global: {g.n_modes} {g.method} modes of {g.n_frames} frames on "
                   f"{len(g.atom_indices)} atoms, var% {_var(g)}")
    p = data.projection
    if p is not None:
        noise = (p["model_id"] == 0).mean()
        typer.echo(f"  projections: {len(p['frame'])} frames on global modes 1-3"
                   + (f", {100 * noise:.1f}% in no deposited cluster" if noise else ""))
    ens = {id(e.model_group): e for e in system.ensembles}
    for group, model in system._all_models():
        e = ens.get(id(group))
        n_atoms = sum(1 for _ in model.get_atoms())
        line = f"  model {model._id}: {model.name}, {n_atoms} atoms"
        if e is not None:
            line += f", {e.num_models} frames, precision {e.precision:.2f} A"
        ms = data.mode_sets.get(int(model._id))
        if ms is not None:
            line += (f", {ms.n_modes} {ms.method} modes on {len(ms.atom_indices)} atoms, "
                     f"var% {_var(ms)}")
        typer.echo(line)
        if group.details:
            typer.echo(f"      {group.details}")


@app.command()
def landscape(
    path: Path,
    output: Path = typer.Option(None, help="Image file (default: <file>_landscape.png)."),
    reference: str = typer.Option(
        None, help="Also project the models of a PDB entry (ID or mmCIF file) onto the PCs."),
) -> None:
    """Plot the frames on global PCs 1-3, coloured by cluster, from the file alone."""
    from .plot import plot_landscape

    out = output or path.with_name(path.stem + "_landscape.png")
    plot_landscape(path, out, reference=reference)
    typer.echo(f"wrote {out}")


if __name__ == "__main__":
    app()
