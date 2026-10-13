import shutil
from pathlib import Path

import numpy as np
import pytest

from mdcif.modes import cluster_modes, cumulative_overlap, overlap_matrix, project
from mdcif.provenance import Provenance, provenance_from_dict, provenance_to_dict
from mdcif.trajectory import apply_transform, element_masses, guess_format, load_universe, superpose


def random_rotation(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q *= np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def test_superpose_recovers_rigid_motion():
    rng = np.random.default_rng(1)
    ref = rng.normal(size=(30, 3)) * 5
    frames = np.array([ref @ random_rotation(rng) + rng.normal(size=3) * 10 for _ in range(4)])
    aligned, rot, centres, ref_c = superpose(frames, ref)
    assert np.allclose(aligned, ref, atol=1e-8)
    # the per-frame transform moves other atoms of that frame consistently
    assert np.allclose(apply_transform(frames[2], rot[2], centres[2], ref_c), ref, atol=1e-8)


def test_guess_format_from_contents(tmp_path):
    datafiles = pytest.importorskip("MDAnalysisTests.datafiles")
    for source, name, fmt in [(datafiles.NCDF, "meld.nc.001", "NC"),
                              (datafiles.DCD, "run.v2.001", "DCD"),
                              (datafiles.XTC, "#traj.xtc.1#", "XTC"),
                              (datafiles.TRR, "x.trr.002", "TRR")]:
        path = tmp_path / name
        path.write_bytes(Path(source).read_bytes()[:64])
        assert guess_format(path) == fmt, name


def test_guess_format_from_name(tmp_path):
    pytest.importorskip("MDAnalysis")
    for name, fmt in [("traj.mdcrd.001", "MDCRD"), ("run.v2.001", None), ("prod.1.001", None),
                      ("meld.001", None), ("traj.dcd", None)]:
        path = tmp_path / name
        path.write_text("title\n")                    # no binary signature
        assert guess_format(path) == fmt, name


def test_load_universe_with_numeric_suffix(adk_files, tmp_path):
    topology, dcd = adk_files
    copy = tmp_path / "adk.dcd.001"
    shutil.copy(dcd, copy)
    assert load_universe(topology, copy).trajectory.n_frames == 98
    unknown = tmp_path / "run.v2.001"
    unknown.write_text("title\n")
    with pytest.raises(ValueError, match="trajectory format"):
        load_universe(topology, unknown)


def test_cluster_modes_finds_dominant_direction():
    rng = np.random.default_rng(2)
    base = rng.normal(size=(20, 3)) * 5
    direction = rng.normal(size=(20, 3))
    direction /= np.linalg.norm(direction)
    amp = rng.normal(size=500) * 3.0
    coords = base + amp[:, None, None] * direction + rng.normal(size=(500, 20, 3)) * 0.05
    ms = cluster_modes(coords, n_modes=3)
    assert ms.eigenvectors.shape == (3, 20, 3)
    assert np.allclose(np.linalg.norm(ms.eigenvectors.reshape(3, -1), axis=1), 1.0)
    assert np.all(np.diff(ms.eigenvalues) <= 0)
    assert abs(np.sum(ms.eigenvectors[0] * direction)) > 0.99
    assert ms.variance_fraction[0] > 0.9
    assert overlap_matrix(ms, ms)[0, 0] == pytest.approx(1.0)
    assert cumulative_overlap(ms, ms) == pytest.approx(np.ones(3))
    # projections recover the amplitudes along the dominant direction (up to sign)
    assert ms.origin.shape == (20, 3) and ms.total_variance == pytest.approx(
        ms.eigenvalues[0] / ms.variance_fraction[0])
    p = project(coords, ms)[:, 0]
    assert abs(np.corrcoef(p, amp)[0, 1]) > 0.99
    assert np.var(p, ddof=1) == pytest.approx(ms.eigenvalues[0])


def test_mass_weighted_vectors_are_cartesian_unit_vectors():
    rng = np.random.default_rng(3)
    coords = rng.normal(size=(200, 10, 3))
    masses = rng.uniform(1, 30, size=10)
    ms = cluster_modes(coords, n_modes=2, masses=masses)
    assert ms.method == "quasi-harmonic" and ms.mass_weighted
    assert np.allclose(np.linalg.norm(ms.eigenvectors.reshape(2, -1), axis=1), 1.0)
    assert np.allclose(ms.origin, coords.mean(axis=0))
    # eigenvalues in amu A^2; the Cartesian variance along each unit vector in A^2
    assert np.allclose(ms.cartesian_variance,
                       ms.eigenvalues / ((ms.eigenvectors ** 2).sum(axis=2) @ masses))
    with pytest.raises(ValueError):
        project(coords, ms)


def test_uniform_masses_only_rescale_eigenvalues():
    """CA-only modes: every weight is the same, so quasi-harmonic = PCA up to units."""
    rng = np.random.default_rng(4)
    coords = rng.normal(size=(150, 12, 3)) * rng.uniform(0.5, 3, size=(12, 3))
    pca = cluster_modes(coords, n_modes=3)
    qha = cluster_modes(coords, n_modes=3, masses=np.full(12, 12.011))
    assert np.allclose(qha.eigenvectors, pca.eigenvectors, atol=1e-8)
    assert np.allclose(qha.eigenvalues, 12.011 * pca.eigenvalues)
    assert np.allclose(qha.cartesian_variance, pca.eigenvalues)
    assert np.array_equal(pca.cartesian_variance, pca.eigenvalues)


def test_mass_weighted_cartesian_variance_is_the_cartesian_amplitude():
    """One motion of known size: the masses come out again in cartesian_variance."""
    rng = np.random.default_rng(5)
    mean = rng.normal(size=(8, 3)) * 5
    d = rng.normal(size=24)
    d /= np.linalg.norm(d)
    a = rng.normal(size=600) * 2.0
    coords = mean + (a[:, None] * d).reshape(-1, 8, 3) + rng.normal(size=(600, 8, 3)) * 0.01
    ms = cluster_modes(coords, n_modes=1, masses=rng.uniform(1, 30, size=8))
    assert abs(ms.eigenvectors[0].ravel() @ d) == pytest.approx(1.0, abs=1e-4)
    assert ms.cartesian_variance[0] == pytest.approx(np.var(a, ddof=1), rel=1e-3)
    assert ms.eigenvalues[0] > 2 * ms.cartesian_variance[0]  # amu A^2, not A^2


def test_element_masses_ignore_repartitioned_topology_masses(adk_files):
    import MDAnalysis as mda

    u = mda.Universe(*adk_files)               # PSF: no elements, guessed from names
    ca = u.select_atoms("name CA")
    ca.masses = 9.994                          # as after hydrogen mass repartitioning
    assert np.allclose(element_masses(ca), 12.011)
    bb = u.select_atoms("resid 1 and name N CA C O")
    assert np.allclose(element_masses(bb), [14.007, 12.011, 12.011, 15.999])


@pytest.mark.parametrize("n_frames", [40, 400])  # fewer / more frames than 3N = 150
def test_cluster_modes_match_full_eigendecomposition(n_frames):
    rng = np.random.default_rng(7)
    basis = np.linalg.qr(rng.normal(size=(150, 5)))[0]
    X = rng.normal(size=(n_frames, 5)) * [5.0, 4.0, 3.0, 2.0, 1.5] @ basis.T
    X += rng.normal(size=X.shape) * 0.1
    ms = cluster_modes(X.reshape(n_frames, 50, 3), n_modes=5)
    evals, evecs = np.linalg.eigh(np.cov(X, rowvar=False))
    assert ms.eigenvalues == pytest.approx(evals[::-1][:5])
    assert ms.total_variance == pytest.approx(evals.sum())
    v = ms.eigenvectors.reshape(5, -1)
    assert np.allclose(np.abs(np.sum(v * evecs[:, ::-1][:, :5].T, axis=1)), 1.0)
    # sign convention: the largest component of each vector is positive
    assert np.all(v[np.arange(5), np.abs(v).argmax(axis=1)] > 0)


def test_hdbscan_separates_states_in_pc_space():
    from mdcif.clustering import clusters_from_labels, default_min_cluster_size, hdbscan_labels

    rng = np.random.default_rng(4)
    base = rng.normal(size=(15, 3)) * 5
    shift = rng.normal(size=(15, 3)) * 2.0
    a = base + rng.normal(size=(300, 15, 3)) * 0.2
    b = base + shift + rng.normal(size=(100, 15, 3)) * 0.2
    coords = np.concatenate([a, b])
    ms = cluster_modes(coords, n_modes=3)
    labels, params = hdbscan_labels(project(coords, ms))
    c = clusters_from_labels(coords, labels, method="HDBSCAN", feature="PCs 1-3")
    assert c.n_clusters == 2 and c.ihm_method == "Density based threshold-clustering"
    assert list(c.populations) == [300, 100]  # relabelled by population
    assert np.all(c.labels[:300] == 0) and np.all(c.labels[300:] == 1)
    assert params["min_cluster_size"] == 25            # 1% of 400 frames is below the floor
    assert default_min_cluster_size(10_000) == 100


def test_fast_medoid_matches_exact():
    from mdcif.clustering import medoid

    rng = np.random.default_rng(5)
    # two unequal blobs: the medoid is not simply the frame closest to the mean
    X = np.concatenate([rng.normal(size=(1800, 30)), rng.normal(size=(900, 30)) + 1.5])
    exact = medoid(X, 10, exact_max=len(X))
    assert medoid(X, 10, exact_max=100) == exact


def test_anm_hessian_has_rigid_body_null_space():
    from mdcif.modes import anm_hessian, anm_modes

    rng = np.random.default_rng(6)
    x = np.cumsum(rng.normal(size=(40, 3)) * 2.2, axis=0)  # chain-like, all within cutoff
    h = anm_hessian(x, cutoff=50.0)
    assert np.allclose(h, h.T)
    evals = np.linalg.eigvalsh(h)
    assert np.sum(np.abs(evals) < 1e-8 * evals.max()) == 6 and evals.min() > -1e-8
    ms = anm_modes(x, n_modes=4, cutoff=50.0)
    assert np.all(np.diff(ms.eigenvalues) <= 0)
    # a pure translation is a zero mode, so it is orthogonal to the ANM modes
    t = np.tile([1.0, 0, 0], 40) / np.sqrt(40)
    assert np.allclose(ms.eigenvectors.reshape(4, -1) @ t, 0, atol=1e-8)


def test_anm_modes_keep_soft_modes_of_extended_chain():
    from mdcif.modes import anm_hessian, anm_modes

    # stiff, extended chain: its softest real modes are ~1e-7 of the largest eigenvalue
    rng = np.random.default_rng(0)
    d, x = np.array([1.0, 0, 0]), [np.zeros(3)]
    for _ in range(149):
        d = d + rng.normal(size=3) * 0.8
        d /= np.linalg.norm(d)
        x.append(x[-1] + 3.8 * d)
    x = np.array(x)
    evals = np.linalg.eigvalsh(anm_hessian(x))
    assert evals[6] < 1e-6 * evals[-1]  # below the old cut
    ms = anm_modes(x, n_modes=3)
    assert ms.eigenvalues == pytest.approx(1 / evals[6:9], rel=1e-6)
    assert ms.variance_fraction == pytest.approx(ms.eigenvalues / np.sum(1 / evals[6:]))


def test_provenance_dict_round_trip():
    prov = provenance_from_dict({
        "dataset": {"title": "t", "authors": ["A, B."], "trajectory_accession": "MD-X"},
        "system": {"name": "n", "source_structure": "1ABC", "reference_PDB_id": "2ABC"},
        "simulation": {"engine": "GROMACS", "temperature_K": 300, "unknown_key": 1},
        "modes": {"n_modes": 5},
    })
    assert prov.simulation.temperature_K == 300 and prov.modes.n_modes == 5
    assert prov.reference_PDB_id == "2ABC"
    assert provenance_from_dict(provenance_to_dict(prov)) == prov


def test_hdbscan_build_stores_global_modes_and_projections(adk_mdcif):
    from mdcif.reader import read_mdcif
    from mdcif.validate import check_modes

    out, summary = adk_mdcif
    data = read_mdcif(out)
    assert summary.clusters.method == "HDBSCAN"
    models = [m for _, m in data.system._all_models()]
    assert len(models) == summary.n_models == summary.clusters.n_clusters
    g = data.global_modes
    assert g is not None and g.scope == "global" and g.n_modes == 3 and g.model_id == 1
    assert np.allclose(g.eigenvectors, summary.global_modes.eigenvectors, atol=1e-5)
    assert np.allclose(g.origin, summary.global_modes.origin, atol=1e-3)
    p = data.projection
    assert len(p["frame"]) == summary.n_frames == 98
    # the stored cluster of every frame matches the clustering (0 = noise)
    labels = summary.clusters.labels
    assert np.array_equal(p["model_id"], np.where(labels >= 0, labels + 1, 0))
    # projecting each deposited model with the stored origin and vectors puts it among the
    # frames of its own cluster
    from mdcif.plot import model_projections

    mp = model_projections(data)
    for m in range(1, len(models) + 1):
        own = p["coords"][p["model_id"] == m]
        d = np.linalg.norm(p["coords"] - mp[m], axis=1)
        assert np.min(np.linalg.norm(own - mp[m], axis=1)) == pytest.approx(d.min(), abs=0.05)
    assert check_modes(out) == []


@pytest.mark.parametrize("suffix", [".cif", ".bcif"])
def test_build_and_read_round_trip(adk_files, tmp_path, suffix):
    from mdcif.pipeline import build_mdcif
    from mdcif.reader import read_mdcif
    from mdcif.validate import check_modes

    out = tmp_path / f"adk{suffix}"
    prov = Provenance(name="AdK test", source_structure="4AKE")
    # labels from "any other tool": closed half / open half of the transition
    labels = (np.arange(98) >= 60).astype(int)
    summary = build_mdcif(*adk_files, prov, out, labels=labels, method="by hand", n_modes=3)
    assert summary.clusters.n_clusters == 2 and list(summary.clusters.populations) == [60, 38]
    assert summary.clusters.populations.sum() == summary.n_frames == 98
    assert summary.clusters.method == "by hand" and summary.clusters.ihm_method == "Other"

    data = read_mdcif(out)
    models = [m for _, m in data.system._all_models()]
    assert len(models) == 2
    assert all(sum(1 for _ in m.get_atoms()) == summary.n_deposited_atoms for m in models)
    assert sorted(data.mode_sets) == [1, 2]
    for c, model in enumerate(models):
        written, read = summary.mode_sets[c], data.mode_sets[int(model._id)]
        assert np.allclose(read.eigenvalues, written.eigenvalues, rtol=1e-5)
        assert np.allclose(read.eigenvectors, written.eigenvectors, atol=1e-5)
        # mode vectors point at the CA atoms of their own model
        names = [a.atom_id for a in model.get_atoms()]
        assert {names[i] for i in read.atom_indices} == {"CA"}
    assert data.global_modes is not None and data.projection is not None
    assert check_modes(out) == []


def test_mass_weighted_build_round_trip(adk_mdcif_mw):
    from mdcif.reader import read_mdcif
    from mdcif.validate import validate_mdcif

    out, summary = adk_mdcif_mw
    data = read_mdcif(out)
    assert not data.global_modes.mass_weighted   # the clustering axes stay plain PCA
    for c, (_, ms) in enumerate(sorted(data.mode_sets.items())):
        written = summary.mode_sets[c]
        assert ms.mass_weighted and ms.method == "quasi-harmonic"
        assert np.allclose(ms.cartesian_variance, written.cartesian_variance, rtol=1e-5)
        # CA atoms only, all 12.011 amu: the masses only change the units
        assert np.allclose(ms.cartesian_variance, ms.eigenvalues / 12.011, rtol=1e-5)
        assert "mass-weighted" in data.mode_set_info[c + 2]["details"]
    validate_mdcif(out)


@pytest.mark.parametrize("suffix", [".cif", ".bcif"])
def test_validate_against_dictionaries(adk_files, tmp_path, suffix):
    """Full validation: PDBx + IHM + mdcif_ext dictionaries (downloaded once) and mode checks."""
    from mdcif.pipeline import build_mdcif
    from mdcif.validate import validate_mdcif

    out = tmp_path / f"adk{suffix}"
    build_mdcif(*adk_files, Provenance(name="AdK test"), out, anm=False, min_cluster_size=5)
    validate_mdcif(out)  # raises MdcifValidationError listing every problem


def test_build_from_several_replicas(adk_files, tmp_path):
    from mdcif.pipeline import build_mdcif
    from mdcif.reader import read_mdcif

    top, dcd = adk_files
    out = tmp_path / "two_replicas.cif"
    summary = build_mdcif(top, [dcd, dcd], Provenance(name="AdK x2"), out, min_cluster_size=5)
    assert summary.n_frames == 196
    data = read_mdcif(out)
    assert all("replica" in m.name for _, m in data.system._all_models())


def test_frame_range_and_projection_subsample(adk_files, tmp_path):
    from mdcif.pipeline import build_mdcif
    from mdcif.reader import read_mdcif

    out = tmp_path / "range.cif"
    summary = build_mdcif(*adk_files, Provenance(name="AdK range"), out, start=10, stop=90,
                          step=2, max_projection_points=25, max_models=2, anm=True,
                          min_cluster_size=5)
    assert summary.n_frames == 40 and summary.n_models <= 2
    overlaps = [ov for ov in summary.anm_overlap if ov is not None]
    assert overlaps and all(len(ov) == 3 and np.all((ov > 0) & (ov <= 1 + 1e-9))
                            for ov in overlaps)
    assert "ANM cumulative overlap" in summary.report()
    data = read_mdcif(out)
    p = data.projection
    assert len(p["frame"]) == 25 and p["frame"][0] == 10 and p["frame"][-1] == 88
    assert set(np.unique(p["model_id"])) <= {0, 1, 2}
    assert len([m for _, m in data.system._all_models()]) == summary.n_models


def test_reference_entry_is_reused(adk_files, adk_reference, tmp_path):
    from mdcif.pipeline import build_mdcif
    from mdcif.reader import read_mdcif
    from mdcif.reference import load_reference
    from mdcif.trajectory import rmsd_to

    out = tmp_path / "with_reference.cif"
    summary = build_mdcif(*adk_files, Provenance(name="AdK ref"), out, anm=False,
                          reference=adk_reference, min_cluster_size=5)
    data = read_mdcif(out)
    ref = load_reference(adk_reference)
    # the MD chain (segid 4AKE) is mapped onto the reference chain and entity
    assert [a.id for a in data.system.asym_units] == [a.id for a in ref.system.asym_units]
    assert len(data.system.entities) == 1
    # model 1 is optimally superposed on the reference (the others move rigidly with it)
    ref_ca = np.array([[a.x, a.y, a.z] for a in ref.first_model().get_atoms()
                       if a.atom_id == "CA"])
    model = [m for _, m in data.system._all_models()][0]
    ca = np.array([[a.x, a.y, a.z] for a in model.get_atoms() if a.atom_id == "CA"])
    in_place = np.sqrt(((ca - ref_ca) ** 2).sum(axis=1).mean())
    assert in_place == pytest.approx(rmsd_to(ca[None], ref_ca)[0], abs=1e-3)
    assert summary.reference_rmsd[0] == pytest.approx(in_place, abs=1e-3)
    # the entry's citations are kept, none of them as this file's primary citation
    assert len(data.system.citations) == len(ref.system.citations)
    assert not any(c.is_primary for c in data.system.citations)
    # the validation step records the RMSDs
    steps = [s for p in data.system._all_protocols() for a in p.analyses for s in a.steps]
    assert any(s.type == "validation" and "RMSD" in (s.details or "") for s in steps)


def test_landscape_plot(adk_mdcif, adk_reference, tmp_path):
    pytest.importorskip("matplotlib")
    from mdcif.plot import plot_landscape

    out = plot_landscape(adk_mdcif[0], tmp_path / "landscape.png", reference=adk_reference)
    assert out.stat().st_size > 10_000


def test_reference_sequences_are_repaired_or_dropped():
    import ihm
    import ihm.reference

    from mdcif.reference import _repair_references

    full = "MEKEKTYKLILQQ"  # database sequence; the entity is residues 6-11
    entity = ihm.Entity("TYKLIL")
    # python-ihm pads with db_begin - 1 gaps even when the stored sequence starts at residue 1
    padded = ihm.reference.UniProtSequence("X_STRSG", "P00001", "-" * 5 + full)
    wrong = ihm.reference.UniProtSequence("Y_STRSG", "P00002", "AAAAAAAAAAAA")
    for ref in (padded, wrong):
        ref.alignments.append(ihm.reference.Alignment(db_begin=6, db_end=11))
    entity.references.extend([padded, wrong])
    with pytest.warns(UserWarning, match="dropped sequence reference"):
        _repair_references(entity)
    assert entity.references == [padded] and padded.sequence == full


def test_any_clustering_can_be_used(adk_files, tmp_path):
    from sklearn.cluster import AgglomerativeClustering

    from mdcif.pipeline import cluster, load_trajectory, save_mdcif
    from mdcif.reader import read_mdcif

    traj = load_trajectory(*adk_files)
    assert traj.projections.shape == (98, 3) and traj.n_frames == 98
    c = cluster(traj, AgglomerativeClustering(n_clusters=3))
    assert c.n_clusters == 3 and c.method == "AgglomerativeClustering"
    assert c.ihm_method == "Hierarchical" and c.software[0] == "sklearn"
    assert c.params["n_clusters"] == 3 and c.params["n_pcs"] == 3
    c2 = cluster(traj, labels=np.where(np.arange(98) < 49, 0, -1), method="half")
    assert c2.n_clusters == 1 and c2.noise_fraction == pytest.approx(0.5)
    with pytest.raises(ValueError, match="one label per analysed frame"):
        cluster(traj, labels=np.zeros(10, dtype=int))
    cart = cluster(traj, AgglomerativeClustering(n_clusters=2), features="cartesian")
    assert cart.ihm_feature == "RMSD" and "n_pcs" not in cart.params
    with pytest.raises(ValueError):
        cluster(traj, features="tica")
    out = tmp_path / "agglo.cif"
    save_mdcif(out, traj, c, Provenance(name="AdK agglomerative"))
    ens = read_mdcif(out).system.ensembles
    assert [e.clustering_method for e in ens] == ["Hierarchical"] * 3


def test_saving_twice_does_not_move_the_trajectory(adk_files, adk_reference, tmp_path):
    from mdcif.pipeline import cluster, load_trajectory, save_mdcif

    traj = load_trajectory(*adk_files)
    before = traj.global_modes.eigenvectors.copy()
    c = cluster(traj, min_cluster_size=5)
    a = save_mdcif(tmp_path / "a.cif", traj, c, Provenance(name="a"), reference=adk_reference)
    b = save_mdcif(tmp_path / "b.bcif", traj, c, Provenance(name="b"), reference=adk_reference)
    assert np.array_equal(traj.global_modes.eigenvectors, before)
    assert np.allclose(a.global_modes.eigenvectors, b.global_modes.eigenvectors)
    assert a.reference_rmsd == pytest.approx(b.reference_rmsd)


def test_cli_build_with_labels(adk_files, tmp_path):
    from typer.testing import CliRunner

    from mdcif.cli import app
    from mdcif.reader import read_mdcif

    prov = tmp_path / "prov.yaml"
    prov.write_text("system:\n  name: AdK cli\nclustering:\n  software: my-msm-tool\n"
                    "  method: PCCA+\n")
    labels = tmp_path / "labels.txt"
    np.savetxt(labels, (np.arange(98) // 33).astype(int), fmt="%d")
    out = tmp_path / "cli.cif"
    result = CliRunner().invoke(app, ["build", *map(str, adk_files), "--provenance", str(prov),
                                      "--labels", str(labels), "--output", str(out)])
    assert result.exit_code == 0, result.output
    assert "PCCA+" in result.output and "3 clusters" in result.output
    sw = [s.name for s in read_mdcif(out).system.software]
    assert "my-msm-tool" in sw
