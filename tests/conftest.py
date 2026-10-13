import pytest


@pytest.fixture(scope="session")
def adk_files():
    """(topology, trajectory) of the AdK closed -> open (DIMS) run shipped with MDAnalysisTests."""
    datafiles = pytest.importorskip("MDAnalysisTests.datafiles")
    return datafiles.PSF, datafiles.DCD


@pytest.fixture(scope="session")
def adk_universe(adk_files):
    """AdK closed -> open transition (DIMS) trajectory shipped with MDAnalysisTests."""
    mda = pytest.importorskip("MDAnalysis")
    return mda.Universe(*adk_files)


@pytest.fixture(scope="session")
def adk_reference(tmp_path_factory):
    """PDBx/mmCIF reference entry for AdK (open form, 4AKE) made from adk_open.pdb with gemmi."""
    datafiles = pytest.importorskip("MDAnalysisTests.datafiles")
    gemmi = pytest.importorskip("gemmi")
    st = gemmi.read_structure(datafiles.PDB_small)
    st.setup_entities()
    st.remove_ligands_and_waters()
    st.remove_hydrogens()
    st[0][0].name = "A"
    for res in st[0]["A"]:  # adk_open.pdb uses CHARMM names; a PDB entry says HIS
        if res.name in ("HSD", "HSE", "HSP"):
            res.name = "HIS"
    st.setup_entities()
    for ent in st.entities:
        ent.full_sequence = [r.name for r in st[0]["A"]]
    st.assign_label_seq_id()
    path = tmp_path_factory.mktemp("reference") / "adk_open.cif"
    st.make_mmcif_document().write_file(str(path))
    return path


@pytest.fixture(scope="session")
def adk_mdcif(adk_files, tmp_path_factory):
    """mdcif file built from AdK with the default HDBSCAN clustering (and its BuildSummary)."""
    from mdcif.pipeline import build_mdcif
    from mdcif.provenance import Provenance

    out = tmp_path_factory.mktemp("mdcif") / "adk_hdbscan.cif"
    # 98 frames: too few for the default (at least 25 frames per cluster) to find any
    summary = build_mdcif(*adk_files, Provenance(name="AdK test"), out, anm=False,
                          min_cluster_size=5)
    return out, summary


@pytest.fixture(scope="session")
def adk_mdcif_mw(adk_files, tmp_path_factory):
    """mdcif file from AdK with mass-weighted (quasi-harmonic) cluster modes, 2 clusters by hand."""
    import numpy as np

    from mdcif.pipeline import build_mdcif
    from mdcif.provenance import Provenance

    out = tmp_path_factory.mktemp("mdcif") / "adk_qha.cif"
    labels = (np.arange(98) >= 60).astype(int)  # closed half / open half of the transition
    summary = build_mdcif(*adk_files, Provenance(name="AdK test"), out, labels=labels,
                          method="by hand", mass_weighted=True)
    return out, summary
