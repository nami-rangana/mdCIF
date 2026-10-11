import importlib

import pytest


def test_package_imports():
    import mdcif

    assert mdcif.__version__


@pytest.mark.parametrize(
    "name",
    ["numpy", "scipy", "MDAnalysis", "sklearn.cluster", "ihm", "gemmi", "yaml", "typer"],
)
def test_dependency_imports(name):
    importlib.import_module(name)


def test_blas_works():
    # Catches the conda-forge MKL / Intel-OpenMP crash (needs MKL_THREADING_LAYER=TBB).
    import numpy as np

    a = np.arange(9.0).reshape(3, 3)
    assert np.allclose(a @ np.eye(3), a)


def test_adk_fixture(adk_universe):
    assert adk_universe.trajectory.n_frames > 1
    assert len(adk_universe.select_atoms("name CA")) == 214
