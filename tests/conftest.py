"""
pytest configuration and shared fixtures for PyMOL testing.
Enforces OSMesa headless environment and provides synthetic test assets.
"""

import os
import pytest

# Ensure headless environment variable is set before any PyMOL imports occur
os.environ["PYOPENGL_PLATFORM"] = "osmesa"

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
ASSETS_DIR = os.path.join(ROOT_DIR, "assets")


@pytest.fixture(scope="session")
def pymol_cmd():
    """
    Session-scoped PyMOL initialization fixture.
    Initializes headless PyMOL once and yields the `cmd` module.
    """
    import pymol
    pymol.pymol_argv = ["pymol", "-cq"]
    pymol.finish_launching()
    from pymol import cmd, util

    yield cmd

    try:
        cmd.quit()
    except Exception:
        pass


@pytest.fixture(autouse=True)
def clean_pymol_session(pymol_cmd):
    """
    Autouse fixture that resets the PyMOL session between tests to ensure test isolation.
    """
    pymol_cmd.reinitialize()
    yield
    pymol_cmd.delete("all")


@pytest.fixture
def synthetic_monomer_path():
    return os.path.join(ASSETS_DIR, "synthetic_monomer.pdb")


@pytest.fixture
def synthetic_complex_cif_path():
    return os.path.join(ASSETS_DIR, "synthetic_complex.cif")


@pytest.fixture
def synthetic_ligand_pocket_path():
    return os.path.join(ASSETS_DIR, "synthetic_ligand_pocket.pdb")


@pytest.fixture
def synthetic_nmr_ensemble_path():
    return os.path.join(ASSETS_DIR, "synthetic_nmr_ensemble.pdb")


@pytest.fixture
def synthetic_homolog_path():
    return os.path.join(ASSETS_DIR, "synthetic_homolog.pdb")
