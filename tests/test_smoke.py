"""
test_smoke.py
Smoke tests verifying engine initialization, license compliance, and synthetic assets.
"""

import os
import pytest


def test_license_file_exists():
    """Verify compliance protocol license notice file exists and contains official URL."""
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    license_path = os.path.join(root_dir, ".licenses", "pymol_LICENSE.txt")
    assert os.path.exists(license_path), "License notice file .licenses/pymol_LICENSE.txt missing"

    with open(license_path, "r", encoding="utf-8") as f:
        content = f.read()
    assert "https://www.pymol.org/" in content
    assert "Notice Date" in content or "Acknowledged at" in content


def test_apache_license_file_exists():
    """Verify Apache License 2.0 file exists in root directory."""
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    license_path = os.path.join(root_dir, "LICENSE")
    assert os.path.exists(license_path), "LICENSE file in project root missing"

    with open(license_path, "r", encoding="utf-8") as f:
        content = f.read()
    assert "Apache License" in content
    assert "Version 2.0" in content


def test_osmesa_headless_environment():
    """Verify off-screen mesa rasterization variable is configured."""
    assert os.environ.get("PYOPENGL_PLATFORM") == "osmesa"


def test_pymol_cmd_initialization(pymol_cmd):
    """Verify PyMOL cmd API responds and version is queryable."""
    version_info = pymol_cmd.get_version()
    assert version_info is not None
    assert len(version_info) >= 2
    # PyMOL version string (e.g. '3.2.0a')
    assert isinstance(version_info[0], str)
    assert len(version_info[0]) > 0


def test_synthetic_assets_exist(synthetic_monomer_path, synthetic_complex_cif_path,
                                synthetic_ligand_pocket_path, synthetic_nmr_ensemble_path,
                                synthetic_homolog_path):
    """Verify all synthetic coordinate files are present and non-empty."""
    for path in [synthetic_monomer_path, synthetic_complex_cif_path,
                 synthetic_ligand_pocket_path, synthetic_nmr_ensemble_path,
                 synthetic_homolog_path]:
        assert os.path.exists(path), f"Asset file missing: {path}"
        assert os.path.getsize(path) > 0, f"Asset file is empty: {path}"


def test_readme_image_asset_exists():
    """Verify README workflow diagram image exists and has valid PNG header."""
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    img_path = os.path.join(root_dir, "assets", "skill_workflow_diagram.png")
    assert os.path.exists(img_path), "README workflow diagram missing from assets directory"
    assert os.path.getsize(img_path) > 0, "Workflow diagram image file is empty"
    with open(img_path, "rb") as f:
        header = f.read(8)
    # PNG signature: 89 50 4E 47 0D 0A 1A 0A
    assert header == b"\x89PNG\r\n\x1a\n", "Image file is not a valid PNG format"


def test_skill_directory_and_manifest_structure():
    """Verify skills/pymol-visualize directory, SKILL.md, and references exist."""
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    skill_dir = os.path.join(root_dir, "skills", "pymol-visualize")
    assert os.path.isdir(skill_dir), "skills/pymol-visualize directory missing"

    skill_md = os.path.join(skill_dir, "SKILL.md")
    assert os.path.exists(skill_md), "skills/pymol-visualize/SKILL.md missing"
    with open(skill_md, "r", encoding="utf-8") as f:
        content = f.read()
    assert "name: pymol-visualize" in content, "SKILL.md missing 'name: pymol-visualize'"

    for ref_file in ["PYMOL_REFERENCE.md", "RECIPES.md", "citation.bib"]:
        ref_path = os.path.join(skill_dir, "references", ref_file)
        assert os.path.exists(ref_path), f"skills/pymol-visualize/references/{ref_file} missing"


def test_pytest_ini_exists():
    """Verify pytest.ini exists and configures warning filters."""
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    ini_path = os.path.join(root_dir, "pytest.ini")
    assert os.path.exists(ini_path), "pytest.ini missing from root directory"
    with open(ini_path, "r", encoding="utf-8") as f:
        content = f.read()
    assert "filterwarnings" in content, "pytest.ini missing filterwarnings"


