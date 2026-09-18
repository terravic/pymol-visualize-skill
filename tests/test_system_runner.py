"""
test_system_runner.py
End-to-end system tests verifying CLI interfaces, subprocess isolation,
and error-handling workflows for scripts/pymol_runner.py and scripts/validate_structure.py.
"""

import json
import os
import subprocess
import sys
import pytest

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCRIPTS_DIR = os.path.join(ROOT_DIR, "scripts")
ASSETS_DIR = os.path.join(ROOT_DIR, "assets")


def test_cli_validate_structure_pdb():
    """Verify validate_structure CLI successfully parses PDB file."""
    script = os.path.join(SCRIPTS_DIR, "validate_structure.py")
    pdb_path = os.path.join(ASSETS_DIR, "synthetic_monomer.pdb")

    proc = subprocess.run([sys.executable, script, pdb_path, "--json"], capture_output=True, text=True)
    assert proc.returncode == 0
    data = json.loads(proc.stdout)
    assert data["valid"] is True
    assert data["format"] == "PDB"
    assert data["total_atoms"] == 100
    assert data["chains"] == ["A"]


def test_cli_validate_structure_cif():
    """Verify validate_structure CLI successfully parses mmCIF file with multi-chain complex."""
    script = os.path.join(SCRIPTS_DIR, "validate_structure.py")
    cif_path = os.path.join(ASSETS_DIR, "synthetic_complex.cif")

    proc = subprocess.run([sys.executable, script, cif_path, "--json"], capture_output=True, text=True)
    assert proc.returncode == 0
    data = json.loads(proc.stdout)
    assert data["valid"] is True
    assert data["format"] == "mmCIF"
    assert data["total_atoms"] == 100
    assert data["chains"] == ["A", "B"]


def test_cli_validate_structure_nmr_ensemble():
    """Verify validate_structure CLI detects multi-model NMR ensemble."""
    script = os.path.join(SCRIPTS_DIR, "validate_structure.py")
    nmr_path = os.path.join(ASSETS_DIR, "synthetic_nmr_ensemble.pdb")

    proc = subprocess.run([sys.executable, script, nmr_path, "--json"], capture_output=True, text=True)
    assert proc.returncode == 0
    data = json.loads(proc.stdout)
    assert data["is_multi_state"] is True
    assert data["model_count"] == 2


def test_cli_validate_structure_nonexistent_file():
    """Verify validate_structure CLI exits non-zero on non-existent file."""
    script = os.path.join(SCRIPTS_DIR, "validate_structure.py")
    proc = subprocess.run([sys.executable, script, "nonexistent.pdb"], capture_output=True, text=True)
    assert proc.returncode != 0
    assert "not found" in proc.stderr.lower()


def test_cli_validate_structure_empty_file(tmp_path):
    """Verify validate_structure CLI exits non-zero on empty coordinate file."""
    empty_file = str(tmp_path / "empty.pdb")
    open(empty_file, "w").close()

    script = os.path.join(SCRIPTS_DIR, "validate_structure.py")
    proc = subprocess.run([sys.executable, script, empty_file], capture_output=True, text=True)
    assert proc.returncode != 0


def test_cli_pymol_runner_execution(tmp_path):
    """Verify pymol_runner CLI renders recipe and generates PNG and PSE outputs."""
    script = os.path.join(SCRIPTS_DIR, "pymol_runner.py")
    pdb_path = os.path.join(ASSETS_DIR, "synthetic_monomer.pdb")
    out_png = str(tmp_path / "cli_render.png")
    out_pse = str(tmp_path / "cli_session.pse")

    cmd = [
        sys.executable, script,
        "--recipe", "alphafold_plddt",
        "--structure", pdb_path,
        "--out-png", out_png,
        "--out-pse", out_pse,
        "--json",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, f"STDOUT: {proc.stdout}\nSTDERR: {proc.stderr}"

    data = json.loads(proc.stdout)
    assert data["success"] is True
    assert os.path.exists(out_png)
    assert os.path.exists(out_pse)


def test_cli_pymol_runner_script_generation(tmp_path):
    """Verify pymol_runner generates valid standalone PEP 723 script."""
    script = os.path.join(SCRIPTS_DIR, "pymol_runner.py")
    pdb_path = os.path.join(ASSETS_DIR, "synthetic_monomer.pdb")
    out_script = str(tmp_path / "standalone_secstruct.py")
    out_png = str(tmp_path / "gen_render.png")
    out_pse = str(tmp_path / "gen_session.pse")

    gen_cmd = [
        sys.executable, script,
        "--recipe", "secondary_structure",
        "--structure", pdb_path,
        "--out-png", out_png,
        "--out-pse", out_pse,
        "--generate-script", out_script,
    ]
    proc = subprocess.run(gen_cmd, capture_output=True, text=True)
    assert proc.returncode == 0
    assert os.path.exists(out_script)

    # Check PEP 723 header present
    with open(out_script, "r", encoding="utf-8") as f:
        script_content = f.read()
    assert "# /// script" in script_content
    assert "pymol-open-source-whl" in script_content
    assert "cmd.quit()" in script_content

    # Execute generated standalone script via uv run
    run_cmd = ["uv", "run", out_script]
    run_proc = subprocess.run(run_cmd, capture_output=True, text=True)
    assert run_proc.returncode == 0, f"STDOUT: {run_proc.stdout}\nSTDERR: {run_proc.stderr}"
    assert os.path.exists(out_png)
    assert os.path.exists(out_pse)


def test_cli_pymol_runner_missing_structure(tmp_path):
    """Verify pymol_runner returns non-zero error when structure is absent."""
    script = os.path.join(SCRIPTS_DIR, "pymol_runner.py")
    cmd = [
        sys.executable, script,
        "--recipe", "secondary_structure",
        "--structure", "nonexistent_file.pdb",
        "--out-png", str(tmp_path / "out.png"),
        "--out-pse", str(tmp_path / "out.pse"),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode != 0


def test_cli_validate_structure_gzipped(tmp_path):
    """Verify validate_structure CLI handles gzip-compressed (.pdb.gz) structures."""
    import gzip
    script = os.path.join(SCRIPTS_DIR, "validate_structure.py")
    src_pdb = os.path.join(ASSETS_DIR, "synthetic_monomer.pdb")
    gz_pdb = str(tmp_path / "test_compressed.pdb.gz")

    with open(src_pdb, "rb") as f_in, gzip.open(gz_pdb, "wb") as f_out:
        f_out.writelines(f_in)

    proc = subprocess.run([sys.executable, script, gz_pdb, "--json"], capture_output=True, text=True)
    assert proc.returncode == 0
    data = json.loads(proc.stdout)
    assert data["valid"] is True
    assert data["is_compressed"] is True
    assert data["total_atoms"] == 100


def test_cli_validate_structure_ligand_detection():
    """Verify validate_structure CLI detects small-molecule organic ligands."""
    script = os.path.join(SCRIPTS_DIR, "validate_structure.py")
    pocket_path = os.path.join(ASSETS_DIR, "synthetic_ligand_pocket.pdb")

    proc = subprocess.run([sys.executable, script, pocket_path, "--json"], capture_output=True, text=True)
    assert proc.returncode == 0
    data = json.loads(proc.stdout)
    assert "LIG" in data["detected_ligands"]


def test_cli_fetch_structure_invalid_id():
    """Verify fetch_structure CLI safely handles invalid PDB ID length."""
    script = os.path.join(SCRIPTS_DIR, "fetch_structure.py")
    proc = subprocess.run([sys.executable, script, "--pdb", "TOO_LONG_PDB_ID"], capture_output=True, text=True)
    assert proc.returncode != 0
    assert "Invalid PDB accession code" in proc.stderr or "Invalid PDB accession code" in proc.stdout


def test_cli_pymol_runner_gzipped_input(tmp_path):
    """Verify pymol_runner executes cleanly when given gzip-compressed coordinates."""
    import gzip
    script = os.path.join(SCRIPTS_DIR, "pymol_runner.py")
    src_pdb = os.path.join(ASSETS_DIR, "synthetic_monomer.pdb")
    gz_pdb = str(tmp_path / "monomer.pdb.gz")
    out_png = str(tmp_path / "gz_render.png")
    out_pse = str(tmp_path / "gz_session.pse")

    with open(src_pdb, "rb") as f_in, gzip.open(gz_pdb, "wb") as f_out:
        f_out.writelines(f_in)

    cmd = [
        sys.executable, script,
        "--recipe", "secondary_structure",
        "--structure", gz_pdb,
        "--out-png", out_png,
        "--out-pse", out_pse,
        "--json",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, f"STDOUT: {proc.stdout}\nSTDERR: {proc.stderr}"
    assert os.path.exists(out_png)
    assert os.path.exists(out_pse)


def test_cli_generate_dashboard(tmp_path):
    """Verify generate_dashboard CLI produces a valid interactive HTML dashboard."""
    script = os.path.join(SCRIPTS_DIR, "generate_dashboard.py")
    pdb_path = os.path.join(ASSETS_DIR, "synthetic_monomer.pdb")
    out_html = str(tmp_path / "test_dashboard.html")

    cmd = [
        sys.executable, script,
        "--structure", pdb_path,
        "--recipe", "secondary_structure",
        "--output", out_html,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, f"STDOUT: {proc.stdout}\nSTDERR: {proc.stderr}"
    assert os.path.exists(out_html)

    with open(out_html, "r", encoding="utf-8") as f:
        content = f.read()

    assert "<!DOCTYPE html>" in content
    assert "themeToggle" in content
    assert "3Dmol" in content
    assert "synthetic_monomer.pdb" in content
    assert "100" in content  # Total atoms
    assert "Confidence Distribution" in content


def test_prebuilt_ui_dashboard_exists():
    """Verify that ui/dashboard.html is pre-built with synthetic data."""
    repo_root = os.path.dirname(SCRIPTS_DIR)
    ui_dashboard = os.path.join(repo_root, "ui", "dashboard.html")
    assert os.path.exists(ui_dashboard), "ui/dashboard.html must be pre-built in repository"
    with open(ui_dashboard, "r", encoding="utf-8") as f:
        content = f.read()
    assert "<!DOCTYPE html>" in content
    assert "3Dmol" in content
    assert "synthetic_monomer.pdb" in content


def test_pymol_runner_dashboard_integration(tmp_path):
    """Verify that pymol_runner generates dashboard when --dashboard is passed."""
    script = os.path.join(SCRIPTS_DIR, "pymol_runner.py")
    pdb_path = os.path.join(ASSETS_DIR, "synthetic_monomer.pdb")
    out_png = str(tmp_path / "test_dash_render.png")
    out_pse = str(tmp_path / "test_dash_session.pse")
    out_dash = str(tmp_path / "test_runner_dash.html")

    cmd = [
        sys.executable, script,
        "--recipe", "secondary_structure",
        "--structure", pdb_path,
        "--out-png", out_png,
        "--out-pse", out_pse,
        "--dashboard", out_dash,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, f"STDOUT: {proc.stdout}\nSTDERR: {proc.stderr}"
    assert os.path.exists(out_dash)
    with open(out_dash, "r", encoding="utf-8") as f:
        content = f.read()
    assert "synthetic_monomer.pdb" in content
    assert "mol-container" in content


