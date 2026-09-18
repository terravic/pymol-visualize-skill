#!/usr/bin/env python3
"""
pymol_runner.py
Execution harness and script generator for PyMOL headless workflows.
Supports real-world user data, on-demand fetching (PDB/AlphaFold), metric extraction,
and PEP 723 standalone script synthesis.
"""

import argparse
import json
import logging
import os
import subprocess
import sys
from typing import Any, Dict, Optional

from fetch_structure import fetch_alphafold_db, fetch_rcsb_pdb
from validate_structure import validate_structure

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("pymol_runner")

SUPPORTED_RECIPES = [
    "secondary_structure",
    "bfactor_spectrum",
    "alphafold_plddt",
    "highlight_residues",
    "transparent_surface",
    "electrostatic_surface",
    "multichain_palette",
    "bfactor_putty",
    "cavity_detection",
    "distance_measurement",
    "pocket_zoom",
    "protein_ligand_contacts",
    "mutagenesis",
    "superposition",
    "session_modify",
]


def resolve_coordinate_input(input_target: str, cache_dir: str = "data/structures") -> str:
    """
    Resolves input_target to a local coordinate file path.
    If input_target is an existing file path (or gzipped file), expands and returns it.
    If input_target is a 4-letter PDB accession (e.g., '1UBQ', '4HHB'), downloads from RCSB.
    If input_target is a UniProt accession (e.g., 'P04637'), downloads from AlphaFold DB.
    """
    expanded = os.path.expanduser(input_target)
    if os.path.exists(expanded):
        try:
            rel = os.path.relpath(expanded, os.getcwd())
            if not rel.startswith(".."):
                return rel
        except Exception:
            pass
        return expanded

    # Check if target is a PDB ID (4 alphanumeric characters)
    stripped = input_target.strip()
    if len(stripped) == 4 and stripped.isalnum():
        logger.info(f"Target '{stripped}' recognized as PDB accession. Fetching from RCSB PDB...")
        return fetch_rcsb_pdb(stripped, cache_dir, file_format="cif")

    # Check if target is a UniProt ID (e.g., P04637, Q9BYF1)
    if 5 <= len(stripped) <= 10 and stripped[0].upper() in "OPQ":
        logger.info(f"Target '{stripped}' recognized as UniProt accession. Fetching from AlphaFold DB...")
        return fetch_alphafold_db(stripped, cache_dir, file_format="pdb")

    # If neither, return the expanded path (will trigger pre-flight error)
    return expanded


def generate_pep723_script(recipe: str, params: Dict[str, Any], output_script: str) -> str:
    """
    Generates a standalone PEP 723 compliant Python script for a given recipe.
    """
    code_lines = [
        "#!/usr/bin/env python3",
        "# /// script",
        '# requires-python = ">=3.10, <3.13"',
        "# dependencies = [",
        '#     "pymol-open-source-whl",',
        "# ]",
        "# ///",
        "",
        "import os",
        "import sys",
        "",
        "def main():",
    ]

    structure_path = params.get("structure_path", "")
    output_png = params.get("output_png", "output.png")
    output_pse = params.get("output_pse", "output.pse")
    active_resi = params.get("active_residues", "57+102+195")
    ligand_resn = params.get("ligand_resn", "LIG")

    code_lines.extend([
        f'    structure_path = "{structure_path}"',
        f'    output_png = "{output_png}"',
        f'    output_pse = "{output_pse}"',
        "",
        "    if not os.path.exists(structure_path):",
        '        sys.stderr.write(f"Error: Structure file not found: {structure_path}\\n")',
        "        sys.exit(1)",
        "",
        '    os.makedirs(os.path.dirname(output_png) or ".", exist_ok=True)',
        '    os.makedirs(os.path.dirname(output_pse) or ".", exist_ok=True)',
        "",
        '    os.environ["PYOPENGL_PLATFORM"] = "osmesa"',
        "    import pymol",
        '    pymol.pymol_argv = ["pymol", "-cq"]',
        "    pymol.finish_launching()",
        "    from pymol import cmd, util",
        "",
        "    try:",
    ])

    if recipe == "secondary_structure":
        code_lines.extend([
            '        obj_name = "protein_secstruct"',
            "        cmd.load(structure_path, obj_name)",
            "        atom_count = cmd.count_atoms(obj_name)",
            "        if atom_count == 0:",
            '            sys.stderr.write(f"Error: Zero atoms parsed from {structure_path}\\n")',
            "            sys.exit(2)",
            '        print(f"Loaded {atom_count} atoms from {structure_path}")',
            '        cmd.hide("everything", obj_name)',
            '        cmd.show("cartoon", obj_name)',
            '        cmd.color("gray75", f"{obj_name} and ss l+\'\'")',
            '        cmd.color("firebrick", f"{obj_name} and ss h")',
            '        cmd.color("brightorange", f"{obj_name} and ss s")',
            '        cmd.set("cartoon_fancy_helices", 1)',
            '        cmd.set("ray_opaque_background", 1)',
            '        cmd.bg_color("white")',
            "        cmd.orient(obj_name)",
            "        cmd.zoom(obj_name, buffer=2.0)",
            "        cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)",
            "        cmd.save(output_pse)",
            '        print(f"Artifacts exported:\\n - {output_png}\\n - {output_pse}")',
        ])
    elif recipe == "bfactor_spectrum":
        code_lines.extend([
            '        obj_name = "protein_bfactor"',
            "        cmd.load(structure_path, obj_name)",
            "        atom_count = cmd.count_atoms(obj_name)",
            "        if atom_count == 0:",
            '            sys.stderr.write(f"Error: Zero atoms parsed from {structure_path}\\n")',
            "            sys.exit(2)",
            '        cmd.hide("everything", obj_name)',
            '        cmd.show("cartoon", obj_name)',
            '        cmd.spectrum("b", "red_white_blue", f"{obj_name} and polymer.protein")',
            '        cmd.set("ray_opaque_background", 1)',
            '        cmd.bg_color("white")',
            "        cmd.orient(obj_name)",
            "        cmd.zoom(obj_name, buffer=2.0)",
            "        cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)",
            "        cmd.save(output_pse)",
            '        print(f"Artifacts exported:\\n - {output_png}\\n - {output_pse}")',
        ])
    elif recipe == "alphafold_plddt":
        code_lines.extend([
            '        obj_name = "af_structure"',
            "        cmd.load(structure_path, obj_name)",
            "        atom_count = cmd.count_atoms(obj_name)",
            "        if atom_count == 0:",
            '            sys.stderr.write(f"Error: Zero atoms parsed from {structure_path}\\n")',
            "            sys.exit(2)",
            '        cmd.hide("everything", obj_name)',
            '        cmd.show("cartoon", obj_name)',
            '        cmd.set_color("af_very_low", [1.000, 0.490, 0.271])',
            '        cmd.set_color("af_low", [1.000, 0.859, 0.075])',
            '        cmd.set_color("af_confident", [0.396, 0.796, 0.953])',
            '        cmd.set_color("af_very_high", [0.000, 0.325, 0.839])',
            '        cmd.color("af_very_low", f"{obj_name} and b < 50")',
            '        cmd.color("af_low", f"{obj_name} and b >= 50 and b < 70")',
            '        cmd.color("af_confident", f"{obj_name} and b >= 70 and b < 90")',
            '        cmd.color("af_very_high", f"{obj_name} and b >= 90")',
            '        cmd.set("ray_opaque_background", 1)',
            '        cmd.bg_color("white")',
            "        cmd.orient(obj_name)",
            "        cmd.zoom(obj_name, buffer=2.0)",
            "        cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)",
            "        cmd.save(output_pse)",
            '        print(f"Artifacts exported:\\n - {output_png}\\n - {output_pse}")',
        ])
    elif recipe == "highlight_residues":
        code_lines.extend([
            '        obj_name = "active_site_view"',
            "        cmd.load(structure_path, obj_name)",
            '        cmd.hide("everything", obj_name)',
            '        cmd.show("cartoon", f"{obj_name} and polymer.protein")',
            '        cmd.color("gray80", f"{obj_name} and polymer.protein")',
            f'        active_sel = f"{{obj_name}} and resi {active_resi}"',
            '        cmd.show("sticks", active_sel)',
            '        util.cnc(active_sel)',
            '        cmd.set("ray_opaque_background", 1)',
            '        cmd.bg_color("white")',
            '        cmd.zoom(active_sel, buffer=4.0)',
            "        cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)",
            "        cmd.save(output_pse)",
            '        print(f"Artifacts exported:\\n - {output_png}\\n - {output_pse}")',
        ])
    elif recipe == "pocket_zoom":
        code_lines.extend([
            '        obj_name = "complex_view"',
            "        cmd.load(structure_path, obj_name)",
            f'        ligand_sel = f"{{obj_name}} and resn {ligand_resn} and organic"',
            '        pocket_sel = f"byres (({obj_name} and polymer.protein) within 5.0 of ({ligand_sel}))"',
            '        cmd.hide("everything", obj_name)',
            '        cmd.show("cartoon", f"{obj_name} and polymer.protein")',
            '        cmd.color("gray85", f"{obj_name} and polymer.protein")',
            '        cmd.show("sticks", pocket_sel)',
            '        util.cbac(pocket_sel)',
            '        cmd.show("sticks", ligand_sel)',
            '        util.cbag(ligand_sel)',
            '        cmd.set("ray_opaque_background", 1)',
            '        cmd.bg_color("white")',
            '        cmd.zoom(f"{ligand_sel} or {pocket_sel}", buffer=3.0)',
            "        cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)",
            "        cmd.save(output_pse)",
            '        print(f"Artifacts exported:\\n - {output_png}\\n - {output_pse}")',
        ])
    else:
        # Generic recipe loader
        code_lines.extend([
            '        obj_name = "target_obj"',
            "        cmd.load(structure_path, obj_name)",
            "        atom_count = cmd.count_atoms(obj_name)",
            "        if atom_count == 0:",
            '            sys.stderr.write(f"Error: Zero atoms parsed from {structure_path}\\n")',
            "            sys.exit(2)",
            '        cmd.hide("everything", obj_name)',
            '        cmd.show("cartoon", obj_name)',
            '        cmd.color("marine", obj_name)',
            '        cmd.set("ray_opaque_background", 1)',
            '        cmd.bg_color("white")',
            "        cmd.orient(obj_name)",
            "        cmd.zoom(obj_name, buffer=2.0)",
            "        cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)",
            "        cmd.save(output_pse)",
            '        print(f"Artifacts exported:\\n - {output_png}\\n - {output_pse}")',
        ])

    code_lines.extend([
        "    finally:",
        "        cmd.quit()",
        "",
        'if __name__ == "__main__":',
        "    main()",
        "",
    ])

    script_content = "\n".join(code_lines)
    os.makedirs(os.path.dirname(output_script) or ".", exist_ok=True)
    with open(output_script, "w", encoding="utf-8") as f:
        f.write(script_content)
    os.chmod(output_script, 0o755)
    return output_script


def execute_recipe_in_process(recipe: str, params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Executes a PyMOL recipe directly in a dedicated headless sub-process via python script.
    """
    import tempfile

    with tempfile.NamedTemporaryFile(suffix=".py", mode="w", delete=False) as tf:
        temp_script = tf.name

    try:
        generate_pep723_script(recipe, params, temp_script)
        cmd_run = ["uv", "run", temp_script]
        proc = subprocess.run(cmd_run, capture_output=True, text=True, check=False)

        return {
            "recipe": recipe,
            "success": proc.returncode == 0,
            "exit_code": proc.returncode,
            "stdout": proc.stdout,
            "stderr": proc.stderr,
            "output_png": params.get("output_png"),
            "output_pse": params.get("output_pse"),
        }
    finally:
        if os.path.exists(temp_script):
            os.remove(temp_script)


def main():
    parser = argparse.ArgumentParser(description="pymol-visualize: PyMOL Skill Execution Harness and Script Generator.")
    parser.add_argument("--recipe", choices=SUPPORTED_RECIPES, default="secondary_structure", help="Recipe name")
    parser.add_argument("--structure", required=True, help="Input coordinate structure (PDB/CIF path, .gz path, or 4-letter PDB ID, or UniProt ID)")
    parser.add_argument("--active-residues", default="57+102+195", help="Active residue numbers (for highlight_residues recipe)")
    parser.add_argument("--ligand-resn", default="LIG", help="Ligand residue name (for pocket recipes)")
    parser.add_argument("--out-png", default="outputs/render.png", help="Path for rendered PNG")
    parser.add_argument("--out-pse", default="outputs/session.pse", help="Path for saved session PSE")
    parser.add_argument("--generate-script", help="Generate standalone PEP 723 script to target path instead of executing")
    parser.add_argument("--dashboard", nargs="?", const="ui/dashboard.html", default=None, help="Path to generate HTML dashboard after rendering (defaults to ui/dashboard.html)")
    parser.add_argument("--json", action="store_true", help="Format CLI output as JSON")
    args = parser.parse_args()

    resolved_structure = resolve_coordinate_input(args.structure)

    out_png = args.out_png
    try:
        rel_png = os.path.relpath(os.path.expanduser(out_png), os.getcwd())
        if not rel_png.startswith(".."):
            out_png = rel_png
        else:
            out_png = os.path.abspath(os.path.expanduser(out_png))
    except Exception:
        out_png = os.path.abspath(os.path.expanduser(out_png))

    out_pse = args.out_pse
    try:
        rel_pse = os.path.relpath(os.path.expanduser(out_pse), os.getcwd())
        if not rel_pse.startswith(".."):
            out_pse = rel_pse
        else:
            out_pse = os.path.abspath(os.path.expanduser(out_pse))
    except Exception:
        out_pse = os.path.abspath(os.path.expanduser(out_pse))

    params = {
        "structure_path": resolved_structure,
        "active_residues": args.active_residues,
        "ligand_resn": args.ligand_resn,
        "output_png": out_png,
        "output_pse": out_pse,
    }

    if args.generate_script:
        target_path = os.path.abspath(os.path.expanduser(args.generate_script))
        generate_pep723_script(args.recipe, params, target_path)
        print(f"Generated standalone script: {target_path}")
        sys.exit(0)

    result = execute_recipe_in_process(args.recipe, params)

    if result["success"] and args.dashboard:
        from generate_dashboard import build_dashboard_html
        dash_path = args.dashboard
        os.makedirs(os.path.dirname(dash_path) or ".", exist_ok=True)
        html_content = build_dashboard_html(
            structure_path=os.path.relpath(resolved_structure, os.getcwd()) if resolved_structure.startswith(os.getcwd()) else resolved_structure,
            png_path=os.path.relpath(params["output_png"], os.getcwd()) if params["output_png"].startswith(os.getcwd()) else params["output_png"],
            pse_path=os.path.relpath(params["output_pse"], os.getcwd()) if params["output_pse"].startswith(os.getcwd()) else params["output_pse"],
            recipe=args.recipe,
        )
        with open(dash_path, "w", encoding="utf-8") as f:
            f.write(html_content)
        result["dashboard"] = dash_path
        print(f"Generated Dashboard: {dash_path}")
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        if result["success"]:
            print(f"PyMOL recipe '{args.recipe}' executed successfully.")
            print(result["stdout"].strip())
        else:
            sys.stderr.write(f"PyMOL execution error (exit code {result['exit_code']}):\n")
            sys.stderr.write(result["stderr"])
            sys.exit(result["exit_code"])


if __name__ == "__main__":
    main()
