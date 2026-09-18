"""
test_unit_recipes.py
Comprehensive unit test suite for all 16 PyMOL molecular visualization recipes.
All tests execute against 100% synthetic, HCLS-compliant coordinate assets.
"""

import os
import pytest
from pymol import util


def test_recipe_01_secondary_structure(pymol_cmd, synthetic_monomer_path, tmp_path):
    """Recipe 1: Cartoon with secondary structure coloring (helices, sheets, loops)."""
    obj_name = "protein_secstruct"
    pymol_cmd.load(synthetic_monomer_path, obj_name)
    assert pymol_cmd.count_atoms(obj_name) == 100

    pymol_cmd.hide("everything", obj_name)
    pymol_cmd.show("cartoon", obj_name)
    pymol_cmd.color("gray75", f"{obj_name} and ss l+''")
    pymol_cmd.color("firebrick", f"{obj_name} and ss h")
    pymol_cmd.color("brightorange", f"{obj_name} and ss s")

    pymol_cmd.set("cartoon_fancy_helices", 1)
    pymol_cmd.set("ray_opaque_background", 1)
    pymol_cmd.bg_color("white")
    pymol_cmd.orient(obj_name)
    pymol_cmd.zoom(obj_name, buffer=2.0)

    png_path = str(tmp_path / "recipe1.png")
    pse_path = str(tmp_path / "recipe1.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)

    assert os.path.exists(png_path)
    assert os.path.getsize(png_path) > 0
    assert os.path.exists(pse_path)
    assert os.path.getsize(pse_path) > 0


def test_recipe_02_bfactor_continuous(pymol_cmd, synthetic_monomer_path, tmp_path):
    """Recipe 2: B-factor (pLDDT) continuous coloring."""
    obj_name = "protein_bfactor"
    pymol_cmd.load(synthetic_monomer_path, obj_name)
    assert pymol_cmd.count_atoms(obj_name) == 100

    pymol_cmd.hide("everything", obj_name)
    pymol_cmd.show("cartoon", obj_name)
    pymol_cmd.spectrum("b", "red_white_blue", f"{obj_name} and polymer.protein")

    png_path = str(tmp_path / "recipe2.png")
    pse_path = str(tmp_path / "recipe2.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)

    assert os.path.exists(png_path)
    assert os.path.exists(pse_path)


def test_recipe_03_alphafold_plddt_thresholds(pymol_cmd, synthetic_monomer_path, tmp_path):
    """Recipe 3: AlphaFold canonical 4-tier pLDDT coloring."""
    obj_name = "af_structure"
    pymol_cmd.load(synthetic_monomer_path, obj_name)

    pymol_cmd.hide("everything", obj_name)
    pymol_cmd.show("cartoon", obj_name)

    pymol_cmd.set_color("af_very_low", [1.000, 0.490, 0.271])
    pymol_cmd.set_color("af_low", [1.000, 0.859, 0.075])
    pymol_cmd.set_color("af_confident", [0.396, 0.796, 0.953])
    pymol_cmd.set_color("af_very_high", [0.000, 0.325, 0.839])

    pymol_cmd.color("af_very_low", f"{obj_name} and b < 50")
    pymol_cmd.color("af_low", f"{obj_name} and b >= 50 and b < 70")
    pymol_cmd.color("af_confident", f"{obj_name} and b >= 70 and b < 90")
    pymol_cmd.color("af_very_high", f"{obj_name} and b >= 90")

    # Verify counts in each synthetic pLDDT bucket
    assert pymol_cmd.count_atoms(f"{obj_name} and b < 50") > 0
    assert pymol_cmd.count_atoms(f"{obj_name} and b >= 50 and b < 70") > 0
    assert pymol_cmd.count_atoms(f"{obj_name} and b >= 70 and b < 90") > 0
    assert pymol_cmd.count_atoms(f"{obj_name} and b >= 90") > 0

    png_path = str(tmp_path / "recipe3.png")
    pse_path = str(tmp_path / "recipe3.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)
    assert os.path.exists(png_path)


def test_recipe_04_highlight_residues(pymol_cmd, synthetic_ligand_pocket_path, tmp_path):
    """Recipe 4: Muted background cartoon with active site residues displayed as sticks."""
    obj_name = "active_site_view"
    pymol_cmd.load(synthetic_ligand_pocket_path, obj_name)

    pymol_cmd.hide("everything", obj_name)
    pymol_cmd.show("cartoon", f"{obj_name} and polymer.protein")
    pymol_cmd.color("gray80", f"{obj_name} and polymer.protein")

    active_sel = f"{obj_name} and resi 57+195"
    pymol_cmd.show("sticks", active_sel)
    util.cnc(active_sel)

    assert pymol_cmd.count_atoms(active_sel) > 0

    png_path = str(tmp_path / "recipe4.png")
    pse_path = str(tmp_path / "recipe4.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)
    assert os.path.exists(png_path)


def test_recipe_05_transparent_surface(pymol_cmd, synthetic_monomer_path, tmp_path):
    """Recipe 5: Cartoon visible through semi-transparent surface overlay."""
    obj_name = "surface_model"
    pymol_cmd.load(synthetic_monomer_path, obj_name)

    pymol_cmd.hide("everything", obj_name)
    pymol_cmd.show("cartoon", obj_name)
    pymol_cmd.show("surface", obj_name)
    pymol_cmd.set("transparency", 0.35, obj_name)
    pymol_cmd.color("white", f"{obj_name} and surface")

    png_path = str(tmp_path / "recipe5.png")
    pse_path = str(tmp_path / "recipe5.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)
    assert os.path.exists(png_path)


def test_recipe_06_electrostatic_surface(pymol_cmd, synthetic_monomer_path, tmp_path):
    """Recipe 6: Qualitative vacuum electrostatics with fallback."""
    obj_name = "esp_model"
    pymol_cmd.load(synthetic_monomer_path, obj_name)

    pymol_cmd.hide("everything", obj_name)
    pymol_cmd.show("cartoon", obj_name)

    try:
        util.protein_vacuum_esp(obj_name, mode=2, quiet=1)
    except Exception:
        # Fallback for environments without compiled _champ
        pymol_cmd.alter(f"{obj_name}", "q=0.0")
        pymol_cmd.alter(f"{obj_name} and resn LYS+ARG and name NZ+NH1+NH2", "q=1.0")
        pymol_cmd.alter(f"{obj_name} and resn ASP+GLU and name OD1+OD2+OE1+OE2", "q=-0.5")
        pymol_cmd.alter(f"{obj_name} and name N", "q=-0.3")
        pymol_cmd.alter(f"{obj_name} and name O", "q=-0.38")
        pymol_cmd.alter(f"{obj_name} and name C", "q=0.38")

        map_name = f"{obj_name}_e_map"
        pot_name = f"{obj_name}_e_pot"
        pymol_cmd.map_new(map_name, "coulomb_local", 1.0, obj_name, 5.0)
        pymol_cmd.ramp_new(pot_name, map_name, selection=obj_name, zero=1)
        pymol_cmd.show("surface", obj_name)
        pymol_cmd.set("surface_color", pot_name, obj_name)

    pymol_cmd.set("two_sided_lighting", 1)
    png_path = str(tmp_path / "recipe6.png")
    pse_path = str(tmp_path / "recipe6.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)
    assert os.path.exists(png_path)
    assert os.path.exists(pse_path)


def test_recipe_07_multichain_palette(pymol_cmd, synthetic_complex_cif_path, tmp_path):
    """Recipe 7: Multi-chain complex palette with dynamic chain discovery."""
    obj_name = "complex_assembly"
    pymol_cmd.load(synthetic_complex_cif_path, obj_name)

    chains = pymol_cmd.get_chains(obj_name)
    assert "A" in chains and "B" in chains

    palette = ["deepteal", "salmon", "forest", "orange"]
    pymol_cmd.hide("everything", obj_name)
    pymol_cmd.show("cartoon", obj_name)

    for i, chain_id in enumerate(chains):
        color_choice = palette[i % len(palette)]
        pymol_cmd.color(color_choice, f"{obj_name} and chain {chain_id}")

    png_path = str(tmp_path / "recipe7.png")
    pse_path = str(tmp_path / "recipe7.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)
    assert os.path.exists(png_path)


def test_recipe_08_bfactor_putty(pymol_cmd, synthetic_monomer_path, tmp_path):
    """Recipe 8: B-factor putty analysis with scale limits."""
    obj_name = "putty_model"
    pymol_cmd.load(synthetic_monomer_path, obj_name)

    pymol_cmd.hide("everything", obj_name)
    pymol_cmd.show("cartoon", obj_name)
    pymol_cmd.cartoon("putty", obj_name)

    pymol_cmd.set("cartoon_putty_scale_min", 0.3)
    pymol_cmd.set("cartoon_putty_scale_max", 3.0)
    pymol_cmd.spectrum("b", "rainbow", f"{obj_name} and polymer.protein")

    png_path = str(tmp_path / "recipe8.png")
    pse_path = str(tmp_path / "recipe8.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)
    assert os.path.exists(png_path)


def test_recipe_09_cavity_detection(pymol_cmd, synthetic_monomer_path, tmp_path):
    """Recipe 9: Cavity and pocket detection with cavity culling."""
    obj_name = "cavity_target"
    pymol_cmd.load(synthetic_monomer_path, obj_name)

    pymol_cmd.hide("everything", obj_name)
    pymol_cmd.show("cartoon", obj_name)
    pymol_cmd.set("surface_cavity_mode", 1)
    pymol_cmd.set("cavity_cull", 15)
    pymol_cmd.show("surface", obj_name)

    png_path = str(tmp_path / "recipe9.png")
    pse_path = str(tmp_path / "recipe9.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)
    assert os.path.exists(png_path)


def test_recipe_10_batch_rendering(pymol_cmd, synthetic_monomer_path, synthetic_homolog_path, tmp_path):
    """Recipe 10: Multi-structure batch rendering with memory cleanup."""
    batch_files = [synthetic_monomer_path, synthetic_homolog_path]
    rendered = []

    for idx, fpath in enumerate(batch_files):
        obj_name = f"batch_obj_{idx}"
        pymol_cmd.load(fpath, obj_name)
        count = pymol_cmd.count_atoms(obj_name)
        assert count > 0

        pymol_cmd.show("cartoon", obj_name)
        png_out = str(tmp_path / f"batch_{idx}.png")
        pymol_cmd.png(png_out, width=400, height=300, ray=0)
        pymol_cmd.delete(obj_name)
        assert os.path.exists(png_out)
        rendered.append(png_out)

    assert len(rendered) == 2


def test_recipe_11_distance_measurement(pymol_cmd, synthetic_monomer_path, tmp_path):
    """Recipe 11: Distance measurement between residues."""
    obj_name = "dist_target"
    pymol_cmd.load(synthetic_monomer_path, obj_name)

    sel_a = f"{obj_name} and chain A and resi 1 and name CA"
    sel_b = f"{obj_name} and chain A and resi 10 and name CA"

    dist_obj = "dist_1_10"
    pymol_cmd.distance(dist_obj, sel_a, sel_b)
    dist_val = pymol_cmd.get_distance(sel_a, sel_b)
    assert dist_val > 0.0
    assert 10.0 < dist_val < 25.0

    pymol_cmd.set("dash_color", "ruby")
    png_path = str(tmp_path / "recipe11.png")
    pse_path = str(tmp_path / "recipe11.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)
    assert os.path.exists(png_path)


def test_recipe_12_binding_pocket_zoom(pymol_cmd, synthetic_ligand_pocket_path, tmp_path):
    """Recipe 12: Zoom into binding pocket residues within radius of ligand."""
    obj_name = "complex_view"
    pymol_cmd.load(synthetic_ligand_pocket_path, obj_name)

    ligand_sel = f"{obj_name} and resn LIG and organic"
    assert pymol_cmd.count_atoms(ligand_sel) == 6

    pocket_sel = f"byres (({obj_name} and polymer.protein) within 5.0 of ({ligand_sel}))"
    assert pymol_cmd.count_atoms(pocket_sel) > 0

    pymol_cmd.hide("everything", obj_name)
    pymol_cmd.show("cartoon", f"{obj_name} and polymer.protein")
    pymol_cmd.show("sticks", pocket_sel)
    pymol_cmd.show("sticks", ligand_sel)
    util.cbac(pocket_sel)
    util.cbag(ligand_sel)

    pymol_cmd.zoom(f"{ligand_sel} or {pocket_sel}", buffer=3.0)

    png_path = str(tmp_path / "recipe12.png")
    pse_path = str(tmp_path / "recipe12.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)
    assert os.path.exists(png_path)


def test_recipe_13_protein_ligand_contacts(pymol_cmd, synthetic_ligand_pocket_path, tmp_path):
    """Recipe 13: Detailed protein-ligand interactions and polar contacts."""
    obj_name = "interaction_complex"
    pymol_cmd.load(synthetic_ligand_pocket_path, obj_name)

    lig_sel = f"{obj_name} and resn LIG and organic"
    pocket_waters = f"({obj_name} and solvent) within 4.0 of ({lig_sel})"

    polar_obj = "polar_contacts"
    pymol_cmd.distance(polar_obj, lig_sel, f"({obj_name} and polymer.protein) or ({pocket_waters})", cutoff=3.5, mode=2)
    pymol_cmd.set("dash_color", "yellow")
    pymol_cmd.hide("labels", polar_obj)

    png_path = str(tmp_path / "recipe13.png")
    pse_path = str(tmp_path / "recipe13.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)
    assert os.path.exists(png_path)


def test_recipe_14_in_silico_mutagenesis(pymol_cmd, synthetic_monomer_path, tmp_path):
    """Recipe 14: In silico mutagenesis wizard."""
    obj_name = "mutant_model"
    pymol_cmd.load(synthetic_monomer_path, obj_name)

    # Launch mutagenesis wizard
    pymol_cmd.wizard("mutagenesis")
    pymol_cmd.refresh_wizard()
    pymol_cmd.get_wizard().set_mode("ALA")
    pymol_cmd.get_wizard().do_select(f"/{obj_name}//A/2")
    pymol_cmd.get_wizard().apply()
    pymol_cmd.set_wizard(None)

    # Check residue 2 is present
    assert pymol_cmd.count_atoms(f"{obj_name} and chain A and resi 2") > 0

    png_path = str(tmp_path / "recipe14.png")
    pse_path = str(tmp_path / "recipe14.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)
    assert os.path.exists(png_path)


def test_recipe_15_superposition_adaptive_fallback(pymol_cmd, synthetic_monomer_path, synthetic_homolog_path, tmp_path):
    """Recipe 15: Two-structure superposition with sequence alignment and CEAlign fallback."""
    target_obj = "target_struct"
    mobile_obj = "mobile_struct"

    pymol_cmd.load(synthetic_monomer_path, target_obj)
    pymol_cmd.load(synthetic_homolog_path, mobile_obj)

    align_result = pymol_cmd.align(mobile_obj, target_obj, cutoff=2.0, cycles=5)
    rmsd = align_result[0]
    aligned_atoms = align_result[1]
    assert aligned_atoms > 0
    assert rmsd >= 0.0

    # Test that cealign also executes cleanly as fallback mechanism
    ce_result = pymol_cmd.cealign(target=target_obj, mobile=mobile_obj)
    assert "RMSD" in ce_result
    assert ce_result["RMSD"] >= 0.0
    assert ce_result["alignment_length"] > 0

    png_path = str(tmp_path / "recipe15.png")
    pse_path = str(tmp_path / "recipe15.pse")
    pymol_cmd.png(png_path, width=400, height=300, ray=0)
    pymol_cmd.save(pse_path)
    assert os.path.exists(png_path)


def test_recipe_16_session_file_modification(pymol_cmd, synthetic_monomer_path, tmp_path):
    """Recipe 16: Session file loading, modification, and re-saving."""
    initial_pse = str(tmp_path / "initial.pse")
    modified_pse = str(tmp_path / "modified.pse")
    modified_png = str(tmp_path / "modified.png")

    # Step A: Create and save initial session
    obj_name = "test_obj"
    pymol_cmd.load(synthetic_monomer_path, obj_name)
    pymol_cmd.save(initial_pse)
    pymol_cmd.delete("all")

    # Step B: Load session, modify and re-save
    pymol_cmd.load(initial_pse)
    assert pymol_cmd.count_atoms("all") == 100

    pymol_cmd.show("surface", "all")
    pymol_cmd.set("transparency", 0.3)
    pymol_cmd.color("lightteal", "all")

    pymol_cmd.png(modified_png, width=400, height=300, ray=0)
    pymol_cmd.save(modified_pse)

    assert os.path.exists(modified_png)
    assert os.path.exists(modified_pse)
    assert os.path.getsize(modified_pse) > 0
