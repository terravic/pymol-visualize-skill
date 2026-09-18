# PyMOL Molecular Visualization Recipes

This document provides 16 modular, standalone PyMOL recipes for common structural biology and molecular visualization tasks. Each snippet assumes that the headless engine has been initialized and that the PEP 723 header and teardown protocols are active.

```python
# Headless initialization sequence assumed for all recipes:
# import os
# os.environ["PYOPENGL_PLATFORM"] = "osmesa"
# import pymol
# pymol.pymol_argv = ["pymol", "-cq"]
# pymol.finish_launching()
# from pymol import cmd, util
```

---

## Recipe 1: Cartoon with Secondary Structure Coloring

Renders secondary structural elements using distinct, high-contrast colors: helices in red, beta sheets in yellow, and coils/loops in light gray.

```python
def render_secondary_structure(structure_path, output_png, output_pse):
    obj_name = "protein_secstruct"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)

    # Secondary structure assignment styling
    cmd.color("gray75", f"{obj_name} and ss l+''")
    cmd.color("firebrick", f"{obj_name} and ss h")
    cmd.color("brightorange", f"{obj_name} and ss s")

    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.orient(obj_name)
    cmd.zoom(obj_name, buffer=2.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
```

---

## Recipe 2: B-Factor (pLDDT) Continuous Coloring

Applies a continuous color spectrum across crystallographic temperature factors or predicted local distance difference test (pLDDT) scores.

```python
def render_bfactor_continuous(structure_path, output_png, output_pse):
    obj_name = "protein_bfactor"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)

    # Map B-factor values to red-white-blue gradient on protein atoms
    cmd.spectrum("b", "red_white_blue", f"{obj_name} and polymer.protein")

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.orient(obj_name)
    cmd.zoom(obj_name, buffer=2.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
```

---

## Recipe 3: AlphaFold Canonical pLDDT Thresholds

Colors structural models according to the standardized four-bin AlphaFold confidence scale.

```python
def render_alphafold_plddt(structure_path, output_png, output_pse):
    obj_name = "af_structure"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)

    # Define standardized AlphaFold hex colors
    cmd.set_color("af_very_low", [1.000, 0.490, 0.271])   # Orange (< 50)
    cmd.set_color("af_low",      [1.000, 0.859, 0.075])   # Yellow (50 - 70)
    cmd.set_color("af_confident",[0.396, 0.796, 0.953])   # Light Blue (70 - 90)
    cmd.set_color("af_very_high",[0.000, 0.325, 0.839])   # Dark Blue (>= 90)

    # Apply binned selections sequentially
    cmd.color("af_very_low",  f"{obj_name} and b < 50")
    cmd.color("af_low",       f"{obj_name} and b >= 50 and b < 70")
    cmd.color("af_confident", f"{obj_name} and b >= 70 and b < 90")
    cmd.color("af_very_high", f"{obj_name} and b >= 90")

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.orient(obj_name)
    cmd.zoom(obj_name, buffer=2.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
```

---

## Recipe 4: Highlighting Specific Residues

Mutes the structural framework in translucent gray cartoon while highlighting key functional or active site residues in full stick detail with heteroatom coloring.

```python
def render_highlighted_residues(structure_path, active_residues, output_png, output_pse):
    # Example active_residues: "57+102+195"
    obj_name = "active_site_view"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", f"{obj_name} and polymer.protein")
    cmd.color("gray80", f"{obj_name} and polymer.protein")

    active_sel = f"{obj_name} and resi {active_residues}"
    cmd.show("sticks", active_sel)
    util.cnc(active_sel)  # Color non-carbon atoms: N blue, O red, S yellow

    cmd.set("stick_radius", 0.25)
    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.orient(active_sel)
    cmd.zoom(active_sel, buffer=4.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
```

---

## Recipe 5: Transparent Surface Rendering

Displays an outer semi-transparent solvent-accessible surface envelope over an internal ribbon backbone.

```python
def render_transparent_surface(structure_path, output_png, output_pse):
    obj_name = "surface_model"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)
    cmd.color("slate", obj_name)

    cmd.show("surface", obj_name)
    cmd.set("transparency", 0.35, obj_name)
    cmd.set("surface_quality", 1)
    cmd.color("white", f"{obj_name} and surface")

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.orient(obj_name)
    cmd.zoom(obj_name, buffer=2.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
```

---

## Recipe 6: Electrostatic Surface Rendering

Calculates and displays qualitative vacuum electrostatics across the solvent-accessible surface with balanced two-sided lighting. In open-source PyMOL builds lacking the compiled `_champ` module, calculates the Coulombic potential grid directly via `cmd.map_new` and `cmd.ramp_new`.

```python
def render_electrostatic_surface(structure_path, output_png, output_pse):
    obj_name = "esp_model"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)

    try:
        # Compute vacuum electrostatics surface potential via PyMOL utility
        util.protein_vacuum_esp(obj_name, mode=2, quiet=0)
    except Exception:
        # Fallback for environments lacking the optional compiled _champ extension:
        # Assign basic formal and partial charges directly to macromolecule atoms
        cmd.alter(f"{obj_name}", "q=0.0")
        cmd.alter(f"{obj_name} and resn LYS+ARG and name NZ+NH1+NH2", "q=1.0")
        cmd.alter(f"{obj_name} and resn ASP+GLU and name OD1+OD2+OE1+OE2", "q=-0.5")
        cmd.alter(f"{obj_name} and name N", "q=-0.3")
        cmd.alter(f"{obj_name} and name O", "q=-0.38")
        cmd.alter(f"{obj_name} and name C", "q=0.38")

        map_name = f"{obj_name}_e_map"
        pot_name = f"{obj_name}_e_pot"
        cmd.map_new(map_name, "coulomb_local", 1.0, obj_name, 5.0)
        cmd.ramp_new(pot_name, map_name, selection=obj_name, zero=1)
        cmd.show("surface", obj_name)
        cmd.set("surface_color", pot_name, obj_name)

    cmd.set("two_sided_lighting", 1)
    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.orient(obj_name)
    cmd.zoom(obj_name, buffer=2.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
```

---

## Recipe 7: Multi-Chain Complex Palette

Discovers all polypeptide chains in a complex, maps each chain to a distinct high-contrast color palette, and reports residue stoichiometry.

```python
def render_multichain_palette(structure_path, output_png, output_pse):
    obj_name = "complex_assembly"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    chains = cmd.get_chains(obj_name)
    print(f"Discovered {len(chains)} chains: {chains}")

    palette = [
        "deepteal", "salmon", "forest", "orange",
        "marine", "ruby", "warmpink", "slate", "sand"
    ]

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)

    for i, chain_id in enumerate(chains):
        color_choice = palette[i % len(palette)]
        chain_sel = f"{obj_name} and chain {chain_id}"
        cmd.color(color_choice, chain_sel)

        # Count residue count for reporting
        residues = cmd.get_model(f"{chain_sel} and name CA").atom
        print(f"Chain {chain_id}: {len(residues)} amino acid residues, color: {color_choice}")

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.orient(obj_name)
    cmd.zoom(obj_name, buffer=2.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
```

---

## Recipe 8: B-Factor Putty Analysis

Visualizes dynamic flexibility by modulating cartoon tube diameter according to crystallographic B-factors or positional uncertainty.

```python
def render_bfactor_putty(structure_path, output_png, output_pse):
    obj_name = "putty_model"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)
    cmd.cartoon("putty", obj_name)

    # Constrain putty scaling parameters
    cmd.set("cartoon_putty_scale_min", 0.3)
    cmd.set("cartoon_putty_scale_max", 3.0)
    cmd.set("cartoon_putty_radius", 0.35)

    cmd.spectrum("b", "rainbow", f"{obj_name} and polymer.protein")

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.orient(obj_name)
    cmd.zoom(obj_name, buffer=2.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
```

---

## Recipe 9: Cavity and Pocket Detection

Identifies internal solvent-excluded cavities and surface binding pockets, applying cavity culling to filter insignificant structural voids.

```python
def render_cavity_detection(structure_path, output_png, output_pse):
    obj_name = "cavity_target"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)
    cmd.color("gray80", obj_name)

    # Configure cavity mode and culling thresholds
    cmd.set("surface_cavity_mode", 1)     # Cavities and pockets
    cmd.set("surface_cavity_radius", 5)   # Detection radius in grid points
    cmd.set("surface_cavity_cutoff", 4)   # Boundary cutoff
    cmd.set("cavity_cull", 15)            # Discard cavities smaller than 15 grid points

    cmd.show("surface", obj_name)
    cmd.set("transparency", 0.45, obj_name)
    cmd.color("tv_yellow", f"{obj_name} and surface")

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.orient(obj_name)
    cmd.zoom(obj_name, buffer=2.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
```

---

## Recipe 10: Multi-Structure Batch Rendering

Iterates systematically over a directory of coordinate files, checks atom counts, renders each file independently to prevent memory leaks, and deletes objects between iterations.

```python
def render_structure_batch(structure_dir, output_dir):
    import os
    valid_exts = {".pdb", ".cif", ".mmcif"}
    entries = sorted([f for f in os.listdir(structure_dir) if os.path.splitext(f)[1].lower() in valid_exts])

    results = []
    for entry in entries:
        full_path = os.path.join(structure_dir, entry)
        base_name = os.path.splitext(entry)[0]
        obj_name = f"obj_{base_name}"
        png_path = os.path.join(output_dir, f"{base_name}.png")
        pse_path = os.path.join(output_dir, f"{base_name}.pse")

        cmd.load(full_path, obj_name)
        atom_count = cmd.count_atoms(obj_name)
        if atom_count == 0:
            print(f"Skipping empty structure: {entry}")
            cmd.delete(obj_name)
            continue

        cmd.hide("everything", obj_name)
        cmd.show("cartoon", obj_name)
        cmd.color("marine", obj_name)
        cmd.set("ray_opaque_background", 1)
        cmd.bg_color("white")
        cmd.orient(obj_name)
        cmd.zoom(obj_name, buffer=2.0)

        cmd.png(png_path, width=1200, height=900, dpi=150, ray=1)
        cmd.save(pse_path)
        cmd.delete(obj_name)  # Clear object from memory
        results.append((entry, atom_count, png_path))

    return results
```

---

## Recipe 11: Distance Measurement Between Residues

Measures the Euclidean distance between alpha carbons of two specified residues, creates a visible distance object, and extracts the numerical distance value.

```python
def measure_residue_distance(structure_path, chain_a, resi_a, chain_b, resi_b, output_png, output_pse):
    obj_name = "dist_target"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)
    cmd.color("gray75", obj_name)

    sel_a = f"{obj_name} and chain {chain_a} and resi {resi_a} and name CA"
    sel_b = f"{obj_name} and chain {chain_b} and resi {resi_b} and name CA"

    if cmd.count_atoms(sel_a) == 0 or cmd.count_atoms(sel_b) == 0:
        raise ValueError("Target C-alpha atoms not found for specified residue indices.")

    dist_name = f"dist_{resi_a}_{resi_b}"
    cmd.distance(dist_name, sel_a, sel_b)
    dist_value = cmd.get_distance(sel_a, sel_b)
    print(f"Euclidean distance between {chain_a}:{resi_a} and {chain_b}:{resi_b}: {dist_value:.3f} Angstroms")

    cmd.set("dash_color", "ruby")
    cmd.set("dash_width", 3.0)
    cmd.set("dash_gap", 0.3)

    cmd.show("sticks", f"byres ({sel_a} or {sel_b})")
    util.cnc(f"byres ({sel_a} or {sel_b})")

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.zoom(f"byres ({sel_a} or {sel_b})", buffer=5.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
    return dist_value
```

---

## Recipe 12: Zoom into Binding Pocket

Isolates and visualizes all protein residues situated within a defined radius of an organic ligand.

```python
def render_binding_pocket_zoom(structure_path, ligand_resn, radius_angstrom, output_png, output_pse):
    obj_name = "complex_view"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    ligand_sel = f"{obj_name} and resn {ligand_resn} and organic"
    if cmd.count_atoms(ligand_sel) == 0:
        raise ValueError(f"Ligand {ligand_resn} not detected in {structure_path}")

    pocket_sel = f"byres (({obj_name} and polymer.protein) within {radius_angstrom} of ({ligand_sel}))"

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", f"{obj_name} and polymer.protein")
    cmd.color("gray85", f"{obj_name} and polymer.protein")

    # Display pocket residues and ligand
    cmd.show("sticks", pocket_sel)
    util.cbac(pocket_sel)  # Cyan carbon sticks for pocket residues

    cmd.show("sticks", ligand_sel)
    util.cbag(ligand_sel)  # Green carbon sticks for ligand
    cmd.set("stick_radius", 0.3, ligand_sel)

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.orient(f"{ligand_sel} or {pocket_sel}")
    cmd.zoom(f"{ligand_sel} or {pocket_sel}", buffer=3.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
```

---

## Recipe 13: Detailed Protein-Ligand Interactions

Isolates a single ligand instance, computes polar contacts and hydrogen bonds, highlights pocket water molecules, and formats interaction dashes.

```python
def render_ligand_interactions(structure_path, output_png, output_pse):
    obj_name = "interaction_complex"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    # Isolate first ligand instance to avoid multimer ambiguity
    ligand_atoms = cmd.get_model(f"{obj_name} and organic").atom
    if not ligand_atoms:
        raise ValueError("No organic ligand detected in structure.")

    first_lig = ligand_atoms[0]
    lig_sel = f"{obj_name} and chain {first_lig.chain} and resi {first_lig.resi} and resn {first_lig.resn}"
    pocket_residues = f"byres (({obj_name} and polymer.protein) within 4.5 of ({lig_sel}))"
    pocket_waters = f"({obj_name} and solvent) within 3.5 of ({lig_sel})"

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", f"{obj_name} and polymer.protein")
    cmd.color("gray90", f"{obj_name} and polymer.protein")

    # Style pocket residues and ligand
    cmd.show("sticks", pocket_residues)
    util.cbac(pocket_residues)

    cmd.show("sticks", lig_sel)
    util.cbag(lig_sel)
    cmd.set("stick_radius", 0.35, lig_sel)

    # Style solvent waters
    if cmd.count_atoms(pocket_waters) > 0:
        cmd.show("spheres", pocket_waters)
        cmd.set("sphere_scale", 0.25, pocket_waters)
        cmd.color("red", pocket_waters)

    # Calculate polar contacts (mode=2 filters for potential H-bonds)
    polar_obj = "polar_contacts"
    cmd.distance(polar_obj, lig_sel, f"({pocket_residues}) or ({pocket_waters})", cutoff=3.5, mode=2)
    cmd.set("dash_color", "yellow")
    cmd.set("dash_width", 3.0)
    cmd.set("dash_gap", 0.25)
    cmd.hide("labels", polar_obj)

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.zoom(f"{lig_sel} or {pocket_residues}", buffer=3.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
```

---

## Recipe 14: In Silico Mutagenesis

Invokes PyMOL's mutagenesis wizard to substitute an amino acid at a target position, applies the preferred rotamer, and exports the modified structure and session.

```python
def render_in_silico_mutagenesis(structure_path, chain_id, resi_num, target_aa, output_png, output_pse):
    # Example: chain_id="A", resi_num="57", target_aa="ALA"
    obj_name = "mutant_model"
    cmd.load(structure_path, obj_name)
    if cmd.count_atoms(obj_name) == 0:
        raise ValueError(f"No atoms found in {structure_path}")

    # Launch mutagenesis wizard
    cmd.wizard("mutagenesis")
    cmd.refresh_wizard()

    # Set the target replacement amino acid
    cmd.get_wizard().set_mode(target_aa.upper())
    target_sel = f"/{obj_name}//{chain_id}/{resi_num}"
    cmd.get_wizard().do_select(target_sel)

    # Select rotamer 1 and apply substitution
    cmd.get_wizard().apply()
    cmd.set_wizard(None)  # Dismiss wizard cleanly

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)
    cmd.color("gray80", obj_name)

    mutant_resi_sel = f"{obj_name} and chain {chain_id} and resi {resi_num}"
    cmd.show("sticks", mutant_resi_sel)
    util.cnc(mutant_resi_sel)

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.zoom(mutant_resi_sel, buffer=5.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
```

---

## Recipe 15: Two-Structure Superposition with Adaptive Fallback

Superimposes a mobile structure onto a target structure using sequence alignment (`cmd.align`), with an automatic fallback to Combinatorial Extension (`cmd.cealign`) when sequence identity or alignment count is insufficient.

```python
def superimpose_structures(target_path, mobile_path, output_png, output_pse):
    target_obj = "target_struct"
    mobile_obj = "mobile_struct"

    cmd.load(target_path, target_obj)
    cmd.load(mobile_path, mobile_obj)

    if cmd.count_atoms(target_obj) == 0 or cmd.count_atoms(mobile_obj) == 0:
        raise ValueError("Failed to load atomic coordinates for alignment.")

    rmsd = None
    aligned_atoms = 0
    method = "align"

    try:
        # Attempt sequence-based structural alignment
        align_result = cmd.align(mobile_obj, target_obj, cutoff=2.0, cycles=5)
        rmsd = align_result[0]
        aligned_atoms = align_result[1]
        print(f"Sequence alignment succeeded: RMSD={rmsd:.3f} across {aligned_atoms} atoms")

        if aligned_atoms < 15:
            raise RuntimeError("Fewer than 15 atoms aligned; sequence divergence detected.")
    except Exception as exc:
        print(f"Standard align failed or produced insufficient matches ({exc}). Falling back to cealign.")
        method = "cealign"
        ce_result = cmd.cealign(target=target_obj, mobile=mobile_obj)
        rmsd = ce_result["RMSD"]
        aligned_atoms = ce_result["alignment_length"]
        print(f"CEAlign succeeded: RMSD={rmsd:.3f} across {aligned_atoms} aligned residues")

    cmd.hide("everything", "all")
    cmd.show("cartoon", "all")
    cmd.color("marine", target_obj)
    cmd.color("salmon", mobile_obj)

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.orient(target_obj)
    cmd.zoom("all", buffer=2.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)

    return {"method": method, "rmsd": rmsd, "aligned_atoms": aligned_atoms}
```

---

## Recipe 16: Session File Modification

Loads an existing serialized PyMOL session (`.pse`), modifies molecular representations and surface parameters, and saves the updated session.

```python
def modify_existing_session(input_pse, output_png, output_pse):
    # Load session state into current namespace
    cmd.load(input_pse)

    # Restyle all currently active objects
    cmd.show("surface", "polymer.protein")
    cmd.set("transparency", 0.3, "polymer.protein")
    cmd.color("lightteal", "polymer.protein")

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.zoom("all", buffer=2.0)

    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.save(output_pse)
    print(f"Successfully updated and re-saved session to {output_pse}")
```
