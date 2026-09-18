# PyMOL Technical Reference and API Guide

## 1. Headless Engine Setup and Architecture

### 1.1 OSMesa Software Rasterization
In containerized, CI/CD, and remote agent environments, a physical GPU and X11/Wayland display server are generally unavailable. PyMOL supports software-rasterized OpenGL rendering via the Off-Screen Mesa (OSMesa) library.
- Environment Variable: `PYOPENGL_PLATFORM=osmesa`
- The environment variable must be exported before importing any graphical or OpenGL-linked Python modules.

### 1.2 Startup Invocation Flags
PyMOL must be initialized using quiet, non-interactive command-line switches:
```python
import os
os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()
from pymol import cmd, util
```
- `-c`: Disables the graphical user interface, terminal prompt loop, and event listener threads.
- `-q`: Suppresses verbose startup banners and license splash logs.

### 1.3 Rasterization and Export Primitives
- `cmd.draw()`: Performs standard OpenGL hardware/software viewport capture without ray tracing. Fast but lacks advanced shadows, depth cueing, and surface occlusion.
- `cmd.ray(width, height)`: Executes the PyMOL built-in ray tracer on the current view buffer, calculating accurate shadows, lighting reflection models, and depth perspective.
- `cmd.png(filename, width=1600, height=1200, dpi=300, ray=1)`: Combines ray-tracing and raster export directly into a PNG file. When `ray=1` is specified, `cmd.png` internally executes ray tracing before serializing pixels, which is preferred for all publication-quality output.

---

## 2. Selection Syntax Guide

PyMOL selections define subsets of atoms, residues, chains, or geometric centroids.

### 2.1 Atomic and Residue Identifiers
- `chain <id>`: Matches chain identifier (e.g., `chain A`, `chain B`). Selection matches are case-sensitive.
- `resi <index>`: Matches residue sequence numbering.
  - Single residue: `resi 45`
  - Continuous range: `resi 45-60`
  - Discontinuous list: `resi 45+52+88`
  - Combined ranges: `resi 1-10+20-30`
- `resn <name>`: Matches 3-letter or 4-letter residue name (e.g., `resn ALA`, `resn HEM`, `resn ATP`).
- `name <atom_name>`: Matches PDB atom names (e.g., `name CA` for alpha carbons, `name CB` for beta carbons, `name N` for backbone amide nitrogen).
- `ss <type>`: Matches secondary structure assignments:
  - Helices: `ss h`
  - Beta sheets/strands: `ss s`
  - Loops and unstructured coils: `ss l+''`

### 2.2 Structural Classes
- `polymer.protein`: Canonical protein residues excluding solvent, nucleic acids, and small molecules.
- `polymer.nucleic`: DNA/RNA backbone and nucleotide bases.
- `organic`: Small-molecule ligands, prosthetic groups, non-polymeric organic molecules (excluding water).
- `solvent`: Water molecules (e.g., `resn HOH`, `resn WAT`) and common crystallization ions.
- `hetatm`: Non-standard hetero-atoms recorded in PDB HETATM records.
- `all`: Every coordinate point in the current session.

### 2.3 Boolean Operators and Grouping
Operators evaluate with standard precedence (parentheses override):
- `and`: Intersection of sets (e.g., `chain A and polymer.protein and resi 1-50`).
- `or`: Union of sets (e.g., `resi 57 or resi 102 or resi 195`).
- `not`: Complement / negation (e.g., `polymer.protein and not name H*`).
- Parentheses: Complex composite selectors must use parentheses to prevent ambiguity:
  `((chain A and resi 1-20) or (chain B and resi 50-70)) and not solvent`

### 2.4 Proximity and Neighborhood Selectors
- `within <radius> of (<selection>)`: Selects individual atoms located within `<radius>` angstroms of the target selection.
- `byres (<selection>)`: Expands any atom selection to include all atoms belonging to the corresponding complete residues.
- `byres ((selection_a) within 5.0 of (selection_b))`: Standard pattern for identifying whole amino acid residues forming a binding pocket.
- `around <radius>`: Selects atoms in a shell between 0 and `<radius>` angstroms from the target, excluding the target itself.

### 2.5 Selection Identifier Rules
- Selection names created via `cmd.select("name", ...)` must begin with an alphabetic character.
- Permitted characters: letters, numbers, and underscores (`[a-zA-Z0-9_]`).
- Prohibited characters: spaces, hyphens, slashes, punctuation, or PyMOL reserved keywords (`all`, `organic`, `polymer`).
- Invalid: `pocket-site`, `site 1`, `4pocket`.
- Valid: `pocket_site`, `site_1`, `pocket_4`.

---

## 3. Common API Commands

### 3.1 Structure Loading
```python
cmd.load(filename, object_name)
```
Loads coordinate files (`.pdb`, `.cif`, `.mmcif`) into an explicit named object. Always supply an explicit `object_name`.

### 3.2 Representations
```python
cmd.show(representation, selection)
cmd.hide(representation, selection)
```
Standard representations:
- `cartoon`: Ribbon diagram tracing secondary structure elements.
- `sticks`: Cylindrical chemical bonds with optional sphere vertices.
- `spheres`: Space-filling van der Waals spheres.
- `surface`: Solvent-accessible or molecular surface contour.
- `lines`: Wireframe representation for fast broad rendering.
- `ribbon`: Thin spline backbone trace.

### 3.3 Coloring Commands
```python
cmd.color("color_name", selection)
cmd.set_color("custom_color", [r, g, b])
cmd.spectrum(expression="b", palette="blue_white_red", selection=..., minimum=..., maximum=...)
```
Common color names: `cyan`, `salmon`, `slate`, `forest`, `orange`, `marine`, `ruby`, `tv_blue`, `tv_green`, `gray70`, `white`.
`util.cnc(selection)`: Colors non-carbon heteroatoms (nitrogen blue, oxygen red, sulfur yellow).
`util.cbac(selection)`: Colors carbon atoms cyan while applying CPK colors to heteroatoms.
`util.cbag(selection)`: Colors carbon atoms green while applying CPK colors to heteroatoms.

### 3.4 Camera, Alignment, and Geometric Analysis
- `cmd.orient(selection)`: Aligns the principal moments of inertia of the selection along the Cartesian axes and centers the view.
- `cmd.zoom(selection, buffer=2.0, complete=1)`: Frames the selection within the camera viewport with padding.
- `cmd.distance(object_name, selection_a, selection_b, cutoff=3.5, mode=2)`: Creates an explicit distance object showing dashed connector lines. Mode 0 measures all pairs, mode 2 restricts to polar contacts (hydrogen bonds).
- `cmd.get_distance(selection_a, selection_b)`: Returns floating-point Euclidean distance in angstroms between two atom coordinates.
- `cmd.align(mobile, target, cutoff=2.0, cycles=5)`: Performs sequence alignment followed by dynamic programming structural superposition. Returns tuple: `(rmsd, aligned_atoms, cycles, raw_rmsd, initial_atoms)`.
- `cmd.super(mobile, target)`: Structure-based superposition independent of sequence identity, useful for remote homologs.
- `cmd.cealign(target, mobile)`: Combinatorial Extension algorithm for structural alignment of divergent proteins. Note argument order: `(target, mobile)`. Returns dictionary with `"RMSD"` and `"alignment_length"`.

### 3.5 Global and Rendering Settings
- `cmd.set("ray_opaque_background", 1)`: Enforces non-transparent background in ray-traced images.
- `cmd.bg_color("white")`: Sets the viewport background color.
- `cmd.set("transparency", 0.3, selection)`: Applies 30% surface transparency to the target selection.
- `cmd.set("cartoon_fancy_helices", 1)`: Renders alpha-helices with stylized beveled ribbons.
- `cmd.set("antialias", 2)`: Enables 2x supersampling for clean edge anti-aliasing.

---

## 4. Comprehensive Pitfalls and Defensive Rules

### 4.1 Mandatory Process Teardown
PyMOL instances started via `pymol.finish_launching()` maintain open thread pools and memory allocations. Failure to execute `cmd.quit()` causes headless CLI runs to hang indefinitely, blocking agent tool executors. Always enclose PyMOL logic within `try...finally: cmd.quit()`.

### 4.2 Case Sensitivity in Selections
Identifiers in crystallographic and cryo-EM models are strictly case-sensitive. `chain A` and `chain a` designate distinct macromolecular entities. When generating dynamic scripts, inspect available chains via `cmd.get_chains(object_name)` rather than guessing casing.

### 4.3 Offline Isolation and Sandbox Policy
PyMOL contains a built-in `cmd.fetch()` command that queries online RCSB servers. In sandboxed execution environments, outbound internet access may be disabled or restricted. Relying on `cmd.fetch()` creates nondeterministic failures. All scripts must load pre-fetched local files via `cmd.load()`.

### 4.4 Measurement Node Selection Trap
Calling `cmd.distance("dist_obj", ...)` creates a measurement display object, not an atomic selection. Executing `cmd.count_atoms("dist_obj")` or `cmd.color("red", "dist_obj")` results in runtime errors or zero counts. Atoms must be styled directly on the macromolecule objects, and measurement labels/dashes styled via `cmd.set("dash_color", ...)` or `cmd.hide("labels", "dist_obj")`.

### 4.5 Multi-State Ensembles (NMR Models)
Structures derived from NMR or conformational sampling often contain multiple states (`MODEL 1`, `MODEL 2`, ...). Loading an ensemble without specifying `state` causes all states to render concurrently, producing thick, cluttered ribbons and distorted surface calculations.
- Restrict loading or operations to `state=1`: `cmd.show("cartoon", f"{obj_name} and state 1")` or `cmd.load(path, obj_name, state=1)`.

### 4.6 Surface Memory Constraints and Cavity Tuning
Calculating solvent-excluded surfaces on large macromolecular assemblies (greater than 50,000 atoms) can exhaust process memory. When calculating cavities:
- Constrain calculations using `cmd.set("cavity_cull", 10)` to ignore tiny, insignificant voids.
- Avoid unrestricted `surface_cavity_mode = 1` over whole viral capsids or ribosome assemblies without adequate swap space.

### 4.7 Electrostatic Potential and _champ Binary Dependency
The high-level utility `util.protein_vacuum_esp()` internally invokes `chempy.champ.assign` to automatically assign formal Amber charges and patch missing terminal atoms. However, standard Linux open-source PyMOL wheels (`pymol-open-source-whl`) omit the compiled `_champ` C binary, triggering a `NameError: name '_champ' is not defined`.
When executing in open-source headless environments:
- Enclose `util.protein_vacuum_esp()` within an exception-handling block.
- Fall back to assigning formal and partial charges via `cmd.alter()` and computing the Coulomb potential surface directly via `cmd.map_new(..., "coulomb_local", ...)` and `cmd.ramp_new()`.
