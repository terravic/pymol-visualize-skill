# Product Requirement Document (PRD)

## Project Name: PyMOL Molecular Visualization Skill (pymol-visualize)

---

## 1. Document Overview

### 1.1 Document Purpose
This Product Requirement Document (PRD) specifies the functional, technical, and architectural requirements for the pymol-visualize skill. It serves as the definitive reference for the design, execution, validation, and testing of headless molecular visualization workflows.

### 1.2 Status and Target Version
- Document Status: Baseline Approved
- Version: 1.0.0
- Target Environments: Linux x86_64, macOS
- Supported Python Runtimes: Python >= 3.10, < 3.13

---

## 2. Product Objectives and Scope

### 2.1 Core Objectives
1. Headless Execution: Enable automated agents and command-line interfaces to analyze, process, and render 3D macromolecular structures without requiring an X11 or Wayland display server or physical GPU.
2. Deterministic Script Synthesis: Synthesize isolated, standalone PEP 723 Python scripts that execute via uv run to guarantee reproducibility without modifying host environments.
3. Comprehensive Recipe Coverage: Provide 16 standardized biophysical visualization recipes spanning secondary structure cartoon rendering, AlphaFold confidence mapping, active-site inspection, ligand binding pocket analysis, structural superposition, and in silico mutagenesis.
4. Dual-Artifact Delivery: Export both publication-resolution PNG raster images (ray-traced at 300 DPI) and binary PyMOL session files (.pse) for interactive desktop inspection.
5. Interactive Web Dashboard: Provide a self-contained HTML and WebGL evaluation interface (ui/dashboard.html) equipped with 3D model interaction, metric breakdowns, and a 2D software rendering fallback.

### 2.2 Scope Boundaries
- Supported Data: Three-dimensional atomic coordinates in PDB (.pdb, .pdb.gz) and Macromolecular Crystallographic Information File (.cif, .mmcif, .cif.gz, .mmcif.gz) formats.
- Excluded Operations: Molecular dynamics simulations, thermodynamic trajectory integration, de novo protein folding algorithms, and sequence-only alignments lacking 3D atomic coordinates.

---

## 3. User Personas and Target Workflows

### 3.1 Autonomous Agent Harness
- Profile: Automated agent executing tools within sandboxed environments.
- Workflow: Parses user intent, validates input files via scripts/validate_structure.py, identifies the target recipe, executes scripts/pymol_runner.py or synthesizes a PEP 723 script, extracts numerical metrics, and delivers output file paths.

### 3.2 Computational Biologist and Bioinformatician
- Profile: Researcher performing automated structural triage and comparative analysis across coordinate libraries.
- Workflow: Invokes CLI scripts to batch-process assemblies, calculate structural alignments (RMSD), quantify chain stoichiometry, and inspect ligand-binding contacts.

### 3.3 Non-Technical Domain Researcher
- Profile: Scientist or educator requesting structural visualizations via natural language without writing code.
- Workflow: Submits queries (e.g., highlighting catalytic triads, inspecting AlphaFold prediction confidence), receives rendered images with clear explanations, and opens interactive 3D models in a browser.

---

## 4. Functional Requirements

### 4.1 Input Ingestion and Source Resolution (FR-1)
- FR-1.1 Local Path Resolution: Resolve relative file paths, expanded user directories (~), and verify presence on disk prior to script execution.
- FR-1.2 Database Retrieval: Provide automated fetching via scripts/fetch_structure.py:
  - 4-character PDB accessions (e.g., 1UBQ, 4HHB) retrieved from RCSB PDB (https://files.rcsb.org/).
  - UniProt accessions (e.g., P04637) retrieved from the AlphaFold Database (https://alphafold.ebi.ac.uk/).
- FR-1.3 Raw Coordinate Parsing: Accept raw PDB coordinate text pasted into chat, serialize to a temporary file, validate atom records, and pass to the runner.
- FR-1.4 Gzip Transparency: Handle .pdb.gz and .cif.gz archives across validation and PyMOL loading pipelines.

### 4.2 Pre-Flight Coordinate Validation (FR-2)
- FR-2.1 Verification Tool: Implement scripts/validate_structure.py to inspect input coordinate files prior to launching PyMOL.
- FR-2.2 Sanity Audits: Quantify total atoms, polymer atoms, heteroatoms, unique chains, residue counts, and detect ligand residues.
- FR-2.3 Multi-Model Detection: Identify NMR multi-state ensembles and report model counts.
- FR-2.4 Output Modalities: Support human-readable summary logs and machine-readable JSON output (--json).
- FR-2.5 Early Termination: Terminate with non-zero exit codes if files are empty, missing, or contain zero parsed atoms.

### 4.3 Headless Software Rasterization Engine (FR-3)
- FR-3.1 Display Independence: Configure Off-Screen Mesa software rasterization by setting PYOPENGL_PLATFORM=osmesa prior to OpenGL or PyMOL imports.
- FR-3.2 Execution Parameters: Initialize PyMOL with flags -cq (-c disables GUI and command loop; -q enables quiet execution).
- FR-3.3 Teardown Invariant: Enforce mandatory invocation of cmd.quit() in finally blocks to release framebuffer buffers and prevent thread leakage.

### 4.4 Scientific Visualization Recipes (FR-4)
Provide 16 modular recipes documented in skills/pymol-visualize/references/RECIPES.md:
- Recipe 1 (Secondary Structure): Cartoon representation with high-contrast palette (helices red, beta sheets yellow, loops light gray).
- Recipe 2 (B-Factor Spectrum): Continuous gradient mapping across crystallographic temperature factors or pLDDT fields.
- Recipe 3 (AlphaFold Canonical pLDDT): Four-tier categorical thresholding (Very High >= 90 dark blue, Confident 70-90 light blue, Low 50-70 yellow, Very Low < 50 orange).
- Recipe 4 (Active Site Highlights): Muted backbone scaffold with catalytic residues rendered in full stick representation with CPK heteroatom coloring.
- Recipe 5 (Transparent Surface): Semi-transparent solvent-accessible surface envelope superimposed over internal secondary structure.
- Recipe 6 (Electrostatic Surface): Qualitative vacuum electrostatic surface calculation with Coulomb potential map fallback.
- Recipe 7 (Multi-Chain Palette): Dynamic chain discovery and stoichiometric differentiation across hetero-oligomers.
- Recipe 8 (B-Factor Putty): Flexibility representation utilizing variable tube diameter proportional to thermal mobility.
- Recipe 9 (Cavity Detection): Identification and surface mapping of internal binding pockets and solvent-inaccessible cavities.
- Recipe 10 (Batch Rendering): Directory iteration with atom-count verification for batch dataset rendering.
- Recipe 11 (Distance Measurement): Precision Euclidean distance calculation between user-specified atom or residue centroids.
- Recipe 12 (Pocket Zoom): Camera targeting and 5.0 Angstrom neighborhood framing around ligand binding cavities.
- Recipe 13 (Protein-Ligand Contacts): Detection and visualization of polar hydrogen bonds and water-mediated interaction networks.
- Recipe 14 (In Silico Mutagenesis): Point mutation modeling using the PyMOL mutagenesis rotamer selection engine.
- Recipe 15 (Adaptive Superposition): Sequence-based structural alignment (cmd.align) with automated fallback to Combinatorial Extension (cmd.cealign) for divergent topologies.
- Recipe 16 (Session Modification): Loading, modifying, restyling, and serializing pre-existing .pse binary session files.

### 4.5 Dual-Artifact Export (FR-5)
- FR-5.1 High-Resolution Raster Image: Generate 300 DPI ray-traced PNG renders (cmd.png(..., dpi=300, ray=1)) with opaque white backgrounds and anti-aliasing.
- FR-5.2 Binary Session Persistence: Export complete PyMOL sessions (.pse) via cmd.save() preserving representations, color schemes, camera orientation, and atom selections.

### 4.6 Script Generation and Reproducibility (FR-6)
- FR-6.1 PEP 723 Metadata: Prepend inline dependency definitions (requires-python, dependencies = ["pymol-open-source-whl"]) to all synthesized scripts.
- FR-6.2 Standalone Synthesis: Provide CLI capability via python3 scripts/pymol_runner.py --generate-script <path> to emit self-contained executable scripts.
- FR-6.3 Relative Path Integrity: Ensure generated scripts reference coordinate inputs and output destinations using clean relative paths.

### 4.7 Interactive Evaluation Dashboard (FR-7)
- FR-7.1 WebGL 3D Visualization: Render interactive molecular views via bundled 3Dmol.js with mouse orbit, zoom, pan, and picking.
- FR-7.2 Representation and Scheme Switching: Support dynamic switching between cartoon, stick, and sphere representations, as well as AlphaFold pLDDT, secondary structure, and chain palette colorings.
- FR-7.3 2D Software Fallback: Provide an integrated 2D projection engine (Software2DViewer) that activates automatically if WebGL is unavailable or restricted.
- FR-7.4 Theme Management: Include a theme switcher supporting dark and light viewing modes.
- FR-7.5 Diagnostic Metrics Display: Tabulate atom counts, chain stoichiometry, mean B-factors, and AlphaFold confidence distributions.

---

## 5. Non-Functional Requirements

### 5.1 Performance and Latency
- Validation Execution: scripts/validate_structure.py must complete analysis in under 1.0 second for standard monomers (<10,000 atoms).
- Headless Rendering: Ray-traced PNG and PSE generation for monomeric targets must complete in under 15 seconds under OSMesa software rasterization.
- Dashboard Load Time: ui/dashboard.html must achieve complete DOM initialization and 3D rendering in under 2 seconds in standard modern web browsers.

### 5.2 Isolation and Environment Hygiene
- Package Management: All Python dependencies must resolve via uv without requiring global system package installation.
- Transient Subprocesses: Generated PyMOL execution runs in isolated subshell processes to ensure memory deallocation and state isolation.

### 5.3 Reliability and Error Recovery
- Missing Inputs: Emit structured error messages with clear failure descriptions when coordinate paths do not exist.
- Zero-Atom Detection: Abort execution immediately if parsed coordinate structures contain zero atoms.
- Superposition Fallback: Detect low-homology alignment failures and transition to cmd.cealign without process aborts.

### 5.4 Offline Independence
- Self-Contained Testing: The test suite must run purely offline against local synthetic assets without invoking external network requests.
- Bundled Assets: Store vendor JavaScript dependencies (ui/vendor/3Dmol-min.js) locally to support offline evaluation.

---

## 6. Architecture and System Components

```
+-------------------------------------------------------------------------+
|                              User / Agent                               |
+-------------------------------------------------------------------------+
                                     |
                                     v
+-------------------------------------------------------------------------+
|                           Ingestion Layer                               |
|   - scripts/fetch_structure.py (RCSB PDB, AlphaFold DB accessions)      |
|   - Local file path resolution (.pdb, .cif, .gz)                        |
+-------------------------------------------------------------------------+
                                     |
                                     v
+-------------------------------------------------------------------------+
|                          Validation Layer                               |
|   - scripts/validate_structure.py (Atom counts, chains, NMR detection)  |
+-------------------------------------------------------------------------+
                                     |
                                     v
+-------------------------------------------------------------------------+
|                      Execution & Synthesis Layer                        |
|   - scripts/pymol_runner.py (CLI orchestrator & PEP 723 generator)      |
|   - skills/pymol-visualize/references/RECIPES.md (16 scientific recipes)|
+-------------------------------------------------------------------------+
                                     |
                                     v
+-------------------------------------------------------------------------+
|                         Headless PyMOL Engine                           |
|   - Off-Screen Mesa (PYOPENGL_PLATFORM=osmesa)                          |
|   - PyMOL non-interactive runtime (-cq)                                 |
+-------------------------------------------------------------------------+
                                     |
                                     v
+-------------------------------------------------------------------------+
|                            Artifact Layer                               |
|   - outputs/*.png (300 DPI ray-traced publication render)               |
|   - outputs/*.pse (Binary PyMOL session file)                           |
+-------------------------------------------------------------------------+
                                     |
                                     v
+-------------------------------------------------------------------------+
|                            Dashboard Layer                              |
|   - scripts/generate_dashboard.py -> ui/dashboard.html                  |
|   - WebGL 3Dmol.js viewer with Software2DViewer fallback               |
+-------------------------------------------------------------------------+
```

---

## 7. Verification and Testing Framework

### 7.1 Test Suites
- Smoke Tests (tests/test_smoke.py): Validate license notices, Apache 2.0 terms, OSMesa headless environment, synthetic asset presence, and configuration files.
- System Runner Tests (tests/test_system_runner.py): End-to-end testing of CLI entrypoints (validate_structure.py, fetch_structure.py, pymol_runner.py, generate_dashboard.py), gzipped coordinate handling, and standalone script synthesis.
- Unit Recipe Tests (tests/test_unit_recipes.py): Unit tests for all 16 PyMOL recipes verifying atom selections, rendering calls, and artifact generation.

### 7.2 Test Execution Command
```bash
./scripts/run_tests.sh
```

---

## 8. License and Compliance Requirements

- Repository License: Apache License, Version 2.0.
- PyMOL Open-Source Compliance: Automated license acknowledgment protocol ensuring .licenses/pymol_LICENSE.txt records official terms and attribution prior to script generation.
- Third-Party Libraries: Attribution of 3Dmol.js under BSD open-source licensing.
