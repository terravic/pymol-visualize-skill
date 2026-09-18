#!/usr/bin/env python3
"""
generate_synthetic_assets.py
Generates 100% synthetic, HCLS-compliant macromolecular coordinate files (PDB and mmCIF)
for unit testing, smoke testing, and system verification of PyMOL recipes.
Contains zero PHI/PII.
"""

import math
import os

ASSETS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "assets"))


def generate_alpha_helix_coords(num_residues=12, start_resi=1, chain_id="A", b_factors=None):
    """
    Generates standard alpha-helix coordinates.
    Alpha helix parameters: 3.6 residues/turn, 1.5 A rise per residue, ~2.3 A CA radius.
    """
    if b_factors is None:
        b_factors = [85.0] * num_residues

    amino_acids = ["ALA", "LEU", "GLU", "MET", "GLN", "LYS", "ARG", "ILE", "VAL", "TYR", "ASP", "PHE"]
    records = []
    atom_id = 1

    for i in range(num_residues):
        resi = start_resi + i
        resn = amino_acids[i % len(amino_acids)]
        b_fact = b_factors[i % len(b_factors)]
        theta = (2.0 * math.pi / 3.6) * i
        z_ca = 1.5 * i
        r_ca = 2.3

        ca_x = r_ca * math.cos(theta)
        ca_y = r_ca * math.sin(theta)
        ca_z = z_ca

        # Approximate backbone N, C, O relative to CA
        n_x = ca_x - 0.8 * math.cos(theta)
        n_y = ca_y - 0.8 * math.sin(theta)
        n_z = ca_z - 1.0

        c_x = ca_x + 0.8 * math.cos(theta + 0.3)
        c_y = ca_y + 0.8 * math.sin(theta + 0.3)
        c_z = ca_z + 0.8

        o_x = c_x + 0.6 * math.cos(theta + 0.5)
        o_y = c_y + 0.6 * math.sin(theta + 0.5)
        o_z = c_z + 0.3

        cb_x = ca_x + 1.2 * math.cos(theta - 0.4)
        cb_y = ca_y + 1.2 * math.sin(theta - 0.4)
        cb_z = ca_z + 0.2

        atoms = [
            ("N", "N", n_x, n_y, n_z),
            ("CA", "C", ca_x, ca_y, ca_z),
            ("C", "C", c_x, c_y, c_z),
            ("O", "O", o_x, o_y, o_z),
            ("CB", "C", cb_x, cb_y, cb_z),
        ]

        for name, elem, x, y, z in atoms:
            line = (
                f"ATOM  {atom_id:5d} {name:<4s} {resn:3s} {chain_id}{resi:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{b_fact:6.2f}          {elem:>2s}\n"
            )
            records.append(line)
            atom_id += 1

    return records, atom_id


def generate_beta_strand_coords(num_residues=8, start_resi=13, chain_id="A", b_factors=None, start_atom_id=100):
    """
    Generates extended beta-strand coordinates.
    Extended beta strand parameters: ~3.4 A rise per residue.
    """
    if b_factors is None:
        b_factors = [60.0] * num_residues

    amino_acids = ["VAL", "ILE", "TYR", "THR", "PHE", "TRP", "LEU", "SER"]
    records = []
    atom_id = start_atom_id

    for i in range(num_residues):
        resi = start_resi + i
        resn = amino_acids[i % len(amino_acids)]
        b_fact = b_factors[i % len(b_factors)]
        z = 18.0 + 3.4 * i
        y_offset = 1.0 if (i % 2 == 0) else -1.0
        x_offset = 5.0

        ca_x = x_offset
        ca_y = y_offset
        ca_z = z

        n_x = ca_x - 0.5
        n_y = ca_y - 0.8
        n_z = ca_z - 1.2

        c_x = ca_x + 0.5
        c_y = ca_y + 0.8
        c_z = ca_z + 1.2

        o_x = c_x + 0.4
        o_y = c_y + 1.2
        o_z = c_z + 0.4

        cb_x = ca_x + 1.2 * (-1 if (i % 2 == 0) else 1)
        cb_y = ca_y
        cb_z = ca_z

        atoms = [
            ("N", "N", n_x, n_y, n_z),
            ("CA", "C", ca_x, ca_y, ca_z),
            ("C", "C", c_x, c_y, c_z),
            ("O", "O", o_x, o_y, o_z),
            ("CB", "C", cb_x, cb_y, cb_z),
        ]

        for name, elem, x, y, z in atoms:
            line = (
                f"ATOM  {atom_id:5d} {name:<4s} {resn:3s} {chain_id}{resi:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{b_fact:6.2f}          {elem:>2s}\n"
            )
            records.append(line)
            atom_id += 1

    return records, atom_id


def build_synthetic_monomer():
    """
    Builds synthetic monomer PDB containing helix and sheet with diverse pLDDT/B-factors:
    - Residues 1-4: B-factor < 50 (e.g. 42.0)
    - Residues 5-8: B-factor 50-70 (e.g. 65.0)
    - Residues 9-14: B-factor 70-90 (e.g. 78.0)
    - Residues 15-20: B-factor >= 90 (e.g. 96.0)
    """
    # 12 helix residues
    helix_b = [42.0, 45.0, 48.0, 49.0, 55.0, 62.0, 68.0, 69.5, 75.0, 82.0, 88.0, 89.0]
    helix_lines, next_atom_id = generate_alpha_helix_coords(num_residues=12, start_resi=1, chain_id="A", b_factors=helix_b)

    # 8 sheet residues
    sheet_b = [92.0, 94.0, 96.0, 98.0, 95.0, 93.0, 91.0, 90.5]
    sheet_lines, _ = generate_beta_strand_coords(num_residues=8, start_resi=13, chain_id="A", b_factors=sheet_b, start_atom_id=next_atom_id)

    header = (
        "HEADER    SYNTHETIC PROTEIN MONOMER               15-SEP-26   SYN1              \n"
        "TITLE     SYNTHETIC VALIDATION MONOMER FOR PYMOL SKILL                          \n"
        "HELIX    1   1 ALA A    1  PHE A   12  1                                  12    \n"
        "SHEET    1   A 1 VAL A  13  SER A  20  0                                        \n"
    )
    footer = "END\n"

    content = header + "".join(helix_lines) + "".join(sheet_lines) + footer
    path = os.path.join(ASSETS_DIR, "synthetic_monomer.pdb")
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"Wrote {path}")


def build_synthetic_complex_cif():
    """
    Builds synthetic mmCIF complex with Chain A and Chain B.
    """
    helix_a, _ = generate_alpha_helix_coords(num_residues=10, start_resi=1, chain_id="A", b_factors=[85.0] * 10)
    helix_b, _ = generate_alpha_helix_coords(num_residues=10, start_resi=1, chain_id="B", b_factors=[75.0] * 10)

    # Convert PDB ATOM lines to mmCIF atom_site format
    lines = [
        "data_SYN_COMPLEX\n",
        "#\n",
        "_entry.id   SYN_COMPLEX\n",
        "#\n",
        "_audit_conform.dict_name       mmcif_pdbx.dic\n",
        "_audit_conform.dict_version    5.378\n",
        "#\n",
        "loop_\n",
        "_atom_site.group_PDB\n",
        "_atom_site.id\n",
        "_atom_site.type_symbol\n",
        "_atom_site.label_atom_id\n",
        "_atom_site.label_alt_id\n",
        "_atom_site.label_comp_id\n",
        "_atom_site.label_asym_id\n",
        "_atom_site.label_entity_id\n",
        "_atom_site.label_seq_id\n",
        "_atom_site.Cartn_x\n",
        "_atom_site.Cartn_y\n",
        "_atom_site.Cartn_z\n",
        "_atom_site.occupancy\n",
        "_atom_site.B_iso_or_equiv\n",
        "_atom_site.auth_seq_id\n",
        "_atom_site.auth_comp_id\n",
        "_atom_site.auth_asym_id\n",
        "_atom_site.auth_atom_id\n",
    ]

    cif_atom_id = 1
    # Chain A
    for line in helix_a:
        atom_name = line[12:16].strip()
        comp_id = line[17:20].strip()
        seq_id = line[22:26].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        b = float(line[60:66])
        elem = line[76:78].strip()
        cif_line = (
            f"ATOM {cif_atom_id} {elem} {atom_name} . {comp_id} A 1 {seq_id} "
            f"{x:.3f} {y:.3f} {z:.3f} 1.00 {b:.2f} {seq_id} {comp_id} A {atom_name}\n"
        )
        lines.append(cif_line)
        cif_atom_id += 1

    # Chain B with spatial offset (x + 12.0)
    for line in helix_b:
        atom_name = line[12:16].strip()
        comp_id = line[17:20].strip()
        seq_id = line[22:26].strip()
        x = float(line[30:38]) + 12.0
        y = float(line[38:46]) + 2.0
        z = float(line[46:54])
        b = float(line[60:66])
        elem = line[76:78].strip()
        cif_line = (
            f"ATOM {cif_atom_id} {elem} {atom_name} . {comp_id} B 2 {seq_id} "
            f"{x:.3f} {y:.3f} {z:.3f} 1.00 {b:.2f} {seq_id} {comp_id} B {atom_name}\n"
        )
        lines.append(cif_line)
        cif_atom_id += 1

    lines.append("#\n")
    path = os.path.join(ASSETS_DIR, "synthetic_complex.cif")
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"Wrote {path}")


def build_synthetic_ligand_pocket():
    """
    Builds a synthetic protein-ligand complex:
    - Chain A protein residues (His57, Asp102, Ser195-like pocket)
    - Organic ligand LIG (Residue 201) positioned within 2.8 - 4.0 A of pocket residues
    - Solvent water HOH (Residue 301)
    """
    protein_atoms = [
        # Pocket residue 1: HIS 57
        ("ATOM  ", 1, "N", "HIS", "A", 57, 0.000, 0.000, 0.000, 1.00, 20.00, "N"),
        ("ATOM  ", 2, "CA", "HIS", "A", 57, 1.450, 0.000, 0.000, 1.00, 20.00, "C"),
        ("ATOM  ", 3, "C", "HIS", "A", 57, 2.000, 1.420, 0.000, 1.00, 20.00, "C"),
        ("ATOM  ", 4, "O", "HIS", "A", 57, 1.250, 2.380, 0.000, 1.00, 20.00, "O"),
        ("ATOM  ", 5, "CB", "HIS", "A", 57, 2.000, -0.750, 1.220, 1.00, 22.00, "C"),
        ("ATOM  ", 6, "CG", "HIS", "A", 57, 1.500, -2.150, 1.350, 1.00, 25.00, "C"),
        ("ATOM  ", 7, "ND1", "HIS", "A", 57, 2.250, -3.150, 1.900, 1.00, 25.00, "N"),
        ("ATOM  ", 8, "CE1", "HIS", "A", 57, 1.550, -4.250, 1.850, 1.00, 25.00, "C"),
        ("ATOM  ", 9, "NE2", "HIS", "A", 57, 0.350, -4.000, 1.300, 1.00, 25.00, "N"),
        ("ATOM  ", 10, "CD2", "HIS", "A", 57, 0.300, -2.680, 0.980, 1.00, 25.00, "C"),

        # Pocket residue 2: SER 195
        ("ATOM  ", 11, "N", "SER", "A", 195, -1.500, -2.500, 3.500, 1.00, 22.00, "N"),
        ("ATOM  ", 12, "CA", "SER", "A", 195, -1.800, -1.200, 4.100, 1.00, 22.00, "C"),
        ("ATOM  ", 13, "C", "SER", "A", 195, -3.300, -1.000, 4.150, 1.00, 22.00, "C"),
        ("ATOM  ", 14, "O", "SER", "A", 195, -4.050, -1.950, 3.950, 1.00, 22.00, "O"),
        ("ATOM  ", 15, "CB", "SER", "A", 195, -1.150, -1.100, 5.480, 1.00, 24.00, "C"),
        ("ATOM  ", 16, "OG", "SER", "A", 195, -1.380, 0.180, 6.050, 1.00, 26.00, "O"),

        # Pocket residue 3: ASP 102
        ("ATOM  ", 17, "N", "ASP", "A", 102, 4.500, -3.000, 1.000, 1.00, 20.00, "N"),
        ("ATOM  ", 18, "CA", "ASP", "A", 102, 4.800, -4.200, 1.750, 1.00, 20.00, "C"),
        ("ATOM  ", 19, "C", "ASP", "A", 102, 4.200, -5.450, 1.100, 1.00, 20.00, "C"),
        ("ATOM  ", 20, "O", "ASP", "A", 102, 4.800, -6.520, 1.150, 1.00, 20.00, "O"),
        ("ATOM  ", 21, "CB", "ASP", "A", 102, 4.250, -4.050, 3.180, 1.00, 22.00, "C"),
        ("ATOM  ", 22, "CG", "ASP", "A", 102, 4.600, -5.180, 4.120, 1.00, 24.00, "C"),
        ("ATOM  ", 23, "OD1", "ASP", "A", 102, 4.050, -6.280, 3.950, 1.00, 25.00, "O"),
        ("ATOM  ", 24, "OD2", "ASP", "A", 102, 5.420, -4.980, 5.050, 1.00, 25.00, "O"),
    ]

    # Organic Ligand LIG 201 forming a polar hydrogen bond to HIS 57 NE2 and SER 195 OG
    ligand_atoms = [
        ("HETATM", 25, "C1", "LIG", "A", 201, 0.100, -5.300, 3.200, 1.00, 30.00, "C"),
        ("HETATM", 26, "C2", "LIG", "A", 201, 0.450, -6.450, 2.300, 1.00, 30.00, "C"),
        ("HETATM", 27, "O1", "LIG", "A", 201, 0.200, -4.100, 2.450, 1.00, 30.00, "O"), # ~2.8 A to NE2
        ("HETATM", 28, "N1", "LIG", "A", 201, -0.600, -5.450, 4.350, 1.00, 30.00, "N"),
        ("HETATM", 29, "C3", "LIG", "A", 201, -1.050, -4.250, 5.100, 1.00, 30.00, "C"),
        ("HETATM", 30, "O2", "LIG", "A", 201, -0.550, -3.150, 4.750, 1.00, 30.00, "O"), # ~2.9 A to SER OG
    ]

    # Water HOH 301
    water_atoms = [
        ("HETATM", 31, "O", "HOH", "A", 301, 2.500, -5.800, 3.500, 1.00, 35.00, "O"),
    ]

    lines = [
        "HEADER    SYNTHETIC PROTEIN-LIGAND COMPLEX        15-SEP-26   SYN2              \n",
        "TITLE     SYNTHETIC RECEPTOR AND LIGAND FOR BINDING POCKET TESTS                \n",
    ]

    for grp, num, name, resn, chain, resi, x, y, z, occ, b, elem in protein_atoms + ligand_atoms + water_atoms:
        line = (
            f"{grp:<6s}{num:5d} {name:<4s} {resn:3s} {chain}{resi:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  {occ:4.2f}{b:6.2f}          {elem:>2s}\n"
        )
        lines.append(line)

    lines.append("CONECT   25   26   27   28                                                  \n")
    lines.append("CONECT   28   25   29                                                       \n")
    lines.append("CONECT   29   28   30                                                       \n")
    lines.append("END\n")

    path = os.path.join(ASSETS_DIR, "synthetic_ligand_pocket.pdb")
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"Wrote {path}")


def build_synthetic_nmr_ensemble():
    """
    Builds a multi-model PDB (MODEL 1, MODEL 2) to test NMR ensemble handling.
    """
    model1_lines, _ = generate_alpha_helix_coords(num_residues=8, start_resi=1, chain_id="A", b_factors=[25.0] * 8)
    # Model 2 has slightly displaced coordinates
    model2_lines, _ = generate_alpha_helix_coords(num_residues=8, start_resi=1, chain_id="A", b_factors=[35.0] * 8)

    lines = [
        "HEADER    SYNTHETIC NMR ENSEMBLE                  15-SEP-26   SYN3              \n",
        "TITLE     SYNTHETIC MULTI-STATE ENSEMBLE FOR PYMOL TESTING                      \n",
        "MODEL        1\n",
    ]
    lines.extend(model1_lines)
    lines.append("ENDMDL\n")
    lines.append("MODEL        2\n")
    # Shift coordinates slightly for model 2
    for l in model2_lines:
        x = float(l[30:38]) + 0.3
        y = float(l[38:46]) - 0.2
        z = float(l[46:54]) + 0.1
        line = l[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + l[54:]
        lines.append(line)
    lines.append("ENDMDL\n")
    lines.append("END\n")

    path = os.path.join(ASSETS_DIR, "synthetic_nmr_ensemble.pdb")
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"Wrote {path}")


def build_synthetic_homolog():
    """
    Builds a homologous structure with slight conformational shifts for superposition testing.
    """
    # 12 helix residues slightly rotated and translated
    helix_lines, next_atom_id = generate_alpha_helix_coords(
        num_residues=12, start_resi=1, chain_id="A", b_factors=[50.0] * 12
    )
    sheet_lines, _ = generate_beta_strand_coords(
        num_residues=8, start_resi=13, chain_id="A", b_factors=[50.0] * 8, start_atom_id=next_atom_id
    )

    transformed_lines = []
    cos_t = math.cos(0.08)
    sin_t = math.sin(0.08)
    for l in helix_lines + sheet_lines:
        x = float(l[30:38])
        y = float(l[38:46])
        z = float(l[46:54])
        # Slight rotation around z and translation
        nx = x * cos_t - y * sin_t + 0.4
        ny = x * sin_t + y * cos_t - 0.3
        nz = z + 0.2
        line = l[:30] + f"{nx:8.3f}{ny:8.3f}{nz:8.3f}" + l[54:]
        transformed_lines.append(line)

    header = (
        "HEADER    SYNTHETIC PROTEIN HOMOLOG               15-SEP-26   SYN4              \n",
        "TITLE     SYNTHETIC CONFORMATIONAL HOMOLOG FOR ALIGNMENT BENCHMARK              \n",
    )
    footer = "END\n"

    path = os.path.join(ASSETS_DIR, "synthetic_homolog.pdb")
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(header)
        f.writelines(transformed_lines)
        f.write(footer)
    print(f"Wrote {path}")


def main():
    os.makedirs(ASSETS_DIR, exist_ok=True)
    build_synthetic_monomer()
    build_synthetic_complex_cif()
    build_synthetic_ligand_pocket()
    build_synthetic_nmr_ensemble()
    build_synthetic_homolog()
    print("All synthetic assets generated successfully.")


if __name__ == "__main__":
    main()
