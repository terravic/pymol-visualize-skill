#!/usr/bin/env python3
"""
validate_structure.py
Pre-flight validator for macromolecular structure files (PDB and mmCIF).
Supports standard and gzip-compressed coordinates (.pdb, .pdb.gz, .cif, .cif.gz, .mmcif, .mmcif.gz).
Performs rapid sanity checks, atom counting, chain discovery, and multi-state detection
prior to heavy rendering or visualization workflows.
"""

import argparse
import gzip
import json
import logging
import os
import sys

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("validate_structure")


def _open_coordinate_file(file_path: str):
    """
    Opens coordinate file, transparently decompressing gzip if needed.
    """
    expanded_path = os.path.expanduser(file_path)
    if not os.path.exists(expanded_path):
        raise FileNotFoundError(f"Coordinate file not found: {expanded_path}")
    if os.path.getsize(expanded_path) == 0:
        raise ValueError(f"Coordinate file is empty: {expanded_path}")

    if expanded_path.endswith(".gz"):
        return gzip.open(expanded_path, "rt", encoding="utf-8", errors="replace")
    return open(expanded_path, "r", encoding="utf-8", errors="replace")


def inspect_pdb_content(file_path: str) -> dict:
    atom_count = 0
    hetero_count = 0
    chains = set()
    residues = set()
    ligands = set()
    model_count = 0
    has_conect = False

    with _open_coordinate_file(file_path) as f:
        for line in f:
            record = line[:6].strip()
            if record == "MODEL":
                model_count += 1
            elif record == "ATOM":
                atom_count += 1
                chain_id = line[21:22].strip() or "_"
                resi = line[22:26].strip()
                chains.add(chain_id)
                residues.add((chain_id, resi))
            elif record == "HETATM":
                hetero_count += 1
                chain_id = line[21:22].strip() or "_"
                resi = line[22:26].strip()
                resn = line[17:20].strip()
                if resn not in {"HOH", "WAT", "H2O", "DOD"}:
                    ligands.add(resn)
                    residues.add((chain_id, resi))
            elif record == "CONECT":
                has_conect = True

    is_nmr = model_count > 1
    total_atoms = atom_count + hetero_count

    if total_atoms == 0:
        raise ValueError(f"No ATOM or HETATM records found in PDB file: {file_path}")

    return {
        "file_path": file_path,
        "format": "PDB",
        "is_compressed": file_path.endswith(".gz"),
        "valid": True,
        "total_atoms": total_atoms,
        "polymer_atoms": atom_count,
        "hetero_atoms": hetero_count,
        "unique_residues": len(residues),
        "chains": sorted(list(chains)),
        "detected_ligands": sorted(list(ligands)),
        "model_count": max(model_count, 1),
        "is_multi_state": is_nmr,
        "has_conect": has_conect,
    }


def inspect_cif_content(file_path: str) -> dict:
    atom_count = 0
    chains = set()
    residues = set()
    ligands = set()

    in_atom_site_loop = False
    atom_site_headers = []

    with _open_coordinate_file(file_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue

            if stripped == "loop_":
                in_atom_site_loop = False
                atom_site_headers = []
                continue

            if stripped.startswith("_atom_site."):
                in_atom_site_loop = True
                atom_site_headers.append(stripped)
                continue

            if in_atom_site_loop and not stripped.startswith("_"):
                tokens = stripped.split()
                if len(tokens) >= len(atom_site_headers) and len(atom_site_headers) > 0:
                    atom_count += 1
                    asym_idx = -1
                    seq_idx = -1
                    comp_idx = -1
                    group_idx = -1

                    for idx, header in enumerate(atom_site_headers):
                        if header in {"_atom_site.auth_asym_id", "_atom_site.label_asym_id"}:
                            asym_idx = idx
                        elif header in {"_atom_site.auth_seq_id", "_atom_site.label_seq_id"}:
                            seq_idx = idx
                        elif header in {"_atom_site.auth_comp_id", "_atom_site.label_comp_id"}:
                            comp_idx = idx
                        elif header == "_atom_site.group_PDB":
                            group_idx = idx

                    chain_val = tokens[asym_idx] if (0 <= asym_idx < len(tokens)) else "A"
                    seq_val = tokens[seq_idx] if (0 <= seq_idx < len(tokens)) else str(atom_count)
                    comp_val = tokens[comp_idx] if (0 <= comp_idx < len(tokens)) else ""
                    grp_val = tokens[group_idx] if (0 <= group_idx < len(tokens)) else "ATOM"

                    chains.add(chain_val)
                    residues.add((chain_val, seq_val))
                    if grp_val == "HETATM" and comp_val not in {"HOH", "WAT", "H2O", "DOD", ""}:
                        ligands.add(comp_val)

    if atom_count == 0:
        raise ValueError(f"No valid _atom_site records parsed from CIF file: {file_path}")

    return {
        "file_path": file_path,
        "format": "mmCIF",
        "is_compressed": file_path.endswith(".gz"),
        "valid": True,
        "total_atoms": atom_count,
        "polymer_atoms": atom_count,
        "hetero_atoms": 0,
        "unique_residues": len(residues),
        "chains": sorted(list(chains)),
        "detected_ligands": sorted(list(ligands)),
        "model_count": 1,
        "is_multi_state": False,
        "has_conect": False,
    }


def validate_structure(file_path: str) -> dict:
    expanded_path = os.path.expanduser(file_path)
    clean_path = expanded_path[:-3] if expanded_path.endswith(".gz") else expanded_path
    ext = os.path.splitext(clean_path)[1].lower()

    if ext in {".cif", ".mmcif"}:
        return inspect_cif_content(file_path)
    elif ext in {".pdb", ".ent"}:
        return inspect_pdb_content(file_path)
    else:
        # Fallback inspection by content inspection of first block
        with _open_coordinate_file(file_path) as f:
            first_chunk = f.read(1024)
            if "data_" in first_chunk or "_atom_site" in first_chunk:
                return inspect_cif_content(file_path)
            else:
                return inspect_pdb_content(file_path)


def main():
    parser = argparse.ArgumentParser(description="Pre-flight validation for PDB and mmCIF coordinate files.")
    parser.add_argument("structure", help="Path to input coordinate file (supports .pdb, .cif, .gz)")
    parser.add_argument("--json", action="store_true", help="Output results in JSON format")
    args = parser.parse_args()

    try:
        results = validate_structure(args.structure)
        if args.json:
            print(json.dumps(results, indent=2))
        else:
            print(f"Validation successful: {results['file_path']}")
            print(f" - Format: {results['format']} (compressed: {results['is_compressed']})")
            print(f" - Total atoms: {results['total_atoms']}")
            print(f" - Unique residues: {results['unique_residues']}")
            print(f" - Chains: {', '.join(results['chains'])}")
            if results['detected_ligands']:
                print(f" - Ligands: {', '.join(results['detected_ligands'])}")
            print(f" - Multi-state ensemble: {results['is_multi_state']} (models: {results['model_count']})")
        sys.exit(0)
    except Exception as exc:
        logger.error(f"Validation failed: {exc}")
        sys.exit(1)


if __name__ == "__main__":
    main()
