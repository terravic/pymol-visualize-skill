#!/usr/bin/env python3
"""
fetch_structure.py
Retrieves macromolecular coordinate files from public structural biology repositories
(RCSB Protein Data Bank and AlphaFold DB) when users provide accession identifiers.
"""

import argparse
import os
import sys
import urllib.error
import urllib.request

from validate_structure import validate_structure


def fetch_rcsb_pdb(pdb_id: str, output_dir: str, file_format: str = "cif") -> str:
    """
    Downloads coordinate file from RCSB PDB (https://files.rcsb.org/).
    file_format: 'cif' or 'pdb'.
    """
    pdb_id = pdb_id.strip().upper()
    if len(pdb_id) != 4:
        raise ValueError(f"Invalid PDB accession code: '{pdb_id}'. PDB IDs must be 4 characters.")

    os.makedirs(output_dir, exist_ok=True)
    filename = f"{pdb_id.lower()}.{file_format}"
    dest_path = os.path.join(output_dir, filename)

    if file_format == "cif":
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    elif file_format == "pdb":
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    else:
        raise ValueError(f"Unsupported format: {file_format}. Use 'cif' or 'pdb'.")

    print(f"Fetching {pdb_id} from {url}...")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PyMOL-Visualize-Skill/1.0"})
        with urllib.request.urlopen(req, timeout=30) as response, open(dest_path, "wb") as out_file:
            out_file.write(response.read())
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"HTTP error {exc.code} fetching {pdb_id} from {url}: {exc.reason}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Network error fetching {pdb_id}: {exc.reason}") from exc

    # Immediate pre-flight validation
    val_results = validate_structure(dest_path)
    print(f"Downloaded and validated {dest_path}: {val_results['total_atoms']} atoms, {len(val_results['chains'])} chains.")
    return dest_path


def fetch_alphafold_db(uniprot_id: str, output_dir: str, file_format: str = "pdb") -> str:
    """
    Downloads predicted model from AlphaFold Database (EBI).
    uniprot_id: e.g. 'P04637', 'P69905'.
    file_format: 'pdb' or 'cif'.
    """
    uniprot_id = uniprot_id.strip().upper()
    os.makedirs(output_dir, exist_ok=True)

    filename = f"AF-{uniprot_id}-F1-model_v4.{file_format}"
    dest_path = os.path.join(output_dir, filename)
    url = f"https://alphafold.ebi.ac.uk/files/{filename}"

    print(f"Fetching AlphaFold model for {uniprot_id} from {url}...")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PyMOL-Visualize-Skill/1.0"})
        with urllib.request.urlopen(req, timeout=30) as response, open(dest_path, "wb") as out_file:
            out_file.write(response.read())
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"HTTP error {exc.code} fetching AlphaFold model for {uniprot_id}: {exc.reason}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Network error fetching AlphaFold model: {exc.reason}") from exc

    val_results = validate_structure(dest_path)
    print(f"Downloaded and validated {dest_path}: {val_results['total_atoms']} atoms, {len(val_results['chains'])} chains.")
    return dest_path


def main():
    parser = argparse.ArgumentParser(description="Fetch real structure coordinates from RCSB PDB or AlphaFold DB.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pdb", help="4-character PDB accession ID (e.g. 1UBQ, 4HHB)")
    group.add_argument("--uniprot", help="UniProt accession ID for AlphaFold prediction (e.g. P04637)")

    parser.add_argument("--output-dir", default="data/structures", help="Output directory to save coordinates")
    parser.add_argument("--format", choices=["cif", "pdb"], default="cif", help="Coordinate format (cif or pdb)")

    args = parser.parse_args()

    try:
        if args.pdb:
            dest = fetch_rcsb_pdb(args.pdb, args.output_dir, args.format)
        else:
            dest = fetch_alphafold_db(args.uniprot, args.output_dir, args.format)
        print(f"Target structure ready: {dest}")
        sys.exit(0)
    except Exception as exc:
        sys.stderr.write(f"Error: {exc}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
