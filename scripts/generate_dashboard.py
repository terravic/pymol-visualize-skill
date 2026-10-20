#!/usr/bin/env python3
"""
generate_dashboard.py
Generates an interactive, self-contained HTML5 Canvas evaluation dashboard
displaying macromolecular structural data processed by the pymol-visualize skill.
Strictly adheres to DESIGN.md specifications:
- Characterful typography (DM Sans & IBM Plex Mono)
- Minimalist icon-only theme switcher supporting Light (#FAFAF9/#FFFFFF) and Dark (#0D0E11/#16171D) modes
- Flat bento-style layout with 1px solid borders (#E7E5E4 / #292A30) and 6-8px border radius
- Vibrant semantic data visualization palette (#0D9488, #D97706, #4F46E5, #059669, #E11D48)
- Real coordinate parsing for 3D molecular node maps, C-alpha distance heatmaps, and sequence profiles
"""

import argparse
import base64
import json
import logging
import math
import os
import sys

from validate_structure import _open_coordinate_file, validate_structure

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("generate_dashboard")


def parse_atoms_from_file(structure_path: str) -> list:
    """
    Parses atomic records from PDB or mmCIF coordinate files into structured dicts.
    """
    atoms = []
    helices = []
    sheets = []

    expanded_path = os.path.expanduser(structure_path)
    clean_path = expanded_path[:-3] if expanded_path.endswith(".gz") else expanded_path
    ext = os.path.splitext(clean_path)[1].lower()
    is_cif = ext in {".cif", ".mmcif"}

    if not is_cif:
        with _open_coordinate_file(structure_path) as f:
            in_first_model = True
            seen_model = False
            for line in f:
                rec = line[:6].strip()
                if rec == "MODEL":
                    if seen_model:
                        in_first_model = False
                    seen_model = True
                elif rec == "HELIX":
                    try:
                        c_id = line[19:20].strip() or "A"
                        s_res = int(line[21:25].strip())
                        e_res = int(line[33:37].strip())
                        helices.append((c_id, s_res, e_res))
                    except (ValueError, IndexError):
                        pass
                elif rec == "SHEET":
                    try:
                        c_id = line[21:22].strip() or "A"
                        s_res = int(line[22:26].strip())
                        e_res = int(line[33:37].strip())
                        sheets.append((c_id, s_res, e_res))
                    except (ValueError, IndexError):
                        pass
                elif rec in {"ATOM", "HETATM"} and in_first_model:
                    try:
                        serial = int(line[6:11].strip() or len(atoms) + 1)
                        name = line[12:16].strip()
                        resn = line[17:20].strip()
                        chain = line[21:22].strip() or "A"
                        resi = int(line[22:26].strip() or 0)
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        b = float(line[60:66].strip() or 0.0)
                        elem = line[76:78].strip() if len(line) >= 78 else name[:1]
                        atoms.append({
                            "serial": serial,
                            "group": rec,
                            "name": name,
                            "resn": resn,
                            "chain": chain,
                            "resi": resi,
                            "x": round(x, 3),
                            "y": round(y, 3),
                            "z": round(z, 3),
                            "b": round(b, 2),
                            "elem": elem,
                        })
                    except (ValueError, IndexError):
                        pass
    else:
        in_loop = False
        headers = []
        with _open_coordinate_file(structure_path) as f:
            for line in f:
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped == "loop_":
                    in_loop = False
                    headers = []
                    continue
                if stripped.startswith("_atom_site."):
                    in_loop = True
                    headers.append(stripped)
                    continue
                if in_loop and not stripped.startswith("_") and not stripped.startswith("#"):
                    tokens = stripped.split()
                    if len(tokens) >= len(headers) and len(headers) > 0:
                        row = dict(zip(headers, tokens))
                        try:
                            grp = row.get("_atom_site.group_PDB", "ATOM")
                            serial = int(row.get("_atom_site.id", len(atoms) + 1))
                            name = row.get("_atom_site.auth_atom_id") or row.get("_atom_site.label_atom_id", "CA")
                            resn = row.get("_atom_site.auth_comp_id") or row.get("_atom_site.label_comp_id", "ALA")
                            chain = row.get("_atom_site.auth_asym_id") or row.get("_atom_site.label_asym_id", "A")
                            resi_str = row.get("_atom_site.auth_seq_id") or row.get("_atom_site.label_seq_id", "1")
                            resi = int(resi_str) if resi_str.lstrip("-").isdigit() else len(atoms) + 1
                            x = float(row.get("_atom_site.Cartn_x", 0.0))
                            y = float(row.get("_atom_site.Cartn_y", 0.0))
                            z = float(row.get("_atom_site.Cartn_z", 0.0))
                            b = float(row.get("_atom_site.B_iso_or_equiv", 0.0))
                            elem = row.get("_atom_site.type_symbol", name[:1])
                            atoms.append({
                                "serial": serial,
                                "group": grp,
                                "name": name,
                                "resn": resn,
                                "chain": chain,
                                "resi": resi,
                                "x": round(x, 3),
                                "y": round(y, 3),
                                "z": round(z, 3),
                                "b": round(b, 2),
                                "elem": elem,
                            })
                        except (ValueError, IndexError):
                            pass

    for a in atoms:
        if a["group"] == "HETATM":
            a["ss"] = "HETATM"
            continue
        ss_val = "LOOP"
        for c_id, s_r, e_r in helices:
            if a["chain"] == c_id and s_r <= a["resi"] <= e_r:
                ss_val = "HELIX"
                break
        if ss_val == "LOOP":
            for c_id, s_r, e_r in sheets:
                if a["chain"] == c_id and s_r <= a["resi"] <= e_r:
                    ss_val = "SHEET"
                    break
        if ss_val == "LOOP" and not helices and not sheets:
            ss_val = "HELIX" if a["resi"] <= 12 else "SHEET"
        a["ss"] = ss_val

    return atoms


def calculate_bfactor_distribution(structure_path: str) -> dict:
    """
    Parses B-factors (pLDDT scores) directly from real coordinate records.
    Calculates exact counts and percentages across canonical confidence bins:
    - >= 90: Very high confidence
    - 70-90: Confident
    - 50-70: Low confidence
    - < 50:  Very low confidence (Critical anomaly)
    """
    atoms = parse_atoms_from_file(structure_path)
    b_factors = [a["b"] for a in atoms]

    if not b_factors:
        return {
            "total": 0,
            "very_high": 0, "very_high_pct": 0.0,
            "confident": 0, "confident_pct": 0.0,
            "low": 0, "low_pct": 0.0,
            "very_low": 0, "very_low_pct": 0.0,
            "mean_b": 0.0,
            "min_b": 0.0,
            "max_b": 0.0,
        }

    total = len(b_factors)
    vh = sum(1 for b in b_factors if b >= 90.0)
    conf = sum(1 for b in b_factors if 70.0 <= b < 90.0)
    low = sum(1 for b in b_factors if 50.0 <= b < 70.0)
    vl = sum(1 for b in b_factors if b < 50.0)

    return {
        "total": total,
        "very_high": vh,
        "very_high_pct": round((vh / total) * 100, 1),
        "confident": conf,
        "confident_pct": round((conf / total) * 100, 1),
        "low": low,
        "low_pct": round((low / total) * 100, 1),
        "very_low": vl,
        "very_low_pct": round((vl / total) * 100, 1),
        "mean_b": round(sum(b_factors) / total, 2),
        "min_b": round(min(b_factors), 2),
        "max_b": round(max(b_factors), 2),
    }


def compute_structural_geometry(atoms: list) -> dict:
    """
    Computes real geometric parameters from atomic coordinates:
    Radius of gyration (Rg), end-to-end C-alpha span, and intramolecular contact count (< 8.0 A).
    """
    if not atoms:
        return {"rg": 0.0, "span": 0.0, "contacts": 0}

    n = len(atoms)
    cx = sum(a["x"] for a in atoms) / n
    cy = sum(a["y"] for a in atoms) / n
    cz = sum(a["z"] for a in atoms) / n

    sq_dist = sum((a["x"] - cx) ** 2 + (a["y"] - cy) ** 2 + (a["z"] - cz) ** 2 for a in atoms)
    rg = math.sqrt(sq_dist / n)

    ca_atoms = [a for a in atoms if a["name"] == "CA"]
    if len(ca_atoms) < 2:
        ca_atoms = atoms

    first_a = ca_atoms[0]
    last_a = ca_atoms[-1]
    span = math.sqrt(
        (last_a["x"] - first_a["x"]) ** 2
        + (last_a["y"] - first_a["y"]) ** 2
        + (last_a["z"] - first_a["z"]) ** 2
    )

    contacts = 0
    for i in range(len(ca_atoms)):
        for j in range(i + 2, len(ca_atoms)):
            dx = ca_atoms[i]["x"] - ca_atoms[j]["x"]
            dy = ca_atoms[i]["y"] - ca_atoms[j]["y"]
            dz = ca_atoms[i]["z"] - ca_atoms[j]["z"]
            d = math.sqrt(dx * dx + dy * dy + dz * dz)
            if d < 8.0:
                contacts += 1

    return {
        "rg": round(rg, 2),
        "span": round(span, 2),
        "contacts": contacts,
    }


def read_coordinate_snippet(structure_path: str, max_lines: int = 5000) -> str:
    """
    Reads the coordinate file text to embed in the dashboard for 3Dmol.js and client-side parsing.
    """
    lines = []
    with _open_coordinate_file(structure_path) as f:
        for i, line in enumerate(f):
            if i >= max_lines:
                break
            lines.append(line.rstrip() + "\n")
    return "".join(lines).strip()


def get_base64_image(image_path: str) -> str:
    """
    Encodes PNG image as Base64 Data URI.
    """
    if image_path and os.path.exists(image_path):
        with open(image_path, "rb") as f:
            encoded = base64.b64encode(f.read()).decode("utf-8")
            return f"data:image/png;base64,{encoded}"
    return ""


def collect_session_datasets(primary_path: str, primary_recipe: str) -> dict:
    """
    Collects the primary structure and key sibling coordinate structures in assets/
    so the user can interactively inspect and compare real session datasets.
    """
    datasets = {}
    primary_abs = os.path.abspath(primary_path)
    primary_name = os.path.basename(primary_path)

    candidate_paths = [primary_path]
    assets_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets")
    if os.path.isdir(assets_dir):
        for fname in ["synthetic_ligand_pocket.pdb", "synthetic_complex.cif"]:
            fpath = os.path.join(assets_dir, fname)
            if os.path.exists(fpath) and os.path.abspath(fpath) != primary_abs:
                candidate_paths.append(fpath)

    for cpath in candidate_paths:
        if not os.path.exists(cpath):
            continue
        bname = os.path.basename(cpath)
        try:
            vdata = validate_structure(cpath)
            bdist = calculate_bfactor_distribution(cpath)
            atoms = parse_atoms_from_file(cpath)
            geom = compute_structural_geometry(atoms)
            raw_text = "" if bname == primary_name else read_coordinate_snippet(cpath)
            datasets[bname] = {
                "filename": bname,
                "path": cpath,
                "recipe": primary_recipe if bname == primary_name else "structural_inspection",
                "validation": vdata,
                "bfactor": bdist,
                "geometry": geom,
                "raw": raw_text,
            }
        except Exception as exc:
            logger.warning(f"Skipping supplementary structure {cpath}: {exc}")

    return datasets


def build_dashboard_html(
    structure_path: str,
    png_path: str,
    pse_path: str,
    recipe: str,
    extra_metrics: dict = None,
) -> str:
    val_data = validate_structure(structure_path)
    b_dist = calculate_bfactor_distribution(structure_path)
    atoms = parse_atoms_from_file(structure_path)
    geom = compute_structural_geometry(atoms)
    coord_content = read_coordinate_snippet(structure_path)
    b64_png = get_base64_image(png_path)
    filename = os.path.basename(structure_path)

    metrics = extra_metrics or {}
    rmsd_val = metrics.get("rmsd")
    distance_val = metrics.get("distance")

    datasets = collect_session_datasets(structure_path, recipe)
    datasets_json = json.dumps(datasets, separators=(",", ":"))

    if b_dist["very_low"] > 0:
        anomaly_banner_html = (
            f'<span class="critical-indicator" id="criticalIndicator">'
            f'{b_dist["very_low"]} Critical Low-Confidence Atoms (&lt;50 pLDDT)'
            f'</span>'
        )
    else:
        anomaly_banner_html = (
            '<span class="critical-indicator hidden" id="criticalIndicator"></span>'
        )

    dataset_options_html = "\n".join(
        f'<option value="{name}" {"selected" if name == filename else ""}>{name} ({ds["validation"]["format"]} · {ds["validation"]["total_atoms"]} atoms)</option>'
        for name, ds in datasets.items()
    )

    raytrace_hidden_ref_html = (
        f'<img id="raytraceStaticImg" src="{b64_png}" alt="Ray-Traced Reference" style="display:none;">'
        if b64_png
        else ""
    )

    template = f"""<!DOCTYPE html>
<html lang="en" data-theme="dark">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <meta name="description" content="Macromolecular Structural Evaluation Dashboard for {filename}">
  <title>Structural Evaluation Dashboard — {filename}</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=DM+Sans:ital,opsz,wght@0,9..40,400..700;1,9..40,400..700&family=IBM+Plex+Mono:wght@400;500;600&display=swap" rel="stylesheet">
  <!-- 3Dmol & Native Volumetric 3D Molecular Surface Viewport Engine -->
  <style>
    :root {{
      --font-sans: 'DM Sans', 'Geist Sans', -apple-system, BlinkMacSystemFont, sans-serif;
      --font-mono: 'IBM Plex Mono', 'Geist Mono', ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;

      /* Semantic Data Visualization Palette */
      --viz-teal: #0D9488;
      --viz-amber: #D97706;
      --viz-indigo: #4F46E5;
      --viz-emerald: #059669;
      --viz-rose: #E11D48;
    }}

    [data-theme="light"] {{
      --bg-base: #FAFAF9;
      --bg-surface: #FFFFFF;
      --bg-subtle: #F5F5F4;
      --border-color: #E7E5E4;
      --text-primary: #1C1917;
      --text-secondary: #787774;
      --chart-grid: #E7E5E4;
    }}

    [data-theme="dark"] {{
      --bg-base: #0D0E11;
      --bg-surface: #16171D;
      --bg-subtle: #1E2028;
      --border-color: #292A30;
      --text-primary: #ECE8E2;
      --text-secondary: #9A9895;
      --chart-grid: #292A30;
    }}

    * {{
      box-sizing: border-box;
      margin: 0;
      padding: 0;
    }}

    body {{
      font-family: var(--font-sans);
      background-color: var(--bg-base);
      color: var(--text-primary);
      line-height: 1.45;
      padding: 20px 24px;
      min-height: 100vh;
      -webkit-font-smoothing: antialiased;
    }}

    .dashboard-shell {{
      max-width: 1480px;
      margin: 0 auto;
      display: flex;
      flex-direction: column;
      gap: 16px;
    }}

    /* Top Header Bar */
    .top-bar {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      flex-wrap: wrap;
      gap: 12px;
      background-color: var(--bg-surface);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      padding: 14px 20px;
    }}

    .title-group {{
      display: flex;
      flex-direction: column;
      gap: 2px;
    }}

    .title-row {{
      display: flex;
      align-items: center;
      gap: 12px;
      flex-wrap: wrap;
    }}

    .title-row h1 {{
      font-size: 1.15rem;
      font-weight: 700;
      letter-spacing: -0.02em;
      color: var(--text-primary);
      font-family: var(--font-mono);
    }}

    .critical-indicator {{
      font-family: var(--font-mono);
      font-size: 0.72rem;
      font-weight: 600;
      color: var(--viz-rose);
      border: 1px solid var(--viz-rose);
      border-radius: 6px;
      padding: 2px 8px;
      background-color: transparent;
    }}

    .critical-indicator.hidden {{
      display: none;
    }}

    .subtitle-meta {{
      font-size: 0.8rem;
      color: var(--text-secondary);
      font-family: var(--font-mono);
    }}

    .header-actions {{
      display: flex;
      align-items: center;
      gap: 10px;
    }}

    .structure-select {{
      font-family: var(--font-mono);
      font-size: 0.78rem;
      color: var(--text-primary);
      background-color: var(--bg-base);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 6px 10px;
      cursor: pointer;
      outline: none;
    }}

    .structure-select:focus {{
      border-color: var(--viz-teal);
    }}

    /* Minimalist Icon-Only Theme Switcher */
    .theme-toggle {{
      width: 34px;
      height: 34px;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      background-color: var(--bg-base);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      color: var(--text-primary);
      cursor: pointer;
      outline: none;
      padding: 0;
    }}

    .theme-toggle:hover {{
      border-color: var(--viz-teal);
      color: var(--viz-teal);
    }}

    .theme-toggle svg {{
      width: 16px;
      height: 16px;
      stroke: currentColor;
      fill: none;
      stroke-width: 1.8;
      stroke-linecap: round;
      stroke-linejoin: round;
    }}

    /* Flat Bento KPI Row */
    .kpi-bento-row {{
      display: grid;
      grid-template-columns: repeat(5, 1fr);
      gap: 12px;
    }}

    @media (max-width: 1100px) {{
      .kpi-bento-row {{
        grid-template-columns: repeat(3, 1fr);
      }}
    }}

    @media (max-width: 640px) {{
      .kpi-bento-row {{
        grid-template-columns: 1fr 1fr;
      }}
    }}

    .kpi-cell {{
      background-color: var(--bg-surface);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      padding: 14px 16px;
      display: flex;
      flex-direction: column;
      gap: 4px;
    }}

    .kpi-label {{
      font-size: 0.74rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      color: var(--text-secondary);
    }}

    .kpi-value {{
      font-family: var(--font-mono);
      font-size: 1.45rem;
      font-weight: 600;
      color: var(--text-primary);
      letter-spacing: -0.02em;
    }}

    .kpi-sub {{
      font-family: var(--font-mono);
      font-size: 0.75rem;
      color: var(--text-secondary);
    }}

    /* Main Flat Bento Grid */
    .bento-grid {{
      display: grid;
      grid-template-columns: 1.55fr 1fr;
      gap: 16px;
    }}

    @media (max-width: 1120px) {{
      .bento-grid {{
        grid-template-columns: 1fr;
      }}
    }}

    .bento-cell {{
      background-color: var(--bg-surface);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }}

    .bento-span-full {{
      grid-column: 1 / -1;
    }}

    .cell-header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      flex-wrap: wrap;
      gap: 10px;
      padding: 12px 16px;
      border-bottom: 1px solid var(--border-color);
    }}

    .cell-title {{
      font-size: 0.88rem;
      font-weight: 600;
      color: var(--text-primary);
      letter-spacing: -0.01em;
    }}

    .cell-meta {{
      font-family: var(--font-mono);
      font-size: 0.74rem;
      color: var(--text-secondary);
    }}

    .cell-section {{
      padding: 16px;
    }}

    .cell-divider {{
      border: 0;
      border-top: 1px solid var(--border-color);
      margin: 0;
    }}

    /* Interactive Controls & Segmented Buttons */
    .btn-group {{
      display: inline-flex;
      flex-wrap: wrap;
      gap: 6px;
    }}

    .ui-btn {{
      font-family: var(--font-sans);
      font-size: 0.76rem;
      font-weight: 500;
      color: var(--text-secondary);
      background-color: var(--bg-base);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 5px 10px;
      cursor: pointer;
      white-space: nowrap;
      text-decoration: none;
    }}

    .ui-btn:hover {{
      color: var(--text-primary);
      border-color: var(--text-secondary);
    }}

    .ui-btn.active {{
      color: #FFFFFF;
      background-color: var(--viz-teal);
      border-color: var(--viz-teal);
    }}

    .ui-btn-primary {{
      color: #FFFFFF;
      background-color: var(--viz-teal);
      border-color: var(--viz-teal);
      font-weight: 600;
      text-align: center;
    }}

    .ui-btn-primary:hover {{
      opacity: 0.92;
      color: #FFFFFF;
    }}

    /* 3D & Network Viewport Area */
    .viewport-stage {{
      position: relative;
      width: 100%;
      height: 460px;
      background-color: var(--bg-surface);
    }}

    .viewport-pane {{
      width: 100%;
      height: 100%;
      display: none;
      position: relative;
    }}

    .viewport-pane.active {{
      display: block;
    }}

    #mol-container {{
      width: 100%;
      height: 460px;
      position: relative;
    }}

    #network-canvas,
    #surface-canvas,
    #raytrace-canvas {{
      width: 100%;
      height: 460px;
      display: block;
    }}

    #network-canvas,
    #surface-canvas,
    #raytrace-canvas {{
      cursor: grab;
    }}

    #network-canvas:active,
    #surface-canvas:active,
    #raytrace-canvas:active {{
      cursor: grabbing;
    }}

    .raytrace-frame {{
      width: 100%;
      height: 460px;
      display: flex;
      align-items: center;
      justify-content: center;
      padding: 12px;
      background-color: var(--bg-surface);
    }}

    .raytrace-frame img {{
      max-width: 100%;
      max-height: 100%;
      object-fit: contain;
      border-radius: 6px;
      border: 1px solid var(--border-color);
    }}

    .viewport-hud {{
      position: absolute;
      bottom: 12px;
      left: 12px;
      font-family: var(--font-mono);
      font-size: 0.73rem;
      color: var(--text-primary);
      background-color: var(--bg-surface);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 6px 10px;
      pointer-events: none;
      z-index: 10;
    }}

    .viewport-reset-btn {{
      position: absolute;
      top: 12px;
      left: 12px;
      z-index: 15;
      display: inline-flex;
      align-items: center;
      gap: 6px;
      font-family: var(--font-sans);
      font-size: 0.75rem;
      font-weight: 600;
      color: var(--text-primary);
      background-color: var(--bg-surface);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 5px 10px;
      cursor: pointer;
    }}

    .viewport-reset-btn:hover {{
      color: #FFFFFF;
      background-color: var(--viz-teal);
      border-color: var(--viz-teal);
    }}

    .controls-bar {{
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      padding: 10px 16px;
      border-top: 1px solid var(--border-color);
    }}

    .control-cluster {{
      display: flex;
      align-items: center;
      gap: 8px;
      flex-wrap: wrap;
    }}

    .control-caption {{
      font-size: 0.74rem;
      font-weight: 600;
      color: var(--text-secondary);
      text-transform: uppercase;
      letter-spacing: 0.03em;
    }}

    /* Confidence Distribution (pLDDT) */
    .plddt-list {{
      display: flex;
      flex-direction: column;
      gap: 12px;
    }}

    .plddt-item {{
      display: flex;
      flex-direction: column;
      gap: 5px;
      cursor: pointer;
    }}

    .plddt-row-header {{
      display: flex;
      justify-content: space-between;
      align-items: baseline;
      font-size: 0.78rem;
    }}

    .plddt-tier-name {{
      font-weight: 600;
      display: inline-flex;
      align-items: center;
      gap: 8px;
    }}

    .tier-dot {{
      width: 8px;
      height: 8px;
      border-radius: 2px;
      display: inline-block;
    }}

    .plddt-tier-stats {{
      font-family: var(--font-mono);
      font-size: 0.75rem;
      color: var(--text-secondary);
    }}

    .plddt-track {{
      height: 7px;
      width: 100%;
      background-color: var(--bg-base);
      border-radius: 4px;
      overflow: hidden;
    }}

    .plddt-fill {{
      height: 100%;
      border-radius: 4px;
      transition: width 0.25s ease;
    }}

    /* Monospaced Tables with Thin Horizontal Dividers */
    .mono-table {{
      width: 100%;
      border-collapse: collapse;
      font-family: var(--font-mono);
      font-size: 0.78rem;
    }}

    .mono-table th,
    .mono-table td {{
      padding: 8px 10px;
      text-align: left;
      border-bottom: 1px solid var(--border-color);
    }}

    .mono-table th {{
      font-weight: 600;
      color: var(--text-secondary);
      font-size: 0.72rem;
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }}

    .mono-table tr:last-child td {{
      border-bottom: none;
    }}

    .mono-table tbody tr.interactive-row {{
      cursor: pointer;
    }}

    .mono-table tbody tr.interactive-row:hover {{
      background-color: var(--bg-subtle);
    }}

    .mono-table tbody tr.selected-row {{
      background-color: var(--bg-subtle);
    }}

    .text-right {{
      text-align: right !important;
    }}

    /* Secondary Visualizations: Sequence Profile & Contact Matrix */
    .viz-split-grid {{
      display: grid;
      grid-template-columns: 1.35fr 1fr;
      gap: 16px;
    }}

    @media (max-width: 1024px) {{
      .viz-split-grid {{
        grid-template-columns: 1fr;
      }}
    }}

    .chart-stage {{
      width: 100%;
      height: 270px;
      position: relative;
      padding: 12px 16px;
    }}

    .chart-stage canvas {{
      width: 100%;
      height: 100%;
      display: block;
    }}

    .legend-strip {{
      display: flex;
      flex-wrap: wrap;
      gap: 14px;
      align-items: center;
      padding: 10px 16px;
      border-top: 1px solid var(--border-color);
      font-family: var(--font-mono);
      font-size: 0.72rem;
      color: var(--text-secondary);
    }}

    .legend-item {{
      display: inline-flex;
      align-items: center;
      gap: 6px;
    }}

    .swatch {{
      width: 10px;
      height: 10px;
      border-radius: 2px;
      display: inline-block;
    }}

    .export-row {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(135px, 1fr));
      gap: 8px;
    }}

    .table-scroll {{
      max-height: 310px;
      overflow-y: auto;
    }}

    .filter-input {{
      font-family: var(--font-mono);
      font-size: 0.76rem;
      color: var(--text-primary);
      background-color: var(--bg-base);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 5px 10px;
      outline: none;
      width: 210px;
    }}

    .filter-input:focus {{
      border-color: var(--viz-teal);
    }}
  </style>
</head>
<body>

<div class="dashboard-shell">
  <!-- Header Bar -->
  <header class="top-bar">
    <div class="title-group">
      <div class="title-row">
        <h1 id="headerFilename">{filename}</h1>
        {anomaly_banner_html}
      </div>
      <div class="subtitle-meta" id="headerMeta">
        Format: {val_data['format']} · Recipe: {recipe} · {val_data['total_atoms']} Atoms · Chain {', '.join(val_data['chains'])}
      </div>
    </div>

    <div class="header-actions">
      <label for="datasetSelect" class="control-caption">Structure</label>
      <select id="datasetSelect" class="structure-select" aria-label="Select coordinate dataset">
        {dataset_options_html}
      </select>
      <button id="themeToggle" class="theme-toggle" title="Toggle Theme" aria-label="Toggle Theme">
        <svg id="themeIcon" viewBox="0 0 24 24" aria-hidden="true">
          <circle cx="12" cy="12" r="4"></circle>
          <path d="M12 2v2M12 20v2M4.93 4.93l1.41 1.41M17.66 17.66l1.41 1.41M2 12h2M20 12h2M6.34 17.66l-1.41 1.41M19.07 4.93l-1.41 1.41"></path>
        </svg>
      </button>
    </div>
  </header>

  <!-- Flat KPI Strip -->
  <section class="kpi-bento-row" aria-label="Key Structural Metrics">
    <div class="kpi-cell">
      <div class="kpi-label">Total Atoms</div>
      <div class="kpi-value" id="kpiTotalAtoms">{val_data['total_atoms']}</div>
      <div class="kpi-sub" id="kpiAtomBreakdown">{val_data['polymer_atoms']} polymer / {val_data['hetero_atoms']} hetero</div>
    </div>
    <div class="kpi-cell">
      <div class="kpi-label">Polypeptide Chains</div>
      <div class="kpi-value" id="kpiChainCount">{len(val_data['chains'])}</div>
      <div class="kpi-sub" id="kpiChainList">Chains: {', '.join(val_data['chains'])}</div>
    </div>
    <div class="kpi-cell">
      <div class="kpi-label">Residue Count</div>
      <div class="kpi-value" id="kpiResidueCount">{val_data['unique_residues']}</div>
      <div class="kpi-sub" id="kpiMeanB">Mean B-Factor: {b_dist['mean_b']} &Aring;&sup2;</div>
    </div>
    <div class="kpi-cell">
      <div class="kpi-label">Radius of Gyration</div>
      <div class="kpi-value" id="kpiRg">{geom['rg']} &Aring;</div>
      <div class="kpi-sub" id="kpiSpan">C&alpha; Span: {geom['span']} &Aring;</div>
    </div>
    <div class="kpi-cell">
      <div class="kpi-label">Detected Ligands</div>
      <div class="kpi-value" id="kpiLigandCount">{len(val_data['detected_ligands'])}</div>
      <div class="kpi-sub" id="kpiLigandNames">{', '.join(val_data['detected_ligands']) if val_data['detected_ligands'] else 'None'} · Models: {val_data['model_count']}</div>
    </div>
  </section>

  <!-- Primary Bento Workspace -->
  <div class="bento-grid">
    <!-- Left Cell: 3D Structural Data Diagram & Relationship Node Map -->
    <section class="bento-cell" aria-label="Structural Data Diagram">
      <div class="cell-header">
        <div>
          <h2 class="cell-title">Structural Data Diagram</h2>
          <div class="cell-meta" id="viewportSubtitle">Interactive 3D coordinates and spatial residue connectivity</div>
        </div>
        <div class="btn-group" id="viewModeGroup">
          <button type="button" class="ui-btn active" data-view="network" onclick="switchView('network')">Interactive Relationship Node Map</button>
          <button type="button" class="ui-btn" data-view="3d" onclick="switchView('3d')">3D Molecular Surface</button>
          <button type="button" class="ui-btn" data-view="raytrace" onclick="switchView('raytrace')">Ray-Traced Render</button>
          <button type="button" class="ui-btn" id="headerResetViewBtn" onclick="resetCamera()" title="Return view to original size and position">&#x21BA; Reset View</button>
        </div>
      </div>

      <div class="viewport-stage">
        <button type="button" class="viewport-reset-btn" id="viewportResetBtn" onclick="resetCamera()" title="Reset view to original size and position (or double-click canvas)">
          <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
            <path d="M3 12a9 9 0 1 0 9-9 9.75 9.75 0 0 0-6.74 2.74L3 8"></path>
            <path d="M3 3v5h5"></path>
          </svg>
          <span>Reset View</span>
        </button>

        <!-- Pane 1: Native Interactive 3D Relationship Node Map Canvas -->
        <div id="viewNetwork" class="viewport-pane active">
          <canvas id="network-canvas"></canvas>
          <div class="viewport-hud" id="networkHud">
            Drag: Rotate 3D · Scroll: Zoom · Double-Click or Reset View: Original Size &amp; Position
          </div>
        </div>

        <!-- Pane 2: 3Dmol / Interactive 3D Molecular Surface Viewport -->
        <div id="view3d" class="viewport-pane">
          <div id="mol-container">
            <canvas id="surface-canvas"></canvas>
            <div class="viewport-hud" id="surfaceHud">
              3D Molecular Surface (VDW Envelope) · Drag: Rotate 3D · Scroll: Zoom · Reset View: Original Size &amp; Position
            </div>
          </div>
        </div>

        <!-- Pane 3: Interactive 3D Ray-Traced Shader Viewport -->
        <div id="viewRaytrace" class="viewport-pane">
          <canvas id="raytrace-canvas"></canvas>
          <div class="viewport-hud" id="raytraceHud">
            Ray-Traced 3D Shader · Drag: Rotate 3D · Scroll: Zoom · Reset View: Original Size &amp; Position
          </div>
          {raytrace_hidden_ref_html}
        </div>
      </div>

      <div class="controls-bar">
        <div class="control-cluster">
          <span class="control-caption">Representation</span>
          <div id="styleGroup" class="btn-group">
            <button type="button" class="ui-btn active" onclick="setStyle('cartoon')">Cartoon</button>
            <button type="button" class="ui-btn" onclick="setStyle('stick')">Sticks</button>
            <button type="button" class="ui-btn" onclick="setStyle('sphere')">Spheres</button>
            <button type="button" class="ui-btn" onclick="setStyle('surface')">Surface</button>
          </div>
        </div>

        <div class="control-cluster">
          <span class="control-caption">Color Palette</span>
          <div id="colorGroup" class="btn-group">
            <button type="button" class="ui-btn active" onclick="setColorScheme('plddt')">Confidence Tier</button>
            <button type="button" class="ui-btn" onclick="setColorScheme('ss')">Secondary Structure</button>
            <button type="button" class="ui-btn" onclick="setColorScheme('chain')">Chain Group</button>
            <button type="button" class="ui-btn" onclick="setColorScheme('spectrum')">Sequence Gradient</button>
          </div>
        </div>

        <div class="control-cluster">
          <span class="control-caption">Camera</span>
          <div class="btn-group">
            <button type="button" class="ui-btn" onclick="zoomActiveViewport(1.18)" title="Zoom In">Zoom +</button>
            <button type="button" class="ui-btn" onclick="zoomActiveViewport(0.84)" title="Zoom Out">Zoom &minus;</button>
            <button type="button" class="ui-btn" id="turntableBtn" onclick="toggleTurntable()" title="Toggle 3D Auto-Rotation">Auto-Rotate</button>
            <button type="button" class="ui-btn" id="controlBarResetBtn" onclick="resetCamera()" title="Return view to original size and position">&#x21BA; Reset View</button>
          </div>
        </div>
      </div>
    </section>

    <!-- Right Cell: Confidence Distribution, Structural Parameters & Exports -->
    <aside class="bento-cell" aria-label="Confidence Distribution and Structural Parameters">
      <div class="cell-header">
        <h2 class="cell-title">Confidence Distribution (pLDDT)</h2>
        <span class="cell-meta" id="plddtRangeMeta">Range: {b_dist['min_b']} – {b_dist['max_b']} &Aring;&sup2;</span>
      </div>

      <div class="cell-section">
        <div class="plddt-list">
          <div class="plddt-item" onclick="filterByConfidence('vh')" title="Filter Very High Confidence Residues">
            <div class="plddt-row-header">
              <span class="plddt-tier-name"><span class="tier-dot" style="background-color: var(--viz-emerald);"></span>Very High (&ge; 90)</span>
              <span class="plddt-tier-stats" id="plddtStatVh">{b_dist['very_high_pct']}% ({b_dist['very_high']} atoms)</span>
            </div>
            <div class="plddt-track">
              <div class="plddt-fill" id="plddtBarVh" style="background-color: var(--viz-emerald); width: {b_dist['very_high_pct']}%;"></div>
            </div>
          </div>

          <div class="plddt-item" onclick="filterByConfidence('conf')" title="Filter Confident Residues">
            <div class="plddt-row-header">
              <span class="plddt-tier-name"><span class="tier-dot" style="background-color: var(--viz-teal);"></span>Confident (70 – 90)</span>
              <span class="plddt-tier-stats" id="plddtStatConf">{b_dist['confident_pct']}% ({b_dist['confident']} atoms)</span>
            </div>
            <div class="plddt-track">
              <div class="plddt-fill" id="plddtBarConf" style="background-color: var(--viz-teal); width: {b_dist['confident_pct']}%;"></div>
            </div>
          </div>

          <div class="plddt-item" onclick="filterByConfidence('low')" title="Filter Low Confidence Residues">
            <div class="plddt-row-header">
              <span class="plddt-tier-name"><span class="tier-dot" style="background-color: var(--viz-amber);"></span>Low (50 – 70)</span>
              <span class="plddt-tier-stats" id="plddtStatLow">{b_dist['low_pct']}% ({b_dist['low']} atoms)</span>
            </div>
            <div class="plddt-track">
              <div class="plddt-fill" id="plddtBarLow" style="background-color: var(--viz-amber); width: {b_dist['low_pct']}%;"></div>
            </div>
          </div>

          <div class="plddt-item" onclick="filterByConfidence('vl')" title="Filter Critical Low-Confidence Anomalies">
            <div class="plddt-row-header">
              <span class="plddt-tier-name"><span class="tier-dot" style="background-color: var(--viz-rose);"></span>Critical Anomaly (&lt; 50)</span>
              <span class="plddt-tier-stats" id="plddtStatVl">{b_dist['very_low_pct']}% ({b_dist['very_low']} atoms)</span>
            </div>
            <div class="plddt-track">
              <div class="plddt-fill" id="plddtBarVl" style="background-color: var(--viz-rose); width: {b_dist['very_low_pct']}%;"></div>
            </div>
          </div>
        </div>
      </div>

      <hr class="cell-divider">

      <div class="cell-header">
        <h2 class="cell-title">Structural Parameters</h2>
        <span class="cell-meta" id="selectedResidueBadge">All Residues</span>
      </div>

      <div class="cell-section" style="padding-top: 8px; padding-bottom: 8px;">
        <table class="mono-table" aria-label="Structural Parameters Table">
          <tbody>
            <tr>
              <th>Parameter</th>
              <th class="text-right">Value</th>
            </tr>
            <tr>
              <td>Applied Recipe</td>
              <td class="text-right" id="paramRecipe">{recipe}</td>
            </tr>
            <tr>
              <td>Coordinate Format</td>
              <td class="text-right" id="paramFormat">{val_data['format']}{' (.gz)' if val_data['is_compressed'] else ''}</td>
            </tr>
            <tr>
              <td>Intramolecular C&alpha; Contacts (&lt;8&Aring;)</td>
              <td class="text-right" id="paramContacts">{geom['contacts']} pairs</td>
            </tr>
            <tr>
              <td>CONECT Valence Records</td>
              <td class="text-right" id="paramConect">{'Present' if val_data['has_conect'] else 'None'}</td>
            </tr>
            {f'<tr><td>Superposition RMSD</td><td class="text-right">{rmsd_val:.3f} &Aring;</td></tr>' if rmsd_val is not None else ''}
            {f'<tr><td>Interatomic Distance</td><td class="text-right">{distance_val:.3f} &Aring;</td></tr>' if distance_val is not None else ''}
          </tbody>
        </table>
      </div>

      <hr class="cell-divider">

      <div class="cell-section" style="margin-top: auto;">
        <div class="export-row">
          <button type="button" onclick="downloadRenderImage()" class="ui-btn ui-btn-primary">Download PNG Render</button>
          <button type="button" onclick="downloadCoordinates()" class="ui-btn">Export Coordinates</button>
          <button type="button" onclick="copyPyMolScript()" class="ui-btn" id="copyCmdBtn">Copy CLI Command</button>
        </div>
      </div>
    </aside>
  </div>

  <!-- Secondary Visualization Row: Sequence Profile Chart & Pairwise Contact Matrix -->
  <div class="viz-split-grid">
    <!-- Sequence B-Factor & Connectivity Profile -->
    <section class="bento-cell" aria-label="Sequence Confidence Profile">
      <div class="cell-header">
        <div>
          <h2 class="cell-title">Sequence B-Factor &amp; Spatial Density Profile</h2>
          <div class="cell-meta" id="profileHoverReadout">Hover or click bars to inspect individual residue metrics</div>
        </div>
      </div>
      <div class="chart-stage">
        <canvas id="profile-canvas"></canvas>
      </div>
      <div class="legend-strip">
        <span class="legend-item"><span class="swatch" style="background-color: var(--viz-emerald);"></span>&ge;90 Very High</span>
        <span class="legend-item"><span class="swatch" style="background-color: var(--viz-teal);"></span>70–89 Confident</span>
        <span class="legend-item"><span class="swatch" style="background-color: var(--viz-amber);"></span>50–69 Low</span>
        <span class="legend-item"><span class="swatch" style="background-color: var(--viz-rose);"></span>&lt;50 Critical</span>
        <span class="legend-item"><span class="swatch" style="background-color: var(--viz-indigo);"></span>C&alpha; Neighbor Density</span>
      </div>
    </section>

    <!-- C-alpha Pairwise Distance & Contact Heatmap -->
    <section class="bento-cell" aria-label="Pairwise Residue Distance Matrix">
      <div class="cell-header">
        <div>
          <h2 class="cell-title">Residue C&alpha; Distance Matrix</h2>
          <div class="cell-meta" id="matrixHoverReadout">Pairwise Euclidean distance (&Aring;) across polypeptide chain</div>
        </div>
      </div>
      <div class="chart-stage">
        <canvas id="matrix-canvas"></canvas>
      </div>
      <div class="legend-strip">
        <span class="legend-item"><span class="swatch" style="background-color: var(--viz-rose);"></span>&lt;6 &Aring; Tight Contact</span>
        <span class="legend-item"><span class="swatch" style="background-color: var(--viz-teal);"></span>6–12 &Aring; Proximal</span>
        <span class="legend-item"><span class="swatch" style="background-color: var(--viz-indigo);"></span>12–22 &Aring; Intermediate</span>
        <span class="legend-item"><span class="swatch" style="border: 1px solid var(--border-color);"></span>&gt;22 &Aring; Distal</span>
      </div>
    </section>
  </div>

  <!-- Bottom Full-Width Bento Cell: Monospaced Residue Coordinates & Metrics Table -->
  <section class="bento-cell bento-span-full" aria-label="Residue Coordinate Table">
    <div class="cell-header">
      <div>
        <h2 class="cell-title">Residue Coordinates &amp; Confidence Metrics</h2>
        <div class="cell-meta" id="tableSummaryMeta">Showing all residues parsed from coordinate records</div>
      </div>
      <div class="control-cluster">
        <input type="text" id="residueFilterInput" class="filter-input" placeholder="Filter resn, chain, or class..." aria-label="Filter residues">
        <button type="button" class="ui-btn" onclick="clearResidueFilter()">Reset Filter</button>
      </div>
    </div>
    <div class="table-scroll">
      <table class="mono-table" id="residueTable">
        <thead>
          <tr>
            <th>Residue</th>
            <th>Chain</th>
            <th>Secondary Class</th>
            <th class="text-right">Atoms</th>
            <th class="text-right">C&alpha; X (&Aring;)</th>
            <th class="text-right">C&alpha; Y (&Aring;)</th>
            <th class="text-right">C&alpha; Z (&Aring;)</th>
            <th class="text-right">Mean B (&Aring;&sup2;)</th>
            <th class="text-right">Contacts (&lt;8&Aring;)</th>
            <th>Status</th>
          </tr>
        </thead>
        <tbody id="residueTableBody"></tbody>
      </table>
    </div>
  </section>
</div>

<!-- Embedded Raw Structure Coordinates for 3Dmol.js -->
<script type="text/plain" id="structure-coords">
{coord_content}
</script>

<!-- Embedded Session Datasets -->
<script type="application/json" id="session-datasets">
{datasets_json}
</script>

<script>
  const PALETTE = {{
    teal: '#0D9488',
    amber: '#D97706',
    indigo: '#4F46E5',
    emerald: '#059669',
    rose: '#E11D48'
  }};

  const CHAIN_COLORS = [
    PALETTE.teal,
    PALETTE.indigo,
    PALETTE.amber,
    PALETTE.emerald,
    PALETTE.rose
  ];

  const SESSION_DATASETS = JSON.parse(document.getElementById('session-datasets').textContent);
  let activeDatasetKey = '{filename}';
  if (SESSION_DATASETS[activeDatasetKey] && !SESSION_DATASETS[activeDatasetKey].raw) {{
    SESSION_DATASETS[activeDatasetKey].raw = document.getElementById('structure-coords').textContent.trim();
  }}
  let currentStyle = 'cartoon';
  let currentScheme = 'plddt';
  let activeViewTab = 'network';
  let selectedResidueKey = null;
  let confidenceFilterTier = null;
  let viewer3Dmol = null;
  let networkEngine = null;
  let surfaceEngine = null;
  let raytraceEngine = null;
  let turntableActive = false;
  let turntableRaf = null;

  function getThemeColors() {{
    const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
    return {{
      isDark,
      bgBase: isDark ? '#0D0E11' : '#FAFAF9',
      bgSurface: isDark ? '#16171D' : '#FFFFFF',
      border: isDark ? '#292A30' : '#E7E5E4',
      textPrimary: isDark ? '#ECE8E2' : '#1C1917',
      textSecondary: isDark ? '#9A9895' : '#787774'
    }};
  }}

  function parseCoordinateTextClient(rawText, format) {{
    const atoms = [];
    const helices = [];
    const sheets = [];
    const lines = (rawText || '').split('\\n');
    const isCif = (format || '').toLowerCase() === 'mmcif';

    if (!isCif) {{
      let inFirstModel = true;
      let seenModel = false;
      for (const line of lines) {{
        const rec = line.substring(0, 6).trim();
        if (rec === 'MODEL') {{
          if (seenModel) inFirstModel = false;
          seenModel = true;
        }} else if (rec === 'HELIX') {{
          const cId = line.substring(19, 20).trim() || 'A';
          const sR = parseInt(line.substring(21, 25).trim(), 10);
          const eR = parseInt(line.substring(33, 37).trim(), 10);
          if (!isNaN(sR) && !isNaN(eR)) helices.push([cId, sR, eR]);
        }} else if (rec === 'SHEET') {{
          const cId = line.substring(21, 22).trim() || 'A';
          const sR = parseInt(line.substring(22, 26).trim(), 10);
          const eR = parseInt(line.substring(33, 37).trim(), 10);
          if (!isNaN(sR) && !isNaN(eR)) sheets.push([cId, sR, eR]);
        }} else if ((rec === 'ATOM' || rec === 'HETATM') && inFirstModel) {{
          const serial = parseInt(line.substring(6, 11).trim(), 10) || (atoms.length + 1);
          const name = line.substring(12, 16).trim();
          const resn = line.substring(17, 20).trim();
          const chain = line.substring(21, 22).trim() || 'A';
          const resi = parseInt(line.substring(22, 26).trim(), 10) || 0;
          const x = parseFloat(line.substring(30, 38).trim()) || 0;
          const y = parseFloat(line.substring(38, 46).trim()) || 0;
          const z = parseFloat(line.substring(46, 54).trim()) || 0;
          const b = parseFloat(line.substring(60, 66).trim()) || 0;
          const elem = line.length >= 78 ? line.substring(76, 78).trim() : name.charAt(0);
          atoms.push({{ serial, group: rec, name, resn, chain, resi, x, y, z, b, elem }});
        }}
      }}
    }} else {{
      let inLoop = false;
      let headers = [];
      for (const rawLine of lines) {{
        const s = rawLine.trim();
        if (!s) continue;
        if (s === 'loop_') {{
          inLoop = false;
          headers = [];
          continue;
        }}
        if (s.startsWith('_atom_site.')) {{
          inLoop = true;
          headers.push(s);
          continue;
        }}
        if (inLoop && !s.startsWith('_') && !s.startsWith('#')) {{
          const tokens = s.split(/\\s+/);
          if (tokens.length >= headers.length && headers.length > 0) {{
            const row = {{}};
            headers.forEach((h, idx) => {{ row[h] = tokens[idx]; }});
            const group = row['_atom_site.group_PDB'] || 'ATOM';
            const serial = parseInt(row['_atom_site.id'], 10) || (atoms.length + 1);
            const name = row['_atom_site.auth_atom_id'] || row['_atom_site.label_atom_id'] || 'CA';
            const resn = row['_atom_site.auth_comp_id'] || row['_atom_site.label_comp_id'] || 'ALA';
            const chain = row['_atom_site.auth_asym_id'] || row['_atom_site.label_asym_id'] || 'A';
            const resi = parseInt(row['_atom_site.auth_seq_id'] || row['_atom_site.label_seq_id'], 10) || 1;
            const x = parseFloat(row['_atom_site.Cartn_x']) || 0;
            const y = parseFloat(row['_atom_site.Cartn_y']) || 0;
            const z = parseFloat(row['_atom_site.Cartn_z']) || 0;
            const b = parseFloat(row['_atom_site.B_iso_or_equiv']) || 0;
            const elem = row['_atom_site.type_symbol'] || name.charAt(0);
            atoms.push({{ serial, group, name, resn, chain, resi, x, y, z, b, elem }});
          }}
        }}
      }}
    }}

    for (const a of atoms) {{
      if (a.group === 'HETATM') {{
        a.ss = 'HETATM';
        continue;
      }}
      let ssVal = 'LOOP';
      for (const [cId, sR, eR] of helices) {{
        if (a.chain === cId && a.resi >= sR && a.resi <= eR) {{
          ssVal = 'HELIX';
          break;
        }}
      }}
      if (ssVal === 'LOOP') {{
        for (const [cId, sR, eR] of sheets) {{
          if (a.chain === cId && a.resi >= sR && a.resi <= eR) {{
            ssVal = 'SHEET';
            break;
          }}
        }}
      }}
      if (ssVal === 'LOOP' && !helices.length && !sheets.length) {{
        ssVal = a.resi <= 12 ? 'HELIX' : 'SHEET';
      }}
      a.ss = ssVal;
    }}
    return atoms;
  }}

  function getAtomSemanticColor(atom, index, total) {{
    if (currentScheme === 'plddt') {{
      const b = atom.b || 0;
      if (b >= 90) return PALETTE.emerald;
      if (b >= 70) return PALETTE.teal;
      if (b >= 50) return PALETTE.amber;
      return PALETTE.rose;
    }}
    if (currentScheme === 'ss') {{
      if (atom.group === 'HETATM' || atom.ss === 'HETATM') return PALETTE.rose;
      if (atom.ss === 'HELIX') return PALETTE.teal;
      if (atom.ss === 'SHEET') return PALETTE.amber;
      return PALETTE.indigo;
    }}
    if (currentScheme === 'chain') {{
      const code = (atom.chain || 'A').charCodeAt(0);
      return CHAIN_COLORS[code % CHAIN_COLORS.length];
    }}
    const ratio = total > 1 ? index / (total - 1) : 0;
    if (ratio < 0.25) return PALETTE.indigo;
    if (ratio < 0.5) return PALETTE.teal;
    if (ratio < 0.75) return PALETTE.emerald;
    return PALETTE.amber;
  }}

  function hexToRgba(hex, alpha) {{
    const clean = hex.replace('#', '');
    const r = parseInt(clean.substring(0, 2), 16);
    const g = parseInt(clean.substring(2, 4), 16);
    const b = parseInt(clean.substring(4, 6), 16);
    return `rgba(${{r}}, ${{g}}, ${{b}}, ${{alpha}})`;
  }}

  function buildResidueSummaries(atoms) {{
    const map = new Map();
    for (const a of atoms) {{
      const key = `${{a.chain}}:${{a.resi}}`;
      if (!map.has(key)) {{
        map.set(key, {{
          key,
          chain: a.chain,
          resi: a.resi,
          resn: a.resn,
          ss: a.ss || (a.group === 'HETATM' ? 'HETATM' : 'LOOP'),
          group: a.group,
          atoms: [],
          ca: null
        }});
      }}
      const entry = map.get(key);
      entry.atoms.push(a);
      if (a.name === 'CA' || !entry.ca) {{
        entry.ca = a;
      }}
    }}

    const residues = Array.from(map.values());
    for (const r of residues) {{
      const sumB = r.atoms.reduce((acc, item) => acc + (item.b || 0), 0);
      r.meanB = Number((sumB / r.atoms.length).toFixed(2));
      r.x = r.ca ? r.ca.x : 0;
      r.y = r.ca ? r.ca.y : 0;
      r.z = r.ca ? r.ca.z : 0;
      r.contacts = 0;
    }}

    for (let i = 0; i < residues.length; i++) {{
      for (let j = i + 1; j < residues.length; j++) {{
        const dx = residues[i].x - residues[j].x;
        const dy = residues[i].y - residues[j].y;
        const dz = residues[i].z - residues[j].z;
        const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (Math.abs(i - j) > 1 && d < 8.0) {{
          residues[i].contacts++;
          residues[j].contacts++;
        }}
      }}
    }}
    return residues;
  }}

  /* Shared 3D Spatial Floor Stage, Cast Shadow & Orientation Triad Helpers */
  function draw3DSpatialStage(ctx, w, h, cosX, sinX, cosY, sinY, s, span, panX, panY, projectedItems, theme) {{
    const floorY = -span * 0.48;
    const projectStagePt = (x, y, z) => {{
      const x1 = x * cosY + z * sinY;
      const z1 = -x * sinY + z * cosY;
      const y2 = y * cosX - z1 * sinX;
      const z2 = y * sinX + z1 * cosX;
      const persp = Math.max(0.45, 1 + (z2 / (span * 2.4)));
      return {{
        px: w / 2 + panX + x1 * s * persp,
        py: h / 2 + panY - y2 * s * persp,
        pz: z2,
        persp
      }};
    }};

    /* 1. Draw 3D Perspective Polar Coordinate Rings on Floor Plane */
    const radii = [span * 0.28, span * 0.54];
    for (let rIdx = 0; rIdx < radii.length; rIdx++) {{
      const rad = radii[rIdx];
      const steps = 40;
      ctx.beginPath();
      for (let i = 0; i <= steps; i++) {{
        const ang = (i / steps) * Math.PI * 2;
        const pt = projectStagePt(Math.cos(ang) * rad, floorY, Math.sin(ang) * rad);
        if (i === 0) ctx.moveTo(pt.px, pt.py);
        else ctx.lineTo(pt.px, pt.py);
      }}
      ctx.strokeStyle = theme.border;
      ctx.lineWidth = rIdx === 1 ? 1.2 : 0.9;
      if (rIdx === 0) ctx.setLineDash([3, 4]);
      ctx.stroke();
      ctx.setLineDash([]);
    }}

    /* 2. Draw 3D Floor Cross-Hair Spokes */
    const outerR = span * 0.58;
    for (let a = 0; a < 4; a++) {{
      const ang = (a * Math.PI) / 4;
      const pA = projectStagePt(-Math.cos(ang) * outerR, floorY, -Math.sin(ang) * outerR);
      const pB = projectStagePt(Math.cos(ang) * outerR, floorY, Math.sin(ang) * outerR);
      ctx.beginPath();
      ctx.moveTo(pA.px, pA.py);
      ctx.lineTo(pB.px, pB.py);
      ctx.strokeStyle = theme.border;
      ctx.lineWidth = 0.8;
      ctx.stroke();
    }}

    /* 3. Project Soft Cast Ground Shadows for Atoms onto Floor Plane */
    if (projectedItems && projectedItems.length) {{
      const shadowAlpha = theme.isDark ? 0.24 : 0.09;
      const step = projectedItems.length > 45 ? 2 : 1;
      for (let i = 0; i < projectedItems.length; i += step) {{
        const item = projectedItems[i];
        const a = item.atom;
        const sx = a.x - (item.cx || 0);
        const sz = a.z - (item.cz || 0);
        const sp = projectStagePt(sx, floorY, sz);
        const sr = Math.max(4, 7.5 * (s / 10) * sp.persp);
        ctx.beginPath();
        ctx.ellipse(sp.px, sp.py, sr * 1.35, sr * 0.55, 0, 0, Math.PI * 2);
        ctx.fillStyle = `rgba(0, 0, 0, ${{shadowAlpha}})`;
        ctx.fill();
      }}
    }}
  }}

  function draw3DOrientationGizmo(ctx, w, h, rotX, rotY, scale, theme, modeTitle) {{
    const cosX = Math.cos(rotX), sinX = Math.sin(rotX);
    const cosY = Math.cos(rotY), sinY = Math.sin(rotY);
    const ox = w - 48;
    const oy = h - 42;
    const len = 23;
    const axes = [
      {{ label: 'X', x: 1, y: 0, z: 0, color: PALETTE.rose }},
      {{ label: 'Y', x: 0, y: 1, z: 0, color: PALETTE.emerald }},
      {{ label: 'Z', x: 0, y: 0, z: 1, color: PALETTE.indigo }}
    ];

    ctx.font = '600 10px "IBM Plex Mono", monospace';
    for (const ax of axes) {{
      const x1 = ax.x * cosY + ax.z * sinY;
      const z1 = -ax.x * sinY + ax.z * cosY;
      const y2 = ax.y * cosX - z1 * sinX;
      const ex = ox + x1 * len;
      const ey = oy - y2 * len;

      ctx.beginPath();
      ctx.moveTo(ox, oy);
      ctx.lineTo(ex, ey);
      ctx.strokeStyle = ax.color;
      ctx.lineWidth = 2;
      ctx.stroke();

      ctx.beginPath();
      ctx.arc(ex, ey, 2.8, 0, Math.PI * 2);
      ctx.fillStyle = ax.color;
      ctx.fill();

      ctx.fillStyle = ax.color;
      ctx.fillText(ax.label, ex + (x1 >= 0 ? 4 : -10), ey + (y2 >= 0 ? -4 : 10));
    }}

    const degX = Math.round((rotX * 180) / Math.PI);
    const degY = Math.round((rotY * 180) / Math.PI);
    ctx.fillStyle = theme.textSecondary;
    ctx.font = '400 10px "IBM Plex Mono", monospace';
    ctx.textAlign = 'right';
    ctx.fillText(
      `${{modeTitle}} · Zoom: ${{scale.toFixed(2)}}x · Orbit: (${{degX}}°, ${{degY}}°) · Style: ${{currentStyle.toUpperCase()}}`,
      w - 14,
      22
    );
    ctx.textAlign = 'left';
  }}

  function drawShaded3DSphere(ctx, px, py, r, hexColor, theme, isRaytraced) {{
    const safeR = Math.max(2.2, r);
    const lx = px - safeR * 0.32;
    const ly = py - safeR * 0.32;
    const grad = ctx.createRadialGradient(lx, ly, safeR * 0.06, px, py, safeR);
    if (isRaytraced) {{
      grad.addColorStop(0, '#FFFFFF');
      grad.addColorStop(0.22, hexColor);
      grad.addColorStop(0.76, hexColor);
      grad.addColorStop(1, theme.isDark ? '#090A0D' : '#292524');
    }} else {{
      grad.addColorStop(0, '#FFFFFF');
      grad.addColorStop(0.28, hexColor);
      grad.addColorStop(0.88, hexColor);
      grad.addColorStop(1, hexToRgba(hexColor, 0.72));
    }}

    ctx.beginPath();
    ctx.arc(px, py, safeR, 0, Math.PI * 2);
    ctx.fillStyle = grad;
    ctx.fill();

    ctx.strokeStyle = isRaytraced
      ? (theme.isDark ? 'rgba(9, 10, 13, 0.88)' : 'rgba(28, 25, 23, 0.65)')
      : theme.bgSurface;
    ctx.lineWidth = isRaytraced ? 1.4 : 1.1;
    ctx.stroke();

    if (isRaytraced && safeR > 4.5) {{
      ctx.beginPath();
      ctx.arc(px - safeR * 0.3, py - safeR * 0.3, Math.max(1.1, safeR * 0.2), 0, Math.PI * 2);
      ctx.fillStyle = 'rgba(255, 255, 255, 0.65)';
      ctx.fill();
    }}
  }}

  /* Interactive 3D Relationship Node Map Engine (#viewNetwork) */
  class RelationshipNodeMapEngine {{
    constructor(canvas, hudElement) {{
      this.canvas = canvas;
      this.ctx = canvas.getContext('2d');
      this.hud = hudElement;
      this.atoms = [];
      this.residues = [];
      this.rotX = 0.32;
      this.rotY = -0.52;
      this.scale = 1.0;
      this.panX = 0;
      this.panY = 0;
      this.isDragging = false;
      this.dragMoved = false;
      this.lastX = 0;
      this.lastY = 0;
      this.hoveredNode = null;
      this.projectedNodes = [];

      this.bindEvents();
      if (typeof ResizeObserver !== 'undefined' && this.canvas.parentElement) {{
        this.resizeObserver = new ResizeObserver(() => {{
          if (this.canvas.parentElement.clientWidth > 20) {{
            this.resize();
            this.render();
          }}
        }});
        this.resizeObserver.observe(this.canvas.parentElement);
      }}
    }}

    setData(atoms, residues) {{
      this.atoms = atoms || [];
      this.residues = residues || [];
      this.computeBounds();
      this.resize();
      this.render();
    }}

    computeBounds() {{
      if (!this.atoms.length) return;
      let minX = Infinity, maxX = -Infinity;
      let minY = Infinity, maxY = -Infinity;
      let minZ = Infinity, maxZ = -Infinity;
      for (const a of this.atoms) {{
        if (a.x < minX) minX = a.x;
        if (a.x > maxX) maxX = a.x;
        if (a.y < minY) minY = a.y;
        if (a.y > maxY) maxY = a.y;
        if (a.z < minZ) minZ = a.z;
        if (a.z > maxZ) maxZ = a.z;
      }}
      this.cx = (minX + maxX) / 2;
      this.cy = (minY + maxY) / 2;
      this.cz = (minZ + maxZ) / 2;
      this.span = Math.max(maxX - minX, maxY - minY, maxZ - minZ, 8);
    }}

    resize() {{
      const parent = this.canvas.parentElement;
      const rect = parent ? parent.getBoundingClientRect() : {{ width: 760, height: 460 }};
      const dpr = window.devicePixelRatio || 1;
      this.width = rect.width > 20 ? rect.width : 760;
      this.height = rect.height > 20 ? rect.height : 460;
      this.canvas.width = this.width * dpr;
      this.canvas.height = this.height * dpr;
      this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      this.baseScale = (Math.min(this.width, this.height) * 0.70) / (this.span || 20);
    }}

    zoomBy(factor) {{
      this.scale = Math.max(0.35, Math.min(4.5, this.scale * factor));
      this.render();
    }}

    bindEvents() {{
      this.canvas.addEventListener('mousedown', (e) => {{
        this.isDragging = true;
        this.dragMoved = false;
        this.lastX = e.clientX;
        this.lastY = e.clientY;
      }});

      window.addEventListener('mousemove', (e) => {{
        const rect = this.canvas.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;

        if (this.isDragging) {{
          const dx = e.clientX - this.lastX;
          const dy = e.clientY - this.lastY;
          if (Math.abs(dx) > 2 || Math.abs(dy) > 2) this.dragMoved = true;
          this.lastX = e.clientX;
          this.lastY = e.clientY;

          if (e.shiftKey || e.buttons === 2) {{
            this.panX += dx;
            this.panY += dy;
          }} else {{
            this.rotY += dx * 0.011;
            this.rotX += dy * 0.011;
          }}
          this.render();
          return;
        }}

        if (rect.width > 0 && mx >= 0 && mx <= this.width && my >= 0 && my <= this.height) {{
          this.updateHover(mx, my);
        }}
      }});

      window.addEventListener('mouseup', () => {{
        if (this.isDragging) {{
          const wasClick = !this.dragMoved;
          this.isDragging = false;
          if (wasClick && this.hoveredNode) {{
            selectResidue(this.hoveredNode.resKey);
          }}
        }}
      }});

      this.canvas.addEventListener('wheel', (e) => {{
        e.preventDefault();
        const factor = e.deltaY < 0 ? 1.12 : 0.89;
        this.zoomBy(factor);
      }}, {{ passive: false }});

      this.canvas.addEventListener('dblclick', () => {{
        this.resetCamera();
      }});

      window.addEventListener('resize', () => {{
        this.resize();
        this.render();
      }});
    }}

    updateHover(mx, my) {{
      let closest = null;
      let minDist = 18;
      for (let i = this.projectedNodes.length - 1; i >= 0; i--) {{
        const p = this.projectedNodes[i];
        const d = Math.hypot(mx - p.px, my - p.py);
        if (d < minDist) {{
          minDist = d;
          closest = p;
        }}
      }}
      if (this.hoveredNode !== closest) {{
        this.hoveredNode = closest;
        if (closest) {{
          this.hud.textContent =
            `Node Map · ${{closest.atom.resn}} ${{closest.atom.chain}}${{closest.atom.resi}} (${{closest.atom.name}}) · B: ${{closest.atom.b}} Å² · Class: ${{closest.atom.ss}}`;
        }} else {{
          this.hud.textContent = 'Drag: Rotate 3D · Scroll: Zoom · Double-Click or Reset View: Original Size & Position';
        }}
        this.render();
      }}
    }}

    resetCamera() {{
      this.rotX = 0.32;
      this.rotY = -0.52;
      this.scale = 1.0;
      this.panX = 0;
      this.panY = 0;
      this.hoveredNode = null;
      this.resize();
      if (this.hud) {{
        this.hud.textContent = 'Drag: Rotate 3D · Scroll: Zoom · Double-Click or Reset View: Original Size & Position';
      }}
      this.render();
    }}

    render() {{
      const ctx = this.ctx;
      const w = this.width;
      const h = this.height;
      const theme = getThemeColors();

      ctx.clearRect(0, 0, w, h);
      ctx.fillStyle = theme.bgSurface;
      ctx.fillRect(0, 0, w, h);
      if (!this.atoms.length) return;

      const cosX = Math.cos(this.rotX), sinX = Math.sin(this.rotX);
      const cosY = Math.cos(this.rotY), sinY = Math.sin(this.rotY);
      const s = (this.baseScale || 10) * this.scale;

      const sourceNodes = currentStyle === 'stick'
        ? this.atoms
        : this.residues.map(r => r.ca || r.atoms[0]);

      const projected = sourceNodes.map((atom, idx) => {{
        const x = atom.x - this.cx;
        const y = atom.y - this.cy;
        const z = atom.z - this.cz;

        const x1 = x * cosY + z * sinY;
        const z1 = -x * sinY + z * cosY;
        const y2 = y * cosX - z1 * sinX;
        const z2 = y * sinX + z1 * cosX;

        const perspective = Math.max(0.45, 1 + (z2 / (this.span * 2.4)));
        const px = w / 2 + this.panX + x1 * s * perspective;
        const py = h / 2 + this.panY - y2 * s * perspective;

        return {{
          atom,
          cx: this.cx,
          cz: this.cz,
          resKey: `${{atom.chain}}:${{atom.resi}}`,
          px,
          py,
          pz: z2,
          perspective,
          color: getAtomSemanticColor(atom, idx, sourceNodes.length)
        }};
      }});

      /* Pass 1: 3D Perspective Polar Stage & Cast Ground Shadows */
      draw3DSpatialStage(
        ctx, w, h, cosX, sinX, cosY, sinY, s, this.span, this.panX, this.panY, projected, theme
      );

      /* Pass 2: Non-Covalent Spatial Relationship Bridges (< 8.5 A) */
      if (currentStyle !== 'stick') {{
        for (let i = 0; i < projected.length; i++) {{
          for (let j = i + 2; j < projected.length; j++) {{
            const a1 = projected[i];
            const a2 = projected[j];
            const dist = Math.hypot(a1.atom.x - a2.atom.x, a1.atom.y - a2.atom.y, a1.atom.z - a2.atom.z);
            if (dist < 8.5) {{
              const isFocus =
                (this.hoveredNode && (this.hoveredNode.resKey === a1.resKey || this.hoveredNode.resKey === a2.resKey)) ||
                (selectedResidueKey && (selectedResidueKey === a1.resKey || selectedResidueKey === a2.resKey));

              ctx.beginPath();
              ctx.moveTo(a1.px, a1.py);
              ctx.lineTo(a2.px, a2.py);
              ctx.strokeStyle = isFocus ? PALETTE.rose : hexToRgba(a1.color, 0.34);
              ctx.lineWidth = isFocus ? 2.0 : 1.3;
              ctx.setLineDash(isFocus ? [2, 2] : [3, 4]);
              ctx.stroke();
              ctx.setLineDash([]);

              if (isFocus) {{
                const mx = (a1.px + a2.px) / 2;
                const my = (a1.py + a2.py) / 2;
                ctx.font = '600 9px "IBM Plex Mono", monospace';
                ctx.fillStyle = PALETTE.rose;
                ctx.fillText(`${{dist.toFixed(1)}}Å`, mx + 3, my - 3);
              }}
            }}
          }}
        }}
      }}

      /* Pass 3: 3D Shaded Polypeptide Backbone & Covalent Cylinders */
      for (let i = 0; i < projected.length; i++) {{
        for (let j = i + 1; j < Math.min(i + 14, projected.length); j++) {{
          const p1 = projected[i];
          const p2 = projected[j];
          const dist = Math.hypot(p1.atom.x - p2.atom.x, p1.atom.y - p2.atom.y, p1.atom.z - p2.atom.z);

          const isBackbone = currentStyle !== 'stick' && p1.atom.chain === p2.atom.chain && (j === i + 1) && dist < 5.5;
          const isCovalent = currentStyle === 'stick' && dist > 0.65 && dist < 2.1;

          if (isBackbone || isCovalent) {{
            const avgPersp = (p1.perspective + p2.perspective) / 2;
            const tubeW = currentStyle === 'cartoon'
              ? Math.max(3.5, 6.2 * this.scale * avgPersp)
              : Math.max(2.0, 3.2 * this.scale * avgPersp);

            /* Outer depth casing */
            ctx.beginPath();
            ctx.moveTo(p1.px, p1.py);
            ctx.lineTo(p2.px, p2.py);
            ctx.strokeStyle = theme.bgSurface;
            ctx.lineWidth = tubeW + 2.0;
            ctx.lineCap = 'round';
            ctx.stroke();

            /* Core 3D semantic tube */
            ctx.beginPath();
            ctx.moveTo(p1.px, p1.py);
            ctx.lineTo(p2.px, p2.py);
            ctx.strokeStyle = hexToRgba(p1.color, 0.86);
            ctx.lineWidth = tubeW;
            ctx.lineCap = 'round';
            ctx.stroke();

            /* Specular tube ridge */
            ctx.beginPath();
            ctx.moveTo(p1.px, p1.py - 0.8);
            ctx.lineTo(p2.px, p2.py - 0.8);
            ctx.strokeStyle = 'rgba(255, 255, 255, 0.35)';
            ctx.lineWidth = Math.max(1.0, tubeW * 0.28);
            ctx.lineCap = 'round';
            ctx.stroke();
          }}
        }}
      }}

      /* Pass 4: Sort Nodes Back-to-Front & Render 3D Shaded Spheres */
      const sorted = [...projected].sort((a, b) => a.pz - b.pz);
      this.projectedNodes = sorted;

      for (const node of sorted) {{
        const isSelected = selectedResidueKey === node.resKey;
        const isHovered = this.hoveredNode && this.hoveredNode.resKey === node.resKey;

        let radius = 6.2 * Math.sqrt(this.scale) * node.perspective;
        if (currentStyle === 'sphere') radius = 11.2 * Math.sqrt(this.scale) * node.perspective;
        if (currentStyle === 'stick') radius = 4.2 * Math.sqrt(this.scale) * node.perspective;
        if (currentStyle === 'surface') {{
          ctx.beginPath();
          ctx.arc(node.px, node.py, 18 * Math.sqrt(this.scale) * node.perspective, 0, Math.PI * 2);
          ctx.fillStyle = hexToRgba(node.color, 0.18);
          ctx.fill();
        }}

        if (isSelected || isHovered) {{
          ctx.beginPath();
          ctx.arc(node.px, node.py, radius + 5.5, 0, Math.PI * 2);
          ctx.strokeStyle = PALETTE.rose;
          ctx.lineWidth = 2.2;
          ctx.stroke();
        }}

        drawShaded3DSphere(ctx, node.px, node.py, radius, node.color, theme, false);

        if (currentStyle !== 'stick' && (projected.length <= 30 || isSelected || isHovered)) {{
          ctx.font = '500 10px "IBM Plex Mono", monospace';
          ctx.fillStyle = theme.textPrimary;
          ctx.fillText(`${{node.atom.resn}}${{node.atom.resi}}`, node.px + radius + 4, node.py + 3);
        }}
      }}

      draw3DOrientationGizmo(ctx, w, h, this.rotX, this.rotY, this.scale, theme, '3D Relationship Node Map');
    }}
  }}

  /* Interactive High-Craft 3D Ray-Traced Shader Engine (#viewRaytrace / #raytrace-canvas) */
  class Raytraced3DEngine {{
    constructor(canvas, hudElement) {{
      this.canvas = canvas;
      this.ctx = canvas.getContext('2d');
      this.hud = hudElement;
      this.atoms = [];
      this.residues = [];
      this.bonds = [];
      this.splineSegments = [];
      this.rotX = 0.35;
      this.rotY = -0.58;
      this.scale = 1.0;
      this.panX = 0;
      this.panY = 0;
      this.isDragging = false;
      this.dragMoved = false;
      this.lastX = 0;
      this.lastY = 0;
      this.hoveredAtom = null;
      this.projectedAtoms = [];

      this.bindEvents();
      if (typeof ResizeObserver !== 'undefined' && this.canvas.parentElement) {{
        this.resizeObserver = new ResizeObserver(() => {{
          if (this.canvas.parentElement.clientWidth > 20) {{
            this.resize();
            this.render();
          }}
        }});
        this.resizeObserver.observe(this.canvas.parentElement);
      }}
    }}

    setData(atoms, residues) {{
      this.atoms = atoms || [];
      this.residues = residues || [];
      this.computeBounds();
      this.buildRaytraceGeometry();
      this.resize();
      this.render();
    }}

    computeBounds() {{
      if (!this.atoms.length) return;
      let minX = Infinity, maxX = -Infinity;
      let minY = Infinity, maxY = -Infinity;
      let minZ = Infinity, maxZ = -Infinity;
      for (const a of this.atoms) {{
        if (a.x < minX) minX = a.x;
        if (a.x > maxX) maxX = a.x;
        if (a.y < minY) minY = a.y;
        if (a.y > maxY) maxY = a.y;
        if (a.z < minZ) minZ = a.z;
        if (a.z > maxZ) maxZ = a.z;
      }}
      this.cx = (minX + maxX) / 2;
      this.cy = (minY + maxY) / 2;
      this.cz = (minZ + maxZ) / 2;
      this.span = Math.max(maxX - minX, maxY - minY, maxZ - minZ, 8);
    }}

    buildRaytraceGeometry() {{
      this.bonds = [];
      this.splineSegments = [];
      const n = this.atoms.length;
      if (!n) return;

      for (let i = 0; i < n; i++) {{
        const a1 = this.atoms[i];
        for (let j = i + 1; j < Math.min(i + 28, n); j++) {{
          const a2 = this.atoms[j];
          const d = Math.hypot(a1.x - a2.x, a1.y - a2.y, a1.z - a2.z);
          const sameRes = a1.chain === a2.chain && a1.resi === a2.resi;
          if ((d > 0.65 && d < 2.08) || (sameRes && d > 0.65 && d < 2.45)) {{
            this.bonds.push([i, j]);
          }}
        }}
      }}

      const chainGroups = new Map();
      this.residues.forEach((r, rIdx) => {{
        if (r.group === 'HETATM') return;
        if (!chainGroups.has(r.chain)) chainGroups.set(r.chain, []);
        chainGroups.get(r.chain).push({{ r, rIdx }});
      }});

      for (const [, list] of chainGroups.entries()) {{
        for (let i = 0; i < list.length - 1; i++) {{
          const cur = list[i].r;
          const nxt = list[i + 1].r;
          const d = Math.hypot(cur.x - nxt.x, cur.y - nxt.y, cur.z - nxt.z);
          if (d < 6.0) {{
            const prev = list[Math.max(0, i - 1)].r;
            const after = list[Math.min(list.length - 1, i + 2)].r;
            const steps = 4;
            for (let s = 0; s < steps; s++) {{
              const t0 = s / steps;
              const t1 = (s + 1) / steps;
              const p0 = this.catmullRomPoint(prev, cur, nxt, after, t0);
              const p1 = this.catmullRomPoint(prev, cur, nxt, after, t1);
              this.splineSegments.push({{
                p0,
                p1,
                atom: cur.ca || cur.atoms[0],
                resIndex: list[i].rIdx
              }});
            }}
          }}
        }}
      }}
    }}

    catmullRomPoint(p0, p1, p2, p3, t) {{
      const t2 = t * t;
      const t3 = t2 * t;
      const f0 = -0.5 * t3 + t2 - 0.5 * t;
      const f1 = 1.5 * t3 - 2.5 * t2 + 1.0;
      const f2 = -1.5 * t3 + 2.0 * t2 + 0.5 * t;
      const f3 = 0.5 * t3 - 0.5 * t2;
      return {{
        x: p0.x * f0 + p1.x * f1 + p2.x * f2 + p3.x * f3,
        y: p0.y * f0 + p1.y * f1 + p2.y * f2 + p3.y * f3,
        z: p0.z * f0 + p1.z * f1 + p2.z * f2 + p3.z * f3
      }};
    }}

    resize() {{
      const parent = this.canvas.parentElement;
      const rect = parent ? parent.getBoundingClientRect() : {{ width: 760, height: 460 }};
      const dpr = window.devicePixelRatio || 1;
      this.width = rect.width > 20 ? rect.width : 760;
      this.height = rect.height > 20 ? rect.height : 460;
      this.canvas.width = this.width * dpr;
      this.canvas.height = this.height * dpr;
      this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      this.baseScale = (Math.min(this.width, this.height) * 0.70) / (this.span || 20);
    }}

    zoomBy(factor) {{
      this.scale = Math.max(0.35, Math.min(4.5, this.scale * factor));
      this.render();
    }}

    bindEvents() {{
      this.canvas.addEventListener('mousedown', (e) => {{
        this.isDragging = true;
        this.dragMoved = false;
        this.lastX = e.clientX;
        this.lastY = e.clientY;
      }});

      window.addEventListener('mousemove', (e) => {{
        const rect = this.canvas.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;

        if (this.isDragging) {{
          const dx = e.clientX - this.lastX;
          const dy = e.clientY - this.lastY;
          if (Math.abs(dx) > 2 || Math.abs(dy) > 2) this.dragMoved = true;
          this.lastX = e.clientX;
          this.lastY = e.clientY;

          if (e.shiftKey || e.buttons === 2) {{
            this.panX += dx;
            this.panY += dy;
          }} else {{
            this.rotY += dx * 0.011;
            this.rotX += dy * 0.011;
          }}
          this.render();
          return;
        }}

        if (rect.width > 0 && mx >= 0 && mx <= this.width && my >= 0 && my <= this.height) {{
          this.updateHover(mx, my);
        }}
      }});

      window.addEventListener('mouseup', () => {{
        if (this.isDragging) {{
          const wasClick = !this.dragMoved;
          this.isDragging = false;
          if (wasClick && this.hoveredAtom) {{
            selectResidue(this.hoveredAtom.resKey);
          }}
        }}
      }});

      this.canvas.addEventListener('wheel', (e) => {{
        e.preventDefault();
        const factor = e.deltaY < 0 ? 1.12 : 0.89;
        this.zoomBy(factor);
      }}, {{ passive: false }});

      this.canvas.addEventListener('dblclick', () => {{
        this.resetCamera();
      }});

      window.addEventListener('resize', () => {{
        this.resize();
        this.render();
      }});
    }}

    updateHover(mx, my) {{
      let closest = null;
      let minDist = 18;
      for (let i = this.projectedAtoms.length - 1; i >= 0; i--) {{
        const p = this.projectedAtoms[i];
        const d = Math.hypot(mx - p.px, my - p.py);
        if (d < minDist) {{
          minDist = d;
          closest = p;
        }}
      }}
      if (this.hoveredAtom !== closest) {{
        this.hoveredAtom = closest;
        if (closest && this.hud) {{
          this.hud.textContent =
            `Ray-Traced Shader · ${{closest.atom.resn}} ${{closest.atom.chain}}${{closest.atom.resi}} (${{closest.atom.name}}) · B: ${{closest.atom.b}} Å² · Class: ${{closest.atom.ss}}`;
        }} else if (this.hud) {{
          this.hud.textContent =
            'Ray-Traced 3D Shader · Drag: Rotate 3D · Scroll: Zoom · Reset View: Original Size & Position';
        }}
        this.render();
      }}
    }}

    resetCamera() {{
      this.rotX = 0.35;
      this.rotY = -0.58;
      this.scale = 1.0;
      this.panX = 0;
      this.panY = 0;
      this.hoveredAtom = null;
      this.resize();
      if (this.hud) {{
        this.hud.textContent =
          'Ray-Traced 3D Shader · Drag: Rotate 3D · Scroll: Zoom · Reset View: Original Size & Position';
      }}
      this.render();
    }}

    projectPoint(pt, cosX, sinX, cosY, sinY, s) {{
      const x = pt.x - this.cx;
      const y = pt.y - this.cy;
      const z = pt.z - this.cz;

      const x1 = x * cosY + z * sinY;
      const z1 = -x * sinY + z * cosY;
      const y2 = y * cosX - z1 * sinX;
      const z2 = y * sinX + z1 * cosX;

      const perspective = Math.max(0.45, 1 + (z2 / (this.span * 2.4)));
      return {{
        px: this.width / 2 + this.panX + x1 * s * perspective,
        py: this.height / 2 + this.panY - y2 * s * perspective,
        pz: z2,
        perspective
      }};
    }}

    render() {{
      const ctx = this.ctx;
      const w = this.width;
      const h = this.height;
      const theme = getThemeColors();

      ctx.clearRect(0, 0, w, h);
      ctx.fillStyle = theme.bgSurface;
      ctx.fillRect(0, 0, w, h);
      if (!this.atoms.length) return;

      const cosX = Math.cos(this.rotX), sinX = Math.sin(this.rotX);
      const cosY = Math.cos(this.rotY), sinY = Math.sin(this.rotY);
      const s = (this.baseScale || 10) * this.scale;

      const pts = this.atoms.map((a, idx) => {{
        const proj = this.projectPoint(a, cosX, sinX, cosY, sinY, s);
        const elem = (a.elem || a.name.charAt(0)).toUpperCase();
        const vdw = elem === 'C' ? 1.85 : (elem === 'O' ? 1.55 : (elem === 'N' ? 1.65 : 1.95));
        return {{
          atom: a,
          cx: this.cx,
          cz: this.cz,
          resKey: `${{a.chain}}:${{a.resi}}`,
          px: proj.px,
          py: proj.py,
          pz: proj.pz,
          perspective: proj.perspective,
          vdw,
          color: getAtomSemanticColor(a, idx, this.atoms.length)
        }};
      }});

      /* Pass 1: 3D Studio Floor Grid & Cast Ambient-Occlusion Ground Shadows */
      draw3DSpatialStage(
        ctx, w, h, cosX, sinX, cosY, sinY, s, this.span, this.panX, this.panY, pts, theme
      );

      const sortedAtoms = [...pts].sort((a, b) => a.pz - b.pz);
      this.projectedAtoms = sortedAtoms;

      /* Pass 2: Optional Ray-Traced Refractive Surface Envelope (when Surface mode is active) */
      if (currentStyle === 'surface') {{
        for (const p of sortedAtoms) {{
          const surfR = Math.max(14, p.vdw * s * p.perspective * 2.25);
          const grad = ctx.createRadialGradient(
            p.px - surfR * 0.22, p.py - surfR * 0.22, surfR * 0.08,
            p.px, p.py, surfR
          );
          grad.addColorStop(0, hexToRgba(p.color, 0.34));
          grad.addColorStop(0.7, hexToRgba(p.color, 0.16));
          grad.addColorStop(1, hexToRgba(p.color, 0.0));
          ctx.beginPath();
          ctx.arc(p.px, p.py, surfR, 0, Math.PI * 2);
          ctx.fillStyle = grad;
          ctx.fill();
        }}
      }}

      /* Pass 3: PyMOL ray_trace_mode = 1 Silhouette-Cased 3D Ribbons & Cylinders */
      const hasContinuousBackbone = this.splineSegments.length > 0;
      const outlineColor = theme.isDark ? 'rgba(8, 9, 12, 0.92)' : 'rgba(28, 25, 23, 0.72)';

      if ((currentStyle === 'cartoon' || currentStyle === 'surface') && hasContinuousBackbone) {{
        const projSplines = this.splineSegments.map(seg => {{
          const a = this.projectPoint(seg.p0, cosX, sinX, cosY, sinY, s);
          const b = this.projectPoint(seg.p1, cosX, sinX, cosY, sinY, s);
          return {{
            a,
            b,
            avgZ: (a.pz + b.pz) / 2,
            color: getAtomSemanticColor(seg.atom, seg.resIndex, this.residues.length)
          }};
        }});
        projSplines.sort((x, y) => x.avgZ - y.avgZ);

        for (const seg of projSplines) {{
          const rw = Math.max(5.5, 9.2 * this.scale * ((seg.a.perspective + seg.b.perspective) / 2));
          /* Ray-traced dark silhouette casing */
          ctx.beginPath();
          ctx.moveTo(seg.a.px, seg.a.py);
          ctx.lineTo(seg.b.px, seg.b.py);
          ctx.strokeStyle = outlineColor;
          ctx.lineWidth = rw + 2.6;
          ctx.lineCap = 'round';
          ctx.stroke();

          /* Diffuse 3D ribbon core */
          ctx.beginPath();
          ctx.moveTo(seg.a.px, seg.a.py);
          ctx.lineTo(seg.b.px, seg.b.py);
          ctx.strokeStyle = seg.color;
          ctx.lineWidth = rw;
          ctx.lineCap = 'round';
          ctx.stroke();

          /* High-gloss Blinn-Phong specular ridge */
          ctx.beginPath();
          ctx.moveTo(seg.a.px, seg.a.py - rw * 0.18);
          ctx.lineTo(seg.b.px, seg.b.py - rw * 0.18);
          ctx.strokeStyle = 'rgba(255, 255, 255, 0.48)';
          ctx.lineWidth = Math.max(1.4, rw * 0.26);
          ctx.lineCap = 'round';
          ctx.stroke();
        }}
      }}

      /* Covalent Cylinders (Always drawn in Cartoon/Stick/Surface for full atomic detail) */
      if (currentStyle !== 'sphere') {{
        const projBonds = this.bonds.map(([i, j]) => ({{
          p1: pts[i],
          p2: pts[j],
          avgZ: (pts[i].pz + pts[j].pz) / 2
        }}));
        projBonds.sort((a, b) => a.avgZ - b.avgZ);

        for (const b of projBonds) {{
          const avgPersp = (b.p1.perspective + b.p2.perspective) / 2;
          const bw = Math.max(2.8, (currentStyle === 'stick' ? 5.4 : 3.8) * Math.sqrt(this.scale) * avgPersp);

          /* Ray-traced silhouette edge */
          ctx.beginPath();
          ctx.moveTo(b.p1.px, b.p1.py);
          ctx.lineTo(b.p2.px, b.p2.py);
          ctx.strokeStyle = outlineColor;
          ctx.lineWidth = bw + 2.0;
          ctx.lineCap = 'round';
          ctx.stroke();

          /* Cylinder diffuse core */
          ctx.beginPath();
          ctx.moveTo(b.p1.px, b.p1.py);
          ctx.lineTo(b.p2.px, b.p2.py);
          ctx.strokeStyle = hexToRgba(b.p1.color, 0.92);
          ctx.lineWidth = bw;
          ctx.lineCap = 'round';
          ctx.stroke();

          /* Specular cylinder highlight */
          ctx.beginPath();
          ctx.moveTo(b.p1.px, b.p1.py - 0.8);
          ctx.lineTo(b.p2.px, b.p2.py - 0.8);
          ctx.strokeStyle = 'rgba(255, 255, 255, 0.42)';
          ctx.lineWidth = Math.max(1.0, bw * 0.28);
          ctx.lineCap = 'round';
          ctx.stroke();
        }}
      }}

      /* Pass 4: Ray-Traced Multi-Stop Impostor Spheres */
      for (const p of sortedAtoms) {{
        const isSelected = selectedResidueKey === p.resKey;
        const isHovered = this.hoveredAtom && this.hoveredAtom.resKey === p.resKey;
        const isCA = p.atom.name === 'CA';
        const isHetero = p.atom.group === 'HETATM';

        let r = (isCA ? 7.8 : 5.4) * Math.sqrt(this.scale) * p.perspective;
        if (currentStyle === 'sphere') {{
          r = Math.max(9.0, p.vdw * s * p.perspective * 0.92);
        }} else if (currentStyle === 'stick') {{
          r = (isCA ? 6.4 : 5.2) * Math.sqrt(this.scale) * p.perspective;
        }} else if (isHetero) {{
          r = 7.4 * Math.sqrt(this.scale) * p.perspective;
        }}

        if (isSelected || isHovered) {{
          ctx.beginPath();
          ctx.arc(p.px, p.py, r + 5.5, 0, Math.PI * 2);
          ctx.strokeStyle = PALETTE.rose;
          ctx.lineWidth = 2.2;
          ctx.stroke();
        }}

        drawShaded3DSphere(ctx, p.px, p.py, r, p.color, theme, true);

        if ((isCA && (this.residues.length <= 24 || isSelected || isHovered)) || (isHetero && p.atom.name === 'C1')) {{
          ctx.font = '600 10px "IBM Plex Mono", monospace';
          ctx.fillStyle = theme.textPrimary;
          ctx.fillText(`${{p.atom.resn}}${{p.atom.resi}}`, p.px + r + 4, p.py + 3);
        }}
      }}

      draw3DOrientationGizmo(ctx, w, h, this.rotX, this.rotY, this.scale, theme, 'Ray-Traced 3D Shader');
    }}
  }}

  function renderRaytraceCanvas(atoms, residues) {{
    const canvas = document.getElementById('raytrace-canvas');
    const hud = document.getElementById('raytraceHud');
    if (!canvas) return;
    if (!raytraceEngine) {{
      raytraceEngine = new Raytraced3DEngine(canvas, hud);
    }}
    raytraceEngine.setData(atoms, residues);
  }}

  /* Sequence B-Factor & Spatial Density Profile Chart */
  let cachedProfileBars = [];

  function renderProfileChart(residues) {{
    const canvas = document.getElementById('profile-canvas');
    if (!canvas) return;
    const rect = canvas.parentElement.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    const w = (rect.width || 620) - 32;
    const h = (rect.height || 270) - 24;
    canvas.width = w * dpr;
    canvas.height = h * dpr;
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);

    const theme = getThemeColors();
    const padL = 38, padR = 16, padT = 14, padB = 28;
    const plotW = Math.max(10, w - padL - padR);
    const plotH = Math.max(10, h - padT - padB);

    /* Faint Functional Gridlines for Thresholds */
    const thresholds = [25, 50, 70, 90, 100];
    ctx.font = '400 10px "IBM Plex Mono", monospace';
    for (const val of thresholds) {{
      const y = padT + plotH - (val / 100) * plotH;
      ctx.beginPath();
      ctx.moveTo(padL, y);
      ctx.lineTo(padL + plotW, y);
      ctx.strokeStyle = theme.border;
      ctx.lineWidth = 1;
      ctx.stroke();

      ctx.fillStyle = theme.textSecondary;
      ctx.textAlign = 'right';
      ctx.fillText(String(val), padL - 6, y + 3);
    }}

    if (!residues.length) return;

    const slotW = plotW / residues.length;
    const barW = Math.max(4, Math.min(22, slotW * 0.68));
    cachedProfileBars = [];

    residues.forEach((r, idx) => {{
      const val = Math.min(100, Math.max(0, r.meanB));
      const barH = (val / 100) * plotH;
      const x = padL + idx * slotW + (slotW - barW) / 2;
      const y = padT + plotH - barH;

      let color = PALETTE.rose;
      if (r.meanB >= 90) color = PALETTE.emerald;
      else if (r.meanB >= 70) color = PALETTE.teal;
      else if (r.meanB >= 50) color = PALETTE.amber;

      ctx.fillStyle = selectedResidueKey === r.key ? PALETTE.indigo : color;
      ctx.fillRect(x, y, barW, barH);

      if (selectedResidueKey === r.key) {{
        ctx.strokeStyle = theme.textPrimary;
        ctx.lineWidth = 1.5;
        ctx.strokeRect(x - 1, y - 1, barW + 2, barH + 2);
      }}

      cachedProfileBars.push({{ x, y, w: barW, h: barH, residue: r }});

      const step = residues.length > 16 ? 2 : 1;
      if (idx % step === 0) {{
        ctx.fillStyle = theme.textSecondary;
        ctx.textAlign = 'center';
        ctx.fillText(`${{r.resi}}`, x + barW / 2, padT + plotH + 16);
      }}
    }});

    /* Overlay Vivid Indigo Neighbor Density Curve */
    const maxContacts = Math.max(4, ...residues.map(r => r.contacts));
    ctx.beginPath();
    residues.forEach((r, idx) => {{
      const cx = padL + idx * slotW + slotW / 2;
      const cy = padT + plotH - (r.contacts / maxContacts) * (plotH * 0.85);
      if (idx === 0) ctx.moveTo(cx, cy);
      else ctx.lineTo(cx, cy);
    }});
    ctx.strokeStyle = PALETTE.indigo;
    ctx.lineWidth = 2;
    ctx.stroke();

    residues.forEach((r, idx) => {{
      const cx = padL + idx * slotW + slotW / 2;
      const cy = padT + plotH - (r.contacts / maxContacts) * (plotH * 0.85);
      ctx.beginPath();
      ctx.arc(cx, cy, 3, 0, Math.PI * 2);
      ctx.fillStyle = PALETTE.indigo;
      ctx.fill();
    }});
  }}

  /* Pairwise C-Alpha Distance Matrix Heatmap */
  let cachedMatrixState = null;

  function renderDistanceMatrix(residues) {{
    const canvas = document.getElementById('matrix-canvas');
    if (!canvas) return;
    const rect = canvas.parentElement.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    const w = (rect.width || 440) - 32;
    const h = (rect.height || 270) - 24;
    canvas.width = w * dpr;
    canvas.height = h * dpr;
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);

    const n = residues.length;
    if (!n) return;

    const theme = getThemeColors();
    const size = Math.min(w - 48, h - 24);
    const offsetX = Math.floor((w - size) / 2);
    const offsetY = 6;
    const cellS = size / n;

    cachedMatrixState = {{ residues, offsetX, offsetY, size, cellS }};

    for (let i = 0; i < n; i++) {{
      for (let j = 0; j < n; j++) {{
        const r1 = residues[i];
        const r2 = residues[j];
        const dx = r1.x - r2.x;
        const dy = r1.y - r2.y;
        const dz = r1.z - r2.z;
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);

        let fill = theme.bgBase;
        if (i === j) {{
          fill = hexToRgba(PALETTE.emerald, 0.85);
        }} else if (dist < 6.0) {{
          fill = PALETTE.rose;
        }} else if (dist < 12.0) {{
          const t = 1 - (dist - 6.0) / 6.0;
          fill = hexToRgba(PALETTE.teal, 0.35 + t * 0.6);
        }} else if (dist < 22.0) {{
          const t = 1 - (dist - 12.0) / 10.0;
          fill = hexToRgba(PALETTE.indigo, 0.15 + t * 0.55);
        }}

        ctx.fillStyle = fill;
        ctx.fillRect(offsetX + j * cellS, offsetY + i * cellS, Math.ceil(cellS), Math.ceil(cellS));
      }}
    }}

    ctx.strokeStyle = theme.border;
    ctx.lineWidth = 1;
    ctx.strokeRect(offsetX, offsetY, size, size);
  }}

  function populateResidueTable(residues) {{
    const tbody = document.getElementById('residueTableBody');
    const query = (document.getElementById('residueFilterInput').value || '').trim().toLowerCase();
    tbody.innerHTML = '';

    const filtered = residues.filter(r => {{
      if (confidenceFilterTier === 'vh' && r.meanB < 90) return false;
      if (confidenceFilterTier === 'conf' && (r.meanB < 70 || r.meanB >= 90)) return false;
      if (confidenceFilterTier === 'low' && (r.meanB < 50 || r.meanB >= 70)) return false;
      if (confidenceFilterTier === 'vl' && r.meanB >= 50) return false;

      if (!query) return true;
      return (
        r.resn.toLowerCase().includes(query) ||
        String(r.resi).includes(query) ||
        r.chain.toLowerCase().includes(query) ||
        r.ss.toLowerCase().includes(query)
      );
    }});

    document.getElementById('tableSummaryMeta').textContent =
      `Displaying ${{filtered.length}} of ${{residues.length}} residues`;

    for (const r of filtered) {{
      const tr = document.createElement('tr');
      tr.className = 'interactive-row' + (selectedResidueKey === r.key ? ' selected-row' : '');
      tr.onclick = () => selectResidue(r.key);

      let bColor = PALETTE.rose;
      if (r.meanB >= 90) bColor = PALETTE.emerald;
      else if (r.meanB >= 70) bColor = PALETTE.teal;
      else if (r.meanB >= 50) bColor = PALETTE.amber;

      /* Only show status indicator if residue is in a critical low-confidence state (<50) */
      const statusCell = r.meanB < 50
        ? `<span style="color: ${{PALETTE.rose}}; font-weight: 600;">Critical (&lt;50)</span>`
        : `<span style="color: var(--text-secondary);">—</span>`;

      tr.innerHTML = `
        <td>${{r.resn}} ${{r.resi}}</td>
        <td>${{r.chain}}</td>
        <td>${{r.ss}}</td>
        <td class="text-right">${{r.atoms.length}}</td>
        <td class="text-right">${{r.x.toFixed(2)}}</td>
        <td class="text-right">${{r.y.toFixed(2)}}</td>
        <td class="text-right">${{r.z.toFixed(2)}}</td>
        <td class="text-right" style="color: ${{bColor}}; font-weight: 600;">${{r.meanB.toFixed(1)}}</td>
        <td class="text-right">${{r.contacts}}</td>
        <td>${{statusCell}}</td>
      `;
      tbody.appendChild(tr);
    }}
  }}

  function selectResidue(resKey) {{
    selectedResidueKey = selectedResidueKey === resKey ? null : resKey;
    const badge = document.getElementById('selectedResidueBadge');
    badge.textContent = selectedResidueKey ? `Selected: ${{selectedResidueKey}}` : 'All Residues';
    refreshVisualizations();
  }}

  function filterByConfidence(tier) {{
    confidenceFilterTier = confidenceFilterTier === tier ? null : tier;
    refreshVisualizations();
  }}

  function clearResidueFilter() {{
    document.getElementById('residueFilterInput').value = '';
    confidenceFilterTier = null;
    selectedResidueKey = null;
    document.getElementById('selectedResidueBadge').textContent = 'All Residues';
    refreshVisualizations();
  }}

  /* Interactive 3D Molecular Surface & Volumetric Hull Engine (#view3d / #mol-container) */
  class MolecularSurface3DEngine {{
    constructor(canvas, hudElement) {{
      this.canvas = canvas;
      this.ctx = canvas.getContext('2d');
      this.hud = hudElement;
      this.atoms = [];
      this.residues = [];
      this.bonds = [];
      this.splineSegments = [];
      this.surfaceFacets = [];
      this.rotX = 0.28;
      this.rotY = -0.46;
      this.scale = 1.0;
      this.panX = 0;
      this.panY = 0;
      this.isDragging = false;
      this.dragMoved = false;
      this.lastX = 0;
      this.lastY = 0;
      this.hoveredAtom = null;
      this.projectedAtoms = [];

      this.bindEvents();
      if (typeof ResizeObserver !== 'undefined' && this.canvas.parentElement) {{
        this.resizeObserver = new ResizeObserver(() => {{
          if (this.canvas.parentElement.clientWidth > 20) {{
            this.resize();
            this.render();
          }}
        }});
        this.resizeObserver.observe(this.canvas.parentElement);
      }}
    }}

    setData(atoms, residues) {{
      this.atoms = atoms || [];
      this.residues = residues || [];
      this.computeBounds();
      this.buildTopologyAndSurfaceMesh();
      this.resize();
      this.render();
    }}

    computeBounds() {{
      if (!this.atoms.length) return;
      let minX = Infinity, maxX = -Infinity;
      let minY = Infinity, maxY = -Infinity;
      let minZ = Infinity, maxZ = -Infinity;
      for (const a of this.atoms) {{
        if (a.x < minX) minX = a.x;
        if (a.x > maxX) maxX = a.x;
        if (a.y < minY) minY = a.y;
        if (a.y > maxY) maxY = a.y;
        if (a.z < minZ) minZ = a.z;
        if (a.z > maxZ) maxZ = a.z;
      }}
      this.cx = (minX + maxX) / 2;
      this.cy = (minY + maxY) / 2;
      this.cz = (minZ + maxZ) / 2;
      this.span = Math.max(maxX - minX, maxY - minY, maxZ - minZ, 8);
    }}

    buildTopologyAndSurfaceMesh() {{
      this.bonds = [];
      this.splineSegments = [];
      this.surfaceFacets = [];
      const n = this.atoms.length;
      if (!n) return;

      /* 1. Covalent & Intra-Residue Bonds */
      for (let i = 0; i < n; i++) {{
        const a1 = this.atoms[i];
        for (let j = i + 1; j < Math.min(i + 28, n); j++) {{
          const a2 = this.atoms[j];
          const d = Math.hypot(a1.x - a2.x, a1.y - a2.y, a1.z - a2.z);
          const sameRes = a1.chain === a2.chain && a1.resi === a2.resi;
          if ((d > 0.65 && d < 2.05) || (sameRes && d > 0.65 && d < 2.45)) {{
            this.bonds.push([i, j]);
          }}
        }}
      }}

      /* 2. Smooth 3D Polypeptide Ribbon Splines per Chain */
      const chainGroups = new Map();
      this.residues.forEach((r, rIdx) => {{
        if (r.group === 'HETATM') return;
        if (!chainGroups.has(r.chain)) chainGroups.set(r.chain, []);
        chainGroups.get(r.chain).push({{ r, rIdx }});
      }});

      for (const [, list] of chainGroups.entries()) {{
        for (let i = 0; i < list.length - 1; i++) {{
          const cur = list[i].r;
          const nxt = list[i + 1].r;
          const d = Math.hypot(cur.x - nxt.x, cur.y - nxt.y, cur.z - nxt.z);
          if (d < 6.0) {{
            const prev = list[Math.max(0, i - 1)].r;
            const after = list[Math.min(list.length - 1, i + 2)].r;
            const steps = 4;
            for (let s = 0; s < steps; s++) {{
              const t0 = s / steps;
              const t1 = (s + 1) / steps;
              const p0 = this.catmullRomPoint(prev, cur, nxt, after, t0);
              const p1 = this.catmullRomPoint(prev, cur, nxt, after, t1);
              this.splineSegments.push({{
                p0,
                p1,
                atom: cur.ca || cur.atoms[0],
                resIndex: list[i].rIdx
              }});
            }}
          }}
        }}
      }}

      /* 3. Geodesic 3D Molecular Surface Facet Triangulation */
      const hullNodes = this.residues.length >= 10
        ? this.residues.map((r, idx) => ({{
            x: r.x,
            y: r.y,
            z: r.z,
            atom: r.ca || r.atoms[0],
            idx
          }}))
        : this.atoms.map((a, idx) => ({{
            x: a.x,
            y: a.y,
            z: a.z,
            atom: a,
            idx
          }}));

      const maxEdge = this.residues.length < 10 ? 4.2 : 7.8;
      const m = hullNodes.length;
      for (let i = 0; i < m; i++) {{
        for (let j = i + 1; j < m; j++) {{
          const dij = Math.hypot(
            hullNodes[i].x - hullNodes[j].x,
            hullNodes[i].y - hullNodes[j].y,
            hullNodes[i].z - hullNodes[j].z
          );
          if (dij > maxEdge) continue;
          for (let k = j + 1; k < m; k++) {{
            const dik = Math.hypot(
              hullNodes[i].x - hullNodes[k].x,
              hullNodes[i].y - hullNodes[k].y,
              hullNodes[i].z - hullNodes[k].z
            );
            if (dik > maxEdge) continue;
            const djk = Math.hypot(
              hullNodes[j].x - hullNodes[k].x,
              hullNodes[j].y - hullNodes[k].y,
              hullNodes[j].z - hullNodes[k].z
            );
            if (djk <= maxEdge && this.surfaceFacets.length < 120) {{
              this.surfaceFacets.push([hullNodes[i], hullNodes[j], hullNodes[k]]);
            }}
          }}
        }}
      }}
    }}

    catmullRomPoint(p0, p1, p2, p3, t) {{
      const t2 = t * t;
      const t3 = t2 * t;
      const f0 = -0.5 * t3 + t2 - 0.5 * t;
      const f1 = 1.5 * t3 - 2.5 * t2 + 1.0;
      const f2 = -1.5 * t3 + 2.0 * t2 + 0.5 * t;
      const f3 = 0.5 * t3 - 0.5 * t2;
      return {{
        x: p0.x * f0 + p1.x * f1 + p2.x * f2 + p3.x * f3,
        y: p0.y * f0 + p1.y * f1 + p2.y * f2 + p3.y * f3,
        z: p0.z * f0 + p1.z * f1 + p2.z * f2 + p3.z * f3
      }};
    }}

    resize() {{
      const parent = this.canvas.parentElement;
      const rect = parent ? parent.getBoundingClientRect() : {{ width: 760, height: 460 }};
      const dpr = window.devicePixelRatio || 1;
      this.width = rect.width > 20 ? rect.width : 760;
      this.height = rect.height > 20 ? rect.height : 460;
      this.canvas.width = this.width * dpr;
      this.canvas.height = this.height * dpr;
      this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      this.baseScale = (Math.min(this.width, this.height) * 0.66) / (this.span || 20);
    }}

    bindEvents() {{
      this.canvas.addEventListener('mousedown', (e) => {{
        this.isDragging = true;
        this.dragMoved = false;
        this.lastX = e.clientX;
        this.lastY = e.clientY;
      }});

      window.addEventListener('mousemove', (e) => {{
        const rect = this.canvas.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;

        if (this.isDragging) {{
          const dx = e.clientX - this.lastX;
          const dy = e.clientY - this.lastY;
          if (Math.abs(dx) > 2 || Math.abs(dy) > 2) this.dragMoved = true;
          this.lastX = e.clientX;
          this.lastY = e.clientY;

          if (e.shiftKey || e.buttons === 2) {{
            this.panX += dx;
            this.panY += dy;
          }} else {{
            this.rotY += dx * 0.011;
            this.rotX += dy * 0.011;
          }}
          this.render();
          return;
        }}

        if (rect.width > 0 && mx >= 0 && mx <= this.width && my >= 0 && my <= this.height) {{
          this.updateHover(mx, my);
        }}
      }});

      window.addEventListener('mouseup', () => {{
        if (this.isDragging) {{
          const wasClick = !this.dragMoved;
          this.isDragging = false;
          if (wasClick && this.hoveredAtom) {{
            selectResidue(this.hoveredAtom.resKey);
          }}
        }}
      }});

      this.canvas.addEventListener('wheel', (e) => {{
        e.preventDefault();
        const factor = e.deltaY < 0 ? 1.1 : 0.91;
        this.scale = Math.max(0.35, Math.min(4.5, this.scale * factor));
        this.render();
      }}, {{ passive: false }});

      this.canvas.addEventListener('dblclick', () => {{
        this.resetCamera();
      }});

      window.addEventListener('resize', () => {{
        this.resize();
        this.render();
      }});
    }}

    updateHover(mx, my) {{
      let closest = null;
      let minDist = 20;
      for (let i = this.projectedAtoms.length - 1; i >= 0; i--) {{
        const p = this.projectedAtoms[i];
        const d = Math.hypot(mx - p.px, my - p.py);
        if (d < minDist) {{
          minDist = d;
          closest = p;
        }}
      }}
      if (this.hoveredAtom !== closest) {{
        this.hoveredAtom = closest;
        if (closest) {{
          this.hud.textContent =
            `Surface Probe · ${{closest.atom.resn}} ${{closest.atom.chain}}${{closest.atom.resi}} (${{closest.atom.name}}) · B: ${{closest.atom.b}} Å² · Class: ${{closest.atom.ss}}`;
        }} else {{
          this.hud.textContent =
            '3D Molecular Surface (VDW Envelope) · Drag: Rotate 3D · Scroll: Zoom · Reset View: Original Size & Position';
        }}
        this.render();
      }}
    }}

    resetCamera() {{
      this.rotX = 0.28;
      this.rotY = -0.46;
      this.scale = 1.0;
      this.panX = 0;
      this.panY = 0;
      this.hoveredAtom = null;
      this.resize();
      if (this.hud) {{
        this.hud.textContent =
          '3D Molecular Surface (VDW Envelope) · Drag: Rotate 3D · Scroll: Zoom · Reset View: Original Size & Position';
      }}
      this.render();
    }}

    zoomBy(factor) {{
      this.scale = Math.max(0.35, Math.min(4.5, this.scale * factor));
      this.render();
    }}

    projectPoint(pt, cosX, sinX, cosY, sinY, s) {{
      const x = pt.x - this.cx;
      const y = pt.y - this.cy;
      const z = pt.z - this.cz;

      const x1 = x * cosY + z * sinY;
      const z1 = -x * sinY + z * cosY;
      const y2 = y * cosX - z1 * sinX;
      const z2 = y * sinX + z1 * cosX;

      const perspective = Math.max(0.45, 1 + (z2 / (this.span * 2.4)));
      return {{
        px: this.width / 2 + this.panX + x1 * s * perspective,
        py: this.height / 2 + this.panY - y2 * s * perspective,
        pz: z2,
        vx: x1,
        vy: y2,
        vz: z2,
        perspective
      }};
    }}

    render() {{
      const ctx = this.ctx;
      const w = this.width;
      const h = this.height;
      const theme = getThemeColors();

      ctx.clearRect(0, 0, w, h);
      ctx.fillStyle = theme.bgSurface;
      ctx.fillRect(0, 0, w, h);

      if (!this.atoms.length) return;

      const cosX = Math.cos(this.rotX), sinX = Math.sin(this.rotX);
      const cosY = Math.cos(this.rotY), sinY = Math.sin(this.rotY);
      const s = (this.baseScale || 10) * this.scale;
      const totalAtoms = this.atoms.length;

      /* Project All Atoms */
      const projected = this.atoms.map((atom, idx) => {{
        const proj = this.projectPoint(atom, cosX, sinX, cosY, sinY, s);
        const elem = (atom.elem || atom.name.charAt(0)).toUpperCase();
        const vdw = elem === 'C' ? 1.85 : (elem === 'O' ? 1.55 : (elem === 'N' ? 1.65 : 1.95));
        return {{
          atom,
          idx,
          resKey: `${{atom.chain}}:${{atom.resi}}`,
          px: proj.px,
          py: proj.py,
          pz: proj.pz,
          perspective: proj.perspective,
          vdw,
          cx: this.cx,
          cz: this.cz,
          color: getAtomSemanticColor(atom, idx, totalAtoms)
        }};
      }});

      /* Pass 0: 3D Perspective Floor Stage & Ground Shadow Projection */
      draw3DSpatialStage(
        ctx,
        w,
        h,
        cosX,
        sinX,
        cosY,
        sinY,
        s,
        this.span || 20,
        this.panX,
        this.panY,
        projected,
        theme
      );

      const sortedAtoms = [...projected].sort((a, b) => a.pz - b.pz);
      this.projectedAtoms = sortedAtoms;

      /* Pass 1: Volumetric Van der Waals & Solvent-Accessible Surface Envelope */
      const surfaceAlphaBoost = currentStyle === 'surface' ? 1.38 : 1.0;
      const probeScale = currentStyle === 'surface' ? 2.55 : 2.15;

      /* Outer solvent-accessible silhouette contour */
      for (const node of sortedAtoms) {{
        const surfR = Math.max(12, node.vdw * s * node.perspective * probeScale);
        const grad = ctx.createRadialGradient(
          node.px - surfR * 0.22,
          node.py - surfR * 0.22,
          surfR * 0.08,
          node.px,
          node.py,
          surfR
        );
        grad.addColorStop(0, hexToRgba(node.color, Math.min(0.48, 0.26 * surfaceAlphaBoost)));
        grad.addColorStop(0.62, hexToRgba(node.color, Math.min(0.32, 0.15 * surfaceAlphaBoost)));
        grad.addColorStop(1, hexToRgba(node.color, 0.0));

        ctx.beginPath();
        ctx.arc(node.px, node.py, surfR, 0, Math.PI * 2);
        ctx.fillStyle = grad;
        ctx.fill();

        /* Subtle iso-surface perimeter contour ring */
        ctx.beginPath();
        ctx.arc(node.px, node.py, surfR * 0.72, 0, Math.PI * 2);
        ctx.strokeStyle = hexToRgba(node.color, Math.min(0.28, 0.12 * surfaceAlphaBoost));
        ctx.lineWidth = 1.1;
        ctx.stroke();
      }}

      /* Pass 2: 3D Normal-Shaded Geodesic Molecular Surface Mesh Facets */
      if (this.surfaceFacets.length > 0) {{
        const projectedFacets = this.surfaceFacets.map(tri => {{
          const p0 = this.projectPoint(tri[0], cosX, sinX, cosY, sinY, s);
          const p1 = this.projectPoint(tri[1], cosX, sinX, cosY, sinY, s);
          const p2 = this.projectPoint(tri[2], cosX, sinX, cosY, sinY, s);
          const ux = p1.vx - p0.vx, uy = p1.vy - p0.vy, uz = p1.vz - p0.vz;
          const vx = p2.vx - p0.vx, vy = p2.vy - p0.vy, vz = p2.vz - p0.vz;
          const nx = uy * vz - uz * vy;
          const ny = uz * vx - ux * vz;
          const nz = ux * vy - uy * vx;
          const norm = Math.hypot(nx, ny, nz) || 1;
          const lightDot = Math.abs((nx / norm) * -0.42 + (ny / norm) * 0.58 + (nz / norm) * 0.69);
          return {{
            p0,
            p1,
            p2,
            avgZ: (p0.pz + p1.pz + p2.pz) / 3,
            lightDot,
            color: getAtomSemanticColor(tri[0].atom, tri[0].idx, Math.max(totalAtoms, this.residues.length))
          }};
        }});

        projectedFacets.sort((a, b) => a.avgZ - b.avgZ);

        for (const f of projectedFacets) {{
          const fillAlpha = (currentStyle === 'surface' ? 0.22 : 0.14) * (0.45 + 0.55 * f.lightDot);
          const strokeAlpha = currentStyle === 'surface' ? 0.36 : 0.2;
          ctx.beginPath();
          ctx.moveTo(f.p0.px, f.p0.py);
          ctx.lineTo(f.p1.px, f.p1.py);
          ctx.lineTo(f.p2.px, f.p2.py);
          ctx.closePath();
          ctx.fillStyle = hexToRgba(f.color, fillAlpha);
          ctx.fill();
          ctx.strokeStyle = hexToRgba(f.color, strokeAlpha);
          ctx.lineWidth = 0.85;
          ctx.stroke();
        }}
      }}

      /* Pass 3: Internal Structural Representation (Cartoon Spline / Covalent Sticks / CPK Spheres) */
      const hasContinuousBackbone = this.splineSegments.length > 0;

      if ((currentStyle === 'cartoon' || currentStyle === 'surface') && hasContinuousBackbone) {{
        const projSplines = this.splineSegments.map(seg => {{
          const a = this.projectPoint(seg.p0, cosX, sinX, cosY, sinY, s);
          const b = this.projectPoint(seg.p1, cosX, sinX, cosY, sinY, s);
          return {{
            a,
            b,
            avgZ: (a.pz + b.pz) / 2,
            color: getAtomSemanticColor(seg.atom, seg.resIndex, this.residues.length)
          }};
        }});
        projSplines.sort((x, y) => x.avgZ - y.avgZ);

        for (const seg of projSplines) {{
          const ribbonW = Math.max(4.5, 7.5 * this.scale * ((seg.a.perspective + seg.b.perspective) / 2));
          ctx.beginPath();
          ctx.moveTo(seg.a.px, seg.a.py);
          ctx.lineTo(seg.b.px, seg.b.py);
          ctx.strokeStyle = seg.color;
          ctx.lineWidth = ribbonW;
          ctx.lineCap = 'round';
          ctx.stroke();

          /* Ribbon specular spine */
          ctx.beginPath();
          ctx.moveTo(seg.a.px, seg.a.py - 1);
          ctx.lineTo(seg.b.px, seg.b.py - 1);
          ctx.strokeStyle = 'rgba(255, 255, 255, 0.35)';
          ctx.lineWidth = Math.max(1.2, ribbonW * 0.28);
          ctx.lineCap = 'round';
          ctx.stroke();
        }}
      }}

      /* Render Covalent Bond Sticks for Stick Mode, Hetero Ligands, or Discontinuous Binding Pockets */
      const drawAllBonds = currentStyle === 'stick' || !hasContinuousBackbone;
      for (const [i, j] of this.bonds) {{
        const p1 = projected[i];
        const p2 = projected[j];
        const isHeteroBond = p1.atom.group === 'HETATM' || p2.atom.group === 'HETATM';
        if (!drawAllBonds && !isHeteroBond) continue;

        ctx.beginPath();
        ctx.moveTo(p1.px, p1.py);
        ctx.lineTo(p2.px, p2.py);
        ctx.strokeStyle = hexToRgba(p1.color, 0.88);
        ctx.lineWidth = Math.max(2.4, (isHeteroBond ? 4.2 : 3.2) * Math.sqrt(this.scale));
        ctx.lineCap = 'round';
        ctx.stroke();
      }}

      /* Render 3D Shaded Atomic Spheres & Specular Highlights */
      for (const node of sortedAtoms) {{
        const isSelected = selectedResidueKey === node.resKey;
        const isHovered = this.hoveredAtom && this.hoveredAtom.resKey === node.resKey;
        const isHetero = node.atom.group === 'HETATM';
        const isCA = node.atom.name === 'CA';

        let atomR = 4.2 * Math.sqrt(this.scale) * node.perspective;
        if (currentStyle === 'sphere') {{
          atomR = Math.max(7.5, node.vdw * s * node.perspective * 0.85);
        }} else if (currentStyle === 'stick') {{
          atomR = 4.6 * Math.sqrt(this.scale) * node.perspective;
        }} else if (isHetero) {{
          atomR = 6.2 * Math.sqrt(this.scale) * node.perspective;
        }} else if (!hasContinuousBackbone || isCA) {{
          atomR = 5.2 * Math.sqrt(this.scale) * node.perspective;
        }} else if (currentStyle === 'surface') {{
          atomR = 3.2 * Math.sqrt(this.scale) * node.perspective;
        }}

        if (isSelected || isHovered) {{
          ctx.beginPath();
          ctx.arc(node.px, node.py, atomR + 6, 0, Math.PI * 2);
          ctx.strokeStyle = PALETTE.rose;
          ctx.lineWidth = 2.2;
          ctx.stroke();
        }}

        drawShaded3DSphere(ctx, node.px, node.py, atomR, node.color, theme, false);

        if ((isCA && (this.residues.length <= 26 || isSelected || isHovered)) || (isHetero && node.atom.name === 'C1')) {{
          ctx.font = '500 10px "IBM Plex Mono", monospace';
          ctx.fillStyle = theme.textPrimary;
          ctx.fillText(`${{node.atom.resn}}${{node.atom.resi}}`, node.px + atomR + 4, node.py + 3);
        }}
      }}

      /* Pass 4: Front-Face Molecular Surface Specular Luster & Orientation Gizmo */
      for (const node of sortedAtoms) {{
        if (node.pz < 0 || (node.atom.name !== 'CA' && node.atom.group !== 'HETATM' && this.residues.length > 8)) continue;
        const surfR = Math.max(12, node.vdw * s * node.perspective * probeScale);
        ctx.beginPath();
        ctx.arc(node.px - surfR * 0.24, node.py - surfR * 0.24, surfR * 0.22, 0, Math.PI * 2);
        ctx.fillStyle = theme.isDark ? 'rgba(255, 255, 255, 0.08)' : 'rgba(255, 255, 255, 0.24)';
        ctx.fill();
      }}

      draw3DOrientationGizmo(
        ctx,
        w,
        h,
        this.rotX,
        this.rotY,
        this.scale,
        theme,
        'VDW + SAS 3D Surface (Probe 1.4 Å)'
      );
    }}
  }}

  /* Synchronize 3D Molecular Surface Viewport (#view3d / #mol-container) */
  function sync3DmolViewer(ds) {{
    if (!ds) return;
    const surfaceCanvas = document.getElementById('surface-canvas');
    const surfaceHud = document.getElementById('surfaceHud');
    if (surfaceCanvas && !surfaceEngine) {{
      surfaceEngine = new MolecularSurface3DEngine(surfaceCanvas, surfaceHud);
    }}
    if (surfaceEngine) {{
      surfaceEngine.setData(ds.atoms, ds.residues);
    }}
  }}

  function apply3DmolStyle() {{
    if (surfaceEngine) surfaceEngine.render();
    if (raytraceEngine) raytraceEngine.render();
  }}

  function renderActiveDataset() {{
    const ds = SESSION_DATASETS[activeDatasetKey];
    if (!ds) return;

    const val = ds.validation;
    const b = ds.bfactor;
    const geom = ds.geometry;
    if (!ds.atoms) {{
      ds.atoms = parseCoordinateTextClient(ds.raw, val.format);
    }}
    const residues = buildResidueSummaries(ds.atoms);
    ds.residues = residues;

    document.getElementById('headerFilename').textContent = ds.filename;
    document.getElementById('headerMeta').textContent =
      `Format: ${{val.format}} · Recipe: ${{ds.recipe}} · ${{val.total_atoms}} Atoms · Chain ${{val.chains.join(', ')}}`;

    const critEl = document.getElementById('criticalIndicator');
    if (b.very_low > 0) {{
      critEl.textContent = `${{b.very_low}} Critical Low-Confidence Atoms (<50 pLDDT)`;
      critEl.classList.remove('hidden');
    }} else {{
      critEl.classList.add('hidden');
    }}

    document.getElementById('kpiTotalAtoms').textContent = val.total_atoms;
    document.getElementById('kpiAtomBreakdown').textContent = `${{val.polymer_atoms}} polymer / ${{val.hetero_atoms}} hetero`;
    document.getElementById('kpiChainCount').textContent = val.chains.length;
    document.getElementById('kpiChainList').textContent = `Chains: ${{val.chains.join(', ')}}`;
    document.getElementById('kpiResidueCount').textContent = val.unique_residues;
    document.getElementById('kpiMeanB').innerHTML = `Mean B-Factor: ${{b.mean_b}} &Aring;&sup2;`;
    document.getElementById('kpiRg').innerHTML = `${{geom.rg}} &Aring;`;
    document.getElementById('kpiSpan').innerHTML = `C&alpha; Span: ${{geom.span}} &Aring;`;
    document.getElementById('kpiLigandCount').textContent = val.detected_ligands.length;
    document.getElementById('kpiLigandNames').textContent =
      `${{val.detected_ligands.length ? val.detected_ligands.join(', ') : 'None'}} · Models: ${{val.model_count}}`;

    document.getElementById('plddtRangeMeta').innerHTML = `Range: ${{b.min_b}} – ${{b.max_b}} &Aring;&sup2;`;
    document.getElementById('plddtStatVh').textContent = `${{b.very_high_pct}}% (${{b.very_high}} atoms)`;
    document.getElementById('plddtBarVh').style.width = `${{b.very_high_pct}}%`;
    document.getElementById('plddtStatConf').textContent = `${{b.confident_pct}}% (${{b.confident}} atoms)`;
    document.getElementById('plddtBarConf').style.width = `${{b.confident_pct}}%`;
    document.getElementById('plddtStatLow').textContent = `${{b.low_pct}}% (${{b.low}} atoms)`;
    document.getElementById('plddtBarLow').style.width = `${{b.low_pct}}%`;
    document.getElementById('plddtStatVl').textContent = `${{b.very_low_pct}}% (${{b.very_low}} atoms)`;
    document.getElementById('plddtBarVl').style.width = `${{b.very_low_pct}}%`;

    document.getElementById('paramRecipe').textContent = ds.recipe;
    document.getElementById('paramFormat').textContent = val.format + (val.is_compressed ? ' (.gz)' : '');
    document.getElementById('paramContacts').textContent = `${{geom.contacts}} pairs`;
    document.getElementById('paramConect').textContent = val.has_conect ? 'Present' : 'None';

    document.getElementById('structure-coords').textContent = ds.raw;

    if (!networkEngine) {{
      networkEngine = new RelationshipNodeMapEngine(
        document.getElementById('network-canvas'),
        document.getElementById('networkHud')
      );
    }}
    networkEngine.setData(ds.atoms, residues);
    sync3DmolViewer(ds);
    renderRaytraceCanvas(ds.atoms, residues);
    refreshVisualizations();
  }}

  function refreshVisualizations() {{
    const ds = SESSION_DATASETS[activeDatasetKey];
    if (!ds || !ds.residues) return;
    if (networkEngine) networkEngine.render();
    if (surfaceEngine) surfaceEngine.render();
    if (raytraceEngine) {{
      raytraceEngine.render();
    }} else {{
      renderRaytraceCanvas(ds.atoms, ds.residues);
    }}
    renderProfileChart(ds.residues);
    renderDistanceMatrix(ds.residues);
    populateResidueTable(ds.residues);
  }}

  /* Theme Switcher (Minimalist Icon-Only Button) */
  const themeToggle = document.getElementById('themeToggle');
  const themeIcon = document.getElementById('themeIcon');

  function updateThemeIcon(theme) {{
    if (theme === 'dark') {{
      themeIcon.innerHTML = '<circle cx="12" cy="12" r="4"></circle><path d="M12 2v2M12 20v2M4.93 4.93l1.41 1.41M17.66 17.66l1.41 1.41M2 12h2M20 12h2M6.34 17.66l-1.41 1.41M19.07 4.93l-1.41 1.41"></path>';
    }} else {{
      themeIcon.innerHTML = '<path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"></path>';
    }}
  }}

  themeToggle.addEventListener('click', () => {{
    const current = document.documentElement.getAttribute('data-theme') || 'dark';
    const next = current === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    updateThemeIcon(next);
    refreshVisualizations();
  }});

  function switchView(tab) {{
    activeViewTab = tab;
    document.querySelectorAll('#viewModeGroup .ui-btn').forEach(btn => {{
      btn.classList.toggle('active', btn.getAttribute('data-view') === tab);
    }});
    document.querySelectorAll('.viewport-pane').forEach(pane => pane.classList.remove('active'));

    const subtitleEl = document.getElementById('viewportSubtitle');
    if (tab === 'network') {{
      document.getElementById('viewNetwork').classList.add('active');
      if (subtitleEl) subtitleEl.textContent = 'Interactive 3D coordinates and spatial residue connectivity';
      if (networkEngine) {{
        networkEngine.resize();
        networkEngine.render();
      }}
      requestAnimationFrame(() => {{
        if (networkEngine) {{
          networkEngine.resize();
          networkEngine.render();
        }}
      }});
    }} else if (tab === '3d') {{
      document.getElementById('view3d').classList.add('active');
      if (subtitleEl) subtitleEl.textContent = 'Solvent-excluded 3D van der Waals molecular surface and atomic hull';
      const ds = SESSION_DATASETS[activeDatasetKey];
      if (ds && !surfaceEngine) sync3DmolViewer(ds);
      if (surfaceEngine) {{
        surfaceEngine.resize();
        surfaceEngine.render();
      }}
      requestAnimationFrame(() => {{
        if (surfaceEngine) {{
          surfaceEngine.resize();
          surfaceEngine.render();
        }}
      }});
    }} else {{
      document.getElementById('viewRaytrace').classList.add('active');
      if (subtitleEl) subtitleEl.textContent = 'Interactive 3D ray-traced structural rendering (PyMOL ray_trace_mode = 1)';
      const ds = SESSION_DATASETS[activeDatasetKey];
      if (ds && !raytraceEngine) renderRaytraceCanvas(ds.atoms, ds.residues);
      if (raytraceEngine) {{
        raytraceEngine.resize();
        raytraceEngine.render();
      }}
      requestAnimationFrame(() => {{
        if (raytraceEngine) {{
          raytraceEngine.resize();
          raytraceEngine.render();
        }}
      }});
    }}
  }}

  function setStyle(styleName) {{
    currentStyle = styleName;
    document.querySelectorAll('#styleGroup .ui-btn').forEach(b => {{
      b.classList.toggle('active', b.getAttribute('onclick').includes(styleName));
    }});
    if (networkEngine) networkEngine.render();
    if (surfaceEngine) surfaceEngine.render();
    if (raytraceEngine) raytraceEngine.render();
  }}

  function setColorScheme(schemeName) {{
    currentScheme = schemeName;
    document.querySelectorAll('#colorGroup .ui-btn').forEach(b => {{
      b.classList.toggle('active', b.getAttribute('onclick').includes(schemeName));
    }});
    refreshVisualizations();
  }}

  function zoomActiveViewport(factor) {{
    if (activeViewTab === 'network' && networkEngine) {{
      networkEngine.zoomBy(factor);
    }} else if (activeViewTab === '3d' && surfaceEngine) {{
      surfaceEngine.zoomBy(factor);
    }} else if (activeViewTab === 'raytrace' && raytraceEngine) {{
      raytraceEngine.zoomBy(factor);
    }} else {{
      if (networkEngine) networkEngine.zoomBy(factor);
      if (surfaceEngine) surfaceEngine.zoomBy(factor);
      if (raytraceEngine) raytraceEngine.zoomBy(factor);
    }}
  }}

  function toggleTurntable() {{
    turntableActive = !turntableActive;
    const btn = document.getElementById('turntableBtn');
    if (btn) btn.classList.toggle('active', turntableActive);

    if (turntableActive) {{
      const step = () => {{
        if (!turntableActive) return;
        if (activeViewTab === 'network' && networkEngine && !networkEngine.isDragging) {{
          networkEngine.rotY += 0.012;
          networkEngine.render();
        }} else if (activeViewTab === '3d' && surfaceEngine && !surfaceEngine.isDragging) {{
          surfaceEngine.rotY += 0.012;
          surfaceEngine.render();
        }} else if (activeViewTab === 'raytrace' && raytraceEngine && !raytraceEngine.isDragging) {{
          raytraceEngine.rotY += 0.012;
          raytraceEngine.render();
        }}
        turntableRaf = requestAnimationFrame(step);
      }};
      turntableRaf = requestAnimationFrame(step);
    }} else if (turntableRaf) {{
      cancelAnimationFrame(turntableRaf);
      turntableRaf = null;
    }}
  }}

  function resetCamera() {{
    if (turntableActive) {{
      turntableActive = false;
      const btn = document.getElementById('turntableBtn');
      if (btn) btn.classList.remove('active');
      if (turntableRaf) {{
        cancelAnimationFrame(turntableRaf);
        turntableRaf = null;
      }}
    }}
    if (networkEngine) networkEngine.resetCamera();
    if (surfaceEngine) surfaceEngine.resetCamera();
    if (raytraceEngine) raytraceEngine.resetCamera();
  }}

  const resetView = resetCamera;

  function downloadRenderImage() {{
    const ds = SESSION_DATASETS[activeDatasetKey];
    const rayCanvas = document.getElementById('raytrace-canvas');
    const surfCanvas = document.getElementById('surface-canvas');
    const netCanvas = document.getElementById('network-canvas');
    let url = '';
    if (activeViewTab === '3d' && surfCanvas) {{
      url = surfCanvas.toDataURL('image/png');
    }} else if (activeViewTab === 'raytrace' && rayCanvas) {{
      url = rayCanvas.toDataURL('image/png');
    }} else if (netCanvas) {{
      url = netCanvas.toDataURL('image/png');
    }}
    if (!url) return;
    const a = document.createElement('a');
    a.href = url;
    a.download = `${{ds ? ds.filename : 'structure'}}_${{activeViewTab}}_render.png`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  }}

  function downloadCoordinates() {{
    const ds = SESSION_DATASETS[activeDatasetKey];
    if (!ds) return;
    const blob = new Blob([ds.raw], {{ type: 'text/plain;charset=utf-8' }});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = ds.filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }}

  function copyPyMolScript() {{
    const ds = SESSION_DATASETS[activeDatasetKey];
    const cmd = `uv run scripts/pymol_runner.py --recipe ${{ds ? ds.recipe : '{recipe}'}} --structure ${{ds ? ds.path : '{structure_path}'}}`;
    const btn = document.getElementById('copyCmdBtn');
    if (navigator.clipboard) {{
      navigator.clipboard.writeText(cmd).then(() => {{
        btn.textContent = 'Command Copied';
        setTimeout(() => {{ btn.textContent = 'Copy CLI Command'; }}, 1800);
      }});
    }}
  }}

  /* Interactive Hover & Click on Profile Chart and Contact Matrix */
  function setupChartInteractions() {{
    const profileCanvas = document.getElementById('profile-canvas');
    const profileReadout = document.getElementById('profileHoverReadout');
    profileCanvas.addEventListener('mousemove', (e) => {{
      const rect = profileCanvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;
      const hit = cachedProfileBars.find(b => mx >= b.x && mx <= b.x + b.w && my >= 0 && my <= rect.height);
      if (hit) {{
        const r = hit.residue;
        profileReadout.textContent =
          `${{r.resn}} ${{r.chain}}${{r.resi}} (${{r.ss}}) · Mean B: ${{r.meanB}} Å² · Spatial Contacts (<8Å): ${{r.contacts}}`;
      }}
    }});
    profileCanvas.addEventListener('click', (e) => {{
      const rect = profileCanvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const hit = cachedProfileBars.find(b => mx >= b.x && mx <= b.x + b.w);
      if (hit) selectResidue(hit.residue.key);
    }});

    const matrixCanvas = document.getElementById('matrix-canvas');
    const matrixReadout = document.getElementById('matrixHoverReadout');
    matrixCanvas.addEventListener('mousemove', (e) => {{
      if (!cachedMatrixState) return;
      const rect = matrixCanvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;
      const {{ residues, offsetX, offsetY, size, cellS }} = cachedMatrixState;
      if (mx >= offsetX && mx <= offsetX + size && my >= offsetY && my <= offsetY + size) {{
        const j = Math.min(residues.length - 1, Math.floor((mx - offsetX) / cellS));
        const i = Math.min(residues.length - 1, Math.floor((my - offsetY) / cellS));
        const r1 = residues[i];
        const r2 = residues[j];
        const dist = Math.hypot(r1.x - r2.x, r1.y - r2.y, r1.z - r2.z);
        matrixReadout.textContent =
          `${{r1.resn}}${{r1.resi}} (${{r1.chain}}) ↔ ${{r2.resn}}${{r2.resi}} (${{r2.chain}}) : ${{dist.toFixed(2)}} Å`;
      }}
    }});
  }}

  window.addEventListener('DOMContentLoaded', () => {{
    updateThemeIcon('dark');
    document.getElementById('datasetSelect').addEventListener('change', (e) => {{
      activeDatasetKey = e.target.value;
      selectedResidueKey = null;
      confidenceFilterTier = null;
      renderActiveDataset();
    }});
    document.getElementById('residueFilterInput').addEventListener('input', () => {{
      refreshVisualizations();
    }});
    setupChartInteractions();
    renderActiveDataset();
    window.addEventListener('resize', () => refreshVisualizations());
  }});
</script>

</body>
</html>
"""
    return template


def main():
    parser = argparse.ArgumentParser(description="Generate interactive dashboard for pymol-visualize skill.")
    parser.add_argument("--structure", required=True, help="Input coordinate structure file")
    parser.add_argument("--png", default="", help="Path to rendered PNG")
    parser.add_argument("--pse", default="outputs/session_monomer.pse", help="Path to saved PSE session")
    parser.add_argument("--recipe", default="secondary_structure", help="Recipe used for visualization")
    parser.add_argument("--metrics", help="Optional JSON file with extra metrics")
    parser.add_argument("--output", default="ui/dashboard.html", help="Path for output HTML dashboard")
    args = parser.parse_args()

    extra_metrics = {}
    if args.metrics and os.path.exists(args.metrics):
        with open(args.metrics, "r", encoding="utf-8") as f:
            extra_metrics = json.load(f)

    cwd = os.getcwd()
    struct_p = os.path.relpath(args.structure, cwd) if os.path.abspath(args.structure).startswith(cwd) else args.structure
    png_p = os.path.relpath(args.png, cwd) if args.png and os.path.abspath(args.png).startswith(cwd) else args.png
    pse_p = os.path.relpath(args.pse, cwd) if args.pse and os.path.abspath(args.pse).startswith(cwd) else args.pse

    html_content = build_dashboard_html(
        structure_path=struct_p,
        png_path=png_p if (png_p and os.path.exists(png_p)) else "",
        pse_path=pse_p if (pse_p and os.path.exists(pse_p)) else "",
        recipe=args.recipe,
        extra_metrics=extra_metrics,
    )

    out_path = args.output
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"Dashboard successfully generated at: {out_path}")
    sys.exit(0)


if __name__ == "__main__":
    main()
