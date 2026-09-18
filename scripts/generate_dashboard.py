#!/usr/bin/env python3
"""
generate_dashboard.py
Generates an interactive, self-contained web dashboard (HTML/CSS/JS)
displaying macromolecular data processed by the pymol-visualize skill.
Output directory: ui/dashboard.html.
Supports light/dark theme toggle, interactive 3D WebGL viewing via 3Dmol.js,
high-resolution ray-traced render display, confidence metrics, and session export controls.
"""

import argparse
import base64
import json
import logging
import os
import sys

from validate_structure import _open_coordinate_file, validate_structure

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("generate_dashboard")


def calculate_bfactor_distribution(structure_path: str) -> dict:
    """
    Parses B-factors (pLDDT scores) directly from real coordinate records.
    Calculates exact counts and percentages across canonical AlphaFold bins:
    - >= 90: Very high confidence (Dark Blue)
    - 70-90: Confident (Light Blue)
    - 50-70: Low confidence (Yellow)
    - < 50:  Very low confidence (Orange)
    """
    b_factors = []
    with _open_coordinate_file(structure_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    b_val = float(line[60:66].strip())
                    b_factors.append(b_val)
                except (ValueError, IndexError):
                    pass

    if not b_factors:
        return {
            "total": 0,
            "very_high": 0, "very_high_pct": 0.0,
            "confident": 0, "confident_pct": 0.0,
            "low": 0, "low_pct": 0.0,
            "very_low": 0, "very_low_pct": 0.0,
            "mean_b": 0.0,
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
    }


def read_coordinate_snippet(structure_path: str, max_lines: int = 5000) -> str:
    """
    Reads the coordinate file text to embed in the dashboard for 3Dmol.js.
    """
    lines = []
    with _open_coordinate_file(structure_path) as f:
        for i, line in enumerate(f):
            if i >= max_lines:
                break
            lines.append(line)
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


def build_dashboard_html(
    structure_path: str,
    png_path: str,
    pse_path: str,
    recipe: str,
    extra_metrics: dict = None,
) -> str:
    val_data = validate_structure(structure_path)
    b_dist = calculate_bfactor_distribution(structure_path)
    coord_content = read_coordinate_snippet(structure_path)
    b64_png = get_base64_image(png_path)
    filename = os.path.basename(structure_path)
    coord_format = val_data["format"].lower()

    metrics = extra_metrics or {}
    rmsd_val = metrics.get("rmsd")
    distance_val = metrics.get("distance")

    template = f"""<!DOCTYPE html>
<html lang="en" data-theme="dark">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>PyMOL Structural Analytics Dashboard - {filename}</title>
  <!-- Local bundled 3Dmol-min.js with dynamic fallback to cdnjs -->
  <script src="vendor/3Dmol-min.js"></script>
  <script>
    if (typeof $3Dmol === 'undefined') {{
      var s = document.createElement('script');
      s.src = 'https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js';
      document.head.appendChild(s);
    }}
  </script>
  <style>
    :root {{
      --font-sans: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      --font-mono: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
    }}

    [data-theme="dark"] {{
      --bg-base: #090d16;
      --bg-surface: #111827;
      --bg-card: #1f2937;
      --bg-card-hover: #374151;
      --border-color: #374151;
      --border-focus: #38bdf8;
      --text-main: #f9fafb;
      --text-muted: #9ca3af;
      --accent: #38bdf8;
      --accent-hover: #0284c7;
      --accent-badge: rgba(56, 189, 248, 0.15);
      --plddt-vh: #0053d6;
      --plddt-c: #65cbf3;
      --plddt-l: #ffdb13;
      --plddt-vl: #ff7d45;
    }}

    [data-theme="light"] {{
      --bg-base: #f8fafc;
      --bg-surface: #ffffff;
      --bg-card: #f1f5f9;
      --bg-card-hover: #e2e8f0;
      --border-color: #cbd5e1;
      --border-focus: #0284c7;
      --text-main: #0f172a;
      --text-muted: #64748b;
      --accent: #0284c7;
      --accent-hover: #0369a1;
      --accent-badge: rgba(2, 132, 199, 0.12);
      --plddt-vh: #0053d6;
      --plddt-c: #38bdf8;
      --plddt-l: #eab308;
      --plddt-vl: #f97316;
    }}

    * {{
      box-sizing: border-box;
      margin: 0;
      padding: 0;
      transition: background-color 0.2s ease, border-color 0.2s ease, color 0.2s ease;
    }}

    body {{
      font-family: var(--font-sans);
      background-color: var(--bg-base);
      color: var(--text-main);
      line-height: 1.5;
      padding: 24px;
      min-height: 100vh;
    }}

    .dashboard-container {{
      max-width: 1440px;
      margin: 0 auto;
      display: flex;
      flex-direction: column;
      gap: 20px;
    }}

    /* Header Bar */
    header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      background: var(--bg-surface);
      border: 1px solid var(--border-color);
      border-radius: 12px;
      padding: 16px 24px;
    }}

    .header-info h1 {{
      font-size: 1.25rem;
      font-weight: 700;
      color: var(--text-main);
      display: flex;
      align-items: center;
      gap: 12px;
    }}

    .header-badge {{
      font-size: 0.75rem;
      font-family: var(--font-mono);
      background: var(--accent-badge);
      color: var(--accent);
      padding: 2px 8px;
      border-radius: 6px;
      border: 1px solid var(--accent);
      text-transform: uppercase;
    }}

    .header-sub {{
      font-size: 0.85rem;
      color: var(--text-muted);
      margin-top: 4px;
    }}

    /* Theme Toggle */
    .theme-toggle-btn {{
      background: var(--bg-card);
      border: 1px solid var(--border-color);
      color: var(--text-main);
      padding: 6px 14px;
      border-radius: 6px;
      cursor: pointer;
      font-size: 0.85rem;
      font-weight: 500;
      outline: none;
    }}

    .theme-toggle-btn:hover {{
      background: var(--bg-card-hover);
      border-color: var(--border-focus);
    }}

    /* Metrics Stat Strip */
    .stats-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 16px;
    }}

    .stat-card {{
      background: var(--bg-surface);
      border: 1px solid var(--border-color);
      border-radius: 10px;
      padding: 14px 18px;
    }}

    .stat-label {{
      font-size: 0.75rem;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      color: var(--text-muted);
    }}

    .stat-value {{
      font-size: 1.4rem;
      font-weight: 700;
      color: var(--text-main);
      font-family: var(--font-mono);
      margin-top: 4px;
    }}

    .stat-aux {{
      font-size: 0.8rem;
      color: var(--accent);
      margin-top: 2px;
    }}

    /* Main Workspace Grid */
    .main-grid {{
      display: grid;
      grid-template-columns: 1fr 400px;
      gap: 20px;
    }}

    @media (max-width: 1024px) {{
      .main-grid {{
        grid-template-columns: 1fr;
      }}
    }}

    /* Visualizer Container */
    .viewport-panel {{
      background: var(--bg-surface);
      border: 1px solid var(--border-color);
      border-radius: 12px;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }}

    .viewport-toolbar {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      padding: 12px 18px;
      background: var(--bg-surface);
      border-bottom: 1px solid var(--border-color);
    }}

    .tab-group {{
      display: flex;
      gap: 8px;
    }}

    .tab-btn {{
      background: var(--bg-card);
      border: 1px solid var(--border-color);
      color: var(--text-muted);
      padding: 6px 14px;
      border-radius: 6px;
      font-size: 0.82rem;
      font-weight: 600;
      cursor: pointer;
    }}

    .tab-btn.active {{
      background: var(--accent);
      color: #ffffff;
      border-color: var(--accent);
    }}

    .viewport-view {{
      position: relative;
      height: 520px;
      min-height: 520px;
      background: var(--bg-surface);
      display: none;
    }}

    .viewport-view.active {{
      display: block;
    }}

    #mol-container {{
      width: 100%;
      height: 520px;
      min-height: 520px;
      position: relative;
      display: block;
    }}

    .raytrace-img-container {{
      width: 100%;
      height: 100%;
      display: flex;
      align-items: center;
      justify-content: center;
      background: var(--bg-surface);
      padding: 12px;
    }}

    .raytrace-img-container img {{
      max-width: 100%;
      max-height: 100%;
      object-fit: contain;
      border-radius: 6px;
      border: 1px solid var(--border-color);
    }}

    /* Quick Interactive Controls Under Viewport */
    .controls-strip {{
      padding: 12px 18px;
      background: var(--bg-surface);
      border-top: 1px solid var(--border-color);
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 16px;
    }}

    .control-item {{
      display: inline-flex;
      align-items: center;
      gap: 8px;
      flex-shrink: 0;
      white-space: nowrap;
    }}

    .control-label {{
      font-size: 0.8rem;
      font-weight: 600;
      color: var(--text-muted);
      white-space: nowrap;
    }}

    .btn-group {{
      display: inline-flex;
      gap: 6px;
      flex-shrink: 0;
      white-space: nowrap;
    }}

    .ctrl-btn {{
      background: var(--bg-card);
      border: 1px solid var(--border-color);
      color: var(--text-main);
      padding: 5px 12px;
      border-radius: 6px;
      font-size: 0.78rem;
      font-weight: 500;
      cursor: pointer;
      white-space: nowrap;
    }}

    .ctrl-btn:hover {{
      background: var(--bg-card-hover);
    }}

    .ctrl-btn.active {{
      border-color: var(--accent);
      color: var(--accent);
    }}

    /* Sidebar Details */
    .sidebar-panel {{
      display: flex;
      flex-direction: column;
      gap: 16px;
    }}

    .panel-card {{
      background: var(--bg-surface);
      border: 1px solid var(--border-color);
      border-radius: 12px;
      padding: 18px;
    }}

    .panel-card h3 {{
      font-size: 0.95rem;
      font-weight: 700;
      margin-bottom: 12px;
      display: flex;
      align-items: center;
      justify-content: space-between;
    }}

    /* B-factor / pLDDT Progress Stack */
    .plddt-stack {{
      display: flex;
      flex-direction: column;
      gap: 10px;
    }}

    .plddt-row {{
      display: flex;
      flex-direction: column;
      gap: 4px;
    }}

    .plddt-meta {{
      display: flex;
      justify-content: space-between;
      font-size: 0.78rem;
      font-weight: 500;
    }}

    .plddt-bar-wrap {{
      height: 8px;
      background: var(--bg-card);
      border-radius: 4px;
      overflow: hidden;
    }}

    .plddt-bar-fill {{
      height: 100%;
      border-radius: 4px;
    }}

    .bar-vh {{ background-color: var(--plddt-vh); width: {b_dist['very_high_pct']}%; }}
    .bar-c  {{ background-color: var(--plddt-c);  width: {b_dist['confident_pct']}%; }}
    .bar-l  {{ background-color: var(--plddt-l);  width: {b_dist['low_pct']}%; }}
    .bar-vl {{ background-color: var(--plddt-vl); width: {b_dist['very_low_pct']}%; }}

    /* Download and Actions */
    .action-buttons {{
      display: flex;
      flex-direction: column;
      gap: 10px;
    }}

    .primary-action-btn {{
      display: flex;
      align-items: center;
      justify-content: center;
      gap: 8px;
      background: var(--accent);
      color: #ffffff;
      border: none;
      border-radius: 8px;
      padding: 10px;
      font-size: 0.85rem;
      font-weight: 600;
      cursor: pointer;
      text-decoration: none;
    }}

    .primary-action-btn:hover {{
      background: var(--accent-hover);
    }}

    .secondary-action-btn {{
      display: flex;
      align-items: center;
      justify-content: center;
      gap: 8px;
      background: var(--bg-card);
      color: var(--text-main);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      padding: 9px;
      font-size: 0.82rem;
      font-weight: 500;
      cursor: pointer;
      text-decoration: none;
    }}

    .secondary-action-btn:hover {{
      background: var(--bg-card-hover);
      border-color: var(--border-focus);
    }}

    /* Table */
    .data-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.82rem;
    }}

    .data-table th, .data-table td {{
      padding: 6px 8px;
      border-bottom: 1px solid var(--border-color);
      text-align: left;
    }}

    .data-table th {{
      color: var(--text-muted);
      font-weight: 600;
    }}

    .data-table td.mono {{
      font-family: var(--font-mono);
    }}
  </style>
</head>
<body>

<div class="dashboard-container">
  <!-- Header Bar -->
  <header>
    <div class="header-info">
      <h1>
        {filename}
        <span class="header-badge">{val_data['format']}</span>
        <span class="header-badge">{recipe}</span>
      </h1>
      <div class="header-sub">
        Processed via Headless PyMOL Engine - Validated {val_data['total_atoms']} atoms across chain(s) {', '.join(val_data['chains'])}
      </div>
    </div>
    <!-- Theme Toggle -->
    <button id="themeToggle" class="theme-toggle-btn" title="Toggle Theme" aria-label="Toggle Theme">Theme</button>
  </header>

  <!-- Summary Metrics Cards -->
  <div class="stats-grid">
    <div class="stat-card">
      <div class="stat-label">Total Atoms</div>
      <div class="stat-value">{val_data['total_atoms']}</div>
      <div class="stat-aux">{val_data['polymer_atoms']} polymer / {val_data['hetero_atoms']} hetero</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">Polypeptide Chains</div>
      <div class="stat-value">{len(val_data['chains'])}</div>
      <div class="stat-aux">Chains: {', '.join(val_data['chains'])}</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">Residue Count</div>
      <div class="stat-value">{val_data['unique_residues']}</div>
      <div class="stat-aux">Mean B-Factor: {b_dist['mean_b']} &Aring;&sup2;</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">Ensemble State</div>
      <div class="stat-value">{'NMR' if val_data['is_multi_state'] else 'Single Model'}</div>
      <div class="stat-aux">Models: {val_data['model_count']}</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">Detected Ligands</div>
      <div class="stat-value">{len(val_data['detected_ligands'])}</div>
      <div class="stat-aux">{', '.join(val_data['detected_ligands']) if val_data['detected_ligands'] else 'None'}</div>
    </div>
  </div>

  <!-- Main Viewport Grid -->
  <div class="main-grid">
    <div class="viewport-panel">
      <div class="viewport-toolbar">
        <div class="tab-group">
          <button class="tab-btn active" onclick="switchView('3d')">Interactive 3D View</button>
          <button class="tab-btn" onclick="switchView('raytrace')">Ray-Traced Render</button>
        </div>
        <div class="control-label">Preview</div>
      </div>

      <!-- Tab 1: Interactive 3Dmol.js WebGL Viewport -->
      <div id="view3d" class="viewport-view active">
        <div id="mol-container"></div>
      </div>

      <!-- Tab 2: High-Resolution Ray-Traced Output -->
      <div id="viewRaytrace" class="viewport-view">
        <div class="raytrace-img-container">
          {f'<img src="{b64_png}" alt="Ray-Traced Render" style="max-height:480px;">' if b64_png else '<p style="color:var(--text-muted);">Ray-traced render preview not available</p>'}
        </div>
      </div>

      <!-- Interactive 3D Styling Toolbar -->
      <div class="controls-strip">
        <div class="control-item">
          <span class="control-label">Style:</span>
          <div id="styleGroup" class="btn-group">
            <button class="ctrl-btn active" onclick="setStyle('cartoon')">Cartoon</button>
            <button class="ctrl-btn" onclick="setStyle('stick')">Sticks</button>
            <button class="ctrl-btn" onclick="setStyle('sphere')">Spheres</button>
            <button class="ctrl-btn" onclick="setStyle('surface')">Surface</button>
          </div>
        </div>

        <div class="control-item">
          <span class="control-label">Color Scheme:</span>
          <div id="colorGroup" class="btn-group">
            <button class="ctrl-btn active" onclick="setColorScheme('plddt')">AlphaFold pLDDT</button>
            <button class="ctrl-btn" onclick="setColorScheme('ss')">Secondary Structure</button>
            <button class="ctrl-btn" onclick="setColorScheme('chain')">Chain Palette</button>
            <button class="ctrl-btn" onclick="setColorScheme('spectrum')">Rainbow</button>
          </div>
        </div>

        <button class="ctrl-btn" style="margin-left: auto;" onclick="resetCamera()">Reset Camera</button>
      </div>
    </div>

    <!-- Sidebar Analytical Summary -->
    <div class="sidebar-panel">
      <!-- Confidence / B-factor Distribution -->
      <div class="panel-card">
        <h3>Confidence Distribution (pLDDT)</h3>
        <div class="plddt-stack">
          <div class="plddt-row">
            <div class="plddt-meta">
              <span style="color:var(--plddt-vh); font-weight:600;">Very High (&ge; 90)</span>
              <span>{b_dist['very_high_pct']}% ({b_dist['very_high']})</span>
            </div>
            <div class="plddt-bar-wrap">
              <div class="plddt-bar-fill bar-vh"></div>
            </div>
          </div>
          <div class="plddt-row">
            <div class="plddt-meta">
              <span style="color:var(--plddt-c); font-weight:600;">Confident (70 - 90)</span>
              <span>{b_dist['confident_pct']}% ({b_dist['confident']})</span>
            </div>
            <div class="plddt-bar-wrap">
              <div class="plddt-bar-fill bar-c"></div>
            </div>
          </div>
          <div class="plddt-row">
            <div class="plddt-meta">
              <span style="color:var(--plddt-l); font-weight:600;">Low (50 - 70)</span>
              <span>{b_dist['low_pct']}% ({b_dist['low']})</span>
            </div>
            <div class="plddt-bar-wrap">
              <div class="plddt-bar-fill bar-l"></div>
            </div>
          </div>
          <div class="plddt-row">
            <div class="plddt-meta">
              <span style="color:var(--plddt-vl); font-weight:600;">Very Low (&lt; 50)</span>
              <span>{b_dist['very_low_pct']}% ({b_dist['very_low']})</span>
            </div>
            <div class="plddt-bar-wrap">
              <div class="plddt-bar-fill bar-vl"></div>
            </div>
          </div>
        </div>
      </div>

      <!-- Additional Biophysical Metrics -->
      <div class="panel-card">
        <h3>Biophysical Measurements</h3>
        <table class="data-table">
          <tbody>
            <tr>
              <th>Metric</th>
              <th>Value</th>
            </tr>
            <tr>
              <td>Recipe Applied</td>
              <td class="mono">{recipe}</td>
            </tr>
            {f'<tr><td>Superposition RMSD</td><td class="mono">{rmsd_val:.3f} &Aring;</td></tr>' if rmsd_val is not None else ''}
            {f'<tr><td>Interatomic Distance</td><td class="mono">{distance_val:.3f} &Aring;</td></tr>' if distance_val is not None else ''}
            <tr>
              <td>Compressed on Disk</td>
              <td>{'Yes (.gz)' if val_data['is_compressed'] else 'No'}</td>
            </tr>
            <tr>
              <td>Has CONECT Records</td>
              <td>{'Yes' if val_data['has_conect'] else 'No'}</td>
            </tr>
          </tbody>
        </table>
      </div>

      <!-- Session Export and Downloads -->
      <div class="panel-card">
        <h3>Artifact Exports</h3>
        <div class="action-buttons">
          {f'<a href="{b64_png}" download="{filename}_render.png" class="primary-action-btn">Download Ray-Traced Image</a>' if b64_png else ''}
          {f'<button onclick="downloadSession()" class="secondary-action-btn">Download PyMOL Session (.pse)</button>' if pse_path else ''}
          <button onclick="copyPyMolScript()" class="secondary-action-btn">Copy PyMOL Script</button>
        </div>
      </div>
    </div>
  </div>
</div>

<!-- Embedded Raw Structure Coordinates for 3Dmol.js -->
<script type="text/plain" id="structure-coords">
{coord_content}
</script>

<script>
  let viewer = null;
  let fallback2D = null;
  let currentStyle = 'cartoon';
  let currentScheme = 'plddt';

  // Interactive 2D Software Renderer (Fallback when WebGL is restricted/unavailable)
  class Software2DViewer {{
    constructor(container, coordData, isDark) {{
      this.container = container;
      this.isDark = isDark;
      this.atoms = this.parseCoordinates(coordData);
      this.rotX = 0.35;
      this.rotY = -0.45;
      this.scale = 1.0;
      this.panX = 0;
      this.panY = 0;
      this.isDragging = false;
      this.lastMouseX = 0;
      this.lastMouseY = 0;

      this.setupDOM();
      this.calculateBounds();
      this.setupEvents();
      this.render();
    }}

    parseCoordinates(data) {{
      const atoms = [];
      const lines = data.split("\\n");
      for (let line of lines) {{
        line = line.trim();
        if (line.startsWith("ATOM") || line.startsWith("HETATM")) {{
          const name = line.substring(12, 16).trim();
          const resn = line.substring(17, 20).trim();
          const chain = line.substring(21, 22).trim() || "A";
          const resi = parseInt(line.substring(22, 26)) || 0;
          const x = parseFloat(line.substring(30, 38)) || 0;
          const y = parseFloat(line.substring(38, 46)) || 0;
          const z = parseFloat(line.substring(46, 54)) || 0;
          const b = parseFloat(line.substring(60, 66)) || 0;
          const elem = line.length >= 78 ? line.substring(76, 78).trim() : (name.startsWith("H") ? "H" : name[0]);
          atoms.push({{ name, resn, chain, resi, x, y, z, b, elem }});
        }}
      }}
      return atoms;
    }}

    calculateBounds() {{
      if (!this.atoms.length) return;
      let minX = Infinity, maxX = -Infinity;
      let minY = Infinity, maxY = -Infinity;
      let minZ = Infinity, maxZ = -Infinity;
      for (const a of this.atoms) {{
        minX = Math.min(minX, a.x); maxX = Math.max(maxX, a.x);
        minY = Math.min(minY, a.y); maxY = Math.max(maxY, a.y);
        minZ = Math.min(minZ, a.z); maxZ = Math.max(maxZ, a.z);
      }}
      this.cx = (minX + maxX) / 2;
      this.cy = (minY + maxY) / 2;
      this.cz = (minZ + maxZ) / 2;
      const span = Math.max(maxX - minX, maxY - minY, maxZ - minZ) || 1;
      this.baseScale = (Math.min(this.width, this.height) * 0.7) / span;
    }}

    setupDOM() {{
      this.container.innerHTML = '';
      this.container.style.position = 'relative';

      const badge = document.createElement('div');
      badge.style.cssText = 'position:absolute; top:12px; left:12px; z-index:10; background:rgba(0,0,0,0.65); border:1px solid var(--border-color); color:var(--accent); font-size:0.75rem; padding:4px 10px; border-radius:6px; font-weight:600; font-family:var(--font-mono); pointer-events:none;';
      badge.innerText = 'Interactive 2D Mode';
      this.container.appendChild(badge);

      this.canvas = document.createElement('canvas');
      this.canvas.style.cssText = 'width:100%; height:100%; display:block; cursor:grab;';
      this.container.appendChild(this.canvas);
      this.ctx = this.canvas.getContext('2d');
      this.resize();
    }}

    resize() {{
      const rect = this.container.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      this.width = rect.width || 700;
      this.height = rect.height || 520;
      this.canvas.width = this.width * dpr;
      this.canvas.height = this.height * dpr;
      this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }}

    setupEvents() {{
      const c = this.canvas;
      c.addEventListener('mousedown', (e) => {{
        this.isDragging = true;
        this.lastMouseX = e.clientX;
        this.lastMouseY = e.clientY;
        c.style.cursor = 'grabbing';
      }});

      window.addEventListener('mousemove', (e) => {{
        if (!this.isDragging) return;
        const dx = e.clientX - this.lastMouseX;
        const dy = e.clientY - this.lastMouseY;
        this.lastMouseX = e.clientX;
        this.lastMouseY = e.clientY;

        if (e.shiftKey || e.buttons === 2) {{
          this.panX += dx;
          this.panY += dy;
        }} else {{
          this.rotY += dx * 0.012;
          this.rotX += dy * 0.012;
        }}
        this.render();
      }});

      window.addEventListener('mouseup', () => {{
        if (this.isDragging) {{
          this.isDragging = false;
          c.style.cursor = 'grab';
        }}
      }});

      c.addEventListener('wheel', (e) => {{
        e.preventDefault();
        const factor = e.deltaY < 0 ? 1.1 : 0.9;
        this.scale *= factor;
        this.render();
      }}, {{ passive: false }});

      window.addEventListener('resize', () => {{
        this.resize();
        this.render();
      }});
    }}

    getColor(atom) {{
      if (currentScheme === 'plddt') {{
        const b = atom.b || 0;
        if (b >= 90) return '#0053d6';
        if (b >= 70) return '#65cbf3';
        if (b >= 50) return '#ffdb13';
        return '#ff7d45';
      }}
      if (currentScheme === 'ss') {{
        if (atom.name === 'CA' && atom.resi <= 12) return '#dc2626';
        if (atom.name === 'CA' && atom.resi >= 13) return '#f59e0b';
        return '#9ca3af';
      }}
      if (currentScheme === 'chain') {{
        const palette = ['#0d9488', '#f43f5e', '#16a34a', '#ea580c', '#3b82f6', '#d946ef'];
        return palette[(atom.chain.charCodeAt(0) || 65) % palette.length];
      }}
      const maxRes = Math.max(1, ...this.atoms.map(a => a.resi));
      const hue = ((atom.resi / maxRes) * 260) % 360;
      return 'hsl(' + hue + ', 80%, 50%)';
    }}

    render() {{
      const ctx = this.ctx;
      const w = this.width;
      const h = this.height;
      const isDark = document.documentElement.getAttribute('data-theme') === 'dark';

      ctx.fillStyle = isDark ? '#090d16' : '#f8fafc';
      ctx.fillRect(0, 0, w, h);

      if (!this.atoms.length) return;

      const cosX = Math.cos(this.rotX), sinX = Math.sin(this.rotX);
      const cosY = Math.cos(this.rotY), sinY = Math.sin(this.rotY);
      const totalScale = this.baseScale * this.scale;

      const transformed = [];
      for (const a of this.atoms) {{
        let x = a.x - this.cx;
        let y = a.y - this.cy;
        let z = a.z - this.cz;

        let x1 = x * cosY + z * sinY;
        let z1 = -x * sinY + z * cosY;

        let y2 = y * cosX - z1 * sinX;
        let z2 = y * sinX + z1 * cosX;

        const px = w / 2 + this.panX + x1 * totalScale;
        const py = h / 2 + this.panY - y2 * totalScale;

        transformed.push({{
          atom: a,
          px, py, pz: z2,
          color: this.getColor(a)
        }});
      }}

      const caAtoms = transformed.filter(t => t.atom.name === 'CA');
      if (currentStyle === 'cartoon' || currentStyle === 'surface') {{
        if (caAtoms.length > 1) {{
          for (let i = 0; i < caAtoms.length - 1; i++) {{
            const p1 = caAtoms[i];
            const p2 = caAtoms[i+1];
            if (p2.atom.resi === p1.atom.resi + 1 && p1.atom.chain === p2.atom.chain) {{
              ctx.beginPath();
              ctx.moveTo(p1.px, p1.py);
              ctx.lineTo(p2.px, p2.py);
              ctx.strokeStyle = p1.color;
              ctx.lineWidth = Math.max(3, 7 * (totalScale / 18));
              ctx.lineCap = 'round';
              ctx.stroke();
            }}
          }}
        }}
      }}

      if (currentStyle === 'stick') {{
        for (let i = 0; i < transformed.length; i++) {{
          for (let j = i + 1; j < Math.min(i + 15, transformed.length); j++) {{
            const a1 = transformed[i];
            const a2 = transformed[j];
            const dx = a1.atom.x - a2.atom.x;
            const dy = a1.atom.y - a2.atom.y;
            const dz = a1.atom.z - a2.atom.z;
            const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
            if (dist > 0.8 && dist < 1.9) {{
              ctx.beginPath();
              ctx.moveTo(a1.px, a1.py);
              ctx.lineTo(a2.px, a2.py);
              ctx.strokeStyle = a1.color;
              ctx.lineWidth = Math.max(1.5, 3.5 * (totalScale / 18));
              ctx.stroke();
            }}
          }}
        }}
      }}

      const renderAtoms = currentStyle === 'stick' ? transformed : caAtoms;
      renderAtoms.sort((a, b) => a.pz - b.pz);

      for (const t of renderAtoms) {{
        let r = 4;
        if (currentStyle === 'sphere') {{
          r = Math.max(5, 12 * (totalScale / 18));
        }} else if (currentStyle === 'cartoon') {{
          r = Math.max(3, 6 * (totalScale / 18));
        }} else if (currentStyle === 'surface') {{
          r = Math.max(6, 16 * (totalScale / 18));
        }}

        ctx.beginPath();
        ctx.arc(t.px, t.py, r, 0, Math.PI * 2);
        ctx.fillStyle = t.color;
        if (currentStyle === 'surface') ctx.globalAlpha = 0.6;
        ctx.fill();
        ctx.globalAlpha = 1.0;
        ctx.strokeStyle = isDark ? '#000000' : '#ffffff';
        ctx.lineWidth = 1;
        ctx.stroke();
      }}
    }}

    resetCamera() {{
      this.rotX = 0.35;
      this.rotY = -0.45;
      this.scale = 1.0;
      this.panX = 0;
      this.panY = 0;
      this.render();
    }}
  }}

  // Light/Dark Theme Handling
  const themeToggle = document.getElementById('themeToggle');

  function updateThemeButton(theme) {{
    themeToggle.textContent = theme === 'dark' ? 'Theme: Dark' : 'Theme: Light';
  }}

  themeToggle.addEventListener('click', () => {{
    const current = document.documentElement.getAttribute('data-theme') || 'dark';
    const next = current === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    updateThemeButton(next);
    if (viewer) {{
      viewer.setBackgroundColor(next === 'dark' ? 0x090d16 : 0xf8fafc);
      viewer.render();
    }}
    if (fallback2D) {{
      fallback2D.render();
    }}
  }});

  // Tab View Switcher
  function switchView(tab) {{
    document.querySelectorAll('.tab-btn').forEach(btn => btn.classList.remove('active'));
    document.querySelectorAll('.viewport-view').forEach(view => view.classList.remove('active'));

    if (tab === '3d') {{
      const btn = document.querySelector("button[onclick*='3d']");
      if (btn) btn.classList.add('active');
      const v = document.getElementById('view3d');
      if (v) v.classList.add('active');
      if (viewer) {{
        viewer.resize();
        viewer.render();
      }}
      if (fallback2D) {{
        fallback2D.resize();
        fallback2D.render();
      }}
    }} else {{
      const btn = document.querySelector("button[onclick*='raytrace']");
      if (btn) btn.classList.add('active');
      const v = document.getElementById('viewRaytrace');
      if (v) v.classList.add('active');
    }}
  }}

  // 3Dmol.js Initialization with automatic 2D fallback
  function start3DViewer() {{
    const container = document.getElementById('mol-container');
    const coordElem = document.getElementById('structure-coords');
    if (!container || !coordElem) return;

    const coordData = coordElem.textContent.trim();
    if (!coordData) {{
      console.warn("No coordinate data found. Switching to raytrace view.");
      switchView('raytrace');
      return;
    }}

    const isDark = document.documentElement.getAttribute('data-theme') === 'dark';

    try {{
      if (typeof $3Dmol === 'undefined') {{
        throw new Error("3Dmol script not loaded");
      }}
      const config = {{ backgroundColor: isDark ? '#090d16' : '#f8fafc' }};
      container.innerHTML = '';
      viewer = $3Dmol.createViewer(container, config);

      const fmt = '{coord_format}';
      const model = viewer.addModel(coordData, fmt === 'mmcif' ? 'cif' : 'pdb');
      if (!model) {{
        throw new Error("3Dmol failed to parse coordinate model.");
      }}

      applyCurrentStyle();
      viewer.zoomTo();
      viewer.render();
    }} catch (err) {{
      console.warn("WebGL unavailable or restricted; activating interactive 2D engine:", err);
      viewer = null;
      fallback2D = new Software2DViewer(container, coordData, isDark);
    }}
  }}

  function applyCurrentStyle() {{
    if (!viewer) return;
    viewer.removeAllSurfaces();

    const plddtFunc = function(atom) {{
      const b = atom.b || 0;
      if (b >= 90) return '#0053d6';
      if (b >= 70) return '#65cbf3';
      if (b >= 50) return '#ffdb13';
      return '#ff7d45';
    }};

    const ssFunc = function(atom) {{
      if (atom.ss === 'h' || atom.ss === 'h1' || atom.ss === 'h2') return '#dc2626';
      if (atom.ss === 's') return '#f59e0b';
      return '#9ca3af';
    }};

    const chainPalette = ['#0d9488', '#f43f5e', '#16a34a', '#ea580c', '#3b82f6', '#d946ef'];
    const chainFunc = function(atom) {{
      const charCode = (atom.chain || 'A').charCodeAt(0);
      return chainPalette[charCode % chainPalette.length];
    }};

    let styleObj = {{}};
    if (currentStyle === 'cartoon') {{
      if (currentScheme === 'plddt') {{
        styleObj = {{ cartoon: {{ colorfunc: plddtFunc }} }};
      }} else if (currentScheme === 'ss') {{
        styleObj = {{ cartoon: {{ colorfunc: ssFunc }} }};
      }} else if (currentScheme === 'chain') {{
        styleObj = {{ cartoon: {{ colorfunc: chainFunc }} }};
      }} else {{
        styleObj = {{ cartoon: {{ color: 'spectrum' }} }};
      }}
    }} else if (currentStyle === 'stick') {{
      if (currentScheme === 'plddt') {{
        styleObj = {{ stick: {{ radius: 0.25, colorfunc: plddtFunc }} }};
      }} else if (currentScheme === 'ss') {{
        styleObj = {{ stick: {{ radius: 0.25, colorfunc: ssFunc }} }};
      }} else if (currentScheme === 'chain') {{
        styleObj = {{ stick: {{ radius: 0.25, colorfunc: chainFunc }} }};
      }} else {{
        styleObj = {{ stick: {{ radius: 0.25, color: 'spectrum' }} }};
      }}
    }} else if (currentStyle === 'sphere') {{
      if (currentScheme === 'plddt') {{
        styleObj = {{ sphere: {{ scale: 0.35, colorfunc: plddtFunc }} }};
      }} else if (currentScheme === 'ss') {{
        styleObj = {{ sphere: {{ scale: 0.35, colorfunc: ssFunc }} }};
      }} else if (currentScheme === 'chain') {{
        styleObj = {{ sphere: {{ scale: 0.35, colorfunc: chainFunc }} }};
      }} else {{
        styleObj = {{ sphere: {{ scale: 0.35, color: 'spectrum' }} }};
      }}
    }} else if (currentStyle === 'surface') {{
      let surfOpts = {{ opacity: 0.45 }};
      if (currentScheme === 'plddt') {{
        styleObj = {{ cartoon: {{ opacity: 0.7, colorfunc: plddtFunc }} }};
        surfOpts.colorfunc = plddtFunc;
      }} else if (currentScheme === 'ss') {{
        styleObj = {{ cartoon: {{ opacity: 0.7, colorfunc: ssFunc }} }};
        surfOpts.colorfunc = ssFunc;
      }} else if (currentScheme === 'chain') {{
        styleObj = {{ cartoon: {{ opacity: 0.7, colorfunc: chainFunc }} }};
        surfOpts.colorfunc = chainFunc;
      }} else {{
        styleObj = {{ cartoon: {{ opacity: 0.7, color: 'spectrum' }} }};
        surfOpts.color = 'spectrum';
      }}
      viewer.addSurface($3Dmol.SurfaceType.VDW, surfOpts);
    }}

    viewer.setStyle({{}}, styleObj);
    viewer.render();
  }}

  function setStyle(styleName) {{
    currentStyle = styleName;
    document.querySelectorAll("#styleGroup .ctrl-btn").forEach(b => {{
      b.classList.toggle('active', b.getAttribute('onclick').includes(styleName));
    }});
    if (viewer) applyCurrentStyle();
    if (fallback2D) fallback2D.render();
  }}

  function setColorScheme(schemeName) {{
    currentScheme = schemeName;
    document.querySelectorAll("#colorGroup .ctrl-btn").forEach(b => {{
      b.classList.toggle('active', b.getAttribute('onclick').includes(schemeName));
    }});
    if (viewer) applyCurrentStyle();
    if (fallback2D) fallback2D.render();
  }}

  function resetCamera() {{
    if (viewer) {{
      viewer.zoomTo();
      viewer.render();
    }}
    if (fallback2D) {{
      fallback2D.resetCamera();
    }}
  }}

  function downloadSession() {{
    alert("PyMOL session file is preserved at: {pse_path}\\nOpen directly inside desktop PyMOL.");
  }}

  function copyPyMolScript() {{
    const cmd = "uv run scripts/pymol_runner.py --recipe {recipe} --structure {structure_path}";
    navigator.clipboard.writeText(cmd).then(() => {{
      alert("Reproducibility command copied to clipboard:\\n" + cmd);
    }});
  }}

  // Immediate and polling viewer initializer
  function initViewer() {{
    if (typeof $3Dmol !== 'undefined') {{
      start3DViewer();
    }} else {{
      let attempts = 0;
      const maxAttempts = 20; // 2 seconds
      const interval = setInterval(() => {{
        attempts++;
        if (typeof $3Dmol !== 'undefined') {{
          clearInterval(interval);
          start3DViewer();
        }} else if (attempts >= maxAttempts) {{
          clearInterval(interval);
          start3DViewer(); // Launches 2D engine fallback
        }}
      }}, 100);
    }}
  }}

  window.addEventListener('DOMContentLoaded', () => {{
    updateThemeButton('dark');
    initViewer();
  }});
</script>

</body>
</html>
"""
    return template


def main():
    parser = argparse.ArgumentParser(description="Generate interactive dashboard for pymol-visualize skill.")
    parser.add_argument("--structure", required=True, help="Input coordinate structure file")
    parser.add_argument("--png", default="outputs/render.png", help="Path to rendered PNG")
    parser.add_argument("--pse", default="outputs/session.pse", help="Path to saved PSE session")
    parser.add_argument("--recipe", default="secondary_structure", help="Recipe used for visualization")
    parser.add_argument("--metrics", help="Optional JSON file with extra metrics")
    parser.add_argument("--output", default="ui/dashboard.html", help="Path for output HTML dashboard")
    args = parser.parse_args()

    extra_metrics = {}
    if args.metrics and os.path.exists(args.metrics):
        with open(args.metrics, "r", encoding="utf-8") as f:
            extra_metrics = json.load(f)

    # Normalize paths to relative paths originating strictly from project root
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
