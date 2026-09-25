import math
import numpy as np
import pandas as pd
import html
import json
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors3D

def extract_poses(pdbqt_file_path):
    """
    Extracts individual 3D ligand poses from an AutoDock Vina output PDBQT file.
    Returns a list of pose coordinate strings.
    """
    poses = []
    current_pose = []
    try:
        with open(pdbqt_file_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('MODEL'):
                    current_pose = []
                elif line.startswith('ENDMDL'):
                    if current_pose:
                        poses.append("".join(current_pose))
                else:
                    if not line.startswith(('ENDROOT', 'TORSDOF', 'ROOT')):
                        current_pose.append(line)
        # If no MODEL tag was found but file contains coordinates (single pose)
        if not poses and current_pose:
            poses.append("".join(current_pose))
        return poses
    except Exception as e:
        print(f"Error extracting poses from {pdbqt_file_path}: {e}")
        return []

def calc_interactions(ligand_lines, receptor_pdbqt_path, cutoff=4.0):
    """
    Calculates non-covalent intermolecular interactions (H-bonds, Polar, Hydrophobic)
    between the ligand atoms and receptor binding site residues within cutoff distance (default 4.0 Å).
    """
    if isinstance(ligand_lines, str):
        ligand_lines = ligand_lines.split('\n')

    # 1. Parse Receptor Coordinates & Residues
    receptor_atoms = []
    try:
        with open(receptor_pdbqt_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        res_name = line[17:20].strip()
                        res_num = line[22:26].strip()
                        atom_name = line[12:16].strip()
                        element = line[76:78].strip() or (atom_name[0] if atom_name else 'C')
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        receptor_atoms.append({
                            "res": f"{res_name}-{res_num}",
                            "atom": atom_name,
                            "element": element,
                            "xyz": [x, y, z]
                        })
                    except ValueError:
                        pass
    except Exception as e:
        print(f"Error reading receptor {receptor_pdbqt_path}: {e}")
        return pd.DataFrame()

    # 2. Parse Ligand Coordinates
    ligand_atoms = []
    for idx, line in enumerate(ligand_lines):
        if line.startswith(("ATOM", "HETATM")):
            try:
                atom_name = line[12:16].strip()
                element = line[76:78].strip() or (atom_name[0] if atom_name else 'C')
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                ligand_atoms.append({
                    "id": f"Atom {idx+1} ({element})",
                    "element": element,
                    "xyz": [x, y, z]
                })
            except ValueError:
                pass

    # 3. Pairwise Non-Covalent Contact Detection with Granular Classification
    interactions = []
    seen = set()

    for lig in ligand_atoms:
        lx, ly, lz = lig["xyz"]
        for rec in receptor_atoms:
            rx, ry, rz = rec["xyz"]
            dist = math.sqrt((lx - rx)**2 + (ly - ry)**2 + (lz - rz)**2)
            
            if dist <= cutoff:
                key = (rec["res"], lig["id"])
                if key not in seen:
                    res_code = rec["res"].split("-")[0].upper()
                    
                    # Detailed Classification:
                    if dist <= 3.3 and lig["element"] in ["O", "N", "F", "S", "HD"] and rec["element"] in ["O", "N", "S", "HD"]:
                        bond_type = "Hydrogen Bond"
                        bond_color = "#FF3B30"
                    elif res_code in ["ASP", "GLU", "LYS", "ARG", "HIS"] and dist <= 3.8 and lig["element"] in ["O", "N"]:
                        bond_type = "Salt Bridge / Electrostatic"
                        bond_color = "#FFD60A"
                    elif res_code in ["PHE", "TYR", "TRP", "HIS"] and dist <= 4.0:
                        bond_type = "π-π / Aromatic Contact"
                        bond_color = "#BF5AF2"
                    elif res_code in ["LEU", "ILE", "VAL", "ALA", "PRO", "MET", "CYS"]:
                        bond_type = "Hydrophobic Aliphatic"
                        bond_color = "#64D2FF"
                    else:
                        bond_type = "Van der Waals Contact"
                        bond_color = "#30D158"
                        
                    interactions.append({
                        "Receptor Residue": rec["res"],
                        "Receptor Atom": rec["atom"],
                        "Ligand Atom": lig["id"],
                        "Distance (Å)": round(dist, 2),
                        "Interaction Type": bond_type,
                        "color": bond_color,
                        "Color": bond_color,
                        "lig_xyz": lig["xyz"],
                        "Ligand XYZ": lig["xyz"],
                        "rec_xyz": rec["xyz"],
                        "Receptor XYZ": rec["xyz"]
                    })
                    seen.add(key)

    df_inter = pd.DataFrame(interactions)
    if not df_inter.empty:
        df_inter = df_inter.sort_values(by="Distance (Å)").reset_index(drop=True)
    return df_inter

def calc_advanced_ligand_efficiency(affinity_kcal, smiles, ki_molar):
    """
    Computes rigorous medicinal chemistry efficiency metrics:
    - Heavy Atom Count (N_heavy)
    - Ligand Efficiency (LE = -ΔG / N_heavy)
    - Binding Lipophilicity Efficiency (LipE = pKi - cLogP)
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"n_heavy": 0, "le": 0.0, "lipe": 0.0, "pki": 0.0, "clogp": 0.0}
        
    n_heavy = mol.GetNumHeavyAtoms()
    clogp = Descriptors.MolLogP(mol)
    
    le = (-affinity_kcal / n_heavy) if n_heavy > 0 else 0.0
    pki = -math.log10(max(ki_molar, 1e-15))
    lipe = pki - clogp
    
    return {
        "n_heavy": n_heavy,
        "clogp": round(clogp, 2),
        "pki": round(pki, 2),
        "le": round(le, 3),
        "lipe": round(lipe, 2)
    }

def _extract_pocket_sidechains_pdb(receptor_data: str, ligand_data: str, cutoff: float = 5.5) -> str:
    """
    Extracts receptor residues within `cutoff` Å of the ligand heavy atoms
    so the 3D Cinema Studio can render active-site side-chains as explicit sticks
    and wrap a focused electrostatic pocket surface around the cavity.
    """
    if not receptor_data or not ligand_data:
        return ""
    lig_pts = []
    for line in ligand_data.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                lig_pts.append((x, y, z))
            except Exception:
                pass
    if not lig_pts:
        return ""

    lig_arr = np.array(lig_pts, dtype=float)
    active_res_keys = set()
    rec_lines = receptor_data.splitlines()

    for line in rec_lines:
        if line.startswith("ATOM"):
            try:
                rx = float(line[30:38].strip())
                ry = float(line[38:46].strip())
                rz = float(line[46:54].strip())
                dists = np.linalg.norm(lig_arr - np.array([rx, ry, rz]), axis=1)
                if float(np.min(dists)) <= cutoff:
                    res_key = (line[17:20].strip(), line[21:22].strip(), line[22:26].strip())
                    active_res_keys.add(res_key)
            except Exception:
                pass

    pocket_pdb_lines = []
    for line in rec_lines:
        if line.startswith("ATOM"):
            res_key = (line[17:20].strip(), line[21:22].strip(), line[22:26].strip())
            if res_key in active_res_keys:
                clean_line = line[:66].ljust(76) + line[12:14].strip().rjust(2)
                pocket_pdb_lines.append(clean_line)
    return "\n".join(pocket_pdb_lines)


def build_3dmol_html(container_id, receptor_data, ligand_data, interactions_df=None, receptor_style='cartoon', ligand_style='stick', show_surface=False, height=550):
    """
    Constructs a VMD / ChimeraX / PyMOL-grade 3D Docking Complex Cinema Studio
    featuring:
      - Auto-zoomed Active-Site Pocket Close-Up vs Whole Protein view
      - Explicit Pocket Side-Chain Sticks (< 5.5 Å) with residue callouts
      - Translucent Electrostatic Pocket Cavity Cloud
      - Non-Covalent H-Bond & Hydrophobic Force Vectors with Å distances
      - Interactive Atom-to-Atom Caliper Distance Tool
      - Continuous 360° Turntable Cinema Orbit
      - 1-Click 4K PNG Snapshot & 60-FPS HD WebM Video Recorder
    """
    pocket_pdb_str = _extract_pocket_sidechains_pdb(receptor_data, ligand_data, cutoff=5.5)

    interactions_list = []
    if interactions_df is not None and not interactions_df.empty:
        for _, row in interactions_df.iterrows():
            rx, ry, rz = row.get("Receptor XYZ", row.get("rec_xyz", [0, 0, 0]))
            lx, ly, lz = row.get("Ligand XYZ", row.get("lig_xyz", [0, 0, 0]))
            color = row.get("Color", row.get("color", "#FF3366"))
            res_label = str(row.get("Receptor Residue", "RES"))
            dist = str(row.get("Distance (Å)", "2.9"))
            interactions_list.append({
                "rx": float(rx), "ry": float(ry), "rz": float(rz),
                "lx": float(lx), "ly": float(ly), "lz": float(lz),
                "color": color, "label": f"{res_label} ({dist} Å)"
            })

    rec_json = json.dumps(receptor_data or "")
    lig_json = json.dumps(ligand_data or "")
    pocket_json = json.dumps(pocket_pdb_str or "")
    inter_json = json.dumps(interactions_list)
    init_surf_bool = "true" if show_surface else "false"

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        * {{ box-sizing: border-box; user-select: none; }}
        body, html {{
            margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden;
            background: radial-gradient(circle at 50% 42%, #131B2E 0%, #090D16 100%);
            font-family: -apple-system, BlinkMacSystemFont, "SF Pro Display", "Inter", sans-serif;
        }}
        #studio-wrapper {{
            width: 100%; height: {height}px; position: relative;
            border-radius: 14px; border: 1px solid rgba(255,255,255,0.14);
            overflow: hidden; box-shadow: 0 18px 42px rgba(0,0,0,0.55);
        }}
        #top-bar {{
            position: absolute; top: 10px; left: 10px; right: 10px; z-index: 20;
            display: flex; flex-wrap: wrap; gap: 6px; align-items: center; justify-content: space-between;
            background: rgba(11, 15, 25, 0.86); backdrop-filter: blur(14px);
            padding: 7px 12px; border-radius: 10px; border: 1px solid rgba(255,255,255,0.12);
        }}
        .btn-group {{ display: flex; flex-wrap: wrap; gap: 5px; align-items: center; }}
        .c-chip {{
            background: rgba(255,255,255,0.06); color: #CBD5E1;
            border: 1px solid rgba(255,255,255,0.14); border-radius: 7px;
            padding: 4px 9px; font-size: 10.5px; font-weight: 600; cursor: pointer;
            transition: all 0.15s ease; display: inline-flex; align-items: center; gap: 4px;
        }}
        .c-chip:hover {{ background: rgba(30, 41, 59, 0.95); color: #FFF; border-color: #64D2FF; }}
        .c-chip.active {{
            background: rgba(10, 132, 255, 0.25); color: #64D2FF;
            border-color: rgba(100, 210, 255, 0.55);
        }}
        .c-chip.rec-btn {{
            background: rgba(255, 69, 58, 0.22); color: #FF6961;
            border-color: rgba(255, 69, 58, 0.5); font-weight: 700;
        }}
        .c-chip.rec-btn.recording {{
            background: #FF3B30; color: #FFF; animation: pulseRec 1s infinite;
        }}
        @keyframes pulseRec {{
            0% {{ box-shadow: 0 0 0 0 rgba(255,59,48,0.7); }}
            70% {{ box-shadow: 0 0 0 8px rgba(255,59,48,0); }}
            100% {{ box-shadow: 0 0 0 0 rgba(255,59,48,0); }}
        }}
        #bottom-hud {{
            position: absolute; bottom: 10px; left: 10px; right: 10px; z-index: 20;
            display: flex; justify-content: space-between; align-items: center;
            background: rgba(11, 15, 25, 0.88); backdrop-filter: blur(12px);
            padding: 6px 14px; border-radius: 9px; border: 1px solid rgba(255,255,255,0.12);
            font-size: 11px; color: #E2E8F0;
        }}
        .dot {{ display: inline-block; width: 8px; height: 8px; border-radius: 50%; margin-right: 4px; }}
        #caliper-status {{ color: #FFD60A; font-family: "JetBrains Mono", monospace; font-weight: 600; }}
    </style>
</head>
<body>
    <div id="studio-wrapper">
        <div id="top-bar">
            <div class="btn-group">
                <button id="btn-pocket-zoom" class="c-chip active" title="Auto-Zoom into Active Binding Cavity">🎯 Pocket Close-Up</button>
                <button id="btn-full-zoom" class="c-chip" title="View Full Protein Backbone">🌐 Whole Protein</button>
                <button id="btn-sidechains" class="c-chip active" title="Show Active-Site Side-Chain Sticks (<5.5 Å)">🧬 Side-Chains</button>
                <button id="btn-vectors" class="c-chip active" title="Show H-Bond & Hydrophobic Force Vectors">⚡ H-Bonds ({len(interactions_list)})</button>
                <button id="btn-surface" class="c-chip" title="Show Electrostatic Pocket Cavity Surface">🧊 Cavity Cloud</button>
                <button id="btn-caliper" class="c-chip" title="Click any 2 atoms to measure distance in Å">📏 Caliper Ruler</button>
            </div>
            <div class="btn-group">
                <button id="btn-orbit" class="c-chip active" title="Continuous 360° Cinema Turntable">🌀 360° Orbit</button>
                <button id="btn-snap" class="c-chip" title="Download High-Res PNG Snapshot">📸 Snapshot</button>
                <button id="btn-record" class="c-chip rec-btn" title="Record 60-FPS HD WebM Turntable Movie">🎥 Record HD Video</button>
            </div>
        </div>
        <div id="{container_id}" style="width: 100%; height: 100%;"></div>
        <div id="bottom-hud">
            <div style="display:flex; gap:14px; align-items:center;">
                <span><span class="dot" style="background:#30D158;"></span><b>Docked Ligand (Stick+Sphere)</b></span>
                <span><span class="dot" style="background:#94A3B8;"></span><b>Pocket Side-Chains (&lt;5.5 Å)</b></span>
                <span><span class="dot" style="background:#FF3366;"></span><b>Polar H-Bond</b></span>
                <span><span class="dot" style="background:#00D2FF;"></span><b>Hydrophobic Contact</b></span>
            </div>
            <div id="caliper-status">PyMOL / ChimeraX Cinema Ready</div>
        </div>
    </div>
    <script>
        (function() {{
            var receptorStr = {rec_json};
            var ligandStr = {lig_json};
            var pocketStr = {pocket_json};
            var interactions = {inter_json};
            var showSidechains = true;
            var showVectors = true;
            var showSurface = {init_surf_bool};
            var caliperMode = false;
            var caliperAtom1 = null;
            var autoOrbit = true;
            var orbitTimer = null;
            var surfObj = null;
            var vectorShapes = [];
            var vectorLabels = [];
            var caliperShapes = [];
            var mediaRecorder = null;
            var recordedChunks = [];
            var viewer = null;

            function drawVectors() {{
                for (var i = 0; i < vectorShapes.length; i++) viewer.removeShape(vectorShapes[i]);
                for (var j = 0; j < vectorLabels.length; j++) viewer.removeLabel(vectorLabels[j]);
                vectorShapes = [];
                vectorLabels = [];
                if (!showVectors) return;
                for (var k = 0; k < interactions.length; k++) {{
                    var it = interactions[k];
                    var s = viewer.addCylinder({{
                        start: {{x: it.rx, y: it.ry, z: it.rz}},
                        end: {{x: it.lx, y: it.ly, z: it.lz}},
                        radius: 0.09, dashed: true, color: it.color,
                        fromCap: 1, toCap: 1
                    }});
                    vectorShapes.push(s);
                    var lbl = viewer.addLabel(it.label, {{
                        position: {{x: (it.rx + it.lx)/2.0, y: (it.ry + it.ly)/2.0, z: (it.rz + it.lz)/2.0}},
                        backgroundColor: 'rgba(11, 15, 25, 0.88)',
                        borderColor: it.color, borderThickness: 1,
                        fontColor: '#FFFFFF', fontSize: 10.5, inFront: true
                    }});
                    vectorLabels.push(lbl);
                }}
            }}

            function updateSurface() {{
                if (surfObj !== null) {{
                    try {{ viewer.removeSurface(surfObj); }} catch(e) {{}}
                    surfObj = null;
                }}
                if (showSurface) {{
                    var targetModel = (pocketStr && pocketStr.trim().length > 10) ? 2 : 0;
                    surfObj = viewer.addSurface($3Dmol.SurfaceType.VDW, {{
                        opacity: 0.34, colorscheme: 'whiteCarbon'
                    }}, {{model: targetModel}});
                }}
            }}

            function applyStyles() {{
                if (viewer.getModel(0)) {{
                    viewer.setStyle({{model: 0}}, {{{receptor_style}: {{color: 'spectrum', opacity: 0.55}}}});
                }}
                if (viewer.getModel(1)) {{
                    viewer.setStyle({{model: 1}}, {{
                        stick: {{colorscheme: 'greenCarbon', radius: 0.22}},
                        sphere: {{colorscheme: 'greenCarbon', scale: 0.27}}
                    }});
                }}
                if (viewer.getModel(2)) {{
                    if (showSidechains) {{
                        viewer.setStyle({{model: 2}}, {{
                            stick: {{colorscheme: 'whiteCarbon', radius: 0.13, opacity: 0.92}}
                        }});
                    }} else {{
                        viewer.setStyle({{model: 2}}, {{}});
                    }}
                }}
            }}

            var waitTimer = setInterval(function() {{
                if (typeof $3Dmol !== 'undefined') {{
                    clearInterval(waitTimer);
                    var el = document.getElementById("{container_id}");
                    viewer = $3Dmol.createViewer(el, {{defaultcolors: $3Dmol.rasmolElementColors, antialias: true}});
                    viewer.setBackgroundColor(0x0B0F19, 1.0);

                    if (receptorStr && receptorStr.trim().length > 0) viewer.addModel(receptorStr, "pdb");
                    if (ligandStr && ligandStr.trim().length > 0) {{
                        var fmt = ligandStr.indexOf("$$$$") !== -1 ? "sdf" : "pdb";
                        viewer.addModel(ligandStr, fmt);
                    }}
                    if (pocketStr && pocketStr.trim().length > 0) viewer.addModel(pocketStr, "pdb");

                    applyStyles();
                    drawVectors();
                    if (showSurface) {{
                        document.getElementById('btn-surface').classList.add('active');
                        updateSurface();
                    }}

                    if (viewer.getModel(1)) {{
                        viewer.zoomTo({{model: 1}});
                        viewer.zoom(0.75);
                    }} else {{
                        viewer.zoomTo();
                    }}
                    viewer.render();

                    // Continuous 360 Turntable Orbit
                    orbitTimer = setInterval(function() {{
                        if (autoOrbit && viewer) {{
                            viewer.rotate(0.45, "y");
                            viewer.render();
                        }}
                    }}, 40);

                    // Caliper click handler
                    viewer.setClickable({{}}, true, function(atom) {{
                        if (!caliperMode || !atom) return;
                        if (!caliperAtom1) {{
                            caliperAtom1 = atom;
                            document.getElementById('caliper-status').innerText = '📏 First atom selected (' + (atom.resn || 'LIG') + ':' + atom.atom + '). Click second atom...';
                        }} else {{
                            var dx = atom.x - caliperAtom1.x;
                            var dy = atom.y - caliperAtom1.y;
                            var dz = atom.z - caliperAtom1.z;
                            var dist = Math.sqrt(dx*dx + dy*dy + dz*dz).toFixed(2);
                            var cLine = viewer.addCylinder({{
                                start: {{x: caliperAtom1.x, y: caliperAtom1.y, z: caliperAtom1.z}},
                                end: {{x: atom.x, y: atom.y, z: atom.z}},
                                radius: 0.08, dashed: true, color: '#FFD60A'
                            }});
                            var cLbl = viewer.addLabel('📏 ' + dist + ' Å', {{
                                position: {{x: (caliperAtom1.x+atom.x)/2, y: (caliperAtom1.y+atom.y)/2, z: (caliperAtom1.z+atom.z)/2}},
                                backgroundColor: '#FFD60A', fontColor: '#000000', fontSize: 11, inFront: true
                            }});
                            caliperShapes.push(cLine, cLbl);
                            document.getElementById('caliper-status').innerText = '📏 Measured: ' + dist + ' Å (' + (caliperAtom1.resn||'A') + ':' + caliperAtom1.atom + ' ↔ ' + (atom.resn||'B') + ':' + atom.atom + ')';
                            caliperAtom1 = null;
                            viewer.render();
                        }}
                    }});

                    // Button Controls
                    document.getElementById('btn-pocket-zoom').onclick = function() {{
                        this.classList.add('active');
                        document.getElementById('btn-full-zoom').classList.remove('active');
                        if (viewer.getModel(1)) {{ viewer.zoomTo({{model: 1}}); viewer.zoom(0.75); viewer.render(); }}
                    }};
                    document.getElementById('btn-full-zoom').onclick = function() {{
                        this.classList.add('active');
                        document.getElementById('btn-pocket-zoom').classList.remove('active');
                        viewer.zoomTo({{model: 0}}); viewer.render();
                    }};
                    document.getElementById('btn-sidechains').onclick = function() {{
                        showSidechains = !showSidechains;
                        this.classList.toggle('active', showSidechains);
                        applyStyles(); viewer.render();
                    }};
                    document.getElementById('btn-vectors').onclick = function() {{
                        showVectors = !showVectors;
                        this.classList.toggle('active', showVectors);
                        drawVectors(); viewer.render();
                    }};
                    document.getElementById('btn-surface').onclick = function() {{
                        showSurface = !showSurface;
                        this.classList.toggle('active', showSurface);
                        updateSurface(); viewer.render();
                    }};
                    document.getElementById('btn-caliper').onclick = function() {{
                        caliperMode = !caliperMode;
                        caliperAtom1 = null;
                        this.classList.toggle('active', caliperMode);
                        document.getElementById('caliper-status').innerText = caliperMode ? '📏 Caliper Active: Click any 2 atoms to measure distance' : 'PyMOL / ChimeraX Cinema Ready';
                    }};
                    document.getElementById('btn-orbit').onclick = function() {{
                        autoOrbit = !autoOrbit;
                        this.classList.toggle('active', autoOrbit);
                    }};
                    document.getElementById('btn-snap').onclick = function() {{
                        var canvas = document.querySelector('#{container_id} canvas');
                        if (!canvas) return;
                        var a = document.createElement('a');
                        a.href = canvas.toDataURL('image/png');
                        a.download = 'EthnoDock_3D_Docking_Complex_4K.png';
                        a.click();
                    }};
                    document.getElementById('btn-record').onclick = function() {{
                        var btn = this;
                        var canvas = document.querySelector('#{container_id} canvas');
                        if (!canvas || typeof canvas.captureStream !== 'function') return;
                        if (mediaRecorder && mediaRecorder.state === 'recording') {{
                            mediaRecorder.stop();
                            return;
                        }}
                        recordedChunks = [];
                        autoOrbit = true;
                        document.getElementById('btn-orbit').classList.add('active');
                        var stream = canvas.captureStream(60);
                        var mime = MediaRecorder.isTypeSupported('video/webm;codecs=vp9') ? 'video/webm;codecs=vp9' : 'video/webm';
                        mediaRecorder = new MediaRecorder(stream, {{mimeType: mime, videoBitsPerSecond: 6000000}});
                        mediaRecorder.ondataavailable = function(e) {{ if (e.data && e.data.size > 0) recordedChunks.push(e.data); }};
                        mediaRecorder.onstop = function() {{
                            var blob = new Blob(recordedChunks, {{type: 'video/webm'}});
                            var url = URL.createObjectURL(blob);
                            var a = document.createElement('a');
                            a.href = url;
                            a.download = 'EthnoDock_Docking_Complex_Cinema_60FPS.webm';
                            a.click();
                            btn.classList.remove('recording');
                            btn.innerText = '🎥 Record HD Video';
                        }};
                        btn.classList.add('recording');
                        btn.innerText = '⏺️ Recording (5s)...';
                        mediaRecorder.start();
                        setTimeout(function() {{
                            if (mediaRecorder && mediaRecorder.state === 'recording') mediaRecorder.stop();
                        }}, 5000);
                    }};
                }}
            }}, 80);
        }})();
    </script>
</body>
</html>"""
    return html_content


def build_dual_pose_comparison_3dmol_html(container_id, receptor_data, parent_ligand_data, var_ligand_data, var_interactions_df=None, parent_name="Parent Phytochemical", var_name="Optimized Derivative", height=550):
    """
    Constructs a VMD / ChimeraX / PyMOL-grade Dual-Pose Bioisosteric Superposition & Pharmacophore Studio:
      - Superimposes Natural Parent (Gold Stick) and Bioisosteric Derivative (Cyan Stick)
      - Renders 3D Pharmacophore Volume Envelopes (Translucent Gold vs Electric Cyan Clouds)
        showing steric/electronic sub-pocket burial gains
      - Displays Active-Site Pocket Side-Chains (< 5.5 Å)
      - Interactive Distance Caliper, 360° Orbit, 4K Snapshot & 60-FPS HD WebM Video Recorder
    """
    pocket_pdb_str = _extract_pocket_sidechains_pdb(receptor_data, var_ligand_data or parent_ligand_data, cutoff=5.5)

    interactions_list = []
    if var_interactions_df is not None and not var_interactions_df.empty:
        for _, row in var_interactions_df.iterrows():
            rx, ry, rz = row.get("Receptor XYZ", row.get("rec_xyz", [0, 0, 0]))
            lx, ly, lz = row.get("Ligand XYZ", row.get("lig_xyz", [0, 0, 0]))
            color = row.get("Color", row.get("color", "#00D2FF"))
            res_label = str(row.get("Receptor Residue", "RES"))
            dist = str(row.get("Distance (Å)", "2.8"))
            interactions_list.append({
                "rx": float(rx), "ry": float(ry), "rz": float(rz),
                "lx": float(lx), "ly": float(ly), "lz": float(lz),
                "color": color, "label": f"{res_label} ({dist} Å)"
            })

    rec_json = json.dumps(receptor_data or "")
    parent_json = json.dumps(parent_ligand_data or "")
    var_json = json.dumps(var_ligand_data or "")
    pocket_json = json.dumps(pocket_pdb_str or "")
    inter_json = json.dumps(interactions_list)

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        * {{ box-sizing: border-box; user-select: none; }}
        body, html {{
            margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden;
            background: radial-gradient(circle at 50% 42%, #131B2E 0%, #090D16 100%);
            font-family: -apple-system, BlinkMacSystemFont, "SF Pro Display", "Inter", sans-serif;
        }}
        #viewer-wrapper {{
            width: 100%; height: {height}px; position: relative;
            border-radius: 14px; border: 1px solid rgba(255,255,255,0.14);
            overflow: hidden; box-shadow: 0 18px 42px rgba(0,0,0,0.55);
        }}
        #ctrl-bar {{
            position: absolute; top: 10px; left: 10px; right: 10px; z-index: 20;
            display: flex; flex-wrap: wrap; gap: 6px; align-items: center; justify-content: space-between;
            background: rgba(11, 15, 25, 0.88); backdrop-filter: blur(14px);
            padding: 7px 12px; border-radius: 10px; border: 1px solid rgba(255,255,255,0.12);
        }}
        .btn-group {{ display: flex; flex-wrap: wrap; gap: 5px; align-items: center; }}
        .c-btn {{
            background: rgba(255,255,255,0.06); color: #CBD5E1;
            border: 1px solid rgba(255,255,255,0.14); border-radius: 7px;
            padding: 4px 9px; font-size: 10.5px; font-weight: 600; cursor: pointer;
            transition: all 0.15s ease; display: inline-flex; align-items: center; gap: 4px;
        }}
        .c-btn:hover {{ background: rgba(30, 41, 59, 0.95); color: #FFF; border-color: #64D2FF; }}
        .c-btn.active {{
            background: rgba(10, 132, 255, 0.25); color: #64D2FF;
            border-color: rgba(100, 210, 255, 0.55);
        }}
        .c-btn.rec-btn {{
            background: rgba(255, 69, 58, 0.22); color: #FF6961;
            border-color: rgba(255, 69, 58, 0.5); font-weight: 700;
        }}
        .c-btn.rec-btn.recording {{
            background: #FF3B30; color: #FFF; animation: pulseRec 1s infinite;
        }}
        @keyframes pulseRec {{
            0% {{ box-shadow: 0 0 0 0 rgba(255,59,48,0.7); }}
            70% {{ box-shadow: 0 0 0 8px rgba(255,59,48,0); }}
            100% {{ box-shadow: 0 0 0 0 rgba(255,59,48,0); }}
        }}
        #legend {{
            position: absolute; bottom: 10px; left: 10px; right: 10px; z-index: 15;
            background: rgba(11, 15, 25, 0.9); backdrop-filter: blur(12px);
            padding: 6px 14px; border-radius: 9px; border: 1px solid rgba(255,255,255,0.12);
            font-size: 11px; color: #E2E8F0; display: flex; justify-content: space-between; align-items: center;
        }}
        .dot {{ display: inline-block; width: 9px; height: 9px; border-radius: 50%; margin-right: 5px; }}
    </style>
</head>
<body>
    <div id="viewer-wrapper">
        <div id="ctrl-bar">
            <div class="btn-group">
                <button class="c-btn active" id="btn-show-parent" title="Toggle Natural Parent Compound (Gold)">🟡 Parent</button>
                <button class="c-btn active" id="btn-show-var" title="Toggle Bioisosteric Derivative (Cyan)">🔵 Derivative</button>
                <button class="c-btn active" id="btn-pharm-clouds" title="Show 3D Steric/Electronic Pharmacophore Envelopes">🌌 Pharmacophore Clouds</button>
                <button class="c-btn active" id="btn-sidechains" title="Show Active-Site Side-Chains (<5.5 Å)">🧬 Pocket Sticks</button>
                <button class="c-btn" id="btn-show-surf" title="Show Receptor Cavity Surface">🧊 Cavity Wall</button>
            </div>
            <div class="btn-group">
                <button class="c-btn active" id="btn-spin" title="Continuous 360° Cinema Orbit">🌀 360° Orbit</button>
                <button class="c-btn" id="btn-reset" title="Zoom to Active Pocket">🎯 Pocket Focus</button>
                <button class="c-btn" id="btn-snap" title="Download 4K PNG">📸 Snapshot</button>
                <button class="c-btn rec-btn" id="btn-record" title="Record 60-FPS HD WebM Video">🎥 Record HD Video</button>
            </div>
        </div>
        <div id="{container_id}" style="width: 100%; height: 100%;"></div>
        <div id="legend">
            <div style="display:flex; gap:14px; align-items:center;">
                <span><span class="dot" style="background:#FFD60A;"></span><b>Parent:</b> {parent_name} (Gold Envelope)</span>
                <span><span class="dot" style="background:#64D2FF;"></span><b>Derivative:</b> {var_name} (Cyan Envelope)</span>
                <span><span class="dot" style="background:#FF3366;"></span>Polar H-Bond</span>
            </div>
            <span style="color:#64D2FF; font-family:monospace; font-weight:600;">3D Pharmacophore Superposition Studio</span>
        </div>
    </div>
    <script>
        (function() {{
            var receptorStr = {rec_json};
            var parentStr = {parent_json};
            var varStr = {var_json};
            var pocketStr = {pocket_json};
            var interactions = {inter_json};
            var viewer = null;
            var showParent = true;
            var showVar = true;
            var showPharmClouds = true;
            var showSidechains = true;
            var showSurface = false;
            var isSpinning = true;
            var parentSurf = null;
            var varSurf = null;
            var pocketSurf = null;
            var mediaRecorder = null;
            var recordedChunks = [];

            function updatePharmacophoreEnvelopes() {{
                if (parentSurf !== null) {{ try {{ viewer.removeSurface(parentSurf); }} catch(e) {{}} parentSurf = null; }}
                if (varSurf !== null) {{ try {{ viewer.removeSurface(varSurf); }} catch(e) {{}} varSurf = null; }}
                if (!showPharmClouds) return;
                if (showParent && viewer.getModel(1)) {{
                    parentSurf = viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.25, color: '#FFD60A'}}, {{model: 1}});
                }}
                if (showVar && viewer.getModel(2)) {{
                    varSurf = viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.30, color: '#00D2FF'}}, {{model: 2}});
                }}
            }}

            function updateCavityWall() {{
                if (pocketSurf !== null) {{ try {{ viewer.removeSurface(pocketSurf); }} catch(e) {{}} pocketSurf = null; }}
                if (showSurface) {{
                    var mIdx = (pocketStr && pocketStr.trim().length > 10) ? 3 : 0;
                    pocketSurf = viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.35, color: 'white'}}, {{model: mIdx}});
                }}
            }}

            var timer = setInterval(function() {{
                if (typeof $3Dmol !== 'undefined') {{
                    clearInterval(timer);
                    var element = document.getElementById("{container_id}");
                    viewer = $3Dmol.createViewer(element, {{defaultcolors: $3Dmol.rasmolElementColors, antialias: true}});
                    viewer.setBackgroundColor(0x0B0F19, 1.0);

                    if (receptorStr && receptorStr.trim().length > 0) {{
                        viewer.addModel(receptorStr, "pdb");
                        viewer.setStyle({{model: 0}}, {{cartoon: {{color: 'spectrum', opacity: 0.50}}}});
                    }}
                    if (parentStr && parentStr.trim().length > 0) {{
                        var pFmt = parentStr.indexOf("$$$$") !== -1 ? "sdf" : "pdb";
                        viewer.addModel(parentStr, pFmt);
                        viewer.setStyle({{model: 1}}, {{stick: {{colorscheme: 'goldCarbon', radius: 0.17}}}});
                    }}
                    if (varStr && varStr.trim().length > 0) {{
                        var vFmt = varStr.indexOf("$$$$") !== -1 ? "sdf" : "pdb";
                        viewer.addModel(varStr, vFmt);
                        viewer.setStyle({{model: 2}}, {{
                            stick: {{colorscheme: 'cyanCarbon', radius: 0.22}},
                            sphere: {{colorscheme: 'cyanCarbon', scale: 0.25}}
                        }});
                    }}
                    if (pocketStr && pocketStr.trim().length > 0) {{
                        viewer.addModel(pocketStr, "pdb");
                        viewer.setStyle({{model: 3}}, {{stick: {{colorscheme: 'whiteCarbon', radius: 0.12, opacity: 0.88}}}});
                    }}

                    for (var i = 0; i < interactions.length; i++) {{
                        var it = interactions[i];
                        viewer.addCylinder({{
                            start: {{x: it.rx, y: it.ry, z: it.rz}},
                            end: {{x: it.lx, y: it.ly, z: it.lz}},
                            radius: 0.08, dashed: true, color: it.color
                        }});
                        viewer.addLabel(it.label, {{
                            position: {{x: (it.rx+it.lx)/2, y: (it.ry+it.ly)/2, z: (it.rz+it.lz)/2}},
                            backgroundColor: 'rgba(11, 15, 25, 0.88)',
                            fontColor: '#FFF', fontSize: 10.5, inFront: true
                        }});
                    }}

                    updatePharmacophoreEnvelopes();
                    if (viewer.getModel(2)) {{ viewer.zoomTo({{model: 2}}); viewer.zoom(0.75); }}
                    else if (viewer.getModel(1)) {{ viewer.zoomTo({{model: 1}}); viewer.zoom(0.75); }}
                    else viewer.zoomTo();
                    viewer.render();

                    setInterval(function() {{
                        if (isSpinning && viewer) {{
                            viewer.rotate(0.45, "y");
                            viewer.render();
                        }}
                    }}, 40);

                    document.getElementById('btn-show-parent').onclick = function() {{
                        showParent = !showParent;
                        this.classList.toggle('active', showParent);
                        if (viewer.getModel(1)) {{
                            viewer.setStyle({{model: 1}}, showParent ? {{stick: {{colorscheme: 'goldCarbon', radius: 0.17}}}} : {{}});
                            updatePharmacophoreEnvelopes();
                            viewer.render();
                        }}
                    }};
                    document.getElementById('btn-show-var').onclick = function() {{
                        showVar = !showVar;
                        this.classList.toggle('active', showVar);
                        if (viewer.getModel(2)) {{
                            viewer.setStyle({{model: 2}}, showVar ? {{stick: {{colorscheme: 'cyanCarbon', radius: 0.22}}, sphere: {{colorscheme: 'cyanCarbon', scale: 0.25}}}} : {{}});
                            updatePharmacophoreEnvelopes();
                            viewer.render();
                        }}
                    }};
                    document.getElementById('btn-pharm-clouds').onclick = function() {{
                        showPharmClouds = !showPharmClouds;
                        this.classList.toggle('active', showPharmClouds);
                        updatePharmacophoreEnvelopes();
                        viewer.render();
                    }};
                    document.getElementById('btn-sidechains').onclick = function() {{
                        showSidechains = !showSidechains;
                        this.classList.toggle('active', showSidechains);
                        if (viewer.getModel(3)) {{
                            viewer.setStyle({{model: 3}}, showSidechains ? {{stick: {{colorscheme: 'whiteCarbon', radius: 0.12, opacity: 0.88}}}} : {{}});
                            viewer.render();
                        }}
                    }};
                    document.getElementById('btn-show-surf').onclick = function() {{
                        showSurface = !showSurface;
                        this.classList.toggle('active', showSurface);
                        updateCavityWall();
                        viewer.render();
                    }};
                    document.getElementById('btn-spin').onclick = function() {{
                        isSpinning = !isSpinning;
                        this.classList.toggle('active', isSpinning);
                    }};
                    document.getElementById('btn-reset').onclick = function() {{
                        if (viewer.getModel(2)) {{ viewer.zoomTo({{model: 2}}); viewer.zoom(0.75); }}
                        else viewer.zoomTo();
                        viewer.render();
                    }};
                    document.getElementById('btn-snap').onclick = function() {{
                        var canvas = document.querySelector('#{container_id} canvas');
                        if (!canvas) return;
                        var a = document.createElement('a');
                        a.href = canvas.toDataURL('image/png');
                        a.download = 'EthnoDock_DualPose_Pharmacophore_4K.png';
                        a.click();
                    }};
                    document.getElementById('btn-record').onclick = function() {{
                        var btn = this;
                        var canvas = document.querySelector('#{container_id} canvas');
                        if (!canvas || typeof canvas.captureStream !== 'function') return;
                        if (mediaRecorder && mediaRecorder.state === 'recording') {{ mediaRecorder.stop(); return; }}
                        recordedChunks = [];
                        isSpinning = true;
                        document.getElementById('btn-spin').classList.add('active');
                        var stream = canvas.captureStream(60);
                        var mime = MediaRecorder.isTypeSupported('video/webm;codecs=vp9') ? 'video/webm;codecs=vp9' : 'video/webm';
                        mediaRecorder = new MediaRecorder(stream, {{mimeType: mime, videoBitsPerSecond: 6000000}});
                        mediaRecorder.ondataavailable = function(e) {{ if (e.data && e.data.size > 0) recordedChunks.push(e.data); }};
                        mediaRecorder.onstop = function() {{
                            var blob = new Blob(recordedChunks, {{type: 'video/webm'}});
                            var url = URL.createObjectURL(blob);
                            var a = document.createElement('a');
                            a.href = url;
                            a.download = 'EthnoDock_Bioisostere_Pharmacophore_Cinema_60FPS.webm';
                            a.click();
                            btn.classList.remove('recording');
                            btn.innerText = '🎥 Record HD Video';
                        }};
                        btn.classList.add('recording');
                        btn.innerText = '⏺️ Recording (5s)...';
                        mediaRecorder.start();
                        setTimeout(function() {{
                            if (mediaRecorder && mediaRecorder.state === 'recording') mediaRecorder.stop();
                        }}, 5000);
                    }};
                }}
            }}, 80);
        }})();
    </script>
</body>
</html>"""
    return html_content



def generate_3d_conformer_analysis(smiles):
    """
    Generates an energy-minimized 3D conformer using RDKit ETKDGv3 and MMFF94/UFF forcefields,
    and extracts geometric and stereochemical descriptors.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    res = AllChem.EmbedMolecule(mol_h, params)
    if res != 0:
        res = AllChem.EmbedMolecule(mol_h)
    try:
        AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol_h, maxIters=500)
        except Exception:
            pass
            
    mol_block = Chem.MolToMolBlock(mol_h)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    try:
        vol = round(float(AllChem.ComputeMolVolume(mol_h)), 1)
    except Exception:
        vol = 0.0
    try:
        rg = round(float(Descriptors3D.RadiusOfGyration(mol_h)), 2)
    except Exception:
        rg = 0.0
    try:
        asph = round(float(Descriptors3D.Asphericity(mol_h)), 3)
    except Exception:
        asph = 0.0
        
    return {
        "mol_block": mol_block,
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "total_atoms": mol_h.GetNumAtoms(),
        "chiral_centers": chiral_centers,
        "chiral_count": len(chiral_centers),
        "volume_a3": vol,
        "radius_of_gyration": rg,
        "asphericity": asph
    }

def build_standalone_ligand_3d_html(container_id, mol_block, style='ball_and_stick', show_surface=False, auto_spin=True, height=200, bg_color='0x000000', colorscheme='cyanCarbon'):
    """
    Constructs an interactive 3D WebGL viewer using 3Dmol.js specifically for
    standalone chemical ligand conformer observation and stereochemical analysis.
    """
    mol_json = json.dumps(mol_block)
    auto_spin_js = "true" if auto_spin else "false"
    
    if style == 'stick':
        style_js = f"viewer.setStyle({{}}, {{stick: {{radius: 0.22, colorscheme: '{colorscheme}'}}}});"
    elif style == 'sphere':
        style_js = f"viewer.setStyle({{}}, {{sphere: {{scale: 0.85, colorscheme: '{colorscheme}'}}}});"
    else: # ball_and_stick
        style_js = f"viewer.setStyle({{}}, {{stick: {{radius: 0.16, colorscheme: '{colorscheme}'}}, sphere: {{scale: 0.28, colorscheme: '{colorscheme}'}}}});"
        
    surface_js = "viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.55, color: 'white'});" if show_surface else ""
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body, html {{ margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; background-color: #000000; }}
        #viewer-wrapper {{ width: 100%; height: {height}px; position: relative; border-radius: 10px; }}
        #ctrl-overlay {{ position: absolute; bottom: 8px; right: 8px; z-index: 10; display: flex; gap: 6px; }}
        .ctrl-btn {{ background: rgba(15, 23, 42, 0.85); color: #94A3B8; border: 1px solid rgba(255,255,255,0.15); border-radius: 6px; padding: 3px 8px; font-size: 10px; cursor: pointer; font-family: -apple-system, sans-serif; }}
        .ctrl-btn:hover {{ background: rgba(30, 41, 59, 0.95); color: #FFF; }}
    </style>
</head>
<body>
    <div id="viewer-wrapper">
        <div id="{container_id}" style="width: 100%; height: 100%;"></div>
        <div id="ctrl-overlay">
            <button class="ctrl-btn" onclick="toggleSpin()">Spin</button>
            <button class="ctrl-btn" onclick="resetView()">Reset</button>
        </div>
    </div>
    <script>
        (function() {{
            var molData = {mol_json};
            var isSpinning = {auto_spin_js};
            var viewerInstance = null;
            var initCount = 0;

            var timer = setInterval(function() {{
                initCount++;
                if (typeof $3Dmol !== 'undefined') {{
                    clearInterval(timer);
                    var element = document.getElementById("{container_id}");
                    var viewer = $3Dmol.createViewer(element, {{defaultcolors: $3Dmol.rasmolElementColors}});
                    viewerInstance = viewer;
                    viewer.setBackgroundColor({bg_color});

                    if (molData && molData.trim().length > 0) {{
                        viewer.addModel(molData, "mol");
                        {style_js}
                        {surface_js}
                    }}

                    viewer.zoomTo();
                    viewer.render();
                    if (isSpinning) {{
                        viewer.spin(true);
                    }}

                    window.toggleSpin = function() {{
                        isSpinning = !isSpinning;
                        viewerInstance.spin(isSpinning);
                    }};
                    window.resetView = function() {{
                        viewerInstance.zoomTo();
                        viewerInstance.render();
                    }};
                }} else if (initCount > 50) {{
                    clearInterval(timer);
                }}
            }}, 100);
        }})();
    </script>
</body>
</html>"""
    return html_content

def build_3d_molecular_kit_html(container_id, mol_block, chiral_indices=None, height=460, colorscheme='cyanCarbon'):
    """
    Constructs a comprehensive, client-side 3D Molecular Analysis Kit (Workbench) using 3Dmol.js.
    Provides live interactive tools:
    - Atom Click Inspector (Element, Index, Coordinates)
    - Interactive 3D Distance Measurement (Click 2 atoms to measure interatomic distance in Å)
    - Chiral Stereocenter Halo Highlighter
    - Van der Waals Molecular Surface Envelope
    - Dynamic Style Switcher (Ball & Stick, Licorice, CPK Space-Filling, Wireframe)
    - Auto-Rotation / Spin & Camera Reset
    """
    mol_json = json.dumps(mol_block)
    chiral_json = json.dumps(chiral_indices or [])

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body, html {{ margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; background: #07090E; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; }}
        #kit-wrapper {{ width: 100%; height: {height}px; position: relative; border-radius: 12px; border: 1px solid rgba(255,255,255,0.1); background: #000; overflow: hidden; }}
        #toolbar {{ position: absolute; top: 10px; left: 10px; right: 10px; z-index: 20; display: flex; flex-wrap: wrap; gap: 6px; align-items: center; background: rgba(15, 19, 28, 0.88); backdrop-filter: blur(12px); padding: 8px 12px; border-radius: 8px; border: 1px solid rgba(255,255,255,0.1); }}
        .kit-btn {{ background: rgba(255,255,255,0.06); color: #E2E8F0; border: 1px solid rgba(255,255,255,0.15); border-radius: 6px; padding: 4px 10px; font-size: 11px; font-weight: 500; cursor: pointer; transition: all 0.15s ease; }}
        .kit-btn:hover {{ background: rgba(255,255,255,0.18); color: #FFF; }}
        .kit-btn.active {{ background: #0A84FF; color: #FFF; border-color: #0A84FF; font-weight: 600; box-shadow: 0 0 10px rgba(10,132,255,0.4); }}
        .kit-btn.active-gold {{ background: #FFD60A; color: #000; border-color: #FFD60A; font-weight: 700; box-shadow: 0 0 10px rgba(255,214,10,0.4); }}
        #info-bar {{ position: absolute; bottom: 8px; left: 10px; right: 10px; z-index: 20; background: rgba(15, 19, 28, 0.88); backdrop-filter: blur(12px); padding: 7px 14px; border-radius: 8px; border: 1px solid rgba(255,255,255,0.1); font-size: 11px; color: #94A3B8; display: flex; justify-content: space-between; align-items: center; }}
        .badge-cyan {{ background: rgba(100,210,255,0.15); color: #64D2FF; padding: 2px 6px; border-radius: 4px; font-weight: 600; }}
    </style>
</head>
<body>
    <div id="kit-wrapper">
        <div id="toolbar">
            <span style="color:#64D2FF; font-weight:700; font-size:11px; margin-right:4px;">🧪 3D Molecular Kit:</span>
            <button class="kit-btn active" id="btn-bns" onclick="setStyle('bns')">Ball & Stick</button>
            <button class="kit-btn" id="btn-stick" onclick="setStyle('stick')">Licorice</button>
            <button class="kit-btn" id="btn-cpk" onclick="setStyle('cpk')">Space-Filling (CPK)</button>
            <button class="kit-btn" id="btn-wire" onclick="setStyle('wire')">Wireframe</button>
            <span style="color:rgba(255,255,255,0.2); margin: 0 2px;">|</span>
            <button class="kit-btn" id="btn-surf" onclick="toggleSurface()">VDW Surface</button>
            <button class="kit-btn" id="btn-chiral" onclick="toggleChiral()">Highlight Chiral Centers</button>
            <button class="kit-btn" id="btn-measure" onclick="toggleMeasure()">📏 Measure Distance (Å)</button>
            <span style="color:rgba(255,255,255,0.2); margin: 0 2px;">|</span>
            <button class="kit-btn" id="btn-spin" onclick="toggleSpin()">Auto-Spin</button>
            <button class="kit-btn" onclick="resetView()">Reset View</button>
        </div>
        <div id="{container_id}" style="width: 100%; height: 100%;"></div>
        <div id="info-bar">
            <span id="status-text">💡 <b>Click any atom</b> to inspect element & coordinates, or click <b>"Measure Distance"</b> to measure bond lengths.</span>
            <span class="badge-cyan">WebGL 3D Molecular Engine</span>
        </div>
    </div>

    <script>
        (function() {{
            var molData = {mol_json};
            var chiralIndices = {chiral_json};
            var viewer = null;
            var currentStyle = 'bns';
            var showSurface = false;
            var showChiral = false;
            var isSpinning = false;
            var measureMode = false;
            var selectedAtoms = [];
            var surfaceObj = null;
            var chiralSpheres = [];
            var measureObjects = [];
            var initCount = 0;

            var timer = setInterval(function() {{
                initCount++;
                if (typeof $3Dmol !== 'undefined') {{
                    clearInterval(timer);
                    var element = document.getElementById("{container_id}");
                    viewer = $3Dmol.createViewer(element, {{defaultcolors: $3Dmol.rasmolElementColors}});
                    viewer.setBackgroundColor(0x000000);

                    if (molData && molData.trim().length > 0) {{
                        viewer.addModel(molData, "mol");
                        applyCurrentStyle();
                    }}

                    viewer.zoomTo();
                    viewer.render();

                    // Clickable atom inspector & measurement tool
                    viewer.setClickable({{}}, true, function(atom, viewerInstance, event, container) {{
                        if (!atom) return;
                        var infoBox = document.getElementById("status-text");

                        if (measureMode) {{
                            selectedAtoms.push(atom);
                            if (selectedAtoms.length === 1) {{
                                infoBox.innerHTML = "🎯 <b>Atom 1 Selected:</b> <span style='color:#64D2FF;'>" + atom.elem + " #" + (atom.serial || atom.index) + "</span>. Click <b>second atom</b> to complete distance measurement.";
                            }} else if (selectedAtoms.length >= 2) {{
                                var a1 = selectedAtoms[0];
                                var a2 = selectedAtoms[1];
                                var dx = a1.x - a2.x, dy = a1.y - a2.y, dz = a1.z - a2.z;
                                var dist = Math.sqrt(dx*dx + dy*dy + dz*dz).toFixed(3);

                                var cyl = viewer.addCylinder({{
                                    start: {{x: a1.x, y: a1.y, z: a1.z}},
                                    end: {{x: a2.x, y: a2.y, z: a2.z}},
                                    radius: 0.1,
                                    dashed: true,
                                    color: '#FFD60A'
                                }});
                                var lbl = viewer.addLabel(dist + " Å", {{
                                    position: {{x: (a1.x + a2.x)/2, y: (a1.y + a2.y)/2, z: (a1.z + a2.z)/2}},
                                    backgroundColor: 'rgba(0,0,0,0.85)',
                                    fontColor: '#FFD60A',
                                    fontSize: 12,
                                    inFront: true
                                }});
                                measureObjects.push(cyl);
                                measureObjects.push(lbl);

                                infoBox.innerHTML = "📏 <b>Interatomic Distance:</b> <span style='color:#64D2FF;'>" + a1.elem + " #" + (a1.serial || a1.index) + "</span> ↔ <span style='color:#64D2FF;'>" + a2.elem + " #" + (a2.serial || a2.index) + "</span> = <b style='color:#FFD60A; font-size:12px;'>" + dist + " Å</b>";
                                selectedAtoms = [];
                                viewer.render();
                            }}
                        }} else {{
                            infoBox.innerHTML = "🔍 <b>Atom #" + (atom.serial || atom.index) + ":</b> <span style='color:#64D2FF; font-weight:700;'>" + atom.elem + "</span> • Coordinates: (" + atom.x.toFixed(2) + ", " + atom.y.toFixed(2) + ", " + atom.z.toFixed(2) + ") • Charge: " + (atom.charge || 0);
                            if (viewer._lastInspectLabel) viewer.removeLabel(viewer._lastInspectLabel);
                            viewer._lastInspectLabel = viewer.addLabel(atom.elem + (atom.serial || atom.index), {{
                                position: {{x: atom.x, y: atom.y, z: atom.z}},
                                backgroundColor: 'rgba(10,132,255,0.85)',
                                fontColor: '#FFF',
                                fontSize: 11,
                                inFront: true
                            }});
                            viewer.render();
                        }}
                    }});
                }} else if (initCount > 50) {{
                    clearInterval(timer);
                }}
            }}, 100);

            function applyCurrentStyle() {{
                if (!viewer) return;
                var model = viewer.getModel();
                if (!model) return;

                if (currentStyle === 'stick') {{
                    model.setStyle({{}}, {{stick: {{radius: 0.22, colorscheme: '{colorscheme}'}}}});
                }} else if (currentStyle === 'cpk') {{
                    model.setStyle({{}}, {{sphere: {{scale: 0.85, colorscheme: '{colorscheme}'}}}});
                }} else if (currentStyle === 'wire') {{
                    model.setStyle({{}}, {{line: {{linewidth: 2, colorscheme: '{colorscheme}'}}}});
                }} else {{ // bns
                    model.setStyle({{}}, {{
                        stick: {{radius: 0.16, colorscheme: '{colorscheme}'}},
                        sphere: {{scale: 0.28, colorscheme: '{colorscheme}'}}
                    }});
                }}
                viewer.render();
            }}

            window.setStyle = function(styleKey) {{
                currentStyle = styleKey;
                ['bns', 'stick', 'cpk', 'wire'].forEach(function(s) {{
                    var el = document.getElementById('btn-' + s);
                    if (el) el.classList.toggle('active', s === styleKey);
                }});
                applyCurrentStyle();
            }};

            window.toggleSurface = function() {{
                showSurface = !showSurface;
                var btn = document.getElementById('btn-surf');
                if (btn) btn.classList.toggle('active', showSurface);

                if (surfaceObj) {{
                    viewer.removeSurface(surfaceObj);
                    surfaceObj = null;
                }}
                if (showSurface) {{
                    surfaceObj = viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.55, color: 'white'}});
                }}
                viewer.render();
            }};

            window.toggleChiral = function() {{
                showChiral = !showChiral;
                var btn = document.getElementById('btn-chiral');
                if (btn) btn.classList.toggle('active-gold', showChiral);

                chiralSpheres.forEach(function(s) {{ viewer.removeShape(s); }});
                chiralSpheres = [];

                if (showChiral && chiralIndices && chiralIndices.length > 0) {{
                    var model = viewer.getModel();
                    var atoms = model.selectedAtoms({{}});
                    chiralIndices.forEach(function(idx) {{
                        var at = atoms.find(function(a) {{ return (a.index === idx || a.serial === (idx + 1)); }});
                        if (at) {{
                            var s = viewer.addSphere({{
                                center: {{x: at.x, y: at.y, z: at.z}},
                                radius: 0.55,
                                color: '#FFD60A',
                                opacity: 0.65
                            }});
                            chiralSpheres.push(s);
                        }}
                    }});
                    document.getElementById("status-text").innerHTML = "🌟 <b>Chiral Stereocenters Highlighted:</b> " + chiralIndices.length + " chiral carbon centers identified and marked with golden halos.";
                }}
                viewer.render();
            }};

            window.toggleMeasure = function() {{
                measureMode = !measureMode;
                selectedAtoms = [];
                var btn = document.getElementById('btn-measure');
                if (btn) btn.classList.toggle('active-gold', measureMode);

                if (measureMode) {{
                    document.getElementById("status-text").innerHTML = "📏 <b>Distance Measurement Active:</b> Click on any <b>first atom</b>, then click on a <b>second atom</b> to measure interatomic distance in Å.";
                }} else {{
                    measureObjects.forEach(function(obj) {{
                        try {{ viewer.removeShape(obj); }} catch(e) {{}}
                        try {{ viewer.removeLabel(obj); }} catch(e) {{}}
                    }});
                    measureObjects = [];
                    document.getElementById("status-text").innerHTML = "💡 Distance measurement cleared. Returned to standard 3D atom inspection mode.";
                    viewer.render();
                }}
            }};

            window.toggleSpin = function() {{
                isSpinning = !isSpinning;
                var btn = document.getElementById('btn-spin');
                if (btn) btn.classList.toggle('active', isSpinning);
                viewer.spin(isSpinning);
            }};

            window.resetView = function() {{
                viewer.zoomTo();
                viewer.render();
            }};
        }})();
    </script>
</body>
</html>"""
    return html_content

def build_redocking_superposition_3dmol_html(container_id, receptor_data, crystal_ligand_data, docked_ligand_data, rmsd_val, ligand_name="Experimental Co-Crystal", tier="EXEMPLARY", height=490):
    """
    Constructs an interactive 3D WebGL superposition viewer in 3Dmol.js
    comparing the native crystallographic co-crystal pose (cyan) directly against
    the blind AutoDock Vina redocked pose (gold) within the target receptor pocket.
    """
    rec_json = json.dumps(receptor_data)
    crys_json = json.dumps(crystal_ligand_data)
    dock_json = json.dumps(docked_ligand_data)
    
    badge_color = "#10B981" if rmsd_val <= 2.0 else "#F59E0B"

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body, html {{ margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; background-color: #0B0E14; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; }}
        #viewer-wrapper {{ width: 100%; height: {height}px; position: relative; border-radius: 12px; border: 1px solid rgba(255,255,255,0.12); }}
        #ctrl-bar {{ position: absolute; top: 10px; left: 10px; right: 10px; z-index: 20; display: flex; flex-wrap: wrap; gap: 6px; align-items: center; background: rgba(15, 23, 42, 0.92); backdrop-filter: blur(10px); padding: 7px 12px; border-radius: 8px; border: 1px solid rgba(255,255,255,0.12); }}
        .c-btn {{ background: rgba(255,255,255,0.08); color: #E2E8F0; border: 1px solid rgba(255,255,255,0.15); border-radius: 6px; padding: 3px 9px; font-size: 11px; font-weight: 500; cursor: pointer; transition: all 0.15s ease; }}
        .c-btn:hover {{ background: rgba(255,255,255,0.2); color: #FFF; }}
        .c-btn.active {{ background: #0A84FF; color: #FFF; border-color: #0A84FF; font-weight: 600; }}
        .c-btn.active-cyan {{ background: #00B4D8; color: #000; border-color: #00B4D8; font-weight: 700; }}
        .c-btn.active-gold {{ background: #FFD166; color: #000; border-color: #FFD166; font-weight: 700; }}
        #legend {{ position: absolute; bottom: 10px; left: 10px; z-index: 10; background: rgba(15, 23, 42, 0.92); backdrop-filter: blur(8px); padding: 8px 14px; border-radius: 8px; border: 1px solid rgba(255,255,255,0.15); font-size: 11px; color: #E2E8F0; display: flex; gap: 14px; align-items: center; }}
        .dot-cyan {{ display: inline-block; width: 10px; height: 10px; background: #00E5FF; border-radius: 50%; margin-right: 5px; box-shadow: 0 0 6px #00E5FF; }}
        .dot-gold {{ display: inline-block; width: 10px; height: 10px; background: #FFD166; border-radius: 50%; margin-right: 5px; box-shadow: 0 0 6px #FFD166; }}
        #hud-badge {{ position: absolute; top: 54px; right: 12px; z-index: 15; background: rgba(6, 78, 59, 0.85); border: 1px solid #10B981; border-radius: 6px; padding: 4px 10px; font-size: 11px; font-weight: 700; color: #A7F3D0; }}
    </style>
</head>
<body>
    <div id="viewer-wrapper">
        <div id="ctrl-bar">
            <span style="color:#00E5FF; font-weight:700; font-size:11px; margin-right:4px;">🔬 Co-Crystal Superposition:</span>
            <button class="c-btn active-cyan" id="btn-show-crys" onclick="toggleCrystal()">Native Crystal (Cyan)</button>
            <button class="c-btn active-gold" id="btn-show-dock" onclick="toggleDocked()">Redocked Vina (Gold)</button>
            <button class="c-btn active" id="btn-show-rec" onclick="toggleReceptor()">Receptor Pocket</button>
            <span style="color:rgba(255,255,255,0.2); margin: 0 2px;">|</span>
            <button class="c-btn" id="btn-spin" onclick="toggleSpin()">Auto-Spin</button>
            <button class="c-btn" onclick="resetView()">Reset View</button>
        </div>
        <div id="hud-badge" style="border-color:{badge_color}; color:{badge_color};">
            Heavy-Atom RMSD: {rmsd_val:.2f} Å &bull; {tier}
        </div>
        <div id="{container_id}" style="width: 100%; height: 100%;"></div>
        <div id="legend">
            <span><span class="dot-cyan"></span> <b>X-Ray Co-Crystal:</b> {ligand_name} (Experimental)</span>
            <span><span class="dot-gold"></span> <b>In-Silico Redocked:</b> Mode 1 Vina (Blind Search)</span>
        </div>
    </div>
    <script>
        (function() {{
            var receptorStr = {rec_json};
            var crystalStr = {crys_json};
            var dockedStr = {dock_json};
            var viewer = null;
            var showCrystal = true;
            var showDocked = true;
            var showReceptor = true;
            var isSpinning = false;
            var mRec = null, mCrys = null, mDock = null;
            var initCount = 0;

            var timer = setInterval(function() {{
                initCount++;
                if (typeof $3Dmol !== 'undefined') {{
                    clearInterval(timer);
                    var el = document.getElementById("{container_id}");
                    if (!el) return;
                    viewer = $3Dmol.createViewer(el, {{ backgroundColor: '#0B0E14' }});

                    // 1. Add Receptor
                    if (receptorStr) {{
                        mRec = viewer.addModel(receptorStr, "pdbqt");
                        mRec.setStyle({{}}, {{ cartoon: {{ color: '#64748B', opacity: 0.65 }} }});
                    }}

                    // 2. Add Native Crystal Ligand (Cyan)
                    if (crystalStr) {{
                        mCrys = viewer.addModel(crystalStr, "pdb");
                        mCrys.setStyle({{}}, {{
                            stick: {{ colorscheme: {{ prop: 'elem', map: {{ C: '#00E5FF', O: '#EF4444', N: '#3B82F6', S: '#EAB308', P: '#F97316' }} }}, radius: 0.22 }},
                            sphere: {{ colorscheme: {{ prop: 'elem', map: {{ C: '#00E5FF', O: '#EF4444', N: '#3B82F6', S: '#EAB308', P: '#F97316' }} }}, radius: 0.35 }}
                        }});
                    }}

                    // 3. Add Docked Ligand (Gold)
                    if (dockedStr) {{
                        mDock = viewer.addModel(dockedStr, "pdbqt");
                        mDock.setStyle({{}}, {{
                            stick: {{ colorscheme: {{ prop: 'elem', map: {{ C: '#FFD166', O: '#EF4444', N: '#3B82F6', S: '#EAB308', P: '#F97316' }} }}, radius: 0.18 }},
                            sphere: {{ colorscheme: {{ prop: 'elem', map: {{ C: '#FFD166', O: '#EF4444', N: '#3B82F6', S: '#EAB308', P: '#F97316' }} }}, radius: 0.28 }}
                        }});
                    }}

                    if (mCrys) {{
                        viewer.zoomTo({{ model: mCrys }});
                    }} else if (mDock) {{
                        viewer.zoomTo({{ model: mDock }});
                    }} else {{
                        viewer.zoomTo();
                    }}
                    viewer.render();
                }} else if (initCount > 50) {{
                    clearInterval(timer);
                }}
            }}, 100);

            window.toggleCrystal = function() {{
                showCrystal = !showCrystal;
                var btn = document.getElementById('btn-show-crys');
                if (btn) btn.classList.toggle('active-cyan', showCrystal);
                if (mCrys) {{
                    mCrys.setStyle({{}}, showCrystal ? {{
                        stick: {{ colorscheme: {{ prop: 'elem', map: {{ C: '#00E5FF', O: '#EF4444', N: '#3B82F6', S: '#EAB308', P: '#F97316' }} }}, radius: 0.22 }},
                        sphere: {{ colorscheme: {{ prop: 'elem', map: {{ C: '#00E5FF', O: '#EF4444', N: '#3B82F6', S: '#EAB308', P: '#F97316' }} }}, radius: 0.35 }}
                    }} : {{}});
                }}
                viewer.render();
            }};

            window.toggleDocked = function() {{
                showDocked = !showDocked;
                var btn = document.getElementById('btn-show-dock');
                if (btn) btn.classList.toggle('active-gold', showDocked);
                if (mDock) {{
                    mDock.setStyle({{}}, showDocked ? {{
                        stick: {{ colorscheme: {{ prop: 'elem', map: {{ C: '#FFD166', O: '#EF4444', N: '#3B82F6', S: '#EAB308', P: '#F97316' }} }}, radius: 0.18 }},
                        sphere: {{ colorscheme: {{ prop: 'elem', map: {{ C: '#FFD166', O: '#EF4444', N: '#3B82F6', S: '#EAB308', P: '#F97316' }} }}, radius: 0.28 }}
                    }} : {{}});
                }}
                viewer.render();
            }};

            window.toggleReceptor = function() {{
                showReceptor = !showReceptor;
                var btn = document.getElementById('btn-show-rec');
                if (btn) btn.classList.toggle('active', showReceptor);
                if (mRec) {{
                    mRec.setStyle({{}}, showReceptor ? {{ cartoon: {{ color: '#64748B', opacity: 0.65 }} }} : {{}});
                }}
                viewer.render();
            }};

            window.toggleSpin = function() {{
                isSpinning = !isSpinning;
                var btn = document.getElementById('btn-spin');
                if (btn) btn.classList.toggle('active', isSpinning);
                viewer.spin(isSpinning);
            }};

            window.resetView = function() {{
                if (mCrys) viewer.zoomTo({{ model: mCrys }});
                else viewer.zoomTo();
                viewer.render();
            }};
        }})();
    </script>
</body>
</html>"""
    return html_content



