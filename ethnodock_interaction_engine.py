import math
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

def build_3dmol_html(container_id, receptor_data, ligand_data, interactions_df=None, receptor_style='cartoon', ligand_style='stick', show_surface=False, height=520):
    """
    Constructs an interactive 3D WebGL viewer using 3Dmol.js with custom color schemes,
    dashed interaction cylinders, and 3D residue annotations.
    """
    # Build JS lines for interaction cylinders and labels
    cylinders_js = []
    if interactions_df is not None and not interactions_df.empty:
        for _, row in interactions_df.iterrows():
            rx, ry, rz = row.get("Receptor XYZ", row.get("rec_xyz", [0, 0, 0]))
            lx, ly, lz = row.get("Ligand XYZ", row.get("lig_xyz", [0, 0, 0]))
            color = row.get("Color", row.get("color", "#FF3366"))
            res_label = row["Receptor Residue"]
            dist = row["Distance (Å)"]
            
            # Dashed cylinder
            cylinders_js.append(
                f"viewer.addCylinder({{start:{{x:{rx}, y:{ry}, z:{rz}}}, end:{{x:{lx}, y:{ly}, z:{lz}}}, radius:0.08, dashed:true, color:'{color}'}});"
            )
            # Label
            cylinders_js.append(
                f"viewer.addLabel('{res_label} ({dist}Å)', {{position: {{x:{rx}, y:{ry}, z:{rz}}}, backgroundColor: 'rgba(15, 23, 42, 0.85)', fontColor: 'white', fontSize: 11, inFront: true}});"
            )

    cylinders_script = "\n                            ".join(cylinders_js)
    surface_js = "viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.65, color: 'white'}, {model: 0});" if show_surface else ""

    # Sanitize data for template
    rec_json = json.dumps(receptor_data)
    lig_json = json.dumps(ligand_data)

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        <style>
            body, html {{ margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; background-color: #0E1117; }}
            #viewer-wrapper {{ width: 100%; height: {height}px; position: relative; border-radius: 10px; border: 1px solid #2D3748; }}
            #legend {{ position: absolute; bottom: 10px; left: 10px; z-index: 10; background: rgba(14, 17, 23, 0.85); padding: 6px 12px; border-radius: 6px; border: 1px solid #4A5568; font-family: monospace; font-size: 11px; color: #E2E8F0; }}
            .dot-red {{ display: inline-block; width: 8px; height: 8px; background: #FF3366; border-radius: 50%; margin-right: 4px; }}
            .dot-blue {{ display: inline-block; width: 8px; height: 8px; background: #00D2FF; border-radius: 50%; margin-right: 4px; margin-left: 8px; }}
        </style>
    </head>
    <body>
        <div id="viewer-wrapper">
            <div id="{container_id}" style="width: 100%; height: 100%;"></div>
            <div id="legend">
                <span class="dot-red"></span> Polar / H-Bond
                <span class="dot-blue"></span> Hydrophobic
            </div>
        </div>
        <script>
            (function() {{
                var receptorStr = {rec_json};
                var ligandStr = {lig_json};
                var initCount = 0;

                var timer = setInterval(function() {{
                    initCount++;
                    if (typeof $3Dmol !== 'undefined') {{
                        clearInterval(timer);
                        var element = document.getElementById("{container_id}");
                        var viewer = $3Dmol.createViewer(element, {{defaultcolors: $3Dmol.rasmolElementColors}});
                        viewer.setBackgroundColor(0x0E1117);

                        // 1. Add Receptor Model (Model 0)
                        if (receptorStr && receptorStr.trim().length > 0) {{
                            viewer.addModel(receptorStr, "pdb");
                            viewer.setStyle({{model: 0}}, {{{receptor_style}: {{color: 'spectrum'}} }});
                            {surface_js}
                        }}

                        // 2. Add Ligand Model (Model 1)
                        if (ligandStr && ligandStr.trim().length > 0) {{
                            var format = ligandStr.indexOf("$$$$") !== -1 ? "sdf" : "pdb";
                            viewer.addModel(ligandStr, format);
                            viewer.setStyle({{model: 1}}, {{{ligand_style}: {{colorscheme: 'greenCarbon'}} }});
                        }}

                        // 3. Add Dashed Interactions & Residue Labels
                        try {{
                            {cylinders_script}
                        }} catch (err) {{
                            console.error("Interaction rendering note:", err);
                        }}

                        viewer.zoomTo();
                        viewer.render();
                    }} else if (initCount > 50) {{
                        clearInterval(timer);
                    }}
                }}, 100);
            }})();
        </script>
    </body>
    </html>
    """
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


