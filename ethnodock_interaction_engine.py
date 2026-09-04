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

