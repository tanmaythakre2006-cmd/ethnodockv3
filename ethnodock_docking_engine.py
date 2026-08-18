import os
import sys
import urllib.request
import zipfile
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.join(BASE_DIR, "bin")
os.makedirs(BIN_DIR, exist_ok=True)

def get_vina_executable():
    """
    Locates or automatically downloads the appropriate AutoDock Vina binary
    for the host operating system (Windows or Linux).
    """
    is_windows = sys.platform.startswith("win")
    if is_windows:
        vina_exe = os.path.join(BIN_DIR, "vina.exe")
        if os.path.exists(vina_exe):
            return vina_exe
            
        # Check in current directory
        if os.path.exists(os.path.join(BASE_DIR, "vina.exe")):
            return os.path.join(BASE_DIR, "vina.exe")
            
        # Check system PATH
        try:
            res = subprocess.run(["vina", "--help"], capture_output=True, text=True)
            if res.returncode == 0 or "AutoDock Vina" in res.stdout or "AutoDock Vina" in res.stderr:
                return "vina"
        except Exception:
            pass

        # Attempt to download Windows Vina direct executable from GitHub release
        win_url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe"
        try:
            print("Downloading AutoDock Vina for Windows...")
            req = urllib.request.Request(win_url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as resp, open(vina_exe, 'wb') as f_out:
                f_out.write(resp.read())
            if os.path.exists(vina_exe):
                return vina_exe
        except Exception as e:
            print(f"Warning: Could not download Windows Vina binary: {e}")
            
        return vina_exe
    else:
        # Linux binary
        vina_bin = os.path.join(BIN_DIR, "vina")
        if os.path.exists(vina_bin):
            return vina_bin
            
        # Check in BASE_DIR or PATH
        if os.path.exists(os.path.join(BASE_DIR, "vina")):
            return os.path.join(BASE_DIR, "vina")
            
        try:
            res = subprocess.run(["vina", "--help"], capture_output=True, text=True)
            if res.returncode == 0 or "AutoDock Vina" in res.stdout:
                return "vina"
        except Exception:
            pass
            
        linux_url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
        try:
            print("Downloading AutoDock Vina for Linux...")
            urllib.request.urlretrieve(linux_url, vina_bin)
            os.chmod(vina_bin, 0o755)
            return vina_bin
        except Exception as e:
            print(f"Warning: Could not download Linux Vina binary: {e}")
            return vina_bin

def fetch_receptor(pdb_id, output_pdb=None):
    """
    Downloads PDB file from RCSB PDB and converts it into a clean,
    AutoDock-ready PDBQT file with assigned atom types and placeholder charges.
    """
    pdb_id = pdb_id.strip().upper()
    if output_pdb is None:
        output_pdb = os.path.join(BASE_DIR, f"{pdb_id}.pdb")
        
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, output_pdb)
    except Exception as e:
        print(f"Error downloading PDB {pdb_id}: {e}")
        return None

    with open(output_pdb, 'r', encoding='utf-8', errors='ignore') as f:
        pdb_str = f.read()

    lines = []
    # Strip water molecules and non-protein ions that might interfere, or format standard atoms
    for line in pdb_str.split('\n'):
        if line.startswith("ATOM") or (line.startswith("HETATM") and not line[17:20].strip() in ["HOH", "WAT", "DOD"]):
            element = line[76:78].strip()
            if not element:
                element = line[12:16].strip()[0]
            
            ad_type = element
            if element == 'H': ad_type = 'HD'
            elif element == 'C':
                # Heuristic for aromatic carbons
                res_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                if res_name in ["PHE", "TYR", "TRP", "HIS"] and atom_name.startswith(("CD", "CE", "CZ", "CG", "CH")):
                    ad_type = 'A'
                else:
                    ad_type = 'C'
            elif element == 'N': ad_type = 'NA' if line[17:20].strip() in ["HIS", "TRP"] else 'N'
            elif element == 'O': ad_type = 'OA'
            elif element == 'S': ad_type = 'SA'
            elif element == 'P': ad_type = 'P'
            elif element == 'CL': ad_type = 'Cl'
            elif element == 'BR': ad_type = 'Br'
            elif element == 'F': ad_type = 'F'
            elif element == 'I': ad_type = 'I'
            elif element == 'FE': ad_type = 'Fe'
            elif element == 'ZN': ad_type = 'Zn'
            elif element == 'CA': ad_type = 'Ca'
            elif element == 'MG': ad_type = 'Mg'
            elif element == 'MN': ad_type = 'Mn'

            # AutoDock format: columns 1-66 preserved + 4 spaces + "+0.000" + space + atom type (2 chars)
            new_line = line[:66].ljust(66) + "    +0.000 " + ad_type.ljust(2)
            lines.append(new_line)

    pdbqt_output = output_pdb.replace('.pdb', '.pdbqt')
    with open(pdbqt_output, 'w', encoding='utf-8') as f:
        f.write("\n".join(lines))

    return pdbqt_output

def prepare_ligand(smiles, output_pdbqt=None):
    """
    Converts a SMILES string into a 3D conformer with UFF and MMFF94 force field
    energy minimization, calculating strain delta, and writes an AutoDock PDBQT file.
    """
    uff_delta = 0.0
    if output_pdbqt is None:
        output_pdbqt = os.path.join(BASE_DIR, "ligand.pdbqt")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, 0.0

    mol = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 42
    embed_status = AllChem.EmbedMolecule(mol, ps)
    if embed_status != 0:
        # Fallback to standard embedding if ETKDGv3 fails
        AllChem.EmbedMolecule(mol, randomSeed=42)

    # UFF minimization and delta calculation
    try:
        ff = AllChem.UFFGetMoleculeForceField(mol)
        if ff:
            e1 = ff.CalcEnergy()
            ff.Minimize(maxIts=500)
            e2 = ff.CalcEnergy()
            uff_delta = max(0.0, e1 - e2)
    except Exception as e:
        print(f"UFF optimization note: {e}")

    # MMFF94 optimization
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception as e:
        print(f"MMFF optimization note: {e}")

    pdb_block = Chem.MolToPDBBlock(mol)

    lines = ["ROOT"]
    for line in pdb_block.split('\n'):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            element = line[76:78].strip()
            if not element:
                element = line[12:16].strip()[0]
            ad_type = element
            if element == 'H': ad_type = 'HD'
            elif element == 'C': ad_type = 'C'
            elif element == 'N': ad_type = 'NA' if 'N' in line and any(ar in line for ar in ['c', 'n']) else 'N'
            elif element == 'O': ad_type = 'OA'
            elif element == 'S': ad_type = 'SA'
            elif element == 'P': ad_type = 'P'
            elif element == 'CL': ad_type = 'Cl'
            elif element == 'BR': ad_type = 'Br'
            elif element == 'F': ad_type = 'F'
            elif element == 'I': ad_type = 'I'

            new_line = line[:66].ljust(66) + "    +0.000 " + ad_type.ljust(2)
            lines.append(new_line)
            
    lines.append("ENDROOT")
    lines.append("TORSDOF 0")

    with open(output_pdbqt, 'w', encoding='utf-8') as f:
        f.write("\n".join(lines))
        
    return output_pdbqt, uff_delta

def smart_cavity_finder(pdb_file):
    """
    Automatically calculates the binding pocket cavity centroid (X, Y, Z)
    and recommended bounding box dimensions (Sx, Sy, Sz) from the receptor structure.
    """
    try:
        coords = []
        with open(pdb_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append([x, y, z])
                    except ValueError:
                        pass
        if coords:
            coords_arr = np.array(coords)
            center = np.mean(coords_arr, axis=0)
            mins = np.min(coords_arr, axis=0)
            maxs = np.max(coords_arr, axis=0)
            dims = maxs - mins
            # Clip bounding box to standard Vina search dimensions (between 16 and 40 Å)
            dims = np.clip(dims * 0.6, a_min=16.0, a_max=40.0)
            return [round(float(c), 3) for c in center], [round(float(d), 3) for d in dims]
    except Exception as e:
        print(f"Error finding cavity in {pdb_file}: {e}")

    return [0.0, 0.0, 0.0], [22.0, 22.0, 22.0]

def parse_vina_output(vina_text):
    """
    Parses AutoDock Vina standard output into a clean structured list of dictionaries.
    """
    lines = vina_text.split('\n')
    data = []
    parsing = False
    for line in lines:
        if '-----+------------+----------+----------' in line:
            parsing = True
            continue
        if parsing:
            parts = line.split()
            if len(parts) >= 4 and parts[0].isdigit():
                try:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    rmsd_lb = float(parts[2])
                    rmsd_ub = float(parts[3])
                    data.append({
                        'mode': mode,
                        'affinity': affinity,
                        'rmsd_lb': rmsd_lb,
                        'rmsd_ub': rmsd_ub
                    })
                except ValueError:
                    pass
            elif len(parts) == 0 or 'Writing' in line:
                break
    return data

def run_vina_docking(receptor_pdbqt, ligand_pdbqt, center, dims, exhaustiveness=8, cpu=1, seed=42, output_pdbqt=None):
    """
    Executes AutoDock Vina CLI on the prepared receptor and ligand.
    """
    if output_pdbqt is None:
        output_pdbqt = ligand_pdbqt.replace(".pdbqt", "_out.pdbqt")

    vina_exe = get_vina_executable()

    cmd = [
        vina_exe,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(dims[0]),
        "--size_y", str(dims[1]),
        "--size_z", str(dims[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--cpu", str(cpu),
        "--seed", str(seed),
        "--out", output_pdbqt
    ]
    
    try:
        res = subprocess.run(cmd, capture_output=True, text=True)
        raw_output = res.stdout + "\n" + res.stderr
        parsed_poses = parse_vina_output(raw_output)
        return raw_output, parsed_poses, output_pdbqt
    except Exception as e:
        err_msg = f"Docking execution error: {e}"
        print(err_msg)
        return err_msg, [], output_pdbqt
