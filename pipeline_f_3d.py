import os
from ethnodoc_models import init_db, get_session, Phytochemical
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("[Pipeline F] Error: RDKit is required for offline 3D coordinate generation.")
    exit(1)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
STRUCT_DIR = os.path.join(BASE_DIR, "assets", "structures")
os.makedirs(STRUCT_DIR, exist_ok=True)

def run_pipeline_f():
    print("[Pipeline F] Starting API-Free 3D Coordinate Generation")
    
    db_path = f"sqlite:///{os.path.join(BASE_DIR, 'ethnodoc.db')}"
    engine = init_db(db_path)
    session = get_session(engine)
    
    compounds = session.query(Phytochemical).all()
    print(f"[Pipeline F] Found {len(compounds)} compounds to process.")
    
    success_count = 0
    
    for c in compounds:
        if not c.smiles:
            print(f" -> Skipping {c.compound_name}: No SMILES string.")
            continue
            
        print(f" -> Generating 3D conformer for: {c.compound_name}...")
        
        # Load molecule from SMILES
        mol = Chem.MolFromSmiles(c.smiles)
        if mol is None:
            print(f"    [!] Error parsing SMILES.")
            continue
            
        # Add hydrogens, crucial for 3D geometry
        mol = Chem.AddHs(mol)
        
        # Embed 3D coordinates using ETKDG (distance geometry)
        try:
            res = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            if res != 0:
                print(f"    [!] Failed to embed molecule.")
                continue
                
            # Optimize the geometry using MMFF94 force field
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Save as SDF
            safe_name = c.compound_name.replace(" ", "_").lower()
            filepath = os.path.join(STRUCT_DIR, f"{safe_name}.sdf")
            writer = Chem.SDWriter(filepath)
            writer.write(mol)
            writer.close()
            
            print(f"    [+] Saved 3D SDF: {filepath}")
            success_count += 1
            
        except Exception as e:
            print(f"    [!] Exception during 3D generation: {str(e)}")

    print(f"\n[Pipeline F] Complete. Generated {success_count} 3D structures offline.")

if __name__ == '__main__':
    run_pipeline_f()
