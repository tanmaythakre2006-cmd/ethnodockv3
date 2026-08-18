import os
import time
from datetime import datetime
from ethnodoc_models import init_db, get_session, Species, SpeciesImage

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ASSETS_DIR = os.path.join(BASE_DIR, "assets", "images")
os.makedirs(ASSETS_DIR, exist_ok=True)

def run_pipeline_c_local():
    print("[Pipeline C] Starting STRICTLY LOCAL Image Acquisition (API-Free Mode)")
    
    db_path = f"sqlite:///{os.path.join(BASE_DIR, 'ethnodoc.db')}"
    engine = init_db(db_path)
    session = get_session(engine)
    
    # Get all current species
    species_list = session.query(Species).all()
    species_map = {sp.species_id: sp for sp in species_list}
    
    print(f"[Pipeline C] Scanning local directory: {ASSETS_DIR}")
    
    # In API-Free mode, the user drops images into assets/images/ 
    # named by species ID (e.g. 'entity_42.jpg' or 'entity_42_alt.png')
    # Or, we just catalog whatever is already there.
    
    files = os.listdir(ASSETS_DIR)
    added_count = 0
    
    for filename in files:
        if not (filename.endswith(".jpg") or filename.endswith(".jpeg") or filename.endswith(".png")):
            continue
            
        filepath = os.path.join(ASSETS_DIR, filename)
        
        # Check if already in DB
        existing = session.query(SpeciesImage).filter_by(filename=filename).first()
        if existing:
            continue
            
        # Infer species ID from filename if possible. 
        # Expected format: "speciesid_timestamp.jpg" or "speciesid.jpg"
        # Since Wikipedia pipeline saved them as: sp.species_id + "_" + timestamp
        parts = filename.split('_')
        possible_sp_id = parts[0]
        if "entity_" in filename:
            # OKF entities look like 'entity_0'
            idx = filename.find('entity_')
            end_idx = filename.find('_', idx + 7)
            if end_idx == -1: end_idx = filename.find('.', idx)
            possible_sp_id = filename[idx:end_idx]
            
        # Validate that the species exists
        if possible_sp_id in species_map:
            new_image = SpeciesImage(
                image_id=f"img_local_{filename}",
                species_id=possible_sp_id,
                filename=filename,
                source_name="Local Asset (API-Free)",
                source_url="file://" + filepath,
                license_info="Local Database",
                acquisition_date=datetime.now().isoformat()
            )
            session.add(new_image)
            added_count += 1
            print(f"  [+] Cataloged Local Image: {filename} -> {possible_sp_id}")
            
    session.commit()
    print(f"\n[Pipeline C] Complete. {added_count} local images mapped API-free.")

if __name__ == '__main__':
    run_pipeline_c_local()
