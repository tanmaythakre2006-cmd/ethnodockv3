import os
import json
from ethnodoc_models import init_db, get_session, Species, HistoricalSource, Phytochemical

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OKF_INDEX_FILE = os.path.join(BASE_DIR, "okf_index.json")

def load_okf_data():
    if not os.path.exists(OKF_INDEX_FILE):
        print("Error: OKF Index not found.")
        return {}
    with open(OKF_INDEX_FILE, "r", encoding="utf-8") as f:
        return json.load(f)

def run_pipeline_a():
    print("[Pipeline A] Starting Data Ingestion: Flat OKF to Relational SQLite")
    
    # Init DB
    db_path = f"sqlite:///{os.path.join(BASE_DIR, 'ethnodoc.db')}"
    engine = init_db(db_path)
    session = get_session(engine)
    
    # Load old flat data
    index_data = load_okf_data()
    entities = {k: v for k, v in index_data.items() if isinstance(v, dict)}
    
    print(f"[Pipeline A] Found {len(entities)} OKF entities to process.")
    
    # Tracks for relationships
    added_sources = {}
    
    species_count = 0
    source_count = 0
    compound_count = 0
    
    # 1. Process Texts (Sources)
    for eid, data in entities.items():
        if data.get("type") == "classical_text":
            source_id = data.get("entity_id")
            if source_id not in added_sources:
                title = data.get("title", "Unknown Title")
                
                new_source = HistoricalSource(
                    source_id=source_id,
                    title=title,
                    provenance="OKF Legacy Ingestion"
                )
                session.add(new_source)
                added_sources[source_id] = new_source
                added_sources[title] = new_source # Also map by title for crossref
                source_count += 1

    session.commit()
    
    # 2. Process Species & link to Sources
    for eid, data in entities.items():
        if data.get("type") == "ethnobotanical_species":
            tax = data.get("taxonomy", {})
            english = data.get("title", "").split("(")[0].strip()
            
            new_species = Species(
                species_id=eid,
                english_name=english,
                scientific_name=tax.get("binomial_name", "Unknown"),
                chinese_name=tax.get("chinese_name", ""),
                taxonomy_genus=tax.get("genus", "")
            )
            
            # Map legacy sources
            sources = data.get("sources", [])
            for src in sources:
                if src in added_sources:
                    new_species.historical_sources.append(added_sources[src])
                else:
                    # Create generic source if missing
                    fallback_source_id = f"src_{hash(src)}"
                    if fallback_source_id not in added_sources:
                        s = HistoricalSource(source_id=fallback_source_id, title=src, provenance="Inferred from OKF")
                        session.add(s)
                        added_sources[fallback_source_id] = s
                        added_sources[src] = s
                    new_species.historical_sources.append(added_sources[src])
                    
            session.add(new_species)
            species_count += 1
            
    # 3. Process Phytochemicals
    for eid, data in entities.items():
        if data.get("type") == "phytochemical_compound":
            title = data.get("title", "").split("(")[0].strip()
            new_comp = Phytochemical(
                compound_id=eid,
                compound_name=title,
                provenance="OKF Legacy Ingestion"
            )
            session.add(new_comp)
            compound_count += 1
            
    session.commit()
    
    print(f"[Pipeline A] Ingestion Complete.")
    print(f" - Species Migrated: {species_count}")
    print(f" - Sources Migrated: {len(added_sources)}")
    print(f" - Compounds Migrated: {compound_count}")
    
if __name__ == '__main__':
    run_pipeline_a()
