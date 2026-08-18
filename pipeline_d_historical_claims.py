import os
import hashlib
from ethnodoc_models import init_db, get_session, Species, HistoricalSource, HistoricalClaim, Preparation

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def run_pipeline_d():
    print("[Pipeline D] Starting Historical Claims & Preparation Alignment")
    
    db_path = f"sqlite:///{os.path.join(BASE_DIR, 'ethnodoc.db')}"
    engine = init_db(db_path)
    session = get_session(engine)
    
    # We will map claims based on the species_source_association table.
    # Since the user mandated NO API scraping and NO invented/fake science,
    # we will populate the claims table precisely with what is verifiable locally.
    
    species_list = session.query(Species).all()
    print(f"[Pipeline D] Scanning {len(species_list)} species for linked historical sources...")
    
    claims_added = 0
    preps_added = 0
    
    for sp in species_list:
        for source in sp.historical_sources:
            # Generate a deterministic ID
            claim_id_str = f"{sp.species_id}_{source.source_id}_claim"
            claim_id = hashlib.md5(claim_id_str.encode()).hexdigest()[:12]
            
            # Check if exists
            existing_claim = session.query(HistoricalClaim).filter_by(claim_id=claim_id).first()
            if not existing_claim:
                new_claim = HistoricalClaim(
                    claim_id=claim_id,
                    species_id=sp.species_id,
                    source_id=source.source_id,
                    claim_text="[Awaiting Local Academic Ingestion]",
                    translation="Pending strictly verified translation from Kanripo/CText local corpus.",
                    disease_condition="Pending Evidence Classification"
                )
                session.add(new_claim)
                claims_added += 1
                
                # Link an empty strict preparation requirement
                prep_id = f"prep_{claim_id}"
                new_prep = Preparation(
                    preparation_id=prep_id,
                    claim_id=claim_id,
                    plant_part="Strictly Undefined (Awaiting Evidence)",
                    processing_method="Do not invent missing preparation details. Awaiting local manual ingest.",
                    provenance="EthnoDoc TCM System Enforcement"
                )
                session.add(new_prep)
                preps_added += 1

    session.commit()
    print(f"[Pipeline D] Structuring Complete.")
    print(f" -> Generated {claims_added} strict evidence placeholders linking Species to Sources.")
    print(f" -> Generated {preps_added} strictly enforced preparation blocks (No fake science).")

if __name__ == '__main__':
    run_pipeline_d()
