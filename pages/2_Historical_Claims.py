import streamlit as st
import os
import sys
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# --- Config ---
st.set_page_config(page_title="Historical Claims | EthnoDoc TCM", page_icon="📜", layout="wide")

st.markdown("""
<style>
    .stApp { background-color: #0E1117; color: #FAFAFA; font-family: 'Inter', sans-serif; }
    .claim-header { font-size: 2rem; font-weight: 300; color: #E91E63; border-bottom: 1px solid #333; padding-bottom: 0.5rem; margin-bottom: 1.5rem; }
    .claim-box { background: #1E2127; padding: 1.5rem; border-radius: 8px; border-left: 4px solid #E91E63; margin-bottom: 1rem; }
    .prep-box { background: #2D3748; padding: 1rem; border-radius: 4px; margin-top: 1rem; border: 1px dashed #4A5568; }
    .provenance-tag { display: inline-block; background: #E91E63; color: white; padding: 2px 8px; border-radius: 12px; font-size: 0.8rem; font-weight: bold; margin-bottom: 0.5rem; }
</style>
""", unsafe_allow_html=True)

# Add base directory to path so we can import ethnodoc_models
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from ethnodoc_models import Species, HistoricalClaim, Preparation

# --- Database ---
@st.cache_resource
def get_session():
    db_path = f"sqlite:///{os.path.join(BASE_DIR, 'ethnodoc.db')}"
    engine = create_engine(db_path, echo=False)
    Session = sessionmaker(bind=engine)
    return Session()

session = get_session()

# --- Main Logic ---
if 'selected_species_id' not in st.session_state:
    st.warning("Please search and select a species from the Landing Page first.")
    st.stop()
    
species_id = st.session_state['selected_species_id']
species = session.query(Species).filter_by(species_id=species_id).first()

st.markdown(f"<div class='claim-header'>Historical Claims: {species.english_name}</div>", unsafe_allow_html=True)
st.caption("EthnoDoc strictly enforces provenance. Unverified claims and undefined preparations are explicitly marked rather than computationally invented.")

claims = session.query(HistoricalClaim).filter_by(species_id=species_id).all()

if not claims:
    st.info("No historical claims structured for this species yet.")
else:
    for claim in claims:
        source_title = claim.source.title if claim.source else "Unknown Source"
        st.markdown(f"""
        <div class='claim-box'>
            <div class='provenance-tag'>{source_title}</div>
            <p><strong>Disease / Condition:</strong> <em>{claim.disease_condition}</em></p>
            <p><strong>Original Text:</strong> {claim.claim_text}</p>
            <p><strong>Translation:</strong> {claim.translation}</p>
        """, unsafe_allow_html=True)
        
        # Display Preparations linked to this claim
        preps = session.query(Preparation).filter_by(claim_id=claim.claim_id).all()
        if preps:
            st.markdown("<strong>Preparation Methodology:</strong>", unsafe_allow_html=True)
            for prep in preps:
                st.markdown(f"""
                <div class='prep-box'>
                    <span><strong>Plant Part:</strong> {prep.plant_part}</span><br/>
                    <span><strong>Method:</strong> {prep.processing_method}</span><br/>
                    <span style='color: #A0AEC0; font-size: 0.85rem;'>Provenance: {prep.provenance}</span>
                </div>
                """, unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)

st.markdown("---")
st.page_link("pages/3_Phytochemistry.py", label="View Molecular Profiles & Bioactives", icon="🧪")

