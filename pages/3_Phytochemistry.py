import streamlit as st
import os
import sys
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import py3Dmol
from stmol import showmol

# --- Config ---
st.set_page_config(page_title="Phytochemistry | EthnoDoc TCM", page_icon="🧪", layout="wide")

st.markdown("""
<style>
    .stApp { background-color: #0E1117; color: #FAFAFA; font-family: 'Inter', sans-serif; }
    .page-title { font-size: 2.5rem; font-weight: 300; color: #00BCD4; border-bottom: 1px solid #333; padding-bottom: 0.5rem; margin-bottom: 1.5rem; }
    .mol-card { background: #1a1c23; padding: 1.5rem; border-radius: 8px; border-left: 4px solid #00BCD4; margin-bottom: 1.5rem; }
    
    /* Strict Evidence Tags */
    .tag-putative { background: #FF9800; color: white; padding: 3px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: bold; margin-bottom: 0.5rem; display: inline-block; }
    .tag-validated { background: #4CAF50; color: white; padding: 3px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: bold; margin-bottom: 0.5rem; display: inline-block; }
    .tag-reported { background: #607D8B; color: white; padding: 3px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: bold; margin-bottom: 0.5rem; display: inline-block; }
    
    .smiles-box { font-family: monospace; background: #0E1117; padding: 0.5rem; border-radius: 4px; color: #A0AEC0; word-wrap: break-word; }
    .metrics { display: flex; gap: 2rem; margin-top: 1rem; }
    .metric-label { font-size: 0.8rem; color: #888; }
    .metric-val { font-size: 1.1rem; font-weight: bold; color: #E2E8F0; }
</style>
""", unsafe_allow_html=True)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from ethnodoc_models import HistoricalClaim, Phytochemical, claim_compound_association

@st.cache_resource
def get_session():
    db_path = f"sqlite:///{os.path.join(BASE_DIR, 'ethnodoc.db')}"
    engine = create_engine(db_path, echo=False)
    Session = sessionmaker(bind=engine)
    return Session()

session = get_session()

# We need the selected claim ID from session state if coming from Historical Claims page,
# OR we can just show all phytochemicals for the selected species if coming from Species Profile.

if 'selected_species_id' not in st.session_state:
    st.warning("Please search and select a species from the Landing Page first.")
    st.stop()
    
species_id = st.session_state['selected_species_id']

st.markdown(f"<div class='page-title'>Phytochemical & Molecular Profiles</div>", unsafe_allow_html=True)
st.caption("EthnoDoc TCM strict evidence compliance active: Distinguishing between reported constituents, putative actives, and experimentally validated compounds.")

# Find all claims for this species, then find all compounds linked to those claims
claims = session.query(HistoricalClaim).filter_by(species_id=species_id).all()
claim_ids = [c.claim_id for c in claims]

if not claim_ids:
    st.info("No verified historical claims found to map phytochemicals against.")
    st.stop()

# Query the association table to get the compounds and their evidence types
stmt = session.query(
    Phytochemical, 
    claim_compound_association.c.evidence_type,
    claim_compound_association.c.evidence_strength
).join(
    claim_compound_association, Phytochemical.compound_id == claim_compound_association.c.compound_id
).filter(
    claim_compound_association.c.claim_id.in_(claim_ids)
).all()

# Group by compound to avoid duplicates if linked to multiple claims
compound_map = {}
for comp, ev_type, ev_strength in stmt:
    if comp.compound_id not in compound_map:
        compound_map[comp.compound_id] = {
            'compound': comp,
            'evidence_type': ev_type,
            'evidence_strength': ev_strength
        }
        
if not compound_map:
    st.info("No offline phytochemical data has been mathematically mapped to this species' claims yet. (See Pipeline E)")
else:
    comps_list = list(compound_map.values())
    tabs = st.tabs([d['compound'].compound_name for d in comps_list])
    
    for idx, data in enumerate(comps_list):
        with tabs[idx]:
            c = data['compound']
            ev_type = data['evidence_type'] or "Reported constituent"
            
            # Determine CSS class based on strict spec (Point 13)
            tag_class = "tag-reported"
            if "Putative" in ev_type: tag_class = "tag-putative"
            elif "Experimentally" in ev_type or "validated" in ev_type.lower(): tag_class = "tag-validated"
            
            st.markdown(f"""
            <div class='mol-card'>
                <div><span class='{tag_class}'>{ev_type}</span></div>
                <h3 style='margin-top: 0;'>{c.compound_name}</h3>
                <div class='metrics'>
                    <div><span class='metric-label'>Formula</span><br><span class='metric-val'>{c.molecular_formula or 'Unknown'}</span></div>
                    <div><span class='metric-label'>Mol. Weight</span><br><span class='metric-val'>{c.molecular_weight or 'Unknown'}</span></div>
                    <div><span class='metric-label'>Chemical Class</span><br><span class='metric-val'>{c.chemical_class or 'Unknown'}</span></div>
                </div>
            </div>
            """, unsafe_allow_html=True)
            
            # 3D Structure and SMILES
            col1, col2 = st.columns([1, 1])
            with col1:
                st.markdown("<span class='metric-label'>SMILES String</span>", unsafe_allow_html=True)
                st.markdown(f"<div class='smiles-box'>{c.smiles or 'No SMILES data in offline DB.'}</div>", unsafe_allow_html=True)
                
            with col2:
                st.markdown("<span class='metric-label'>3D Conformation (API-Free)</span>", unsafe_allow_html=True)
                safe_name = c.compound_name.replace(" ", "_").lower()
                sdf_path = os.path.join(BASE_DIR, "assets", "structures", f"{safe_name}.sdf")
                
                if os.path.exists(sdf_path):
                    with open(sdf_path, "r") as f:
                        sdf_data = f.read()
                    
                    # Render with py3Dmol
                    view = py3Dmol.view(width=400, height=300)
                    view.addModel(sdf_data, "sdf")
                    view.setStyle({'stick': {}})
                    view.setBackgroundColor('#0E1117')
                    view.zoomTo()
                    showmol(view, height=300, width=400)
                else:
                    st.info("3D Geometry pending offline calculation.")

st.markdown("---")
st.page_link("pages/4_Computational_Docking.py", label="Run AutoDock Vina Simulations", icon="🧬")
