import streamlit as st
import os
import sys
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# --- Config ---
st.set_page_config(page_title="Species Profile | EthnoDoc TCM", page_icon="🌿", layout="wide")

st.markdown("""
<style>
    .stApp { background-color: #0E1117; color: #FAFAFA; font-family: 'Inter', sans-serif; }
    .species-title { font-size: 2.5rem; font-weight: 300; color: #4CAF50; border-bottom: 1px solid #333; padding-bottom: 1rem; }
    .tax-box { background: #1E2127; padding: 1.5rem; border-radius: 8px; border-left: 4px solid #4CAF50; }
    .gallery-img { border-radius: 8px; object-fit: cover; width: 100%; max-height: 300px; }
    .source-box { background: #1a1c23; padding: 1rem; border-radius: 5px; margin-bottom: 1rem; border: 1px solid #2a2d35; }
</style>
""", unsafe_allow_html=True)

# Add base directory to path so we can import ethnodoc_models
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from ethnodoc_models import Species, SpeciesImage, HistoricalSource

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
    # st.switch_page("app.py") # Streamlit >= 1.30
    st.stop()
    
species_id = st.session_state['selected_species_id']
species = session.query(Species).filter_by(species_id=species_id).first()

if not species:
    st.error("Species not found in database.")
    st.stop()

# --- Header ---
st.markdown(f"<div class='species-title'>{species.english_name} <span style='color: #888;'>({species.chinese_name})</span></div>", unsafe_allow_html=True)

# Create Tabs for cleaner layout
tab1, tab2, tab3, tab4 = st.tabs(["Taxonomy", "TCM Properties", "Visual Evidence", "Historical Evidence"])

with tab1:
    st.markdown("### Taxonomic Classification")
    st.markdown(f"**Genus**: {species.taxonomy_genus}")
    st.markdown(f"**Species**: {species.scientific_name}")

with tab2:
    st.markdown("### TCM Properties")
    st.info("TCM properties (Nature, Flavor, Meridians) are being aggregated from Historical Sources...")
    # Mock data for UI presentation based on the spec
    tcm_cols = st.columns(3)
    tcm_cols[0].markdown(f"**Nature**: Neutral")
    tcm_cols[1].markdown(f"**Flavor**: Sweet / Bitter")
    tcm_cols[2].markdown(f"**Meridians**: Liver, Spleen")

with tab3:
    st.markdown("### Visual Evidence (Local Archive)")
    if not species.images:
        st.info("No verified visual evidence linked to this species.")
    else:
        # Display images in a clean grid
        num_images = len(species.images)
        cols = st.columns(min(num_images, 3))
        for idx, img in enumerate(species.images):
            col = cols[idx % 3]
            img_path = os.path.join(BASE_DIR, "assets", "images", img.filename)
            try:
                col.image(img_path, caption=f"ID: {img.image_id}", use_container_width=True)
                col.markdown(f"<span style='font-size:0.8rem; color:#888;'>Checksum: {img.sha256_checksum[:12]}</span>", unsafe_allow_html=True)
            except Exception:
                col.error(f"Failed to load local asset: {img.filename}")

with tab4:
    st.markdown("### Historical Evidence & Provenance")
    st.markdown("The following historical texts document claims or usage of this species.")
    
    if not species.historical_sources:
        st.info("No historical sources linked to this species.")
    else:
        for source in species.historical_sources:
            st.markdown(f"- 📜 **{source.title}**  <span style='color: #888; font-size: 0.85em;'>({source.provenance})</span>", unsafe_allow_html=True)
            
        st.write("")
        st.page_link("pages/2_Historical_Claims.py", label="View Strict Evidence Claims & Preparation Data", icon="🔍")
