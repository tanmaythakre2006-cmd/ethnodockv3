import os
import sys
import math
import base64
import json
from io import BytesIO
import pandas as pd
import numpy as np
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Draw

# Add root directory to sys.path
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_DIR)

import ethnodock_docking_engine as dock_eng
import ethnodock_interaction_engine as inter_eng
import ethnodock_bioisostere_engine as bio_eng
import ethnodock_admet_engine as admet_eng
import ethnodock_dossier_engine as dossier_eng
import ethnodock_paozhi_engine as paozhi_eng

# --- Page Configuration ---
st.set_page_config(
    page_title="EthnoDock Pro • In-Silico Pharmacognosy",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# --- Apple macOS / iOS Premium Glassmorphic Design System ---
st.markdown(
    """
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap');

    /* Base Reset & Apple Typography */
    html, body, [class*="css"], .stApp {
        font-family: -apple-system, BlinkMacSystemFont, "SF Pro Display", "SF Pro Text", "Inter", sans-serif !important;
        background: #08090C !important;
        color: #F5F5F7 !important;
        letter-spacing: -0.015em;
    }

    /* Ambient Hero Glow */
    .hero-container {
        position: relative;
        background: radial-gradient(circle at 50% 0%, rgba(48, 209, 88, 0.15) 0%, rgba(10, 132, 255, 0.08) 35%, rgba(8, 9, 12, 0.95) 75%);
        border: 1px solid rgba(255, 255, 255, 0.08);
        border-radius: 24px;
        padding: 50px 40px;
        text-align: center;
        margin-bottom: 30px;
        box-shadow: 0 30px 60px rgba(0, 0, 0, 0.5);
    }
    
    .hero-pill {
        display: inline-flex;
        align-items: center;
        gap: 8px;
        background: rgba(255, 255, 255, 0.06);
        border: 1px solid rgba(255, 255, 255, 0.12);
        padding: 6px 16px;
        border-radius: 30px;
        font-size: 0.8rem;
        font-weight: 600;
        color: #30D158;
        letter-spacing: 0.8px;
        text-transform: uppercase;
        margin-bottom: 18px;
    }

    .hero-title {
        font-size: 3.1rem;
        font-weight: 800;
        letter-spacing: -0.035em;
        line-height: 1.15;
        margin: 0 0 16px 0;
        background: linear-gradient(180deg, #FFFFFF 0%, #A1A1A6 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }

    .hero-subtitle {
        font-size: 1.15rem;
        color: #A1A1A6;
        max-width: 760px;
        margin: 0 auto 24px auto;
        line-height: 1.6;
        font-weight: 400;
    }

    /* Stat Bar */
    .hero-stats-grid {
        display: grid;
        grid-template-columns: repeat(4, 1fr);
        gap: 16px;
        max-width: 880px;
        margin: 0 auto;
        padding-top: 24px;
        border-top: 1px solid rgba(255, 255, 255, 0.08);
    }
    .hero-stat-box {
        text-align: center;
    }
    .hero-stat-num {
        font-size: 1.6rem;
        font-weight: 700;
        color: #FFFFFF;
    }
    .hero-stat-lbl {
        font-size: 0.72rem;
        color: #86868B;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        margin-top: 2px;
    }

    /* Command Search Card */
    .search-card-container {
        background: rgba(255, 255, 255, 0.04);
        backdrop-filter: blur(35px);
        -webkit-backdrop-filter: blur(35px);
        border: 1px solid rgba(255, 255, 255, 0.1);
        border-radius: 20px;
        padding: 28px;
        margin-bottom: 35px;
        box-shadow: 0 25px 50px rgba(0, 0, 0, 0.45);
    }

    /* Glass Cards */
    .apple-card {
        background: rgba(255, 255, 255, 0.035);
        backdrop-filter: blur(30px);
        -webkit-backdrop-filter: blur(30px);
        border: 1px solid rgba(255, 255, 255, 0.08);
        border-radius: 18px;
        padding: 24px;
        margin-bottom: 20px;
        box-shadow: 0 20px 40px rgba(0, 0, 0, 0.4);
    }

    .apple-card-compact {
        background: rgba(255, 255, 255, 0.025);
        border: 1px solid rgba(255, 255, 255, 0.06);
        border-radius: 14px;
        padding: 16px 20px;
        margin-bottom: 14px;
    }

    /* Mission Grid */
    .mission-grid {
        display: grid;
        grid-template-columns: repeat(3, 1fr);
        gap: 20px;
        margin-bottom: 35px;
    }
    .mission-box {
        background: rgba(255, 255, 255, 0.025);
        border: 1px solid rgba(255, 255, 255, 0.06);
        border-radius: 16px;
        padding: 26px;
        transition: all 0.2s ease;
    }
    .mission-box:hover {
        background: rgba(255, 255, 255, 0.04);
        border-color: rgba(255, 255, 255, 0.12);
        transform: translateY(-3px);
    }
    .mission-icon {
        font-size: 2rem;
        margin-bottom: 12px;
        display: inline-block;
    }
    .mission-title {
        font-size: 1.1rem;
        font-weight: 600;
        color: #FFFFFF;
        margin-bottom: 8px;
    }
    .mission-desc {
        font-size: 0.88rem;
        color: #86868B;
        line-height: 1.55;
    }

    /* Apple Pill Badges */
    .apple-badge {
        display: inline-flex;
        align-items: center;
        gap: 6px;
        background: rgba(48, 209, 88, 0.12);
        color: #30D158;
        border: 1px solid rgba(48, 209, 88, 0.28);
        padding: 4px 12px;
        border-radius: 20px;
        font-size: 0.75rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.6px;
    }
    .apple-badge-gold {
        background: rgba(255, 214, 10, 0.12);
        color: #FFD60A;
        border: 1px solid rgba(255, 214, 10, 0.28);
    }
    .apple-badge-blue {
        background: rgba(10, 132, 255, 0.12);
        color: #64D2FF;
        border: 1px solid rgba(10, 132, 255, 0.28);
    }
    .apple-badge-purple {
        background: rgba(191, 90, 242, 0.12);
        color: #BF5AF2;
        border: 1px solid rgba(191, 90, 242, 0.28);
    }

    /* Apple Stat Widgets */
    .apple-stat-box {
        background: rgba(255, 255, 255, 0.03);
        border: 1px solid rgba(255, 255, 255, 0.07);
        border-radius: 14px;
        padding: 16px;
        text-align: center;
        transition: all 0.2s ease;
    }
    .apple-stat-box:hover {
        background: rgba(255, 255, 255, 0.05);
        border-color: rgba(255, 255, 255, 0.12);
    }
    .apple-stat-lbl {
        font-size: 0.72rem;
        color: #86868B;
        text-transform: uppercase;
        letter-spacing: 0.8px;
        font-weight: 500;
    }
    .apple-stat-val {
        font-size: 1.6rem;
        font-weight: 700;
        color: #FFFFFF;
        margin-top: 4px;
        letter-spacing: -0.02em;
    }

    /* macOS Window Frame */
    .macos-window {
        border-radius: 14px;
        background: #0E1118;
        border: 1px solid rgba(255, 255, 255, 0.1);
        overflow: hidden;
        box-shadow: 0 25px 50px rgba(0,0,0,0.6);
        margin: 15px 0;
    }
    .macos-bar {
        background: #181C26;
        padding: 10px 14px;
        display: flex;
        align-items: center;
        gap: 8px;
        border-bottom: 1px solid rgba(255, 255, 255, 0.06);
    }
    .macos-dot {
        width: 11px;
        height: 11px;
        border-radius: 50%;
        display: inline-block;
    }
    .dot-red { background: #FF5F56; }
    .dot-yellow { background: #FFBD2E; }
    .dot-green { background: #27C93F; }
    .macos-title {
        color: #86868B;
        font-size: 12px;
        font-weight: 500;
        margin-left: 8px;
    }

    /* Primary Apple Buttons */
    .stButton>button {
        background: #0A84FF !important;
        color: #FFFFFF !important;
        border: none !important;
        border-radius: 22px !important;
        font-weight: 600 !important;
        font-size: 0.88rem !important;
        padding: 0.55rem 1.4rem !important;
        box-shadow: 0 4px 14px rgba(10, 132, 255, 0.3) !important;
        transition: all 0.2s cubic-bezier(0.16, 1, 0.3, 1) !important;
    }
    .stButton>button:hover {
        background: #0071E3 !important;
        transform: scale(1.02) !important;
        box-shadow: 0 6px 20px rgba(10, 132, 255, 0.5) !important;
    }

    /* Segmented Tab Headers */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
        background: rgba(255, 255, 255, 0.03);
        padding: 6px;
        border-radius: 16px;
        border: 1px solid rgba(255, 255, 255, 0.08);
        margin-bottom: 24px;
    }
    .stTabs [data-baseweb="tab"] {
        border-radius: 12px;
        color: #86868B;
        font-weight: 600;
        padding: 8px 20px;
        border: none;
    }
    .stTabs [aria-selected="true"] {
        background: #0A84FF !important;
        color: #FFFFFF !important;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# --- State Clearing Helper ---
def clear_session_docking_state():
    keys_to_clear = [k for k in st.session_state.keys() if any(k.startswith(p) for p in [
        'setup_done_', 'center_', 'dims_', 'rec_pdbqt_', 'docking_done_', 'docking_data_',
        'pdb_id_', 'smiles_', 'uff_delta_', 'interactions_df_', 'docking_var_done_',
        'interactions_var_df_', 'docking_var_data_', 'poses_'
    ])]
    for k in keys_to_clear:
        del st.session_state[k]

# --- Load Master Image Gallery ---
@st.cache_data
def load_image_gallery():
    gallery_path = os.path.join(BASE_DIR, "species_image_gallery.json")
    if os.path.exists(gallery_path):
        with open(gallery_path, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}

gallery_data = load_image_gallery()

def get_species_photo_b64(common_name, botanical_name=""):
    botanical_clean = botanical_name.lower().strip()
    common_clean = common_name.lower().strip()
    
    matched_entry = None
    for k, v in gallery_data.items():
        v_sci = v.get("scientific_name", "").lower()
        v_com = v.get("common_name", "").lower()
        if botanical_clean and (botanical_clean in v_sci or v_sci in botanical_clean):
            matched_entry = v
            break
        if common_clean and (common_clean in v_com or v_com in common_clean):
            matched_entry = v
            break
            
    if matched_entry:
        local_photos = matched_entry.get("local_photos", [])
        for lp in local_photos:
            lp_fixed = os.path.join(BASE_DIR, "species_images", os.path.basename(lp))
            if os.path.exists(lp_fixed):
                try:
                    with open(lp_fixed, "rb") as img_f:
                        ext = "png" if lp_fixed.endswith(".png") else "jpeg"
                        return f"data:image/{ext};base64,{base64.b64encode(img_f.read()).decode()}"
                except Exception:
                    pass
        photos = matched_entry.get("photos", [])
        if photos:
            return photos[0]
            
    return None

# --- Helper: Render 2D Molecular Structure ---
def get_image_base64(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=(380, 280))
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            return base64.b64encode(buffered.getvalue()).decode()
    except Exception:
        pass
    return ""

# --- Load Master TCM Database ---
@st.cache_data
def load_tcm_master():
    csv_path = os.path.join(BASE_DIR, "ethnodock_tcm_master.csv")
    return pd.read_csv(csv_path)

df = load_tcm_master()

# ==========================================
# 🧭 TOP SEGMENTED NAVIGATION
# ==========================================
tab_intro, tab_workbench = st.tabs([
    "🏛️ Project Introduction & Discovery",
    "🔬 In-Silico Molecular Workbench"
])

# ==========================================
# 🌟 TAB 1: DEDICATED INTRODUCTION & SEARCH
# ==========================================
with tab_intro:
    # 1. Hero Introduction
    st.markdown("""
    <div class="hero-container">
        <span class="hero-pill">🌿 EthnoDock Pro • In-Silico Pharmacognosy Initiative</span>
        <h1 class="hero-title">Bridging 2,000 Years of Botanical Wisdom<br>with Western Structural Pharmacology</h1>
        <p class="hero-subtitle">
            EthnoDock is a world-first computational ethnopharmacology suite. We digitize canonical Traditional Chinese Medicine (TCM) pharmacopoeias—from the Han Dynasty's <i>Shennong Bencaojing</i> to the Ming Dynasty's <i>Bencao Gangmu</i>—and systematically validate their active phytochemical mechanisms against verified human macromolecular drug targets using rigorous empirical molecular docking.
        </p>
        
        <div class="hero-stats-grid">
            <div class="hero-stat-box">
                <div class="hero-stat-num" style="color:#30D158;">1,180+</div>
                <div class="hero-stat-lbl">Curated Formulations</div>
            </div>
            <div class="hero-stat-box">
                <div class="hero-stat-num" style="color:#64D2FF;">AutoDock Vina</div>
                <div class="hero-stat-lbl">Empirical Scoring Engine</div>
            </div>
            <div class="hero-stat-box">
                <div class="hero-stat-num" style="color:#FFD60A;">3Dmol.js</div>
                <div class="hero-stat-lbl">WebGL Interaction Studio</div>
            </div>
            <div class="hero-stat-box">
                <div class="hero-stat-num" style="color:#BF5AF2;">Paozhi 炮制</div>
                <div class="hero-stat-lbl">Detoxification Audit</div>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # 2. Main Search & Catalog Explorer
    st.markdown("""
    <div class="search-card-container">
        <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:16px;">
            <div>
                <h3 style="margin:0; font-size:1.3rem; font-weight:700; color:#FFFFFF;">🔍 Pharmacognosy Catalog & Target Search</h3>
                <p style="margin:2px 0 0 0; font-size:0.88rem; color:#86868B;">Select any medicinal plant, therapeutic protein target, or phytochemical to inspect its molecular identity.</p>
            </div>
            <span class="apple-badge apple-badge-blue">Live Query</span>
        </div>
    """, unsafe_allow_html=True)

    col_s1, col_s2 = st.columns([1, 2], gap="large")
    with col_s1:
        search_by = st.selectbox(
            "Search Catalog by Filter:",
            ["Common Name", "Protein Target", "Active Phytochemical"],
            on_change=clear_session_docking_state,
            key="tab1_search_by"
        )
    with col_s2:
        if search_by == "Common Name":
            options = df["Common Name"].unique().tolist()
        elif search_by == "Protein Target":
            options = df["Protein Target"].unique().tolist()
        else:
            options = df["Active Phytochemical"].unique().tolist()

        selected_option = st.selectbox(
            f"Select {search_by}:",
            options,
            on_change=clear_session_docking_state,
            key="tab1_selected_option"
        )

    st.markdown("</div>", unsafe_allow_html=True)

    # 3. Selected Herb Spotlight Card
    selected_data = df[df[search_by] == selected_option]

    if not selected_data.empty:
        row = selected_data.iloc[0]
        plant_photo_b64 = get_species_photo_b64(row['Common Name'], row['Botanical Name'])
        img_b64 = get_image_base64(row['SMILES'])

        col_p1, col_p2, col_p3 = st.columns([1, 1, 1.4], gap="medium")

        with col_p1:
            if plant_photo_b64:
                st.markdown(f"""
                <div class="apple-card" style="padding:14px; text-align:center; height:100%;">
                    <div style="border-radius:12px; overflow:hidden; height:190px; background:#000;">
                        <img src="{plant_photo_b64}" style="width:100%; height:100%; object-fit:cover;"/>
                    </div>
                    <div style="margin-top:12px;">
                        <span class="apple-badge">Verified Field Specimen</span>
                        <div style="font-weight:600; font-size:15px; margin-top:4px;">{row['Common Name']}</div>
                        <div style="font-size:12px; color:#86868B; font-style:italic;">{row['Botanical Name']}</div>
                    </div>
                </div>
                """, unsafe_allow_html=True)

        with col_p2:
            img_tag = f'<img src="data:image/png;base64,{img_b64}" style="width:100%; height:190px; object-fit:contain;"/>' if img_b64 else '<p>2D Topology</p>'
            st.markdown(f"""
            <div class="apple-card" style="padding:14px; text-align:center; height:100%;">
                <div style="background:#000000; border-radius:12px; padding:6px; height:190px; display:flex; align-items:center; justify-content:center;">
                    {img_tag}
                </div>
                <div style="margin-top:12px;">
                    <span class="apple-badge apple-badge-gold">{row['Chemical Class']}</span>
                    <div style="font-weight:600; font-size:15px; margin-top:4px; color:#52B788;">{row['Active Phytochemical']}</div>
                    <div style="font-size:11px; color:#86868B;">PubChem CID: {row['PubChem CID']}</div>
                </div>
            </div>
            """, unsafe_allow_html=True)

        with col_p3:
            st.markdown(f"""
            <div class="apple-card" style="padding:20px; height:100%; display:flex; flex-direction:column; justify-content:space-between;">
                <div>
                    <span class="apple-badge apple-badge-purple">{row['Dynasty']} • {row['Classical Source']}</span>
                    <div style="font-size:14px; color:#F5F5F7; margin-top:8px; font-style:italic; line-height:1.4;">"{row['English Translation']}"</div>
                    <hr style="border:none; border-top:1px solid rgba(255,255,255,0.06); margin:12px 0;">
                    <div style="font-size:12px; color:#86868B; text-transform:uppercase;">Primary Drug Target</div>
                    <div style="font-size:15px; font-weight:700; color:#30D158;">{row['Protein Target']} ({row['Gene Symbol']})</div>
                    <div style="font-size:12px; color:#A1A1A6; margin-top:2px;"><b>RCSB PDB ID:</b> <span style="color:#64D2FF; font-weight:700;">{row['PDB ID']}</span> • <b>UniProt:</b> {row['UniProt ID']}</div>
                </div>
                <div style="margin-top:14px;">
                    <p style="font-size:12px; color:#86868B; margin:0;">👉 Switch to the <b>In-Silico Molecular Workbench</b> tab above to run the live 3D docking simulation for {row['Common Name']}!</p>
                </div>
            </div>
            """, unsafe_allow_html=True)

    # 4. Mission Pillars & Methodology (Who We Are & What We Do)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("""
    <div style="margin-bottom:16px;">
        <span class="apple-badge apple-badge-gold">Core Methodology</span>
        <h2 style="margin:8px 0 4px 0; font-size:1.7rem; font-weight:700;">Our Scientific & Philological Architecture</h2>
        <p style="margin:0; color:#86868B; font-size:0.95rem;">How EthnoDock bridges historical philology with state-of-the-art computational biophysics.</p>
    </div>

    <div class="mission-grid">
        <div class="mission-box">
            <div class="mission-icon">📜</div>
            <div class="mission-title">1. Classical Corpus Digitization</div>
            <div class="mission-desc">We extract verbatim clinical indications and formula structures from canonical Chinese medical texts (*Han, Tang, Song, Ming* dynasties) and cross-reference them with Kew Gardens botanical taxonomy.</div>
        </div>
        <div class="mission-box">
            <div class="mission-icon">⚛️</div>
            <div class="mission-title">2. Empirical Molecular Docking</div>
            <div class="mission-desc">Using AutoDock Vina v1.2.7 and RDKit ETKDGv3 conformer generation, we model binding affinities (&Delta;G in kcal/mol), RMSD conformational hierarchies, and thermodynamic inhibition constants (Ki).</div>
        </div>
        <div class="mission-box">
            <div class="mission-icon">⚗️</div>
            <div class="mission-title">3. Paozhi (炮制) Detoxification</div>
            <div class="mission-desc">We model classical thermal processing alchemy (soaking, wine-steaming, sand-roasting), proving how toxic raw diester alkaloids thermally hydrolyze into safe, non-toxic therapeutic metabolites.</div>
        </div>
    </div>
    """, unsafe_allow_html=True)


# ==========================================
# 🔬 TAB 2: IN-SILICO MOLECULAR WORKBENCH
# ==========================================
with tab_workbench:
    if not selected_data.empty:
        for idx, row in selected_data.iterrows():
            # Check Paozhi Availability
            paozhi_key = paozhi_eng.has_paozhi(row['Common Name'])
            active_compound_name = row['Active Phytochemical']
            active_smiles = row['SMILES']
            active_chemical_class = row['Chemical Class']
            is_processed_state = False

            # ==========================================
            # STAGE 01: BOTANICAL SPECIMEN & METADATA PROFILE
            # ==========================================
            st.markdown("""
            <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:14px;">
                <div style="display:flex; align-items:center; gap:8px;">
                    <span class="apple-badge apple-badge-purple">Stage 01</span>
                    <h3 style="margin:0; font-size:1.25rem; font-weight:600;">Botanical Specimen & Chemical Characterization</h3>
                </div>
            </div>
            """, unsafe_allow_html=True)

            if paozhi_key:
                pz_info = paozhi_eng.get_paozhi_data(paozhi_key)
                col_pz1, col_pz2 = st.columns([1, 2], vertical_alignment="center")
                with col_pz1:
                    paozhi_choice = st.radio(
                        "Pharmaceutical Processing State:",
                        ["Raw Form (生品)", "Processed (炮制品)"],
                        horizontal=True,
                        key=f"pz_switch_tab2_{idx}",
                        on_change=clear_session_docking_state
                    )
                if "Processed" in paozhi_choice:
                    is_processed_state = True
                    active_compound_name = pz_info['processed_compound']
                    active_smiles = pz_info['processed_smiles']
                    active_chemical_class = pz_info['processed_class']
                    with col_pz2:
                        st.markdown(f"""
                        <div style="background:rgba(255, 214, 10, 0.08); border:1px solid rgba(255, 214, 10, 0.25); border-radius:12px; padding:10px 14px; font-size:12px;">
                            <span style="color:#FFD60A; font-weight:600;">⚗️ Classical Paozhi Transformation:</span> {pz_info['reaction_equation']}<br>
                            <span style="color:#A1A1A6;">{pz_info['detox_benefit']}</span>
                        </div>
                        """, unsafe_allow_html=True)
                else:
                    with col_pz2:
                        st.markdown(f"""
                        <div style="background:rgba(255, 69, 58, 0.08); border:1px solid rgba(255, 69, 58, 0.25); border-radius:12px; padding:10px 14px; font-size:12px; color:#FF453A;">
                            ⚠️ <b>Raw Toxicity:</b> {pz_info['raw_toxicity_warning']}
                        </div>
                        """, unsafe_allow_html=True)

            # 3-Column Apple Dashboard Card
            col_img_bot, col_img_mol, col_details = st.columns([1, 1, 1.4], gap="medium")

            # Column 1: Botanical Plant Photo
            with col_img_bot:
                plant_photo_b64 = get_species_photo_b64(row['Common Name'], row['Botanical Name'])
                if plant_photo_b64:
                    st.markdown(f"""
                    <div class="apple-card" style="padding:14px; text-align:center; height:100%;">
                        <div style="border-radius:12px; overflow:hidden; height:200px; background:#000;">
                            <img src="{plant_photo_b64}" style="width:100%; height:100%; object-fit:cover;"/>
                        </div>
                        <div style="margin-top:12px;">
                            <span class="apple-badge" style="font-size:11px;">Verified Specimen</span>
                            <div style="font-weight:600; font-size:14px; margin-top:4px; color:#FFF;">{row['Common Name']}</div>
                            <div style="font-size:12px; color:#86868B; font-style:italic;">{row['Botanical Name']}</div>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)

            # Column 2: 2D Chemical Molecule
            with col_img_mol:
                img_b64 = get_image_base64(active_smiles)
                img_tag = f'<img src="data:image/png;base64,{img_b64}" style="width:100%; height:200px; object-fit:contain;"/>' if img_b64 else '<p style="color:#666;">2D Structure</p>'
                st.markdown(f"""
                <div class="apple-card" style="padding:14px; text-align:center; height:100%;">
                    <div style="background:#000000; border-radius:12px; padding:6px; height:200px; display:flex; align-items:center; justify-content:center;">
                        {img_tag}
                    </div>
                    <div style="margin-top:12px;">
                        <span class="apple-badge apple-badge-gold" style="font-size:11px;">{active_chemical_class}</span>
                        <div style="font-weight:600; font-size:14px; margin-top:4px; color:#52B788;">{active_compound_name}</div>
                        <div style="font-size:11px; color:#86868B;">CID: {row['PubChem CID']} • InChIKey</div>
                    </div>
                </div>
                """, unsafe_allow_html=True)

            # Column 3: Classical & Target Details
            with col_details:
                st.markdown(f"""
                <div class="apple-card" style="padding:18px; height:100%;">
                    <div style="border-bottom:1px solid rgba(255,255,255,0.06); padding-bottom:10px; margin-bottom:10px;">
                        <div style="font-size:11px; color:#86868B; text-transform:uppercase;">Classical Canon • {row['Dynasty']}</div>
                        <div style="font-size:14px; font-weight:600; color:#FFD60A;">{row['Classical Source']}</div>
                        <div style="font-size:13px; color:#F5F5F7; margin-top:4px; font-style:italic;">"{row['English Translation']}"</div>
                    </div>
                    <div>
                        <div style="font-size:11px; color:#86868B; text-transform:uppercase;">Western Target Mechanism</div>
                        <div style="font-size:14px; font-weight:600; color:#30D158;">{row['Protein Target']} ({row['Gene Symbol']})</div>
                        <div style="font-size:12px; color:#A1A1A6; margin-top:2px;"><b>UniProt:</b> {row['UniProt ID']} • <b>PDB:</b> <span style="color:#64D2FF; font-weight:700;">{row['PDB ID']}</span></div>
                        <div style="font-size:11px; color:#86868B; margin-top:4px;"><b>Indication:</b> {row['Western Indication']}</div>
                        <div style="font-size:11px; color:#6E6E73; margin-top:4px; word-break:break-all;"><b>SMILES:</b> <code>{active_smiles}</code></div>
                    </div>
                </div>
                """, unsafe_allow_html=True)

            # ==========================================
            # STAGE 02: RECEPTOR & SMART CAVITY SETUP
            # ==========================================
            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown("""
            <div style="display:flex; align-items:center; gap:8px; margin-bottom:14px;">
                <span class="apple-badge apple-badge-blue">Stage 02</span>
                <h3 style="margin:0; font-size:1.25rem; font-weight:600;">Receptor & Search Cavity Configuration</h3>
            </div>
            """, unsafe_allow_html=True)

            with st.expander(f"⚙️ Target Binding Cavity Controller — {row['Protein Target']} (PDB: {row['PDB ID']})", expanded=True):
                pdb_id = row['PDB ID']
                smiles = active_smiles

                col_btn, col_info = st.columns([1, 2], vertical_alignment="center")
                with col_btn:
                    init_btn = st.button(f"📥 Initialize PDB {pdb_id}", key=f"init_tab2_{idx}", use_container_width=True)
                with col_info:
                    st.caption(f"Downloads verified structure directly from RCSB PDB, calculates pocket centroid, and formats AutoDock atom types.")

                if init_btn:
                    with st.spinner(f"Retrieving structure {pdb_id} from RCSB..."):
                        output_pdb = os.path.join(BASE_DIR, f"{pdb_id}.pdb")
                        receptor_pdbqt = dock_eng.fetch_receptor(pdb_id, output_pdb)
                        if receptor_pdbqt:
                            center, dims = dock_eng.smart_cavity_finder(output_pdb)
                            st.session_state[f'center_{idx}'] = center
                            st.session_state[f'dims_{idx}'] = dims
                            st.session_state[f'rec_pdbqt_{idx}'] = receptor_pdbqt
                            st.session_state[f'setup_done_{idx}'] = True
                            st.success(f"Receptor {pdb_id} initialized!")
                        else:
                            st.error(f"Could not load structure {pdb_id}.")

                if st.session_state.get(f'setup_done_{idx}', False):
                    center = st.session_state[f'center_{idx}']
                    dims = st.session_state[f'dims_{idx}']

                    col_c1, col_c2, col_c3, col_s1, col_s2, col_s3 = st.columns(6)
                    cx = col_c1.number_input("Center X", value=float(center[0]), format="%.2f", key=f"cx_tab2_{idx}")
                    cy = col_c2.number_input("Center Y", value=float(center[1]), format="%.2f", key=f"cy_tab2_{idx}")
                    cz = col_c3.number_input("Center Z", value=float(center[2]), format="%.2f", key=f"cz_tab2_{idx}")
                    sx = col_s1.number_input("Size X", value=float(dims[0]), format="%.2f", key=f"sx_tab2_{idx}")
                    sy = col_s2.number_input("Size Y", value=float(dims[1]), format="%.2f", key=f"sy_tab2_{idx}")
                    sz = col_s3.number_input("Size Z", value=float(dims[2]), format="%.2f", key=f"sz_tab2_{idx}")

                    exhaustiveness = st.slider("Vina Exhaustiveness (Sampling Precision)", min_value=4, max_value=32, value=8, step=4, key=f"exh_tab2_{idx}")

                    # ==========================================
                    # STAGE 03: DOCKING SIMULATION & 3D STUDIO
                    # ==========================================
                    st.markdown("<br>", unsafe_allow_html=True)
                    st.markdown("""
                    <div style="display:flex; align-items:center; gap:8px; margin-bottom:14px;">
                        <span class="apple-badge apple-badge-gold">Stage 03</span>
                        <h3 style="margin:0; font-size:1.25rem; font-weight:600;">In-Silico Docking & 3D WebGL Interaction Studio</h3>
                    </div>
                    """, unsafe_allow_html=True)

                    if st.button(f"🚀 Execute Molecular Simulation for {active_compound_name}", key=f"dock_tab2_{idx}", use_container_width=False):
                        with st.spinner(f"Minimizing conformer and docking into {pdb_id}..."):
                            receptor_pdbqt = st.session_state[f'rec_pdbqt_{idx}']
                            ligand_pdbqt_path = os.path.join(BASE_DIR, "active_ligand.pdbqt")
                            ligand_pdbqt, uff_delta = dock_eng.prepare_ligand(smiles, ligand_pdbqt_path)

                            if receptor_pdbqt and ligand_pdbqt:
                                raw_log, parsed_poses, out_pdbqt = dock_eng.run_vina_docking(
                                    receptor_pdbqt, ligand_pdbqt, [cx, cy, cz], [sx, sy, sz], exhaustiveness=exhaustiveness
                                )
                                if parsed_poses:
                                    st.session_state[f'docking_data_{idx}'] = parsed_poses
                                    st.session_state[f'out_pdbqt_{idx}'] = out_pdbqt
                                    st.session_state[f'docking_done_{idx}'] = True
                                    st.session_state[f'pdb_id_{idx}'] = pdb_id
                                    st.session_state[f'smiles_{idx}'] = smiles
                                    st.session_state[f'uff_delta_{idx}'] = uff_delta
                                    st.success("AutoDock Vina simulation converged!")
                                else:
                                    st.error("Docking failed. Log:")
                                    st.code(raw_log)

            if st.session_state.get(f'docking_done_{idx}', False):
                data = st.session_state[f'docking_data_{idx}']
                out_pdbqt = st.session_state[f'out_pdbqt_{idx}']
                uff_delta = st.session_state.get(f'uff_delta_{idx}', 0.0)
                receptor_pdbqt = st.session_state[f'rec_pdbqt_{idx}']

                table_data = []
                for d in data:
                    aff = d['affinity']
                    rt = 0.001987 * 298.15
                    ki_molar = math.exp(aff / rt)
                    if ki_molar < 1e-6:
                        ki_str = f"{ki_molar * 1e9:.2f} nM"
                    elif ki_molar < 1e-3:
                        ki_str = f"{ki_molar * 1e6:.2f} µM"
                    else:
                        ki_str = f"{ki_molar * 1e3:.2f} mM"

                    table_data.append({
                        "Mode": d['mode'],
                        "Affinity (kcal/mol)": aff,
                        "Estimated Ki": ki_str,
                        "RMSD (l.b.)": d['rmsd_lb'],
                        "RMSD (u.b.)": d['rmsd_ub']
                    })

                # Pose Selector
                col_sel, col_empty = st.columns([2, 1])
                with col_sel:
                    pose_options = [f"Mode {d['Mode']} • Affinity: {d['Affinity (kcal/mol)']} kcal/mol (Ki: {d['Estimated Ki']})" for d in table_data]
                    selected_mode_str = st.selectbox("Inspect Conformation Mode:", pose_options, key=f"pose_select_tab2_{idx}")
                    selected_idx = pose_options.index(selected_mode_str)
                    selected_pose_data = table_data[selected_idx]

                # Extract 3D Pose coordinates
                poses = inter_eng.extract_poses(out_pdbqt)
                if poses and selected_idx < len(poses):
                    selected_pose_str = poses[selected_idx]
                    interactions_df = inter_eng.calc_interactions(selected_pose_str, receptor_pdbqt, cutoff=4.0)
                    st.session_state[f'interactions_df_{idx}'] = interactions_df

                    interacting_res = list(interactions_df["Receptor Residue"].unique()) if not interactions_df.empty else []

                    # Apple Stat Widgets
                    col_w1, col_w2, col_w3, col_w4 = st.columns(4)
                    with col_w1:
                        st.markdown(f"""
                        <div class="apple-stat-box">
                            <div class="apple-stat-lbl">Binding Affinity (ΔG)</div>
                            <div class="apple-stat-val" style="color:#30D158;">{selected_pose_data['Affinity (kcal/mol)']} <span style="font-size:0.9rem;">kcal/mol</span></div>
                        </div>
                        """, unsafe_allow_html=True)
                    with col_w2:
                        st.markdown(f"""
                        <div class="apple-stat-box">
                            <div class="apple-stat-lbl">Estimated Inhibition (Ki)</div>
                            <div class="apple-stat-val" style="color:#64D2FF;">{selected_pose_data['Estimated Ki']}</div>
                        </div>
                        """, unsafe_allow_html=True)
                    with col_w3:
                        st.markdown(f"""
                        <div class="apple-stat-box">
                            <div class="apple-stat-lbl">UFF Minimization (ΔE)</div>
                            <div class="apple-stat-val" style="color:#FFD60A;">{uff_delta:.2f} <span style="font-size:0.9rem;">kcal/mol</span></div>
                        </div>
                        """, unsafe_allow_html=True)
                    with col_w4:
                        st.markdown(f"""
                        <div class="apple-stat-box">
                            <div class="apple-stat-lbl">Pocket Anchors (<4.0Å)</div>
                            <div class="apple-stat-val">{len(interacting_res)} <span style="font-size:0.9rem;">Residues</span></div>
                        </div>
                        """, unsafe_allow_html=True)

                    st.markdown("<br>", unsafe_allow_html=True)

                    # macOS 3D Interaction Studio
                    col_3d, col_table = st.columns([3, 2], gap="large")

                    with col_3d:
                        style_c1, style_c2, style_c3 = st.columns(3)
                        with style_c1:
                            receptor_style = style_c1.selectbox("Receptor View", ['cartoon', 'stick', 'sphere'], index=0, key=f"rec_s_tab2_{idx}")
                        with style_c2:
                            ligand_style = style_c2.selectbox("Ligand View", ['stick', 'sphere', 'cross'], index=0, key=f"lig_s_tab2_{idx}")
                        with style_c3:
                            show_surface = style_c3.checkbox("Pocket Mesh", value=False, key=f"surf_s_tab2_{idx}")

                        with open(receptor_pdbqt, 'r', encoding='utf-8') as rf:
                            receptor_str = rf.read()

                        viewer_html = inter_eng.build_3dmol_html(
                            container_id=f"viewer_macos_tab2_{idx}",
                            receptor_data=receptor_str,
                            ligand_data=selected_pose_str,
                            interactions_df=interactions_df,
                            receptor_style=receptor_style,
                            ligand_style=ligand_style,
                            show_surface=show_surface,
                            height=480
                        )
                        
                        st.markdown("""
                        <div class="macos-window">
                            <div class="macos-bar">
                                <span class="macos-dot dot-red"></span>
                                <span class="macos-dot dot-yellow"></span>
                                <span class="macos-dot dot-green"></span>
                                <span class="macos-title">3D WebGL Studio • Interaction Engine</span>
                            </div>
                        """, unsafe_allow_html=True)
                        components.html(viewer_html, height=490)
                        st.markdown("</div>", unsafe_allow_html=True)

                    with col_table:
                        st.markdown("#### 🕸️ Pocket Contact Network")
                        if not interactions_df.empty:
                            display_df = interactions_df[["Receptor Residue", "Distance (Å)", "Interaction Type"]]
                            st.dataframe(display_df, hide_index=True, use_container_width=True)
                            h_bonds = interactions_df[interactions_df["Interaction Type"] == "Hydrogen / Polar Bond"]
                            st.success(f"**Identified:** {len(h_bonds)} Hydrogen/Polar Bonds (Crimson) and {len(interactions_df)-len(h_bonds)} Hydrophobic Contacts (Cyan).")
                        else:
                            st.info("No close contacts (< 4.0 Å) detected for this conformation.")

                        st.markdown("#### 📋 9-Pose Conformational Hierarchy")
                        st.dataframe(pd.DataFrame(table_data), hide_index=True, use_container_width=True)

                    # ==========================================
                    # STAGE 04: BIOISOSTERE LEAD OPTIMIZATION
                    # ==========================================
                    st.markdown("<br>", unsafe_allow_html=True)
                    st.markdown("""
                    <div style="display:flex; align-items:center; gap:8px; margin-bottom:14px;">
                        <span class="apple-badge apple-badge-purple">Stage 04</span>
                        <h3 style="margin:0; font-size:1.25rem; font-weight:600;">Semi-Synthetic Bioisostere Lead Optimization</h3>
                    </div>
                    """, unsafe_allow_html=True)

                    variants = bio_eng.generate_tcm_derivatives(smiles)

                    if variants:
                        var_labels = [f"{v['name']} — {v['variant_smiles'][:32]}..." for v in variants]
                        selected_var_label = st.selectbox("Select Rational Lead Modification:", var_labels, key=f"var_sel_tab2_{idx}")
                        selected_var_idx = var_labels.index(selected_var_label)
                        chosen_var = variants[selected_var_idx]

                        col_var_info, col_var_dock = st.columns([3, 1], vertical_alignment="center")
                        with col_var_info:
                            st.markdown(f"""
                            <div class="apple-card-compact">
                                <span style="color:#64D2FF; font-weight:600;">Medicinal Chemistry Rationale:</span> {chosen_var['rationale']}<br>
                                <span style="font-size:11px; color:#86868B;"><b>Variant SMILES:</b> <code>{chosen_var['variant_smiles']}</code></span>
                            </div>
                            """, unsafe_allow_html=True)
                        with col_var_dock:
                            dock_var_btn = st.button("⚡ Dock Derivative", key=f"dock_var_tab2_{idx}", use_container_width=True)

                        if dock_var_btn:
                            with st.spinner("Docking semi-synthetic analog..."):
                                var_pdbqt_path = os.path.join(BASE_DIR, "var_ligand.pdbqt")
                                var_ligand_pdbqt, _ = dock_eng.prepare_ligand(chosen_var['variant_smiles'], var_pdbqt_path)
                                if var_ligand_pdbqt:
                                    _, parsed_var_poses, var_out_pdbqt = dock_eng.run_vina_docking(
                                        receptor_pdbqt, var_ligand_pdbqt, [cx, cy, cz], [sx, sy, sz], exhaustiveness=exhaustiveness
                                    )
                                    if parsed_var_poses:
                                        st.session_state[f'docking_var_data_{idx}'] = parsed_var_poses
                                        st.session_state[f'var_out_pdbqt_{idx}'] = var_out_pdbqt
                                        st.session_state[f'docking_var_done_{idx}'] = True

                        if st.session_state.get(f'docking_var_done_{idx}', False):
                            var_data = st.session_state[f'docking_var_data_{idx}']
                            var_best_aff = var_data[0]['affinity']
                            parent_best_aff = selected_pose_data['Affinity (kcal/mol)']
                            delta_aff = var_best_aff - parent_best_aff

                            col_va, col_vb = st.columns([1, 2], vertical_alignment="center")
                            with col_va:
                                st.markdown(f"""
                                <div class="apple-stat-box">
                                    <div class="apple-stat-lbl">Derivative Affinity</div>
                                    <div class="apple-stat-val" style="color:#64D2FF;">{var_best_aff} <span style="font-size:0.9rem;">kcal/mol</span></div>
                                </div>
                                """, unsafe_allow_html=True)
                            with col_vb:
                                if delta_aff < 0:
                                    st.success(f"✨ **Lead Optimization Success:** Derivative enhances binding affinity by **{abs(delta_aff):.2f} kcal/mol** over the natural compound.")
                                else:
                                    st.info("Derivative maintains stable binding compatibility.")

                    # ==========================================
                    # STAGE 05: ADMET & SCIENTIFIC DOSSIER
                    # ==========================================
                    st.markdown("<br>", unsafe_allow_html=True)
                    st.markdown("""
                    <div style="display:flex; align-items:center; gap:8px; margin-bottom:14px;">
                        <span class="apple-badge">Stage 05</span>
                        <h3 style="margin:0; font-size:1.25rem; font-weight:600;">ADMET Pharmacokinetics & Dossier Export</h3>
                    </div>
                    """, unsafe_allow_html=True)

                    orig_adme = admet_eng.get_admet_profile(smiles)
                    adme_data_list = []
                    if orig_adme:
                        state_lbl = "Processed Form (炮制品)" if is_processed_state else "Natural Extract (生品)"
                        orig_adme["Compound Entity"] = f"{active_compound_name} [{state_lbl}]"
                        adme_data_list.append(orig_adme)

                    if variants:
                        for i, v in enumerate(variants[:3]):
                            v_adme = admet_eng.get_admet_profile(v['variant_smiles'])
                            if v_adme:
                                v_adme["Compound Entity"] = f"Derivative {i+1}: {v['name']}"
                                adme_data_list.append(v_adme)

                    if adme_data_list:
                        df_adme = pd.DataFrame(adme_data_list)
                        cols = ['Compound Entity', 'Molecular Weight', 'LogP', 'TPSA (Å²)', 'H-Bond Donors', 'H-Bond Acceptors', 'Rotatable Bonds', 'QED Drug-Likeness', 'Lipinski Violations', 'PAINS Screen']
                        st.dataframe(df_adme[[c for c in cols if c in df_adme.columns]], hide_index=True, use_container_width=True)

                    # Dossier HTML Export
                    poses_html = pd.DataFrame(table_data).to_html(index=False) if 'table_data' in locals() else "<p>None</p>"
                    interactions_clean_html = interactions_df[["Receptor Residue", "Distance (Å)", "Interaction Type"]].to_html(index=False) if ('interactions_df' in locals() and not interactions_df.empty) else "<p>None</p>"

                    variant_dossier_data = None
                    if 'chosen_var' in locals():
                        variant_dossier_data = {
                            "name": chosen_var["name"],
                            "rationale": chosen_var["rationale"],
                            "smiles": chosen_var["variant_smiles"],
                            "affinity": locals().get('var_best_aff', 'N/A')
                        }

                    dossier_html = dossier_eng.generate_tcm_dossier_html(
                        species_name=row['Common Name'],
                        scientific_name=row['Botanical Name'],
                        chinese_name=row['Chinese Name'],
                        source_title=row['Classical Source'],
                        claim_text=row['Ancient Claim'],
                        translation=row['English Translation'],
                        compound_name=active_compound_name,
                        smiles=smiles,
                        target_name=row['Protein Target'],
                        pdb_id=row['PDB ID'],
                        affinity_kcal=selected_pose_data['Affinity (kcal/mol)'],
                        poses_table_html=poses_html,
                        interactions_table_html=interactions_clean_html,
                        admet_dict=orig_adme,
                        variant_info=variant_dossier_data,
                        plant_photo_b64=plant_photo_b64,
                        paozhi_data=(pz_info if is_processed_state else None)
                    )

                    col_dl, col_sig = st.columns([2, 1], vertical_alignment="center")
                    with col_dl:
                        filename = f"EthnoDock_Report_{row['Common Name'].replace(' ', '_')}_{row['PDB ID']}.html"
                        st.download_button(
                            label=f"📥 Download Regulatory Research Dossier (HTML)",
                            data=dossier_html,
                            file_name=filename,
                            mime="text/html",
                            key=f"dl_dossier_tab2_{idx}"
                        )
                    with col_sig:
                        st.markdown("<div style='text-align:right; font-size:12px; color:#86868B;'>EthnoDock Pro • Verified Simulation</div>", unsafe_allow_html=True)
