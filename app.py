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
import ethnodock_reproducibility_engine as repro_eng
import ethnodock_chembl_engine as chembl_eng

# --- Page Configuration ---
st.set_page_config(
    page_title="EthnoDock Pro • Computational Pharmacognosy",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# --- Ultra-Premium Apple / Linear Glassmorphic Design System ---
st.markdown(
    """
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500;600&display=swap');

    /* Global Base Reset */
    html, body, [class*="css"], .stApp {
        font-family: 'Plus Jakarta Sans', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif !important;
        background: #050608 !important;
        color: #F1F5F9 !important;
        letter-spacing: -0.015em;
    }

    /* Ambient Aurora Glow */
    .aurora-hero {
        position: relative;
        background: radial-gradient(circle at 50% -20%, rgba(48, 209, 88, 0.22) 0%, rgba(10, 132, 255, 0.12) 40%, rgba(5, 6, 8, 0.98) 75%);
        border: 1px solid rgba(255, 255, 255, 0.1);
        border-radius: 28px;
        padding: 56px 40px 48px 40px;
        text-align: center;
        margin-bottom: 35px;
        box-shadow: 0 35px 70px rgba(0, 0, 0, 0.6);
        overflow: hidden;
    }

    .aurora-hero::before {
        content: "";
        position: absolute;
        top: 0; left: 20%; right: 20%;
        height: 1px;
        background: linear-gradient(90deg, transparent, rgba(48, 209, 88, 0.6), rgba(100, 210, 255, 0.6), transparent);
    }

    /* Floating Pill Badges */
    .pill-badge {
        display: inline-flex;
        align-items: center;
        gap: 8px;
        background: rgba(255, 255, 255, 0.05);
        backdrop-filter: blur(10px);
        border: 1px solid rgba(255, 255, 255, 0.14);
        padding: 6px 18px;
        border-radius: 30px;
        font-size: 0.78rem;
        font-weight: 700;
        color: #30D158;
        letter-spacing: 1px;
        text-transform: uppercase;
        margin-bottom: 20px;
        box-shadow: 0 4px 15px rgba(48, 209, 88, 0.15);
    }

    .hero-main-title {
        font-size: 3.4rem;
        font-weight: 800;
        letter-spacing: -0.04em;
        line-height: 1.12;
        margin: 0 0 18px 0;
        background: linear-gradient(180deg, #FFFFFF 20%, #94A3B8 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }

    .hero-gradient-text {
        background: linear-gradient(135deg, #30D158 0%, #64D2FF 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }

    .hero-main-desc {
        font-size: 1.18rem;
        color: #94A3B8;
        max-width: 780px;
        margin: 0 auto 28px auto;
        line-height: 1.6;
        font-weight: 400;
    }

    /* Interactive Specimen Gallery Strip */
    .gallery-strip {
        display: grid;
        grid-template-columns: repeat(4, 1fr);
        gap: 16px;
        margin-bottom: 35px;
    }
    .gallery-item {
        background: rgba(255, 255, 255, 0.03);
        border: 1px solid rgba(255, 255, 255, 0.08);
        border-radius: 18px;
        padding: 14px;
        text-align: center;
        transition: all 0.3s cubic-bezier(0.16, 1, 0.3, 1);
        cursor: pointer;
    }
    .gallery-item:hover {
        background: rgba(255, 255, 255, 0.06);
        border-color: rgba(48, 209, 88, 0.4);
        transform: translateY(-4px);
        box-shadow: 0 12px 28px rgba(0, 0, 0, 0.4);
    }
    .gallery-img-box {
        height: 140px;
        border-radius: 12px;
        overflow: hidden;
        background: #000;
        margin-bottom: 10px;
    }
    .gallery-img-box img {
        width: 100%;
        height: 100%;
        object-fit: cover;
        transition: transform 0.4s ease;
    }
    .gallery-item:hover .gallery-img-box img {
        transform: scale(1.06);
    }

    /* Command Center Search Box */
    .search-hub-card {
        background: rgba(255, 255, 255, 0.035);
        backdrop-filter: blur(40px);
        -webkit-backdrop-filter: blur(40px);
        border: 1px solid rgba(255, 255, 255, 0.12);
        border-radius: 24px;
        padding: 32px 36px;
        margin-bottom: 40px;
        box-shadow: 0 25px 60px rgba(0, 0, 0, 0.5);
    }

    /* Graphical Connected Stepper Ribbon */
    .visual-stepper-container {
        display: flex;
        justify-content: space-between;
        align-items: center;
        max-width: 900px;
        margin: 0 auto;
        padding: 24px 10px 10px 10px;
        border-top: 1px solid rgba(255, 255, 255, 0.08);
        position: relative;
    }
    .stepper-node {
        display: flex;
        flex-direction: column;
        align-items: center;
        text-align: center;
        background: rgba(255, 255, 255, 0.03);
        border: 1px solid rgba(255, 255, 255, 0.09);
        border-radius: 18px;
        padding: 16px 20px;
        min-width: 160px;
        transition: all 0.3s cubic-bezier(0.16, 1, 0.3, 1);
        box-shadow: 0 8px 20px rgba(0, 0, 0, 0.3);
    }
    .stepper-node:hover {
        background: rgba(255, 255, 255, 0.06);
        border-color: rgba(48, 209, 88, 0.4);
        transform: translateY(-3px);
        box-shadow: 0 12px 25px rgba(48, 209, 88, 0.15);
    }
    .stepper-icon-circle {
        width: 44px;
        height: 44px;
        border-radius: 50%;
        display: flex;
        align-items: center;
        justify-content: center;
        font-size: 1.4rem;
        margin-bottom: 8px;
        background: rgba(255, 255, 255, 0.05);
        border: 1px solid rgba(255, 255, 255, 0.1);
    }
    /* Executive Scientific Manifesto Card */
    .executive-manifesto-card {
        background: linear-gradient(135deg, rgba(255, 255, 255, 0.04) 0%, rgba(255, 255, 255, 0.015) 100%);
        border: 1px solid rgba(255, 255, 255, 0.1);
        border-radius: 22px;
        padding: 26px 30px;
        margin: 20px auto 10px auto;
        max-width: 940px;
        text-align: left;
        box-shadow: 0 20px 45px rgba(0, 0, 0, 0.4);
    }
    .manifesto-header {
        display: flex;
        justify-content: space-between;
        align-items: center;
        border-bottom: 1px solid rgba(255, 255, 255, 0.07);
        padding-bottom: 12px;
        margin-bottom: 18px;
    }
    .manifesto-dot {
        width: 8px;
        height: 8px;
        border-radius: 50%;
        background: #30D158;
        box-shadow: 0 0 10px #30D158;
        display: inline-block;
    }
    .manifesto-tag {
        font-size: 0.76rem;
        font-weight: 700;
        letter-spacing: 1px;
        color: #FFFFFF;
        text-transform: uppercase;
    }
    .manifesto-id {
        font-family: 'JetBrains Mono', monospace;
        font-size: 0.72rem;
        color: #86868B;
    }
    .manifesto-grid {
        display: grid;
        grid-template-columns: 1fr auto 1fr;
        gap: 20px;
        align-items: stretch;
    }
    .manifesto-col {
        background: rgba(255, 255, 255, 0.02);
        border: 1px solid rgba(255, 255, 255, 0.05);
        border-radius: 14px;
        padding: 16px 18px;
        display: flex;
        flex-direction: column;
        justify-content: space-between;
    }
    .manifesto-col-title {
        font-size: 0.95rem;
        font-weight: 700;
        margin-bottom: 8px;
    }
    .manifesto-p {
        font-size: 0.86rem;
        color: #94A3B8;
        line-height: 1.6;
        margin: 0 0 12px 0;
    }
    .manifesto-sub-badge {
        display: inline-block;
        font-size: 0.72rem;
        color: #FFD60A;
        background: rgba(255, 214, 10, 0.08);
        border: 1px solid rgba(255, 214, 10, 0.2);
        padding: 3px 10px;
        border-radius: 12px;
        width: fit-content;
    }
    .manifesto-divider {
        display: flex;
        align-items: center;
        justify-content: center;
        color: rgba(255, 255, 255, 0.2);
        font-size: 1.4rem;
    }
    .manifesto-footer {
        display: flex;
        justify-content: center;
        flex-wrap: wrap;
        gap: 10px;
        margin-top: 18px;
        padding-top: 14px;
        border-top: 1px solid rgba(255, 255, 255, 0.06);
    }
    .manifesto-pill {
        background: rgba(255, 255, 255, 0.03);
        border: 1px solid rgba(255, 255, 255, 0.08);
        padding: 4px 12px;
        border-radius: 20px;
        font-size: 0.75rem;
        color: #CBD5E1;
        font-weight: 500;
    }

    /* Dual Paradigm Synthesis Cards */
    .synthesis-grid {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 24px;
        margin-bottom: 40px;
    }
    .synthesis-card {
        background: rgba(255, 255, 255, 0.025);
        border: 1px solid rgba(255, 255, 255, 0.08);
        border-radius: 20px;
        padding: 28px;
        transition: all 0.25s ease;
    }
    .synthesis-card:hover {
        background: rgba(255, 255, 255, 0.04);
        border-color: rgba(255, 255, 255, 0.15);
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

    /* Segmented Tabs */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
        background: rgba(255, 255, 255, 0.03);
        padding: 6px;
        border-radius: 18px;
        border: 1px solid rgba(255, 255, 255, 0.08);
        margin-bottom: 28px;
    }
    .stTabs [data-baseweb="tab"] {
        border-radius: 12px;
        color: #86868B;
        font-weight: 600;
        padding: 10px 24px;
        border: none;
        font-size: 0.92rem;
    }
    .stTabs [aria-selected="true"] {
        background: #0A84FF !important;
        color: #FFFFFF !important;
        box-shadow: 0 4px 14px rgba(10, 132, 255, 0.3);
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

# --- Initialize View State ---
if 'current_view' not in st.session_state:
    st.session_state['current_view'] = 'landing'

# ==========================================
# 🌟 VIEW 1: DISCOVERY PORTAL & SEARCH
# ==========================================
if st.session_state['current_view'] == 'landing':
    # 1. Radiant Hero Header with Executive Scientific Architecture
    st.markdown(
        """<div class="aurora-hero">
<span class="pill-badge">🌿 ETHNODOCK PRO • COMPUTATIONAL PHARMACOGNOSY</span>
<h1 class="hero-main-title">Bridging 2,000 Years of Botanical Canon<br>with <span class="hero-gradient-text">Modern Structural Biophysics</span></h1>

<div class="executive-manifesto-card">
<div class="manifesto-header">
<div style="display:flex; align-items:center; gap:8px;">
<span class="manifesto-dot"></span>
<span class="manifesto-tag">Executive Scientific Architecture</span>
</div>
<span class="manifesto-id">PROTOCOL // EDK-TCM-2026</span>
</div>

<div class="manifesto-grid">
<div class="manifesto-col">
<div>
<div class="manifesto-col-title" style="color:#FFD60A;">📜 01. Empirical Dynastic Canon</div>
<p class="manifesto-p">
For millennia, canonical pharmacopoeias—from the Han Dynasty's <i>Shennong Bencaojing</i> (神农本草经) to the Ming Dynasty's <i>Bencao Gangmu</i> (本草纲目)—have codified the clinical efficacy of multi-component botanical remedies through centuries of systematic observation.
</p>
</div>
<div class="manifesto-sub-badge">Kew Taxonomy Verified (POWO)</div>
</div>

<div class="manifesto-divider">
<span>⚛️</span>
</div>

<div class="manifesto-col">
<div>
<div class="manifesto-col-title" style="color:#30D158;">⚡ 02. In-Silico Molecular Pharmacology</div>
<p class="manifesto-p">
<b>EthnoDock Pro</b> translates this ancient literature into atomic-scale structural chemistry. We isolate active metabolites, simulate <b>AutoDock Vina empirical free energy (&Delta;G)</b>, map 3D residue anchor contacts (&lt; 4.0 Å), and audit classical <i>Paozhi</i> (炮制) detoxification pathways.
</p>
</div>
<div class="manifesto-sub-badge" style="border-color:rgba(48,209,88,0.3); color:#30D158;">AutoDock Vina v1.2.7 &bull; RCSB PDB</div>
</div>
</div>

<div class="manifesto-footer">
<span class="manifesto-pill">📜 2,000+ Yrs Codified Canon</span>
<span class="manifesto-pill">🎯 RCSB Macromolecular Pockets</span>
<span class="manifesto-pill">⚡ AutoDock Vina Free Energy (&Delta;G)</span>
<span class="manifesto-pill">⚗️ Paozhi Detoxification Alchemy</span>
</div>
</div>
</div>""",
        unsafe_allow_html=True
    )

    # 2. Curated Specimen Showcase Reel
    st.markdown("""
    <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:14px;">
        <div>
            <h3 style="margin:0; font-size:1.25rem; font-weight:700; color:#FFFFFF;">🌿 Curated Botanical Specimen Gallery</h3>
            <p style="margin:2px 0 0 0; font-size:0.85rem; color:#86868B;">Explore high-resolution field photographs of iconic medicinal species in our verified pharmacopeia.</p>
        </div>
        <span class="apple-badge apple-badge-gold">Kew Verified</span>
    </div>
    """, unsafe_allow_html=True)

    col_g1, col_g2, col_g3, col_g4 = st.columns(4, gap="small")
    showcase_herbs = [
        ("Sweet Wormwood (Qinghao)", "Artemisia annua", "Artemisinin", "SARS-CoV-2 Mpro"),
        ("Baikal Skullcap (Huangqin)", "Scutellaria baicalensis", "Baicalein", "COX-2 Kinase"),
        ("Ginseng (Renshen)", "Panax ginseng", "Ginsenoside Rg1", "Estrogen Receptor"),
        ("Red Sage (Danshen)", "Salvia miltiorrhiza", "Tanshinone IIA", "EGFR Kinase")
    ]

    for col, (cname, bname, comp, target) in zip([col_g1, col_g2, col_g3, col_g4], showcase_herbs):
        with col:
            photo_b64 = get_species_photo_b64(cname, bname)
            img_tag = f'<img src="{photo_b64}"/>' if photo_b64 else '<div style="height:140px; background:#111;"></div>'
            st.markdown(f"""
            <div class="gallery-item">
                <div class="gallery-img-box">
                    {img_tag}
                </div>
                <div style="font-weight:700; font-size:13px; color:#FFF;">{cname}</div>
                <div style="font-size:11px; color:#30D158; margin-top:2px;">{comp} &bull; {target}</div>
            </div>
            """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # 3. Interactive Command Search Center
    st.markdown("""
    <div class="search-hub-card">
        <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:16px;">
            <div>
                <h3 style="margin:0; font-size:1.35rem; font-weight:700; color:#FFFFFF;">🔍 Interactive Molecular Discovery Hub</h3>
                <p style="margin:4px 0 0 0; font-size:0.9rem; color:#86868B;">Select any medicinal plant, therapeutic protein target, or phytochemical to inspect its molecular identity.</p>
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
            key="landing_search_by"
        )
    with col_s2:
        if search_by == "Common Name":
            options = df["Common Name"].unique().tolist()
        elif search_by == "Protein Target":
            options = df["Protein Target"].unique().tolist()
        else:
            options = df["Active Phytochemical"].unique().tolist()

        default_idx = 0
        prev_sel = st.session_state.get('selected_herb_name')
        if prev_sel and prev_sel in options:
            default_idx = options.index(prev_sel)

        selected_option = st.selectbox(
            f"Select {search_by}:",
            options,
            index=default_idx,
            on_change=clear_session_docking_state,
            key="landing_selected_option"
        )

    # Immediately lock persistent selection
    st.session_state['selected_search_by'] = search_by
    st.session_state['selected_herb_name'] = selected_option

    st.markdown("</div>", unsafe_allow_html=True)

    # 4. Selected Herb Spotlight Card
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
                    <div style="border-radius:12px; overflow:hidden; height:200px; background:#000;">
                        <img src="{plant_photo_b64}" style="width:100%; height:100%; object-fit:cover;"/>
                    </div>
                    <div style="margin-top:12px;">
                        <span class="apple-badge">Verified Field Specimen</span>
                        <div style="font-weight:700; font-size:15px; margin-top:4px;">{row['Common Name']}</div>
                        <div style="font-size:12px; color:#86868B; font-style:italic;">{row['Botanical Name']}</div>
                    </div>
                </div>
                """, unsafe_allow_html=True)

        with col_p2:
            img_tag = f'<img src="data:image/png;base64,{img_b64}" style="width:100%; height:200px; object-fit:contain;"/>' if img_b64 else '<p>2D Topology</p>'
            st.markdown(f"""
            <div class="apple-card" style="padding:14px; text-align:center; height:100%;">
                <div style="background:#000000; border-radius:12px; padding:6px; height:200px; display:flex; align-items:center; justify-content:center;">
                    {img_tag}
                </div>
                <div style="margin-top:12px;">
                    <span class="apple-badge apple-badge-gold">{row['Chemical Class']}</span>
                    <div style="font-weight:700; font-size:15px; margin-top:4px; color:#52B788;">{row['Active Phytochemical']}</div>
                    <div style="font-size:11px; color:#86868B;">PubChem CID: {row['PubChem CID']}</div>
                </div>
            </div>
            """, unsafe_allow_html=True)

        with col_p3:
            st.markdown(f"""
            <div class="apple-card" style="padding:22px; height:100%; display:flex; flex-direction:column; justify-content:space-between;">
                <div>
                    <span class="apple-badge apple-badge-purple">{row['Dynasty']} • {row['Classical Source']}</span>
                    <div style="font-size:14px; color:#F5F5F7; margin-top:10px; font-style:italic; line-height:1.45;">"{row['English Translation']}"</div>
                    <hr style="border:none; border-top:1px solid rgba(255,255,255,0.06); margin:14px 0;">
                    <div style="font-size:11px; color:#86868B; text-transform:uppercase; letter-spacing:0.5px;">Validated Human Target</div>
                    <div style="font-size:16px; font-weight:700; color:#30D158; margin-top:2px;">{row['Protein Target']} ({row['Gene Symbol']})</div>
                    <div style="font-size:12px; color:#A1A1A6; margin-top:4px;"><b>RCSB PDB ID:</b> <span style="color:#64D2FF; font-weight:700;">{row['PDB ID']}</span> &bull; <b>UniProt:</b> {row['UniProt ID']}</div>
                </div>
            </div>
            """, unsafe_allow_html=True)

        # 🚀 Prominent Workbench Launch Call-To-Action Button Below Spotlight
        st.markdown("<br>", unsafe_allow_html=True)
        col_cta1, col_cta2, col_cta3 = st.columns([1, 2.2, 1])
        with col_cta2:
            launch_workbench = st.button(
                f"🚀 Launch In-Silico Molecular Workbench for {row['Common Name']} →",
                key="btn_launch_workbench_direct",
                use_container_width=True
            )
            if launch_workbench:
                st.session_state['current_view'] = 'workbench'
                st.session_state['selected_search_by'] = search_by
                st.session_state['selected_herb_name'] = selected_option
                st.rerun()

    # 5. Dual Paradigm Synthesis (Ancient Wisdom vs Modern Biophysics)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("""
    <div style="margin-bottom:18px;">
        <span class="apple-badge apple-badge-gold">Dual-Paradigm Architecture</span>
        <h2 style="margin:8px 0 4px 0; font-size:1.8rem; font-weight:700;">Harmonizing Ancient Canon with Biophysical Rigor</h2>
        <p style="margin:0; color:#86868B; font-size:0.95rem;">How EthnoDock validates traditional pharmacognosy using modern computational chemistry.</p>
    </div>

    <div class="synthesis-grid">
        <div class="synthesis-card" style="border-left: 4px solid #FFD60A;">
            <div style="font-size:1.8rem; margin-bottom:10px;">📜</div>
            <h3 style="margin:0 0 6px 0; color:#FFD60A; font-size:1.15rem;">Ancient Ethnobotanical Canon</h3>
            <p style="margin:0; color:#94A3B8; font-size:0.88rem; line-height:1.55;">
                Extracts verbatim clinical indications, thermal properties (四气五味), and formula synergy principles (君臣佐使) codified over two millennia across the Han, Tang, Song, and Ming dynasties.
            </p>
        </div>
        <div class="synthesis-card" style="border-left: 4px solid #30D158;">
            <div style="font-size:1.8rem; margin-bottom:10px;">⚛️</div>
            <h3 style="margin:0 0 6px 0; color:#30D158; font-size:1.15rem;">Modern Structural Pharmacology</h3>
            <p style="margin:0; color:#94A3B8; font-size:0.88rem; line-height:1.55;">
                Employs AutoDock Vina empirical scoring, RDKit 3D conformer minimization (ETKDGv3/UFF), Euclidean contact calculations (<4.0Å), and Lipinski/PAINS toxicological filters.
            </p>
        </div>
    </div>
    """, unsafe_allow_html=True)


# ==========================================
# 🔬 VIEW 2: IN-SILICO MOLECULAR WORKBENCH
# ==========================================
else:
    # Retrieve current search item securely from persistent session state
    search_by = st.session_state.get('selected_search_by', 'Common Name')
    selected_option = st.session_state.get('selected_herb_name', df['Common Name'].iloc[0])
    selected_data = df[df[search_by] == selected_option]
    if selected_data.empty:
        selected_data = df[df['Common Name'] == selected_option]
    if selected_data.empty:
        selected_data = df.head(1)

    if not selected_data.empty:
        row_active = selected_data.iloc[0]
        # Top Navigation Header with Return Button
        col_back, col_title = st.columns([1, 3], vertical_alignment="center")
        with col_back:
            if st.button("← Return to Discovery Portal", use_container_width=True):
                st.session_state['current_view'] = 'landing'
                st.rerun()
        with col_title:
            st.markdown(
                f"<div style='font-size:1.15rem; font-weight:700; color:#F1F5F9; text-align:right;'>"
                f"🔬 In-Silico Molecular Workbench &bull; <span style='color:#30D158;'>{row_active['Common Name']}</span> × <span style='color:#64D2FF;'>{row_active['Protein Target']} ({row_active['PDB ID']})</span>"
                f"</div>",
                unsafe_allow_html=True
            )
        st.markdown("<hr style='border:none; border-top:1px solid rgba(255,255,255,0.08); margin:10px 0 20px 0;'>", unsafe_allow_html=True)

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

                    col_exh, col_seed = st.columns([1.5, 1])
                    with col_exh:
                        exhaustiveness = st.slider("Vina Exhaustiveness (Sampling Precision)", min_value=4, max_value=32, value=8, step=4, key=f"exh_tab2_{idx}")
                    with col_seed:
                        dock_seed = st.number_input("Deterministic Random Seed (--seed)", value=42, step=1, key=f"seed_tab2_{idx}", help="Ensures 100% bit-for-bit exact peer-reviewed replication.")

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
                                    receptor_pdbqt, ligand_pdbqt, [cx, cy, cz], [sx, sy, sz], exhaustiveness=exhaustiveness, seed=dock_seed
                                )
                                if parsed_poses:
                                    st.session_state[f'docking_data_{idx}'] = parsed_poses
                                    st.session_state[f'out_pdbqt_{idx}'] = out_pdbqt
                                    st.session_state[f'docking_done_{idx}'] = True
                                    st.session_state[f'pdb_id_{idx}'] = pdb_id
                                    st.session_state[f'smiles_{idx}'] = smiles
                                    st.session_state[f'uff_delta_{idx}'] = uff_delta
                                    st.session_state[f'dock_seed_{idx}'] = dock_seed
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

                    # ChEMBL Experimental Wet-Lab Ground Truth Benchmark
                    chembl_data = chembl_eng.get_chembl_ground_truth(active_compound_name)
                    if chembl_data:
                        st.markdown(f"""
                        <div class="apple-card-compact" style="border-left: 4px solid #30D158; background: rgba(48, 209, 88, 0.04); margin-top: 14px; margin-bottom: 18px;">
                            <div style="display: flex; justify-content: space-between; align-items: center;">
                                <span style="font-weight: 700; color: #30D158; font-size: 0.92rem;">🧪 ChEMBL Bioassay Experimental Ground Truth Benchmark</span>
                                <span class="apple-badge apple-badge-gold">{chembl_data['chembl_id']}</span>
                            </div>
                            <div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 14px; margin-top: 10px;">
                                <div>
                                    <div style="font-size: 0.72rem; color: #86868B; text-transform: uppercase;">Experimental Wet-Lab IC50</div>
                                    <div style="font-size: 1.2rem; font-weight: 700; color: #FFFFFF;">{chembl_data['experimental_ic50']}</div>
                                </div>
                                <div>
                                    <div style="font-size: 0.72rem; color: #86868B; text-transform: uppercase;">Predicted In-Silico Ki</div>
                                    <div style="font-size: 1.2rem; font-weight: 700; color: #64D2FF;">{selected_pose_data['Estimated Ki']}</div>
                                </div>
                                <div>
                                    <div style="font-size: 0.72rem; color: #86868B; text-transform: uppercase;">Assay Reference</div>
                                    <div style="font-size: 0.86rem; font-weight: 600; color: #FFD60A;">{chembl_data['pubmed_id']}</div>
                                </div>
                            </div>
                            <div style="font-size: 0.8rem; color: #94A3B8; margin-top: 8px; line-height: 1.45;">
                                <b>Biophysical Validation Note:</b> {chembl_data['correlation_notes']}
                            </div>
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

                    # Prepare raw files for Open-Science Reproducibility Package
                    with open(receptor_pdbqt, 'r', encoding='utf-8') as rf:
                        rec_str = rf.read()
                    
                    lig_str = ""
                    ligand_pdbqt_path = os.path.join(BASE_DIR, "active_ligand.pdbqt")
                    if os.path.exists(ligand_pdbqt_path):
                        with open(ligand_pdbqt_path, 'r', encoding='utf-8') as lf:
                            lig_str = lf.read()
                            
                    vina_out_str = ""
                    if out_pdbqt and os.path.exists(out_pdbqt):
                        with open(out_pdbqt, 'r', encoding='utf-8') as vf:
                            vina_out_str = vf.read()

                    interactions_list = interactions_df.to_dict(orient="records") if ('interactions_df' in locals() and not interactions_df.empty) else []

                    repro_zip_bytes = repro_eng.create_reproducibility_zip_bundle(
                        species_name=row['Common Name'],
                        botanical_name=row['Botanical Name'],
                        classical_source=row['Classical Source'],
                        dynasty=row['Dynasty'],
                        target_name=row['Protein Target'],
                        pdb_id=row['PDB ID'],
                        uniprot_id=row['UniProt ID'],
                        compound_name=active_compound_name,
                        smiles=smiles,
                        receptor_pdbqt_str=rec_str,
                        ligand_pdbqt_str=lig_str,
                        out_pdbqt_str=vina_out_str,
                        center=[cx, cy, cz],
                        size=[sx, sy, sz],
                        exhaustiveness=exhaustiveness,
                        seed=st.session_state.get(f'dock_seed_{idx}', 42),
                        binding_affinity=selected_pose_data['Affinity (kcal/mol)'],
                        interactions_summary=interactions_list
                    )

                    col_dl1, col_dl2, col_sig = st.columns([1.2, 1.2, 1], vertical_alignment="center")
                    with col_dl1:
                        filename = f"EthnoDock_Report_{row['Common Name'].replace(' ', '_')}_{row['PDB ID']}.html"
                        st.download_button(
                            label=f"📄 Download Research Dossier (HTML)",
                            data=dossier_html,
                            file_name=filename,
                            mime="text/html",
                            key=f"dl_dossier_tab2_{idx}",
                            use_container_width=True
                        )
                    with col_dl2:
                        filename_zip = f"EthnoDock_Reproducibility_Package_{row['Common Name'].replace(' ', '_')}_{row['PDB ID']}.zip"
                        st.download_button(
                            label=f"📦 Download Open-Science ZIP Bundle",
                            data=repro_zip_bytes,
                            file_name=filename_zip,
                            mime="application/zip",
                            key=f"dl_zip_tab2_{idx}",
                            use_container_width=True
                        )
                    with col_sig:
                        st.markdown("<div style='text-align:right; font-size:12px; color:#86868B;'>EthnoDock Pro • Verified Simulation & BibTeX</div>", unsafe_allow_html=True)
