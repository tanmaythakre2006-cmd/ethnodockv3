import streamlit as st
import os
import sys
import json
import pandas as pd
import numpy as np
import base64
from io import BytesIO
from rdkit import Chem
from rdkit.Chem import Draw

# Add root directory to path
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)

from ethnodoc_models import Species, HistoricalClaim, Phytochemical, HistoricalSource
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import ethnodock_docking_engine as dock_eng
import ethnodock_interaction_engine as inter_eng
import ethnodock_bioisostere_engine as bio_eng
import ethnodock_admet_engine as admet_eng
import ethnodock_dossier_engine as dossier_eng

# --- Page Config ---
st.set_page_config(
    page_title="Live Docking Lab | EthnoDoc TCM",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Thematic Styling (TCM Jade & Dark Slate Theme) ---
st.markdown("""
<style>
    .stApp { background-color: #0E1117; color: #FAFAFA; font-family: 'Inter', sans-serif; }
    
    .lab-header {
        background: linear-gradient(135deg, #1b4332 0%, #2d6a4f 100%);
        padding: 24px 30px;
        border-radius: 12px;
        border-left: 6px solid #d4af37;
        margin-bottom: 25px;
        box-shadow: 0 4px 15px rgba(0,0,0,0.3);
    }
    .lab-title { font-size: 2.2rem; font-weight: 700; color: #FFFFFF; margin: 0; }
    .lab-subtitle { font-size: 1.05rem; color: #d8f3dc; margin-top: 6px; }
    
    .phase-badge {
        background: #2d6a4f;
        color: #d8f3dc;
        padding: 4px 12px;
        border-radius: 14px;
        font-size: 0.8rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        border: 1px solid #40916c;
        display: inline-block;
        margin-bottom: 8px;
    }
    
    .card-box {
        background: #161B22;
        border: 1px solid #30363D;
        border-radius: 10px;
        padding: 20px;
        margin-bottom: 20px;
    }
    
    .metric-card {
        background: #1E232A;
        border: 1px solid #2D3748;
        border-radius: 8px;
        padding: 14px;
        text-align: center;
    }
    .metric-label { font-size: 0.75rem; color: #A0AEC0; text-transform: uppercase; }
    .metric-value { font-size: 1.5rem; font-weight: 700; color: #52B788; margin-top: 4px; }
</style>
""", unsafe_allow_html=True)

# --- Database Connection ---
@st.cache_resource
def get_db_session():
    db_path = f"sqlite:///{os.path.join(BASE_DIR, 'ethnodoc.db')}"
    engine = create_engine(db_path, echo=False)
    Session = sessionmaker(bind=engine)
    return Session()

session = get_db_session()

# --- Helper: Render 2D Molecule ---
def get_mol_image(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=(320, 260))
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            return base64.b64encode(buffered.getvalue()).decode()
    except Exception:
        pass
    return None

# --- Header Banner ---
st.markdown("""
<div class="lab-header">
    <div class="lab-title">🌿 EthnoDock Live In-Silico Workbench</div>
    <div class="lab-subtitle">Integrated Ethnopharmacology • AutoDock Vina Docking • 3D WebGL Interactions • Bioisostere Lead Optimization</div>
</div>
""", unsafe_allow_html=True)

# Predefined High-Value TCM Disease Targets
CURATED_TARGETS = {
    "SARS-CoV-2 Main Protease (Mpro) [PDB: 6LU7]": "6LU7",
    "Human ACE2 Receptor [PDB: 1R42]": "1R42",
    "Epidermal Growth Factor Receptor (EGFR) [PDB: 1M17]": "1M17",
    "Cyclooxygenase-2 (COX-2) [PDB: 5IKR]": "5IKR",
    "HMG-CoA Reductase (Cardiovascular) [PDB: 1HW9]": "1HW9",
    "Protein Kinase CK2 (Anticancer) [PDB: 2ZJW]": "2ZJW",
    "Penicillin-Binding Protein 2a (PBP2a Antibacterial) [PDB: 1VQQ]": "1VQQ",
    "Matrix Metalloproteinase-9 (MMP-9) [PDB: 1L6J]": "1L6J",
    "Estrogen Receptor Alpha [PDB: 1ERE]": "1ERE",
    "Dopamine D2 Receptor [PDB: 6CM4]": "6CM4",
    "Custom RCSB PDB ID": "CUSTOM"
}

# --- Sidebar Controls & Input Selection ---
st.sidebar.markdown("### 🎛️ Experiment Controls")
input_mode = st.sidebar.radio(
    "Input Mode:",
    ["Select from TCM Knowledge Pool", "Custom SMILES / Molecule"]
)

selected_species = None
selected_compound_name = ""
selected_smiles = ""
historical_claim_obj = None

if input_mode == "Select from TCM Knowledge Pool":
    species_list = session.query(Species).all()
    species_options = {f"{s.english_name or 'Herb'} ({s.scientific_name}) — {s.chinese_name or ''}": s for s in species_list}
    
    species_label = st.sidebar.selectbox("Select Botanical Entity:", list(species_options.keys()))
    selected_species = species_options[species_label]
    
    # Get associated compounds
    claims = session.query(HistoricalClaim).filter_by(species_id=selected_species.species_id).all()
    available_compounds = []
    for c in claims:
        for comp in c.compounds:
            if comp.smiles and comp not in available_compounds:
                available_compounds.append(comp)
                
    if available_compounds:
        comp_dict = {f"{c.compound_name} ({c.chemical_class or 'Phytochemical'})": c for c in available_compounds}
        comp_label = st.sidebar.selectbox("Select Active Phytochemical:", list(comp_dict.keys()))
        comp_obj = comp_dict[comp_label]
        selected_compound_name = comp_obj.compound_name
        selected_smiles = comp_obj.smiles
        if claims:
            historical_claim_obj = claims[0]
    else:
        st.sidebar.warning("No compounds with SMILES found for this species.")
        selected_compound_name = "Baicalein"
        selected_smiles = "O=C1C=C(C2=CC=CC=C2)OC3=C1C(O)=C(O)C(O)=C3"
else:
    selected_compound_name = st.sidebar.text_input("Phytochemical Name:", "Artemisinin")
    selected_smiles = st.sidebar.text_input("Canonical SMILES:", "CC1CCC2C(C(=O)OC3(C2C1C)OOC3(C)C)C")

# --- MAIN WORKBENCH LAYOUT ---
col_left, col_right = st.columns([1, 1], gap="large")

# ==========================================
# PHASE 1: BOTANICAL & CHEMICAL PROFILE
# ==========================================
with col_left:
    st.markdown('<span class="phase-badge">Phase 1</span>', unsafe_allow_html=True)
    st.subheader("🌿 Botanical & Chemical Characterization")
    
    with st.container():
        st.markdown(f"### **{selected_compound_name}**")
        if selected_species:
            st.markdown(f"**Botanical Source:** *{selected_species.scientific_name}* ({selected_species.english_name or 'N/A'}) — **{selected_species.chinese_name or ''}**")
            if historical_claim_obj:
                source_title = historical_claim_obj.source.title if historical_claim_obj.source else "Classical Pharmacopoeia"
                st.info(f"📜 **Classical Reference ({source_title}):**\n\"{historical_claim_obj.claim_text}\"\n\n*Translation:* \"{historical_claim_obj.translation}\"")
        
        st.markdown(f"**SMILES:** `{selected_smiles}`")
        
        # 2D Structure
        img_b64 = get_mol_image(selected_smiles)
        if img_b64:
            st.markdown(
                f'<div style="background:#161B22; padding:15px; border-radius:8px; border:1px solid #30363D; text-align:center;">'
                f'<img src="data:image/png;base64,{img_b64}" style="max-width:280px; border-radius:6px;"/>'
                f'<p style="color:#A0AEC0; font-size:12px; margin-top:8px;">2D Molecular Topology: {selected_compound_name}</p>'
                f'</div>',
                unsafe_allow_html=True
            )

# ==========================================
# PHASE 2: TARGET RECEPTOR & CAVITY SETUP
# ==========================================
with col_right:
    st.markdown('<span class="phase-badge">Phase 2</span>', unsafe_allow_html=True)
    st.subheader("🎯 Target Macromolecule & Grid Cavity")
    
    target_choice = st.selectbox("Select Target Protein Receptor:", list(CURATED_TARGETS.keys()))
    pdb_id = CURATED_TARGETS[target_choice]
    
    if pdb_id == "CUSTOM":
        pdb_id = st.text_input("Enter RCSB PDB ID (4 characters):", "6LU7").strip().upper()
        
    st.markdown(f"**Selected Target Structure:** `{pdb_id}`")
    
    # State tracking for receptor initialization
    init_btn = st.button(f"📥 Initialize Receptor {pdb_id} & Auto-Detect Pocket", key="init_rec")
    
    rec_key = f"rec_ready_{pdb_id}"
    if init_btn:
        with st.spinner(f"Fetching PDB {pdb_id} from RCSB and computing binding cavity..."):
            pdb_out = os.path.join(BASE_DIR, f"{pdb_id}.pdb")
            rec_pdbqt = dock_eng.fetch_receptor(pdb_id, pdb_out)
            if rec_pdbqt and os.path.exists(rec_pdbqt):
                center, dims = dock_eng.smart_cavity_finder(pdb_out)
                st.session_state[f"center_{pdb_id}"] = center
                st.session_state[f"dims_{pdb_id}"] = dims
                st.session_state[f"rec_pdbqt_{pdb_id}"] = rec_pdbqt
                st.session_state[rec_key] = True
                st.success(f"Receptor {pdb_id} initialized with AutoDock atom types!")
            else:
                st.error(f"Failed to fetch or parse PDB structure {pdb_id}. Please check internet connection or PDB ID.")

    if st.session_state.get(rec_key, False):
        center_def = st.session_state.get(f"center_{pdb_id}", [0.0, 0.0, 0.0])
        dims_def = st.session_state.get(f"dims_{pdb_id}", [22.0, 22.0, 22.0])
        
        with st.expander("⚙️ Search Grid Box Coordinates (Auto-Detected)", expanded=True):
            col_c1, col_c2, col_c3 = st.columns(3)
            cx = col_c1.number_input("Center X (Å)", value=float(center_def[0]), format="%.3f")
            cy = col_c2.number_input("Center Y (Å)", value=float(center_def[1]), format="%.3f")
            cz = col_c3.number_input("Center Z (Å)", value=float(center_def[2]), format="%.3f")
            
            col_s1, col_s2, col_s3 = st.columns(3)
            sx = col_s1.number_input("Size X (Å)", value=float(dims_def[0]), format="%.3f")
            sy = col_s2.number_input("Size Y (Å)", value=float(dims_def[1]), format="%.3f")
            sz = col_s3.number_input("Size Z (Å)", value=float(dims_def[2]), format="%.3f")
            
            exhaustiveness = st.slider("Vina Exhaustiveness (Thoroughness)", min_value=4, max_value=32, value=8, step=4)

st.markdown("---")

# ==========================================
# PHASE 3: DOCKING SIMULATION & 3D WebGL
# ==========================================
st.markdown('<span class="phase-badge">Phase 3</span>', unsafe_allow_html=True)
st.header("🧬 Docking Biosimulation & 3D WebGL Interaction Map")

if st.session_state.get(rec_key, False):
    dock_btn = st.button("🚀 Execute In-Silico Docking Simulation", key="dock_run_btn")
    dock_done_key = f"dock_done_{selected_compound_name}_{pdb_id}"
    
    if dock_btn:
        with st.spinner(f"Preparing 3D ligand conformation and executing AutoDock Vina for {selected_compound_name} -> {pdb_id}..."):
            rec_pdbqt = st.session_state[f"rec_pdbqt_{pdb_id}"]
            lig_pdbqt_path = os.path.join(BASE_DIR, "active_ligand.pdbqt")
            lig_pdbqt, uff_delta = dock_eng.prepare_ligand(selected_smiles, lig_pdbqt_path)
            
            if lig_pdbqt:
                center = [cx, cy, cz]
                dims = [sx, sy, sz]
                raw_log, parsed_poses, out_pdbqt = dock_eng.run_vina_docking(
                    rec_pdbqt, lig_pdbqt, center, dims, exhaustiveness=exhaustiveness
                )
                
                if parsed_poses:
                    st.session_state[f"poses_{dock_done_key}"] = parsed_poses
                    st.session_state[f"out_pdbqt_{dock_done_key}"] = out_pdbqt
                    st.session_state[f"uff_{dock_done_key}"] = uff_delta
                    st.session_state[dock_done_key] = True
                    st.success("AutoDock Vina simulation converged successfully!")
                else:
                    st.warning("AutoDock Vina completed with standard output log:")
                    st.code(raw_log)
            else:
                st.error("Failed to generate 3D conformer for ligand SMILES.")

    if st.session_state.get(dock_done_key, False):
        poses = st.session_state.get(f"poses_{dock_done_key}", [])
        out_pdbqt = st.session_state.get(f"out_pdbqt_{dock_done_key}", "")
        uff_delta = st.session_state.get(f"uff_{dock_done_key}", 0.0)
        
        # 1. Conformational Pose Matrix & Metrics
        col_m1, col_m2, col_m3, col_m4 = st.columns(4)
        best_affinity = poses[0]["affinity"] if poses else 0.0
        
        col_m1.markdown(f'<div class="metric-card"><div class="metric-label">Best Binding Affinity</div><div class="metric-value">{best_affinity} kcal/mol</div></div>', unsafe_allow_html=True)
        col_m2.markdown(f'<div class="metric-card"><div class="metric-label">Conformational Modes</div><div class="metric-value">{len(poses)} Poses</div></div>', unsafe_allow_html=True)
        col_m3.markdown(f'<div class="metric-card"><div class="metric-label">UFF Strain Delta (ΔE)</div><div class="metric-value">{uff_delta:.2f} kcal/mol</div></div>', unsafe_allow_html=True)
        col_m4.markdown(f'<div class="metric-card"><div class="metric-label">Target Complex</div><div class="metric-value" style="font-size:1.1rem; color:#A0AEC0;">{pdb_id}</div></div>', unsafe_allow_html=True)
        
        st.markdown("<br>", unsafe_allow_html=True)
        
        # 2. Pose Selector
        pose_options = [f"Mode {p['mode']} (Affinity: {p['affinity']} kcal/mol | RMSD: {p['rmsd_ub']} Å)" for p in poses]
        selected_pose_str = st.selectbox("Select Conformation to Inspect:", pose_options)
        selected_mode_idx = pose_options.index(selected_pose_str)
        selected_pose_data = poses[selected_mode_idx]
        
        # Extract 3D Coordinates for selected pose
        extracted_poses = inter_eng.extract_poses(out_pdbqt)
        if extracted_poses and selected_mode_idx < len(extracted_poses):
            pose_coord_str = extracted_poses[selected_mode_idx]
            rec_pdbqt_path = st.session_state[f"rec_pdbqt_{pdb_id}"]
            
            # Compute Non-covalent Interactions (< 4.0 Å)
            interactions_df = inter_eng.calc_interactions(pose_coord_str, rec_pdbqt_path, cutoff=4.0)
            
            col_v1, col_v2 = st.columns([3, 2])
            
            with col_v1:
                st.markdown("#### 🔮 Interactive 3D Complex (3Dmol.js WebGL)")
                style_c1, style_c2, style_c3 = st.columns(3)
                rec_style = style_c1.selectbox("Receptor Ribbon:", ["cartoon", "stick", "sphere", "line"], index=0)
                lig_style = style_c2.selectbox("Ligand Style:", ["stick", "sphere", "cross"], index=0)
                surf_toggle = style_c3.checkbox("Show Pocket Cavity Mesh", value=False)
                
                with open(rec_pdbqt_path, 'r', encoding='utf-8') as rf:
                    receptor_str = rf.read()
                    
                viewer_html = inter_eng.build_3dmol_html(
                    container_id=f"viewer_{selected_compound_name}_{pdb_id}",
                    receptor_data=receptor_str,
                    ligand_data=pose_coord_str,
                    interactions_df=interactions_df,
                    receptor_style=rec_style,
                    ligand_style=lig_style,
                    show_surface=surf_toggle,
                    height=480
                )
                st.components.v1.html(viewer_html, height=500)
                
            with col_v2:
                st.markdown("#### 🕸️ Intermolecular Contact Residues (< 4.0 Å)")
                if not interactions_df.empty:
                    display_int_df = interactions_df[["Receptor Residue", "Distance (Å)", "Interaction Type"]]
                    st.dataframe(display_int_df, hide_index=True, use_container_width=True)
                    h_bonds = interactions_df[interactions_df["Interaction Type"] == "Hydrogen / Polar Bond"]
                    st.success(f"**Key Anchors:** {len(h_bonds)} Polar/Hydrogen Bonds identified in the binding pocket.")
                else:
                    st.info("No close contacts (< 4.0 Å) detected for this specific orientation.")
else:
    st.info("👆 Please initialize the Target Receptor in Phase 2 to unlock the docking engine.")

st.markdown("---")

# ==========================================
# PHASE 4: BIOISOSTERE LEAD OPTIMIZATION
# ==========================================
st.markdown('<span class="phase-badge">Phase 4</span>', unsafe_allow_html=True)
st.header("🧬 Semi-Synthetic Bioisostere Lead Optimization")
st.caption("Apply rational medicinal chemistry bioisosteric transformations to natural TCM phytochemical scaffolds.")

derivatives = bio_eng.generate_tcm_derivatives(selected_smiles)

if derivatives:
    var_options = [f"{d['name']} — {d['variant_smiles'][:35]}..." for d in derivatives]
    selected_var_str = st.selectbox("Select Rationally Designed Derivative:", var_options)
    selected_var_idx = var_options.index(selected_var_str)
    chosen_variant = derivatives[selected_var_idx]
    
    st.markdown(f"**Medicinal Chemistry Rationale:** {chosen_variant['rationale']}")
    st.markdown(f"**Derivative SMILES:** `{chosen_variant['variant_smiles']}`")
    
    col_var_dock1, col_var_dock2 = st.columns([1, 2])
    with col_var_dock1:
        dock_var_btn = st.button("⚡ Dock Optimized Derivative", key="dock_var_btn")
        
    var_dock_key = f"var_dock_{chosen_variant['variant_smiles']}_{pdb_id}"
    if dock_var_btn and st.session_state.get(rec_key, False):
        with st.spinner("Preparing variant 3D conformer and executing comparative docking..."):
            rec_pdbqt = st.session_state[f"rec_pdbqt_{pdb_id}"]
            var_pdbqt_path = os.path.join(BASE_DIR, "variant_ligand.pdbqt")
            var_pdbqt, _ = dock_eng.prepare_ligand(chosen_variant['variant_smiles'], var_pdbqt_path)
            
            if var_pdbqt:
                center = st.session_state.get(f"center_{pdb_id}", [0.0, 0.0, 0.0])
                dims = st.session_state.get(f"dims_{pdb_id}", [22.0, 22.0, 22.0])
                raw_log, parsed_var_poses, var_out_pdbqt = dock_eng.run_vina_docking(
                    rec_pdbqt, var_pdbqt, center, dims, exhaustiveness=8
                )
                if parsed_var_poses:
                    st.session_state[f"var_poses_{var_dock_key}"] = parsed_var_poses
                    st.session_state[var_dock_key] = True
                    st.success("Derivative docking complete!")
                    
    if st.session_state.get(var_dock_key, False):
        var_poses = st.session_state.get(f"var_poses_{var_dock_key}", [])
        if var_poses:
            var_best_aff = var_poses[0]["affinity"]
            parent_aff = st.session_state.get(f"poses_{dock_done_key}", [{"affinity": 0.0}])[0]["affinity"] if st.session_state.get(dock_done_key, False) else 0.0
            
            delta_affinity = var_best_aff - parent_aff
            col_c1, col_c2 = st.columns(2)
            col_c1.metric("Derivative Affinity", f"{var_best_aff} kcal/mol", delta=f"{-delta_affinity:+.2f} kcal/mol" if parent_aff else None)
            if delta_affinity < 0:
                col_c2.success(f"✨ **Lead Optimization Success:** The derivative enhances binding affinity by {abs(delta_affinity):.2f} kcal/mol over the natural compound!")
            else:
                col_c2.info("The derivative shows comparable binding stability.")
else:
    st.info("No standard bioisosteric transformation rules triggered for this chemical structure.")

st.markdown("---")

# ==========================================
# PHASE 5: ADMET & SCIENTIFIC DOSSIER EXPORT
# ==========================================
st.markdown('<span class="phase-badge">Phase 5</span>', unsafe_allow_html=True)
st.header("🛡️ ADMET Pharmacokinetics & Scientific Dossier Export")

admet_data = admet_eng.get_admet_profile(selected_smiles)

if admet_data:
    col_a1, col_a2, col_a3, col_a4 = st.columns(4)
    col_a1.metric("Molecular Weight", f"{admet_data['Molecular Weight']} Da")
    col_a2.metric("LogP (Lipophilicity)", f"{admet_data['LogP']}")
    col_a3.metric("TPSA (Polar Surface)", f"{admet_data['TPSA (Å²)']} Å²")
    col_a4.metric("QED Drug-Likeness", f"{admet_data['QED Drug-Likeness']}")
    
    st.markdown("#### 🧪 Toxicological & Assay Interference Screening (PAINS)")
    if admet_data["PAINS Clean"]:
        st.success(f"✅ **PAINS Screen Passed:** {admet_data['PAINS Screen']}")
    else:
        st.warning(f"⚠️ **PAINS Alert Identified:** {admet_data['PAINS Screen']}")
        
    st.markdown(f"**Lipinski Rule of 5:** {admet_data['Lipinski Violations']} Violations ({', '.join(admet_data['Lipinski Details'])})")

# Dossier HTML Exporter
st.markdown("#### 📄 Automated Scientific Research Dossier")

species_eng = selected_species.english_name if selected_species else selected_compound_name
species_sci = selected_species.scientific_name if selected_species else "Synthetic/Isolated"
species_chi = selected_species.chinese_name if selected_species else ""
claim_text = historical_claim_obj.claim_text if historical_claim_obj else "Classical TCM Phytochemical Reference"
trans_text = historical_claim_obj.translation if historical_claim_obj else "Validated in-silico screening analysis"
src_title = historical_claim_obj.source.title if (historical_claim_obj and historical_claim_obj.source) else "Bencao Classical Corpus"

# Build Poses HTML table
poses_table_html = "<p>Simulation pending.</p>"
if st.session_state.get(dock_done_key, False):
    poses_list = st.session_state.get(f"poses_{dock_done_key}", [])
    poses_df = pd.DataFrame(poses_list)
    poses_table_html = poses_df.to_html(index=False)

# Build Interactions HTML table
interactions_table_html = "<p>No interaction table calculated.</p>"
if 'interactions_df' in locals() and not interactions_df.empty:
    interactions_table_html = interactions_df[["Receptor Residue", "Distance (Å)", "Interaction Type"]].to_html(index=False)

variant_dossier_info = None
if 'chosen_variant' in locals():
    variant_dossier_info = {
        "name": chosen_variant["name"],
        "rationale": chosen_variant["rationale"],
        "smiles": chosen_variant["variant_smiles"],
        "affinity": locals().get('var_best_aff', 'N/A')
    }

dossier_html = dossier_eng.generate_tcm_dossier_html(
    species_name=species_eng,
    scientific_name=species_sci,
    chinese_name=species_chi,
    source_title=src_title,
    claim_text=claim_text,
    translation=trans_text,
    compound_name=selected_compound_name,
    smiles=selected_smiles,
    target_name=target_choice,
    pdb_id=pdb_id,
    affinity_kcal=locals().get('best_affinity', 'N/A'),
    poses_table_html=poses_table_html,
    interactions_table_html=interactions_table_html,
    admet_dict=admet_data,
    variant_info=variant_dossier_info
)

st.download_button(
    label=f"📥 Download Comprehensive Scientific Dossier ({selected_compound_name}_{pdb_id}.html)",
    data=dossier_html,
    file_name=f"Ethnodock_Dossier_{selected_compound_name}_{pdb_id}.html",
    mime="text/html",
    key="dl_dossier_btn"
)
