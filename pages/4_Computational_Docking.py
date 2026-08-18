import streamlit as st
import os
import sys
import json
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

st.set_page_config(page_title="Computational Docking | EthnoDoc TCM", page_icon="🧬", layout="wide")

st.markdown("""
<style>
    .stApp { background-color: #0E1117; color: #FAFAFA; font-family: 'Inter', sans-serif; }
    .page-title { font-size: 2.5rem; font-weight: 300; color: #52B788; border-bottom: 1px solid #333; padding-bottom: 0.5rem; margin-bottom: 1.5rem; }
    .dock-card { background: #161B22; padding: 1.5rem; border-radius: 8px; border-left: 4px solid #52B788; margin-bottom: 1.5rem; border: 1px solid #30363D; }
    
    .tag-computed { background: #2d6a4f; color: white; padding: 3px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: bold; margin-bottom: 0.5rem; display: inline-block; }
    
    .score-box { background: #1E232A; padding: 1rem; border-radius: 4px; text-align: center; border: 1px solid #2D3748; margin-bottom: 1rem; }
    .score-val { font-size: 2rem; font-weight: bold; color: #52B788; }
    
    .metrics { display: flex; gap: 2rem; margin-top: 1rem; flex-wrap: wrap; }
    .metric-label { font-size: 0.8rem; color: #888; }
    .metric-val { font-size: 1.1rem; font-weight: bold; color: #E2E8F0; }
    
    .lab-cta-box {
        background: linear-gradient(135deg, #1b4332 0%, #2d6a4f 100%);
        padding: 16px 20px;
        border-radius: 8px;
        border: 1px solid #d4af37;
        margin-bottom: 25px;
        display: flex;
        align-items: center;
        justify-content: space-between;
    }
</style>
""", unsafe_allow_html=True)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from ethnodoc_models import Species, HistoricalClaim, Phytochemical, MolecularTarget, DockingExperiment, claim_compound_association
import ethnodock_interaction_engine as inter_eng
import ethnodock_admet_engine as admet_eng
import ethnodock_dossier_engine as dossier_eng

@st.cache_resource
def get_session():
    db_path = f"sqlite:///{os.path.join(BASE_DIR, 'ethnodoc.db')}"
    engine = create_engine(db_path, echo=False)
    Session = sessionmaker(bind=engine)
    return Session()

session = get_session()

if 'selected_species_id' not in st.session_state:
    st.warning("Please search and select a species from the Landing Page first.")
    st.stop()
    
species_id = st.session_state['selected_species_id']
current_species = session.query(Species).filter_by(species_id=species_id).first()

st.markdown("<div class='page-title'>In-Silico AutoDock Vina Docking & Target Profiles</div>", unsafe_allow_html=True)

# CTA to Live Docking Lab
col_cta1, col_cta2 = st.columns([3, 1])
with col_cta1:
    st.markdown("""
    <div style="background:#161B22; padding:12px 18px; border-radius:8px; border-left:4px solid #52B788; border:1px solid #30363D;">
        <span style="font-weight:600; color:#52B788;">⚡ Need on-demand docking or bioisostere design?</span>
        <span style="color:#A0AEC0; font-size:14px; margin-left:8px;">Switch to the 5-Phase Live Docking Lab to dock against any RCSB PDB ID in real-time.</span>
    </div>
    """, unsafe_allow_html=True)
with col_cta2:
    if st.button("🚀 Open Live Docking Lab", use_container_width=True):
        st.switch_page("pages/5_Live_Docking_Lab.py")

st.markdown("<br>", unsafe_allow_html=True)

with st.expander("🎓 Computational Pharmacognosy Guide: Understanding In-Silico Docking"):
    st.markdown("""
    **What is Docking?**
    Molecular docking is finding the optimal binding fit between a TCM natural phytochemical (ligand) and a human/viral disease macromolecule (protein target).
    
    **Interpreting Binding Affinity (ΔG in kcal/mol):**
    - **Stronger binding:** More negative free energy scores (e.g. $-9.5$ kcal/mol represents higher thermodynamic affinity than $-5.0$ kcal/mol).
    - **RMSD:** Evaluates conformational deviation relative to the lowest-energy Mode 1 pose.
    """)

# Get all compounds associated with this species' claims
claims = session.query(HistoricalClaim).filter_by(species_id=species_id).all()
claim_ids = [c.claim_id for c in claims]

if not claim_ids:
    st.info("No claims found for this species.")
    st.stop()

comp_stmt = session.query(Phytochemical).join(
    claim_compound_association, Phytochemical.compound_id == claim_compound_association.c.compound_id
).filter(claim_compound_association.c.claim_id.in_(claim_ids)).all()

# Unique compounds
comps = list({c.compound_id: c for c in comp_stmt}.values())
comp_ids = [c.compound_id for c in comps]

if not comp_ids:
    st.info("No mapped molecules available for docking simulation.")
    st.stop()

# Query docking experiments for these compounds
experiments = session.query(DockingExperiment).filter(DockingExperiment.compound_id.in_(comp_ids)).all()

if not experiments:
    st.info("No offline docking simulations have been indexed for these molecules. Use the Live Docking Lab to compute them live!")
else:
    # Group by compound
    exp_by_comp = {}
    for exp in experiments:
        if exp.compound_id not in exp_by_comp:
            exp_by_comp[exp.compound_id] = []
        exp_by_comp[exp.compound_id].append(exp)
        
    valid_comps = [c for c in comps if c.compound_id in exp_by_comp]
    if valid_comps:
        tabs = st.tabs([f"🌿 {c.compound_name}" for c in valid_comps])
        
        for idx, c in enumerate(valid_comps):
            with tabs[idx]:
                st.markdown(f"### {c.compound_name} Binding Profiles")
                
                # ADMET Summary Box for this compound
                if c.smiles:
                    admet_info = admet_eng.get_admet_profile(c.smiles)
                    if admet_info:
                        col_a, col_b, col_c, col_d = st.columns(4)
                        col_a.metric("Molecular Weight", f"{admet_info['Molecular Weight']} Da")
                        col_b.metric("LogP", f"{admet_info['LogP']}")
                        col_c.metric("TPSA", f"{admet_info['TPSA (Å²)']} Å²")
                        col_d.metric("QED Drug-Likeness", f"{admet_info['QED Drug-Likeness']}")
                
                st.markdown("<br>", unsafe_allow_html=True)
                
                cols = st.columns(2)
                for t_idx, exp in enumerate(exp_by_comp[c.compound_id]):
                    t = exp.target
                    with cols[t_idx % 2]:
                        params = json.loads(exp.parameters_json)
                        grid_sz = f"{params['grid_size'][0]}x{params['grid_size'][1]}x{params['grid_size'][2]}"
                        
                        st.markdown(f"""
                        <div class='dock-card'>
                            <div><span class='tag-computed'>Target: {t.protein_name}</span></div>
                            <h4 style='margin-top: 0; color:#FAFAFA;'>PDB: {t.structure_identifier} | Gene: {t.gene_identifier}</h4>
                            
                            <div class='score-box'>
                                <div class='metric-label'>AutoDock Vina Best Affinity (ΔG)</div>
                                <div class='score-val'>{exp.docking_score} kcal/mol</div>
                            </div>
                        """, unsafe_allow_html=True)
                        
                        # Display 9 Poses Table
                        if exp.poses_json:
                            poses = json.loads(exp.poses_json)
                            md_table = "| Mode | Affinity (kcal/mol) | RMSD l.b. | RMSD u.b. |\n|---|---|---|---|\n"
                            for p in poses:
                                md_table += f"| {p['mode']} | {p['affinity']} | {p['rmsd_lb']} | {p['rmsd_ub']} |\n"
                                
                            st.markdown("**Binding Modes (9 Poses)**")
                            st.markdown(md_table)
                            
                        # 3D Visualizer
                        safe_name = c.compound_name.replace(" ", "_").lower()
                        sdf_path = os.path.join(BASE_DIR, "assets", "structures", f"{safe_name}.sdf")
                        pdb_path = os.path.join(BASE_DIR, "assets", "structures", "targets", f"{t.structure_identifier.lower()}.pdb")
                        
                        st.markdown("**3D Complex Interaction Map (3Dmol.js WebGL)**")
                        if os.path.exists(sdf_path) and os.path.exists(pdb_path):
                            with open(pdb_path, "r", encoding="utf-8") as f: pdb_data = f.read()
                            with open(sdf_path, "r", encoding="utf-8") as f: sdf_data = f.read()
                            
                            viewer_html = inter_eng.build_3dmol_html(
                                container_id=f"viewer_{safe_name}_{t.structure_identifier}",
                                receptor_data=pdb_data,
                                ligand_data=sdf_data,
                                receptor_style='cartoon',
                                ligand_style='stick',
                                height=380
                            )
                            st.components.v1.html(viewer_html, height=400)
                        else:
                            st.info("3D complex structure file ready for real-time visualization.")
                        
                        st.markdown(f"""
                            <div class='metrics'>
                                <div><span class='metric-label'>Software</span><br><span class='metric-val' style='font-size:0.9rem;'>{exp.software_version}</span></div>
                                <div><span class='metric-label'>Grid Size</span><br><span class='metric-val' style='font-size:0.9rem;'>{grid_sz}</span></div>
                                <div><span class='metric-label'>Seed</span><br><span class='metric-val' style='font-size:0.9rem;'>{exp.random_seed}</span></div>
                            </div>
                        </div>
                        """, unsafe_allow_html=True)
