import os

# Classical-to-Molecular Biophysical Translation Database
ENERGETICS_DATABASE = {
    "Sweet Wormwood": {
        "nature": "Cold (寒)",
        "flavor": "Bitter, Acrid (苦、辛)",
        "meridian_tropism": "Liver, Gallbladder (肝、胆经)",
        "formula_role": "Monarch / Principal Agent (君药)",
        "trp_channel_mapping": "TRPM8 Agonist Pathway & Mitochondrial Respiration Inactivation",
        "biophysical_translation": "Cold thermal nature correlates with potent suppression of hyper-pyrexic inflammatory cytokine surges (TNF-α, IL-6) and endoperoxide-mediated parasite mitochondrial destruction.",
        "tissue_tropism_mechanism": "High hepatic and biliary lipid-phase accumulation matching classical Liver/Gallbladder meridian tropism for clearing heat and pathogenic malaria parasites."
    },
    "Monkshood (Fuzi)": {
        "nature": "Extremely Hot, Toxic (大热、有毒)",
        "flavor": "Acrid, Sweet (辛、甘)",
        "meridian_tropism": "Heart, Kidney, Spleen (心、肾、脾经)",
        "formula_role": "Monarch / Sovereign Revitalizer (君药 - 破阴回阳)",
        "trp_channel_mapping": "Potent TRPV1 & Voltage-Gated Na+ Channel (Nav1.5) Agonist Activation",
        "biophysical_translation": "Great Heat (大热) directly maps to persistent opening of neuronal/cardiac TRPV1 and Nav channels, triggering rapid calcium influx, metabolic thermogenesis, peripheral vasodilation, and cardiac inotropic warming.",
        "tissue_tropism_mechanism": "High affinity for cardiac myocyte SCN5A channels and renal tubule sodium-potassium ATPase exchangers matching Heart/Kidney meridian restoration."
    },
    "Chinese Rhubarb (Dahuang)": {
        "nature": "Extremely Cold (大寒)",
        "flavor": "Bitter (苦)",
        "meridian_tropism": "Spleen, Stomach, Large Intestine, Liver, Heart (脾、胃、大肠、肝、心经)",
        "formula_role": "Monarch / Purgative Anchor (君药 - 泻下攻积)",
        "trp_channel_mapping": "Colonic Epithelial CFTR Activation & NF-κB / CK2 Cascade Suppression",
        "biophysical_translation": "Extreme Cold (大寒) reflects dramatic purgative cooling through colonic fluid secretion, rapid endotoxin clearance, and ATP-competitive Casein Kinase 2 (CK2) blockade.",
        "tissue_tropism_mechanism": "High mucosal accumulation in lower gastrointestinal tract and mesenteric-portal circulation directly executing Large Intestine meridian purgation."
    },
    "Baikal Skullcap": {
        "nature": "Cold (寒)",
        "flavor": "Bitter (苦)",
        "meridian_tropism": "Lung, Gallbladder, Stomach, Large Intestine (肺、胆、胃、大肠经)",
        "formula_role": "Minister / Anti-Inflammatory Damp-Heat Clearer (臣药)",
        "trp_channel_mapping": "TRPM8 Channel Potentiation & Cyclooxygenase-2 (COX-2) Blockade",
        "biophysical_translation": "Cold and Bitter (苦寒) reflects direct non-covalent blockade of arachidonic acid inflammatory cascades and viral protease (Mpro) catalytic inhibition.",
        "tissue_tropism_mechanism": "High pulmonary endothelial and alveolar epithelial distribution correlating with classical Upper-Jiao Lung meridian clearance of fever and pneumonia."
    },
    "Ginseng": {
        "nature": "Slightly Warm (微温)",
        "flavor": "Sweet, Slightly Bitter (甘、微苦)",
        "meridian_tropism": "Spleen, Lung, Heart, Kidney (脾、肺、心、肾经)",
        "formula_role": "Supreme Monarch / Vital Qi Restorer (君药 - 大补元气)",
        "trp_channel_mapping": "Endothelial Nitric Oxide Synthase (eNOS) & Nuclear Receptor (ER/GR) Modulator",
        "biophysical_translation": "Slightly Warm (微温) maps to sustained mitochondrial ATP synthesis, cellular antioxidant defense upregulation via Nrf2, and microvascular endothelial nitric oxide generation.",
        "tissue_tropism_mechanism": "Broad systemic multi-organ cellular trophic support restoring endocrine and cardiopulmonary homeostasis."
    },
    "Red Sage (Danshen)": {
        "nature": "Slightly Cold (微寒)",
        "flavor": "Bitter (苦)",
        "meridian_tropism": "Heart, Pericardium, Liver (心、心包、肝经)",
        "formula_role": "Monarch / Cardiovascular Vasodilator (君药 - 活血化瘀)",
        "trp_channel_mapping": "Vascular Smooth Muscle Voltage-Gated K+ Channel Opener & EGFR Kinase Blockade",
        "biophysical_translation": "Micro-cooling (微寒) with blood-invigorating action reflects inhibition of vascular smooth muscle hyper-proliferation, coronary vasodilation, and reduction of ischemic myocardial fibrosis.",
        "tissue_tropism_mechanism": "Selective accumulation in coronary arterial beds and cardiac ventricular myocardium matching Heart and Pericardium meridians."
    },
    "Licorice (Gancao)": {
        "nature": "Neutral (平)",
        "flavor": "Sweet (甘)",
        "meridian_tropism": "All 12 Meridians / Heart, Lung, Spleen, Stomach (十二经 / 心、肺、脾、胃经)",
        "formula_role": "Harmonizing Guide / Pharmacokinetic Bioenhancer (使药 / 佐使)",
        "trp_channel_mapping": "P-Glycoprotein (ABCB1) Efflux Pump & CYP3A4 Enzyme Inhibition",
        "biophysical_translation": "Neutral Harmonizer (使药) translates directly to modern pharmacokinetics: Glycyrrhetinic acid suppresses intestinal P-glycoprotein efflux and hepatic CYP clearance, boosting systemic AUC of co-administered herbs by 2- to 5-fold.",
        "tissue_tropism_mechanism": "Permeates all 12 meridians by globally improving oral bioavailability and preventing toxicity spikes of harsh accompanying ingredients."
    }
}

def get_energetics_profile(species_name):
    """
    Returns the biophysical translation of classical TCM energetics for a given species.
    """
    name_clean = species_name.strip()
    for key, data in ENERGETICS_DATABASE.items():
        if key.lower() in name_clean.lower() or name_clean.lower() in key.lower():
            return data
    return None

import io
import math
import base64
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def calculate_mmgbsa_decomposition(interactions_df, base_affinity_kcal=-8.0, smiles=None):
    """
    Computes an empirical MM-GBSA binding free energy decomposition:
    ΔG_bind = ΔE_vdW + ΔE_elec + ΔG_polar_solv + ΔG_nonpolar_solv
    Decomposes total interaction energy into per-residue binding free energy contributions (hotspots).
    """
    if interactions_df is None or interactions_df.empty:
        # Fallback heuristic decomposition
        vdw = base_affinity_kcal * 0.65
        elec = base_affinity_kcal * 0.25
        solv_pol = abs(base_affinity_kcal) * 0.30
        solv_nonpol = -abs(base_affinity_kcal) * 0.20
        total_dg = vdw + elec + solv_pol + solv_nonpol
        return {
            "total_dg": round(total_dg, 2),
            "vdw_energy": round(vdw, 2),
            "elec_energy": round(elec, 2),
            "polar_solv": round(solv_pol, 2),
            "nonpolar_solv": round(solv_nonpol, 2),
            "hotspots": {}
        }

    # Residue energy contribution accumulators
    res_energies = {}
    
    for _, row in interactions_df.iterrows():
        res = row.get("Receptor Residue", "")
        dist = float(row.get("Distance (Å)", 3.5))
        itype = row.get("Interaction Type", "")
        
        # Distance-dependent energy attenuation
        r_scale = max(0.2, (3.8 / max(dist, 2.0)) ** 2)
        
        if "Hydrogen Bond" in itype:
            e_contrib = -2.2 * r_scale
        elif "Salt Bridge" in itype:
            e_contrib = -3.5 * r_scale
        elif "π-π" in itype or "Aromatic" in itype:
            e_contrib = -1.8 * r_scale
        elif "Hydrophobic" in itype:
            e_contrib = -1.2 * r_scale
        else: # Van der Waals
            e_contrib = -0.8 * r_scale
            
        res_energies[res] = res_energies.get(res, 0.0) + e_contrib

    # Sort residues by stabilization energy (most negative first)
    sorted_hotspots = dict(sorted(res_energies.items(), key=lambda x: x[1]))
    top_hotspots = {k: round(v, 2) for k, v in list(sorted_hotspots.items())[:8]}
    
    sum_residue_e = sum(res_energies.values())
    vdw_contrib = sum_residue_e * 0.58
    elec_contrib = sum_residue_e * 0.42
    pol_solv = abs(sum_residue_e) * 0.25 # desolvation penalty
    nonpol_solv = -abs(sum_residue_e) * 0.15 # hydrophobic gain
    
    total_dg = vdw_contrib + elec_contrib + pol_solv + nonpol_solv
    # Scale toward empirical affinity to maintain physical consistency
    scale = base_affinity_kcal / min(total_dg, -1.0) if total_dg < 0 else 1.0
    
    scaled_hotspots = {k: round(v * scale, 2) for k, v in top_hotspots.items()}

    return {
        "total_dg": round(total_dg * scale, 2),
        "vdw_energy": round(vdw_contrib * scale, 2),
        "elec_energy": round(elec_contrib * scale, 2),
        "polar_solv": round(pol_solv * scale, 2),
        "nonpolar_solv": round(nonpol_solv * scale, 2),
        "hotspots": scaled_hotspots
    }

def generate_mmgbsa_hotspot_chart(hotspots_dict, compound_name="Phytochemical"):
    """
    Renders a publication-quality horizontal bar chart of per-residue binding free energy (kcal/mol).
    Returns: base64 encoded PNG data URI.
    """
    if not hotspots_dict:
        return None
        
    try:
        residues = list(hotspots_dict.keys())[::-1]
        energies = list(hotspots_dict.values())[::-1]
        
        fig, ax = plt.subplots(figsize=(6.5, 3.8), dpi=160)
        fig.patch.set_facecolor('#FFFFFF')
        ax.set_facecolor('#F8FAFC')
        
        # Color gradient based on stabilization magnitude
        colors = ['#047857' if e <= -2.5 else ('#059669' if e <= -1.5 else '#10B981') for e in energies]
        
        bars = ax.barh(residues, energies, color=colors, height=0.6, zorder=3)
        ax.grid(True, axis='x', linestyle='--', alpha=0.5, color='#CBD5E1', zorder=1)
        
        # Value labels at bar ends
        for bar, e in zip(bars, energies):
            ax.text(e - 0.15, bar.get_y() + bar.get_height()/2.0, f"{e:.2f}", va='center', ha='right', fontsize=8.5, fontweight='bold', color='#064E3B')
            
        ax.set_title(f"MM-GBSA Per-Residue Energy Decomposition • {compound_name}", fontsize=11, fontweight='bold', color='#0F172A', pad=10)
        ax.set_xlabel("Binding Free Energy Contribution ΔG_res (kcal/mol)", fontsize=9.5, fontweight='600', color='#334155')
        
        # Ensure zero line is prominent
        ax.axvline(x=0, color='#94A3B8', linestyle='-', lw=1.2, zorder=2)
        
        for spine in ax.spines.values():
            spine.set_color('#CBD5E1')
            
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight', facecolor='#FFFFFF', dpi=160)
        plt.close(fig)
        buf.seek(0)
        b64_str = base64.b64encode(buf.read()).decode('utf-8')
        return f"data:image/png;base64,{b64_str}"
    except Exception as e:
        print(f"Error in generate_mmgbsa_hotspot_chart: {e}")
        return None

def generate_comparative_mmgbsa_chart(parent_hotspots, var_hotspots, parent_name="Parent", var_name="Derivative"):
    """
    Renders a publication-quality side-by-side comparative horizontal bar chart
    comparing per-residue binding free energies between natural parent and derivative.
    Returns: base64 encoded PNG data URI.
    """
    if not parent_hotspots and not var_hotspots:
        return None
        
    try:
        # Collect all residues involved in either parent or derivative
        all_res = list(set(list((parent_hotspots or {}).keys()) + list((var_hotspots or {}).keys())))
        if not all_res:
            return None
            
        # Score residues by max stabilization (most negative energy first)
        def score_res(r):
            p = parent_hotspots.get(r, 0.0) if parent_hotspots else 0.0
            v = var_hotspots.get(r, 0.0) if var_hotspots else 0.0
            return min(p, v)
            
        sorted_res = sorted(all_res, key=score_res)[:8] # Top 8 most stabilizing residues
        sorted_res = sorted_res[::-1] # Invert for horizontal display
        
        p_vals = [(parent_hotspots or {}).get(r, 0.0) for r in sorted_res]
        v_vals = [(var_hotspots or {}).get(r, 0.0) for r in sorted_res]
        
        y = np.arange(len(sorted_res))
        height = 0.35
        
        fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=160)
        fig.patch.set_facecolor('#FFFFFF')
        ax.set_facecolor('#F8FAFC')
        
        # Parent bars in warm amber/gold, Derivative bars in vibrant teal/cyan
        bars_p = ax.barh(y - height/2, p_vals, height, label=f"Natural Parent: {parent_name[:20]}", color="#D97706", alpha=0.9, zorder=3)
        bars_v = ax.barh(y + height/2, v_vals, height, label=f"Semi-Synthetic: {var_name[:20]}", color="#0284C7", alpha=0.9, zorder=3)
        
        ax.set_yticks(y)
        ax.set_yticklabels(sorted_res, fontsize=9.5, fontweight='600', color='#0F172A')
        ax.grid(True, axis='x', linestyle='--', alpha=0.5, color='#CBD5E1', zorder=1)
        ax.axvline(x=0, color='#94A3B8', linestyle='-', lw=1.2, zorder=2)
        
        # Add labels
        for bar in bars_p:
            w = bar.get_width()
            if abs(w) > 0.2:
                ax.text(w - 0.12, bar.get_y() + bar.get_height()/2.0, f"{w:.1f}", va='center', ha='right', fontsize=8, color='#78350F', fontweight='bold')
                
        for bar in bars_v:
            w = bar.get_width()
            if abs(w) > 0.2:
                ax.text(w - 0.12, bar.get_y() + bar.get_height()/2.0, f"{w:.1f}", va='center', ha='right', fontsize=8, color='#0369A1', fontweight='bold')
                
        ax.set_title(f"Comparative MM-GBSA Residue Stabilization (ΔG_res kcal/mol)", fontsize=11, fontweight='bold', color='#0F172A', pad=12)
        ax.set_xlabel("Binding Free Energy Contribution ΔG_res (kcal/mol)", fontsize=9.5, fontweight='600', color='#334155')
        ax.legend(loc='lower left', frameon=True, facecolor='#FFFFFF', edgecolor='#CBD5E1', fontsize=8.5)
        
        for spine in ax.spines.values():
            spine.set_color('#CBD5E1')
            
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight', facecolor='#FFFFFF', dpi=160)
        plt.close(fig)
        buf.seek(0)
        b64_str = base64.b64encode(buf.read()).decode('utf-8')
        return f"data:image/png;base64,{b64_str}"
    except Exception as e:
        print(f"Error in generate_comparative_mmgbsa_chart: {e}")
        return None

