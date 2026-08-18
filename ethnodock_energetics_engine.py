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
