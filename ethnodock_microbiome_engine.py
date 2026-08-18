import os
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, AllChem

# Human Gut Microbiota Biotransformation & In-Vivo Pharmacokinetics Library
MICROBIOME_TRANSFORMATIONS = {
    "Sennoside A": {
        "botanical_source": "Chinese Rhubarb (Rheum officinale / Dahuang)",
        "ingested_scaffold": "Dianthrone Diglucoside (Sennoside A)",
        "ingested_smiles": "O=C1C2=C(C(=O)C3=C1C=CC(=C3)C4=CC=C5C(=C4)C(=O)C6=C(C5=O)C=C(O)C=C6)C=C(O)C=C2",
        "bacterial_enzyme": "Colonic Bacterial β-Glucosidase & Anthrone Reductase (Bacteroides / Bifidobacterium)",
        "metabolic_reaction": "Deglycosylation & Dianthrone Cleavage ──► Free Emodin & Rhein Aglycones",
        "circulating_metabolite": "Emodin (Free Anthraquinone Aglycone)",
        "circulating_smiles": "CC1=CC2=C(C(=C1)O)C(=O)C3=C(C2=O)C=C(C=C3O)O",
        "permeability_raw_papp": "0.18 x 10^-6 cm/s (Impermeable Prodrug)",
        "permeability_active_papp": "18.4 x 10^-6 cm/s (High Human Intestinal Permeability)",
        "in_vivo_pk_note": "Intact sennosides pass through the small intestine unabsorbed to reach the colon, where commensal microflora cleave the sugar bonds to release active therapeutic emodin directly at target mucosa.",
        "phase2_hepatic_fate": "Hepatic Glucuronidation (Emodin-3-O-β-D-glucuronide) with biliary enterohepatic recirculation."
    },
    "Baicalin": {
        "botanical_source": "Baikal Skullcap (Scutellaria baicalensis / Huangqin)",
        "ingested_scaffold": "Flavone 7-O-Glucuronide (Baicalin)",
        "ingested_smiles": "O=C1C=C(C2=CC=CC=C2)OC3=C1C(O)=C(O)C=C3OC4OC(C(O)C(O)C4O)C(O)=O",
        "bacterial_enzyme": "Intestinal Bacterial β-Glucuronidase (E. coli / Clostridium spp.)",
        "metabolic_reaction": "Enzymatic Hydrolysis of 7-O-Glucuronide ──► Free Baicalein Aglycone",
        "circulating_metabolite": "Baicalein (5,6,7-Trihydroxyflavone)",
        "circulating_smiles": "O=C1C=C(C2=CC=CC=C2)OC3=C1C(O)=C(O)C(O)=C3",
        "permeability_raw_papp": "0.42 x 10^-6 cm/s (Poor Oral Bioavailability)",
        "permeability_active_papp": "22.6 x 10^-6 cm/s (Rapid Passive Intestinal Absorption)",
        "in_vivo_pk_note": "Baicalin is hydrolyzed in the cecum into lipophilic Baicalein, which is rapidly absorbed across the intestinal wall into mesenteric circulation.",
        "phase2_hepatic_fate": "UDP-glucuronosyltransferase (UGT1A8/1A9) re-glucuronidation maintaining prolonged systemic anti-inflammatory half-life (T1/2 ~ 11.2 h)."
    },
    "Glycyrrhizin": {
        "botanical_source": "Licorice (Glycyrrhiza uralensis / Gancao)",
        "ingested_scaffold": "Triterpenoid Saponin Diglucuronide (Glycyrrhizin)",
        "ingested_smiles": "CC1(C)C2CCC3(C)C(=O)C=C4C(CCC5(C)C4CCC5(C)C(=O)O)C3(C)C2CCC1OC6OC(C(O)C(O)C6O)C(=O)O",
        "bacterial_enzyme": "Colonic Eubacterium sp. GLH β-Glucuronidase",
        "metabolic_reaction": "Stepwise Cleavage of 2 Glucuronide Sugars ──► 18β-Glycyrrhetinic Acid (18β-GA)",
        "circulating_metabolite": "18β-Glycyrrhetinic Acid (Bioactive Aglycone)",
        "circulating_smiles": "CC1(C)C2CCC3(C)C(=O)C=C4C(CCC5(C)C4CCC5(C)C(=O)O)C3(C)C2CCC1O",
        "permeability_raw_papp": "0.08 x 10^-6 cm/s (Hydrophilic Saponin Membrane Barrier)",
        "permeability_active_papp": "14.2 x 10^-6 cm/s (High Lipophilic Tissue Penetration)",
        "in_vivo_pk_note": "Hydrolysis to 18β-GA unlocks potent 11β-HSD2 enzyme inhibition and P-glycoprotein efflux suppression, explaining why Licorice serves as the classical 'Guide Herb' (使药).",
        "phase2_hepatic_fate": "Sulfate and glucuronide conjugation with high plasma protein binding (> 98%)."
    },
    "Ginsenoside Rg1": {
        "botanical_source": "Ginseng (Panax ginseng / Renshen)",
        "ingested_scaffold": "Protopanaxatriol Bis-Glucoside (Ginsenoside Rg1)",
        "ingested_smiles": "CC(=CCCC(C)(C1CCC2(C)C1C(O)CC3C4(C)CC(O)C(OC5OC(CO)C(O)C(O)C5O)C(C)(C)C4CCC23C)OC6OC(CO)C(O)C(O)C6O)C",
        "bacterial_enzyme": "Intestinal Prevotella & Bacteroides β-Glucosidase",
        "metabolic_reaction": "Sequential Deglycosylation: Rg1 ──► Rh1 ──► Protopanaxatriol (PPT)",
        "circulating_metabolite": "Ginsenoside Rh1 / Protopanaxatriol Aglycone",
        "circulating_smiles": "CC(=CCCC(C)(C1CCC2(C)C1C(O)CC3C4(C)CC(O)C(O)C(C)(C)C4CCC23C)O)C",
        "permeability_raw_papp": "0.35 x 10^-6 cm/s (Bulky Hydrophilic Saponin)",
        "permeability_active_papp": "16.8 x 10^-6 cm/s (High Cellular Membrane Permeability)",
        "in_vivo_pk_note": "Deglycosylated ginsenoside metabolites cross the blood-brain barrier (BBB) and bind intracellular nuclear estrogen receptors with 10-fold higher affinity.",
        "phase2_hepatic_fate": "Phase II fatty acid esterification producing sustained neuroprotective intracellular pools."
    },
    "Hesperidin": {
        "botanical_source": "Tangerine Peel (Citrus reticulata / Chenpi)",
        "ingested_scaffold": "Flavanone 7-O-Rutinoside (Hesperidin)",
        "ingested_smiles": "COC1=C(C=CC(=C1)C2CC(=O)C3=C(C2)C=C(C=C3O)OC4OC(COC5OC(C)C(O)C(O)C5O)C(O)C(O)C4O)O",
        "bacterial_enzyme": "Bifidobacterium & Enterococcus Rutinosidase (α-L-Rhamnosidase + β-D-Glucosidase)",
        "metabolic_reaction": "Rutinoside Cleavage ──► Free Hesperetin Aglycone",
        "circulating_metabolite": "Hesperetin (Bioactive Aglycone)",
        "circulating_smiles": "COC1=C(C=CC(=C1)C2CC(=O)C3=C(C2)C=C(C=C3O)O)O",
        "permeability_raw_papp": "0.22 x 10^-6 cm/s (Impermeable Rutinoside)",
        "permeability_active_papp": "24.5 x 10^-6 cm/s (High Intestinal & Endothelial Permeation)",
        "in_vivo_pk_note": "Free hesperetin directly inhibits vascular endothelial cell adhesion molecule expression (VCAM-1, ICAM-1) to reduce atherosclerosis.",
        "phase2_hepatic_fate": "Glucuronidation at 7-OH and 3'-OH positions with systemic bioavailability."
    }
}

def get_microbiome_data(compound_name):
    """
    Returns human gut microbiota biotransformation profile if compound is a prodrug/glycoside.
    """
    comp_clean = compound_name.strip()
    for name, data in MICROBIOME_TRANSFORMATIONS.items():
        if name.lower() in comp_clean.lower() or comp_clean.lower() in name.lower():
            return data
    return None
