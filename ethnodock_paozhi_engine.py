import os
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

PAOZHI_LIBRARY = {
    "Monkshood (Fuzi)": {
        "botanical": "Aconitum carmichaelii",
        "family": "Ranunculaceae",
        "raw_name": "Raw Fuzi (生附子)",
        "raw_compound": "Aconitine",
        "raw_smiles": "CC(=O)OC1C(C2C3(CC(C2(C1O)O)C4(C3CC(C4O)(OC)OC)N(C)CC)OC)OC(=O)C5=CC=CC=C5",
        "raw_class": "Highly Toxic Diester Diterpenoid Alkaloid",
        "raw_toxicity_warning": "CRITICAL TOXICITY (LD50 ~0.1 mg/kg): Activates voltage-gated Na+ channels causing fatal cardiac arrhythmias and ventricular fibrillation.",
        "paozhi_method": "Water soaking, high-pressure steaming, and baking (水浸煮蒸制 / 黑顺片/白附片)",
        "reaction_type": "Thermal Hydrolysis & C-8 De-esterification",
        "reaction_equation": "Aconitine (Diester) + H2O ──[Steam / 100°C]──► Benzoylaconine (Monoester) + Acetic Acid",
        "processed_name": "Processed Fuzi (熟附子 / 附片)",
        "processed_compound": "Benzoylaconine",
        "processed_smiles": "COC1C2C(C(C1O)O)(C3(CC(C2OC(=O)C4=CC=CC=C4)C5(C3CC(C5O)(OC)OC)N(C)CC)O)O",
        "processed_class": "Detoxified Monoester Diterpenoid Alkaloid",
        "detox_benefit": "Toxicity drops by over 100-fold while preserving central analgesic, anti-inflammatory, and peripheral microcirculation warming properties.",
        "target_protein": "Voltage-Gated Sodium Channel (Nav1.5)",
        "gene_symbol": "SCN5A",
        "uniprot_id": "Q14524",
        "pdb_id": "6LQA"
    },
    "Chinese Rhubarb (Dahuang)": {
        "botanical": "Rheum officinale",
        "family": "Polygonaceae",
        "raw_name": "Raw Dahuang (生大黄)",
        "raw_compound": "Sennoside A",
        "raw_smiles": "O=C1C2=C(C(=O)C3=C1C=CC(=C3)C4=CC=C5C(=C4)C(=O)C6=C(C5=O)C=C(O)C=C6)C=C(O)C=C2",
        "raw_class": "Dianthrone Glycoside (Drastic Purgative)",
        "raw_toxicity_warning": "HIGH IRRITANT (Harsh Purgative): Stimulates violent colonic contractions, severe abdominal cramping, and fluid loss.",
        "paozhi_method": "Rice wine steaming and slow-heat stewing (黄酒蒸制 / 熟大黄/酒大黄)",
        "reaction_type": "Thermal Deglycosylation & Anthraquinone Cleavage",
        "reaction_equation": "Sennoside A (Dianthrone Glycoside) ──[Wine Steam / Heat]──► Free Emodin + Rhein (Aglycones)",
        "processed_name": "Processed Wine Rhubarb (熟大黄 / 酒大黄)",
        "processed_compound": "Emodin",
        "processed_smiles": "CC1=CC2=C(C(=C1)O)C(=O)C3=C(C2=O)C=C(C=C3O)O",
        "processed_class": "Free Anthraquinone Aglycone",
        "detox_benefit": "Eliminates violent bowel irritation and purgative cramping while dramatically increasing anti-inflammatory and vascular endothelial protective properties.",
        "target_protein": "Casein Kinase 2 (CK2)",
        "gene_symbol": "CSNK2A1",
        "uniprot_id": "P68400",
        "pdb_id": "2ZJW"
    },
    "Fleeceflower Root (Heshouwu)": {
        "botanical": "Polygonum multiflorum",
        "family": "Polygonaceae",
        "raw_name": "Raw Heshouwu (生首乌)",
        "raw_compound": "Emodin-8-O-glucoside",
        "raw_smiles": "CC1=CC2=C(C(=C1)O)C(=O)C3=C(C2=O)C=C(C=C3O)OC4C(C(C(C(O4)CO)O)O)O",
        "raw_class": "Anthraquinone Glucoside",
        "raw_toxicity_warning": "MODERATE HEPATOTOXIC BURDEN: Unprocessed anthraquinone conjugates can burden cytochrome CYP450 pathways and cause diarrhea.",
        "paozhi_method": "Black bean decoction steaming 9 times and sun-drying 9 times (黑豆汁九蒸九晒 / 制首乌)",
        "reaction_type": "Polyphenolic Glycoside Conversion & Stilbene Condensation",
        "reaction_equation": "Emodin Glucoside ──[Black Bean Steam 9x]──► 2,3,5,4'-THSG (Active Stilbene Glucoside)",
        "processed_name": "Prepared Heshouwu (制首乌)",
        "processed_compound": "2,3,5,4'-Tetrahydroxystilbene-2-O-glucoside (THSG)",
        "processed_smiles": "C1=CC(=CC=C1C=CC2=C(C(=C(C=C2)O)O)OC3C(C(C(C(O3)CO)O)O)O)O",
        "processed_class": "Bioactive Polyhydroxy Stilbene Glycoside",
        "detox_benefit": "Completely removes laxative/hepatotoxic burden, transforming the root into a premier longevity, lipid-clearing, and hair-darkening tonic.",
        "target_protein": "HMG-CoA Reductase",
        "gene_symbol": "HMGCR",
        "uniprot_id": "P04035",
        "pdb_id": "1HW9"
    },
    "Strychnos Seed (Maqianzi)": {
        "botanical": "Strychnos nux-vomica",
        "family": "Loganiaceae",
        "raw_name": "Raw Maqianzi (生马钱子)",
        "raw_compound": "Strychnine",
        "raw_smiles": "O=C1CC2CN3CCC45C6C3CC2C1C4=CC=CC=C5N6C=O",
        "raw_class": "Indole Alkaloid (Deadly Convulsive Neurotoxin)",
        "raw_toxicity_warning": "CRITICAL LETHAL NEUROTOXICITY (Lethal Dose ~30 mg): Potent competitive glycine receptor antagonist causing violent tetanic convulsions and fatal asphyxiation.",
        "paozhi_method": "Hot sand roasting in sesame oil until swollen and brown (热砂烫制 / 油炸马钱子)",
        "reaction_type": "Thermal Isomerization & N-Oxidation",
        "reaction_equation": "Strychnine ──[Sand Roasting / 230°C]──► Isostrychnine + Strychnine N-Oxide",
        "processed_name": "Processed Roasted Maqianzi (制马钱子)",
        "processed_compound": "Isostrychnine",
        "processed_smiles": "C1=CC=C2C(=C1)N(C=O)C3C24CCC5C6CN(C3)CC=C6CC(=O)C54",
        "processed_class": "Thermally Isomerized Seco-Alkaloid",
        "detox_benefit": "Reduces central nervous system convulsive toxicity by 85% while retaining potent peripheral analgesic and anti-arthritic activity.",
        "target_protein": "Glycine Receptor Alpha-1",
        "gene_symbol": "GLRA1",
        "uniprot_id": "P23415",
        "pdb_id": "5TIN"
    },
    "Licorice Root (Gancao)": {
        "botanical": "Glycyrrhiza uralensis",
        "family": "Fabaceae",
        "raw_name": "Raw Gancao (生甘草)",
        "raw_compound": "Glycyrrhizic Acid",
        "raw_smiles": "CC1(C2CCC3(C(C2(CCC1C4C(C(C(C(O4)C(=O)O)O)O)OC5C(C(C(C(O5)C(=O)O)O)O)O)C)C(=O)C=C6C3(CCC7(C6CC(CC7)(C)C)C(=O)O)C)C)C",
        "raw_class": "Triterpenoid Saponin",
        "raw_toxicity_warning": "MILD PSEUDOALDOSTERONISM: High doses of raw saponins cause sodium retention and potassium loss.",
        "paozhi_method": "Honey-frying and gentle baking (炼蜜炙制 / 炙甘草)",
        "reaction_type": "Flavonoid Deglycosylation & Glycation Saponification",
        "reaction_equation": "Glycyrrhizin + Honey ──[Honey Bake / 120°C]──► Liquiritigenin (Bioactive Flavanone)",
        "processed_name": "Honey-Roasted Gancao (炙甘草)",
        "processed_compound": "Liquiritigenin",
        "processed_smiles": "O=C1CC(C2=CC=C(O)C=C2)OC3=C1C=CC(O)=C3",
        "processed_class": "Flavanone Aglycone",
        "detox_benefit": "Shifts pharmacology from superficial clearing of throat heat to deep tonification of Spleen Qi and stabilization of cardiac arrhythmia.",
        "target_protein": "Estrogen Receptor Beta",
        "gene_symbol": "ESR2",
        "uniprot_id": "Q92731",
        "pdb_id": "1QKM"
    }
}

def has_paozhi(name_str):
    """Checks whether an herb name matches any Paozhi-capable entry."""
    for key in PAOZHI_LIBRARY:
        if key.lower() in name_str.lower() or PAOZHI_LIBRARY[key]["botanical"].lower() in name_str.lower():
            return key
    return None

def get_paozhi_data(key):
    """Retrieves the full Paozhi processing profile."""
    return PAOZHI_LIBRARY.get(key, None)
