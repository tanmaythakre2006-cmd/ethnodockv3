import os
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, FilterCatalog

# 1. Initialize global PAINS catalog for assay interference
_pains_params = FilterCatalog.FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
_pains_catalog = FilterCatalog.FilterCatalog(_pains_params)

# 2. Curated In-Vivo Natural Product Toxicophore & Lethality Catalog
NATURAL_TOXICOPHORES = [
    {
        "name": "Diester Diterpenoid Alkaloid (Aconitine-type Cardiotoxin)",
        "severity": "CRITICAL LETHAL CARDIOTOXICITY",
        "target": "Voltage-Gated Sodium Channels (Nav1.5)",
        "smarts": ["c1ccccc1C(=O)O[#6]", "CC(=O)O[#6]"],
        "smiles": ["CC(=O)OC1C(C2C3(CC(C2(C1O)O)C4(C3CC(C4O)(OC)OC)N(C)CC)OC)OC(=O)C5=CC=CC=C5"],
        "condition": "all",
        "requires_nitrogen": True,
        "description": "Lethal Nav1.5 channel activator (LD50 ~0.1 mg/kg) causing persistent myocardial depolarization, ventricular tachycardia, and fatal cardiac arrest."
    },
    {
        "name": "Strychnos Indole Alkaloid (Strychnine/Brucine Neurotoxin)",
        "severity": "CRITICAL LETHAL NEUROTOXICITY",
        "target": "Spinal Glycine Receptors (GLRA1)",
        "smiles": [
            "O=C1CC2CN3CCC45C6C3CC2C1C4=CC=CC=C5N6C=O",
            "O=C1CC2OCC=C3CN4CCC56C7C4CC3C2C5C1=CC=C7N6",
            "C1=CC=C2C(=C1)N3C(=O)CC4C25CCN6CC=C(CO4)C7CC56C73",
            "COC1=C(C=C2C3=C1N(C(=O)CC4C35CCN6CC=C(CO4)C7CC25C6C7)C)OC"
        ],
        "condition": "any",
        "requires_nitrogen": True,
        "description": "Competitive glycine receptor antagonist causing violent tetanic convulsions, motor neuron spasms, and fatal asphyxiation (Lethal Dose ~30 mg)."
    },
    {
        "name": "Diterpene Triepoxide (Triptolide-type Cytotoxin)",
        "severity": "CRITICAL MULTI-ORGAN CYTOTOXICITY",
        "target": "TFIIH Subunit XPB / Nuclear Pol II",
        "smiles": ["CC(C)C12CC3(C(C1O2)C45C(O4)C6C(O6)C(=O)O5)C7(C3(CC8(C7O8)C)O)C"],
        "condition": "any",
        "requires_nitrogen": False,
        "description": "Irreversibly inhibits XPB helicase causing total transcriptional arrest, acute hepatorenal collapse, and severe bone marrow suppression."
    },
    {
        "name": "Tropane Anticholinergic Alkaloid (Scopolamine/Hyoscyamine)",
        "severity": "CRITICAL ANTICHOLINERGIC NEUROTOXICITY",
        "target": "Muscarinic Acetylcholine Receptors (M1-M5)",
        "smiles": [
            "CN1C2CC(CC1C3C2O3)OC(=O)C(CO)C4=CC=CC=C4",
            "CN1C2CCC1CC(C2)OC(=O)C(CO)c3ccccc3"
        ],
        "condition": "any",
        "requires_nitrogen": True,
        "description": "Potent central and peripheral muscarinic antagonist inducing severe anticholinergic syndrome, hallucinations, hyperthermia, and coma."
    },
    {
        "name": "Aristolochic Acid Core (Nephrotoxic Carcinogen)",
        "severity": "CRITICAL CARCINOGENICITY & NEPHROTOXICITY",
        "target": "Genomic DNA / Renal Proximal Tubules",
        "smiles": [
            "O=C(O)C1=CC2=C(C=C1[N+](=O)[O-])C3=CC(OC)=C4OCOC4=C3C=C2",
            "O=C(O)c1cc2c(cc1[N+](=O)[O-])c3ccc4OCOc4c3cc2"
        ],
        "condition": "any",
        "requires_nitrogen": True,
        "description": "Forms irreversible covalent DNA adducts (dA:T -> T:A) causing rapid bilateral renal interstitial fibrosis and fatal urothelial carcinoma."
    },
    {
        "name": "Pyrrolizidine Alkaloid (Senecionine-type Hepatotoxin)",
        "severity": "CRITICAL HEPATOTOXICITY & GENOTOXICITY",
        "target": "Hepatic Sinusoidal Endothelium / DNA",
        "smarts": ["C1CCN2CC=C(C21)COC(=O)"],
        "condition": "any",
        "requires_nitrogen": True,
        "description": "Metabolized to pyrrolic esters causing irreversible hepatic sinusoidal obstruction syndrome and liver failure."
    }
]

def screen_natural_toxicophores(mol):
    """
    Screens small molecules against lethal and organ-damaging natural product toxicophore patterns
    using both canonical graph comparison and substructure SMARTS filters.
    """
    if not mol:
        return []
    
    alerts = []
    has_n = any(atom.GetSymbol() == "N" for atom in mol.GetAtoms())
    mol_canon = Chem.MolToSmiles(mol)
    
    for tox in NATURAL_TOXICOPHORES:
        if tox.get("requires_nitrogen") and not has_n:
            continue
            
        matched = False
        # 1. Check Canonical SMILES exact/tautomer match
        if "smiles" in tox:
            for s in tox["smiles"]:
                ref_mol = Chem.MolFromSmiles(s)
                if ref_mol:
                    if Chem.MolToSmiles(ref_mol) == mol_canon:
                        matched = True
                        break
                    if mol.HasSubstructMatch(ref_mol):
                        matched = True
                        break
                        
        # 2. Check SMARTS substructure patterns
        if not matched and "smarts" in tox:
            smarts_matches = []
            for s in tox["smarts"]:
                q = Chem.MolFromSmarts(s)
                if q is not None and mol.HasSubstructMatch(q):
                    smarts_matches.append(True)
                else:
                    smarts_matches.append(False)
            if tox["condition"] == "all" and all(smarts_matches) and len(smarts_matches) == len(tox["smarts"]):
                matched = True
            elif tox["condition"] == "any" and any(smarts_matches):
                matched = True
                
        if matched:
            alerts.append(tox)
            
    return alerts

def get_admet_profile(smiles):
    """
    Computes comprehensive ADMET, drug-likeness, toxicological safety parameters,
    Lipinski Rule of Five, Veber oral bioavailability, PAINS, and Natural Product Toxicophore screens.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    aromatic_rings = Lipinski.NumAromaticRings(mol)
    qed_score = QED.qed(mol)

    # 1. Lipinski Violations
    lipinski_violations = []
    if mw > 500:
        lipinski_violations.append("MW > 500 Da")
    if logp > 5.0:
        lipinski_violations.append("LogP > 5.0")
    if hbd > 5:
        lipinski_violations.append("HBD > 5")
    if hba > 10:
        lipinski_violations.append("HBA > 10")

    # 2. Veber Rules (Oral Bioavailability)
    veber_violations = []
    if rotb > 10:
        veber_violations.append("Rotatable Bonds > 10")
    if tpsa > 140:
        veber_violations.append("TPSA > 140 Å²")

    # 3. Assay Interference (PAINS) Screening
    alert_match = _pains_catalog.GetFirstMatch(mol)
    pains_alert = alert_match.GetDescription() if alert_match else "Pass (0 PAINS Alerts)"
    is_pains_clean = alert_match is None

    # 4. In-Vivo Natural Toxicophore & Lethality Screening
    tox_alerts = screen_natural_toxicophores(mol)
    is_toxic_hazardous = len(tox_alerts) > 0
    
    if is_toxic_hazardous:
        top_tox = tox_alerts[0]
        structure_alert_summary = f"CRITICAL HAZARD: {top_tox['name']} ({top_tox['severity']})"
        safety_status = "CRITICAL HAZARD / FAIL"
        safety_badge = "🔴 FAIL"
    elif not is_pains_clean:
        structure_alert_summary = f"PAINS ALERT: {pains_alert}"
        safety_status = "ASSAY INTERFERENCE ALERT"
        safety_badge = "🟡 CAUTION"
    else:
        structure_alert_summary = "PASSED (0 alerts, Verified Clean)"
        safety_status = "PASS (Clean Scaffold)"
        safety_badge = "🟢 PASS"

    # 5. Blood-Brain Barrier (BBB) Estimation Heuristic
    bbb_permeable = (mw < 450) and (1.0 <= logp <= 4.5) and (tpsa < 90) and (hbd <= 3)

    return {
        "SMILES": smiles,
        "Molecular Weight": round(mw, 2),
        "LogP": round(logp, 2),
        "TPSA (Å²)": round(tpsa, 2),
        "H-Bond Donors": hbd,
        "H-Bond Acceptors": hba,
        "Rotatable Bonds": rotb,
        "Aromatic Rings": aromatic_rings,
        "QED Drug-Likeness": round(qed_score, 3),
        "Lipinski Violations": len(lipinski_violations),
        "Lipinski Details": lipinski_violations if lipinski_violations else ["None (Ro5 Compliant)"],
        "Veber Pass": len(veber_violations) == 0,
        "PAINS Screen": pains_alert,
        "PAINS Clean": is_pains_clean,
        "Natural Toxicophore Alerts": tox_alerts,
        "Is Toxicologically Hazardous": is_toxic_hazardous,
        "Structure Alert Screen": structure_alert_summary,
        "Safety Status": safety_status,
        "Safety Badge": safety_badge,
        "Predicted BBB Permeability": "High" if bbb_permeable else "Moderate/Low"
    }
